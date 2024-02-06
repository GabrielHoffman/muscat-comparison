suppressMessages({
    library(dplyr)
    library(jsonlite)
    library(muscat)
    library(scater)
    library(sctransform)
    library(SingleCellExperiment) 
})

# load data & simulation parameters
sce <- readRDS(args$sce)
sim_pars <- fromJSON(args$sim_pars)
set.seed(sim_pars$seed + as.numeric(wcs$i))

# Simulate more cells than needed
# Then downsample later
k_scaling = 2

assignInNamespace( ".check_args_simData", function(u)
    return(list(nk = u$nk, ns = u$ns)), ns="muscat")


# simulate effect size heterogeneity
# with a normal offset added to the logFC 
.nb.replace <- function(cs, d, m, lfc = NULL, f = 1) {
    n_gs <- length(d)
    n_cs <- length(cs)
    if (is.null(lfc)) {
        lfc <- rep(0, n_gs)
    } else {
        lfc[lfc < 0] <- 0
    }
    names(lfc) = names(d)

    # effect size heterogeneity for non-zero effects
    if( any(lfc != 0) ){
        i = lfc != 0
        lfc[i] <- lfc[i] + rnorm(length(lfc[i]), 0, .1) 
        # lfc[lfc < 0] <- 0
    }

    fc <- f * (2 ^ lfc)
    # cat fcs for all cells
    fc <- rep(fc, each = n_cs)

    # cell-level heterogeneity when fc != 1
    # i.e. lfc != 0
    i = (fc != 1)
    s = sqrt(rgamma(1,1, 1000))
    fc[i] = 2^(log2(fc[i]) + rnorm(sum(i), 0, s))
    fc = pmax(1, fc)
    # plot(log(fc))

    ds <- rep(1/d, each = n_cs)
    ms <- c(t(m[, cs])) * fc 
    y <- rnbinom(n_gs * n_cs, size = ds, mu = ms)
    y <- matrix(y, byrow = TRUE, 
        nrow = n_gs, ncol = n_cs, 
        dimnames = list(names(d), cs))
    ms <- split(ms, rep(seq_len(nrow(m)), each = n_cs))
    list(counts = y, means = ms)
}

assignInNamespace(".nb", .nb.replace, ns="muscat")


sim <- simData(sce, 
    paired = FALSE, lfc = .3 ,
    force = TRUE,
    ng = nrow(sce), nc = sim_pars$nc * k_scaling,
    ns = sim_pars$ns, nk = sim_pars$nk,
    p_dd = sim_pars$p_dd, probs = sim_pars$probs)

# filter genes
sim <- sim[rowSums(counts(sim) > 0) >= 10, ]

tab = table(sim$sample_id, sim$cluster_id)

if( k_scaling > 1){

    df_grid = expand.grid(sid = levels(sim$sample_id), 
                            cid = levels(sim$cluster_id))

    keep = lapply( seq(nrow(df_grid)), function(i){

        keep = which( sim$sample_id == df_grid$sid[i] & sim$cluster_id == df_grid$cid[i])

        target = tab[df_grid$sid[i],df_grid$cid[i]]/k_scaling

        # sample cell counts from Negative Binomial 
        # Poisson if theta = Inf
        # additive overdispersion is mu^2/theta
        # variance is 'a' times the Poisson variance 
        # a = 10
        # theta = target / (a-1)
        theta = 3
        ncells = MASS::rnegbin(1, mu=target, theta=theta)
         
        ncells = max(5, ncells)

        if( ncells < length(keep)){
            keep <- sample(keep, ncells)
        }
        keep
    })
    keep = sort(unlist(keep))

    # Subsample
    sim.tmp = sim[,keep]

    # rename cells
    colnames(sim.tmp) = paste0("cell", seq(ncol(sim.tmp)))

    # filter genes
    # sim2 <- sim2[rowSums(counts(sim2) > 0) >= 50, ]

    # set number of cells
    metadata(sim.tmp)$n_cells = table(sim.tmp$sample_id)
    metadata(sim.tmp)$args$nc = sim.tmp$nc

    sim = sim.tmp
}

# subsample genes to get correct number
sim <- sim[sample(nrow(sim), min(nrow(sim), sim_pars$ng)), ]

gi <- metadata(sim)$gene_info 
gi <- dplyr::filter(gi, gene %in% rownames(sim))
metadata(sim)$gene_info <- gi

sim <- computeLibraryFactors(sim)
sim <- logNormCounts(sim)
assays(sim)$cpm <- calculateCPM(sim)

# vst_values <- suppressWarnings(sctransform::vst(counts(sim))$y)
# assays(sim)$vstresiduals <- vst_values

saveRDS(sim, args$sim)