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

sim <- simData(sce, 
    paired = FALSE, lfc = .5 ,
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