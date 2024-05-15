suppressMessages({
    library(dplyr)
    library(jsonlite)
    library(muscat)
    library(scater)
    library(sctransform)
    library(SingleCellExperiment) 
    library(DirichletReg)
})

# load data & simulation parameters
sce <- readRDS(args$sce)
sim_pars <- fromJSON(args$sim_pars)
set.seed(sim_pars$seed + as.numeric(wcs$i))


#--------------------------------
# Remove check on lfc
assignInNamespace( ".check_args_simData", function(u)
    return(list(nk = u$nk, ns = u$ns)), ns="muscat")

# Simulate more cells than needed
# Then downsample later
k_scaling = 5

# sim data
sim <- simData(sce, 
    paired = FALSE, lfc = 0.5,
    ng = nrow(sce), nc = k_scaling*sim_pars$nc,
    ns = sim_pars$ns, nk = sim_pars$nk,
    p_dd = sim_pars$p_dd, probs = sim_pars$probs, force=TRUE)

# count cells for each sample
tab = table(sim$sample_id, sim$cluster_id)

# get the number of cells to sample
alpha = exp(runif(nrow(tab), 0, 3))
prob = gtools::rdirichlet(1, alpha = alpha)
counts1 = rmultinom(1, sim_pars$nc / sim_pars$nk, prob)

alpha = exp(runif(nrow(tab), 0, 3))
prob = gtools::rdirichlet(1, alpha = alpha)
counts2 = rmultinom(1, sim_pars$nc / sim_pars$nk, prob)

# Subsample cells
df_grid = expand.grid(sid = levels(sim$sample_id), 
                        cid = levels(sim$cluster_id))
df_grid$counts = c(counts1, counts2)

keep = sapply( seq(nrow(df_grid)), function(i){

    keep = which( sim$sample_id == df_grid$sid[i] & sim$cluster_id == df_grid$cid[i])

    sample(keep, min(length(keep), df_grid$counts[i]))
})
keep = sort(unlist(keep))

# sim = sim[,keep]
#--------------------------------

# filter after
sim <- sim[rowSums(counts(sim) > 0) >= 10, ]
sim <- sim[sample(nrow(sim), min(nrow(sim), sim_pars$ng)), ]

gi <- metadata(sim)$gene_info 
gi <- dplyr::filter(gi, gene %in% rownames(sim))
metadata(sim)$gene_info <- gi

sim <- computeLibraryFactors(sim)
sim <- logNormCounts(sim)
assays(sim)$cpm <- calculateCPM(sim)
assays(sim)$vstresiduals <- suppressWarnings(
    vst(counts(sim), show_progress = FALSE)$y)

saveRDS(sim, args$sim)