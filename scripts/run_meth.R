suppressMessages({
    library(dplyr)
    library(jsonlite)
    library(purrr)
    library(SingleCellExperiment)
})

sim <- readRDS(args$sim)
meth_pars <- as.list(fromJSON(args$meth_pars))
run_pars <- as.list(fromJSON(args$run_pars))

set.seed(run_pars$seed + as.numeric(wcs$j))

# subset clusters & samples
kids <- levels(sim$cluster_id)
sids <- levels(sim$sample_id)
m <- match(sids, sim$sample_id)
gids <- sim$group_id[m]

if (wcs$k != "x") kids <- sample(kids, wcs$k)
if (wcs$s != "x") sids <- sapply(split(sids, gids), sample, wcs$s)

sim <- .filter_sce(sim, kids, sids)

# subset genes & cells
gs <- rownames(sim)
cs <- colnames(sim)

if (wcs$g != "x") 
    gs <- sample(gs, min(nrow(sim), as.numeric(wcs$g)))

if (wcs$c != "x") {
    # table of cells for each sample
    tab = with(colData(sim), table(sample_id, cluster_id))
    tab = as.data.frame(tab)
    rownames(tab) = with(tab, paste(cluster_id, sample_id, sep='.'))

    # target mean number of cells
    ncells.target = as.numeric(wcs$c)
    nsubj = length(unique(tab$sample_id))

    cs <- split(cs, list(sim$cluster_id, sim$sample_id))
    cs <- unlist(lapply(names(cs), function(id){ 
        u = cs[[id]]

        # total cells in this cluster
        ntotal = sum(tab[tab$cluster_id == tab[id,'cluster_id'],]$Freq)
        # mean number of observed cells per sample in this cluster
        ncells.mean = ntotal / nsubj

        # target size scaled by the fraction of total cells in this subject
        ncells = ncells.target * length(u) / ncells.mean

        sample(u, min(length(u), ncells))
    }))
}

# run method & write results to .rds
source(fun <- args$fun)
fun <- gsub("(.R)", "", basename(fun))
res <- get(fun)(sim[gs, cs], meth_pars)

if (!inherits(res$tbl, "error")) {
    # add metadata
    gi <- metadata(sim)$gene_info %>% 
        dplyr::filter(gene %in% gs) %>% 
        dplyr::mutate_at("cluster_id", as.character) %>% 
        dplyr::select(-"logFC") %>% 
        dplyr::mutate(., sim_lfc = eval(parse(text = ifelse(
            .$sim_mean.B == 0, "0", "log2(sim_mean.B/sim_mean.A)"))))
    
    res$tbl <- left_join(gi, res$tbl, by = c("gene", "cluster_id")) %>%
        {if ("logFC" %in% names(.)) 
            dplyr::rename(., est_lfc = logFC) else .} %>%
        dplyr::mutate(did = wcs$did, sid = wcs$sid, mid = wcs$mid, 
            i = wcs$i, j = wcs$j, g = wcs$g, c = wcs$c, k = wcs$k, s = wcs$s,
            is_de = as.integer(!category %in% c("ee", "ep")))
}

saveRDS(res, args$res)