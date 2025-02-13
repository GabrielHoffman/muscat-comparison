suppressMessages({
    library(limma)
    library(muscat)
    library(scater)
    library(sctransform)
    library(SingleCellExperiment)
})

apply_pb <- function(sce, pars, ds_only = TRUE) {
    t1 <- system.time({
        a <- pars$assay
        if (!ds_only) 
            assay(sce, a) <- switch(a, 
                counts = counts(sce),
                cpm = calculateCPM(counts(sce)),
                logcounts = normalizeCounts(computeLibraryFactors(sce)),
                vstresiduals = vst(counts(sce), show_progress = FALSE)$y)
        pb <- aggregateData(sce, a, fun = pars$fun, scale = pars$scale)
    })[[3]]
    t2 <- system.time({
        if( pars$method == "dreamlet" ){

            library(dreamlet)

            vobj <- processAssays(pb, ~ 1, verbose=FALSE, min.count=3)
            fit <- dreamlet(vobj, ~ group_id, verbose=FALSE )
            tab <- topTable(fit, coef='group_idB', number=Inf)

            tab2 = with(tab, data.frame(gene = ID, cluster_id = assay, logFC, AveExpr, t, p_val=P.Value, B, contrast='B'))

            tab2$p_adj.glb = p.adjust(tab2$p_val, "BH")
            tab2$p_adj.loc = rep(NA, nrow(tab2))

            for( CT in unique(tab2$cluster_id) ){
                idx = which(tab2$cluster_id==CT)
                tab2$p_adj.loc[idx] = p.adjust(tab2$p_val[idx], "BH")
            }
            res = tab2    
        }else if( pars$method == "DESeq2" ){
            library(DESeq2)

            res = lapply(assayNames(pb), function(CT){
                dds = DESeqDataSetFromMatrix( assay(pb, 1), colData(pb), ~ group_id)
                dds = DESeq(dds)
                res = results(dds)

                data.frame(gene = rownames(res),
                    cluster_id = CT,
                    logFC = res$log2FoldChange,
                    AveExpr = res$baseMean, 
                    t = res$stat, 
                    p_val = res$pvalue, 
                    p_adj.loc = p.adjust(res$pvalue, "BH"),
                    contrast='B')
            })
            res = do.call(rbind, res)
            res$p_adj.glb = p.adjust(res$p_val, "BH")

        }else{
            res <- tryCatch(
                do.call(pbDS, c(
                    list(pb = pb, filter = "none", verbose = FALSE),
                    pars[names(pars) %in% names(formals(pbDS))])),
                error = function(e) e)
            if (!inherits(res, "error"))
                res <- dplyr::bind_rows(res$table[[1]])
        }
    })[[3]]
    list(rt = c(t1, t2), tbl = res)
}

