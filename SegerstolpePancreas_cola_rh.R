setwd("/omics/groups/OE0246/internal/guz/cola_hc/examples/SegerstolpePancreas")
library(cola)

process_counts = function(data, column = NULL) {
    mat = assays(data)$counts
    mat = as.matrix(mat)
    s = colSums(mat)
    fa = s/mean(s)
    for(i in 1:ncol(mat)) mat[, i]/fa[i]
    mat = adjust_matrix(log2(mat + 1))

    anno = NULL
    if(!is.null(column)) {
        anno = colData(data)
           anno = as.data.frame(anno)
        anno = anno[, column, drop = FALSE]
    }

    list(mat = mat, anno = anno)
}

library(scRNAseq)
data = readRDS('/omics/groups/OE0246/internal/guz/cola_hc/examples/SegerstolpePancreas/SegerstolpePancreas_data.rds')
lt = process_counts(data, c("cell.type", "disease", "age"))
rh = hierarchical_partition(lt$mat, subset = 500, cores = 4, anno = lt$anno)
saveRDS(rh, file = "SegerstolpePancreas_cola_rh.rds")

cola_report(rh, output = "SegerstolpePancreas_cola_rh_report", title = "cola Report for Hierarchical Partitioning - 'SegerstolpePancreas'")
