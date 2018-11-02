#' @name BuettnerFlorian
#' @title test dataset for SIMLR
#' @description example dataset to test SIMLR from the work by Buettner, Florian, et al.
#' @docType data
#' @usage data(BuettnerFlorian)
#' @format gene expression measurements of individual cells
#' @source Buettner, Florian, et al. "Computational analysis of cell-to-cell heterogeneity in single-cell RNA-sequencing data reveals hidden subpopulations of cells." Nature biotechnology 33.2 (2015): 155-160.
#' @return list of 6:
#'   \enumerate{
#'      \item in_X = input dataset as an (m x n) gene expression measurements of individual cells,
#'      \item n_clust = number of clusters (number of distinct true labels),
#'      \item true_labs = ground true of cluster assignments for each of the n_clust clusters,
#'      \item seed = seed used to compute the results for the example,
#'      \item results = result by SIMLR for the inputs defined as described,
#'      \item nmi = normalized mutual information as a measure of the inferred clusters compared to the true labels
#'   }
NULL


#' @name ZeiselAmit
#' @title test dataset for SIMLR large scale
#' @description example dataset to test SIMLR large scale. This is a reduced version of the dataset from the work by Zeisel, Amit, et al.
#' @docType data
#' @usage data(ZeiselAmit)
#' @format gene expression measurements of individual cells
#' @source Zeisel, Amit, et al. "Cell types in the mouse cortex and hippocampus revealed by single-cell RNA-seq." Science 347.6226 (2015): 1138-1142.
#' @return list of 6:
#'   \enumerate{
#'     \item in_X = input dataset as an (m x n) gene expression measurements of individual cells,
#'     \item n_clust = number of clusters (number of distinct true labels),
#'     \item true_labs = ground true of cluster assignments for each of the n_clust clusters,
#'     \item seed = seed used to compute the results for the example,
#'     \item results = result by SIMLR for the inputs defined as described,
#'     \item nmi = normalized mutual information as a measure of the inferred clusters compared to the true labels
#'   }
NULL


#' @name Buettner
#' @title Data from Buettner et al.
#' @description Data from Buettner et al.
#' @docType data
#' @usage data(Buettner)
#' @format gene expression measurements of individual cells
#' @return list of 6:
#'   \enumerate{
#'     \item in_X = input dataset as an (m x n) gene expression measurements of individual cells,
#'     \item n_clust = number of clusters (number of distinct true labels),
#'     \item true_labs = ground true of cluster assignments for each of the n_clust clusters,
#'     \item seed = seed used to compute the results for the example,
#'     \item results = result by SIMLR for the inputs defined as described,
#'     \item nmi = normalized mutual information as a measure of the inferred clusters compared to the true labels
#'   }
NULL


#' @name Kolodziejczyk
#' @title Data from Kolodziejczyk et al.
#' @description Data from Kolodziejczyk et al.
#' @docType data
#' @usage data(Kolodziejczyk)
#' @format gene expression measurements of individual cells
#' @return list of 6:
#'   \enumerate{
#'     \item in_X = input dataset as an (m x n) gene expression measurements of individual cells,
#'     \item n_clust = number of clusters (number of distinct true labels),
#'     \item true_labs = ground true of cluster assignments for each of the n_clust clusters,
#'     \item seed = seed used to compute the results for the example,
#'     \item results = result by SIMLR for the inputs defined as described,
#'     \item nmi = normalized mutual information as a measure of the inferred clusters compared to the true labels
#'   }
NULL


#' @name Pollen
#' @title Data from Pollen et al.
#' @description Data from Pollen et al.
#' @docType data
#' @usage data(Pollen)
#' @format gene expression measurements of individual cells
#' @return list of 6:
#'   \enumerate{
#'     \item in_X = input dataset as an (m x n) gene expression measurements of individual cells,
#'     \item n_clust = number of clusters (number of distinct true labels),
#'     \item true_labs = ground true of cluster assignments for each of the n_clust clusters,
#'     \item seed = seed used to compute the results for the example,
#'     \item results = result by SIMLR for the inputs defined as described,
#'     \item nmi = normalized mutual information as a measure of the inferred clusters compared to the true labels
#'   }
NULL


#' @name Usoskin
#' @title Data from Usoskin et al.
#' @description Data from Usoskin et al.
#' @docType data
#' @usage data(Usoskin)
#' @format gene expression measurements of individual cells
#' @return list of 6:
#'   \enumerate{
#'     \item in_X = input dataset as an (m x n) gene expression measurements of individual cells,
#'     \item n_clust = number of clusters (number of distinct true labels),
#'     \item true_labs = ground true of cluster assignments for each of the n_clust clusters,
#'     \item seed = seed used to compute the results for the example,
#'     \item results = result by SIMLR for the inputs defined as described,
#'     \item nmi = normalized mutual information as a measure of the inferred clusters compared to the true labels
#'   }
NULL


#' @name Zeisel
#' @title Data from Zeisel et al.
#' @description Data from Zeisel et al.
#' @docType data
#' @usage data(Zeisel)
#' @format gene expression measurements of individual cells
#' @return list of 6:
#'   \enumerate{
#'     \item in_X = input dataset as an (m x n) gene expression measurements of individual cells,
#'     \item n_clust = number of clusters (number of distinct true labels),
#'     \item true_labs = ground true of cluster assignments for each of the n_clust clusters,
#'     \item seed = seed used to compute the results for the example,
#'     \item results = result by SIMLR for the inputs defined as described,
#'     \item nmi = normalized mutual information as a measure of the inferred clusters compared to the true labels
#'   }
NULL


#' @name GliomasReduced
#' @title test dataset for CIMLR
#' @description example dataset to test CIMLR. This is a reduced version of the dataset from the work by The Cancer Genome Atlas Research Network.
#' @docType data
#' @usage data(GliomasReduced)
#' @format multi-omic data of cancer patients
#' @source Cancer Genome Atlas Research Network. "Comprehensive, integrative genomic analysis of diffuse lower-grade gliomas." New England Journal of Medicine 372.26 (2015): 2481-2498.
#' @return list of 1 element:
#'   \enumerate{
#'     \item in_X = input dataset as a list of 4 (reduced) multi-omic data each of which is an (m x n) measurements of cancer patients
#'   }
NULL


#' @name GliomasFull
#' @title test dataset for CIMLR
#' @description This is a dataset from the work by The Cancer Genome Atlas Research Network.
#' @docType data
#' @usage data(GliomasFull)
#' @format multi-omic data of cancer patients
#' @source Cancer Genome Atlas Research Network. "Comprehensive, integrative genomic analysis of diffuse lower-grade gliomas." New England Journal of Medicine 372.26 (2015): 2481-2498.
#' @return list of 1 element:
#'   \enumerate{
#'     \item in_X = input dataset as a list of 4 multi-omic data each of which is an (m x n) measurements of cancer patients
#'   }
NULL
