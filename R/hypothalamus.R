#' MERFISH mouse hypothalamus dataset from Moffitt et al., 2018
#' @description Obtain the MERFISH mouse hypothalamic preoptic region dataset
#' from Moffitt et al., 2018
#' @details The hypothalamus controls essential social behaviors and homeostatic
#' functions. However, the cellular architecture of hypothalamic nuclei, including
#' the molecular identity, spatial organization, and function of distinct cell
#' types, is not well understood.
#'
#' Moffitt et al., 2018, developed an imaging-based cell type identification and
#' mapping method and combined it with single-cell RNA-sequencing to create a
#' molecularly annotated and spatially resolved cell atlas of the mouse hypothalamic
#' preoptic region.  
#'
#' The MERFISH measurements were obtained via combinatorial smFISH imaging for
#' 135 genes (main experiment named \code{"smFISH"}), followed by sequential rounds
#' of non-combinatorial seqFISH for 20 additional genes (stored as an `altExp`
#' named \code{"seqFISH"}). These genes were considered neuronal markers and important
#' for discriminating neuronal cell populations. For behavioral measurements,
#' cFos was added to the set of genes measured with sequential rounds of FISH.
#'   
#' The barcoding scheme contained 140 possible barcodes; 135 of them were used to
#' code the RNAs of the genes assayed via combinatorial smFISH; 5 of these barcodes
#' were left unassigned ("blank"), providing a direct measure of the false-positive
#' rate in MERFISH. Measurements for these 5 blank barcodes is stored in an 
#' `altExp` named \code{"blank"}.
#'
#' @param center.coords logical. Should spatial x- and y-coordinates be centered
#' for each z-layer (bregma slice)? This is useful for making coordinates
#' comparable between bregma slices for visualization and analysis. 
#' Defaults to \code{TRUE}. Use \code{FALSE} to obtain the coordinates as
#' provided in the data release.
#' @return An object of class \code{\linkS4class{SpatialExperiment}}.
#' @references Moffitt et al. (2018) Molecular, spatial, and functional
#' single-cell profiling of the hypothalamic preoptic region. 
#' Science, 362(6416), eaau5324.
#' @source \url{https://doi.org/10.5061/dryad.8t8s248}
#' @examples spe <- MouseHypothalamusMoffitt2018()
#' @export
MouseHypothalamusMoffitt2018 <- function(center.coords = TRUE)
{
    eh <- ExperimentHub::ExperimentHub()
    recs <- AnnotationHub::query(eh, c("MERFISH", "Mofitt2018"))
    
    seg.dat <- .getResource(recs, "_segmentation")

    # (1) expression matrix (assay)
    ind <- grep("Neuron_cluster_ID", colnames(seg.dat))
    genes <- colnames(seg.dat)[(ind + 1):ncol(seg.dat)]
    exprs <- t(as.matrix(seg.dat[,genes]))
    colnames(exprs) <- NULL 
     
    # (2) cell metadata (colData)
    cdat <- seg.dat[,seq_len(ind)]
    colnames(cdat)[c(2:3,5:7)] <- c("sample_id", "sex", "z", "x", "y")
    colnames(cdat) <- tolower(colnames(cdat))
    if(center.coords) cdat <- .centerCoords(cdat)
    rownames(cdat) <- NULL 

    spe <- SpatialExperiment::SpatialExperiment(assays = list(exprs = exprs),
                                                colData = cdat,
                                                spatialCoordsNames = c("x", "y", "z"))
    # (3) store blanks separately as an altExp
    spe <- .separateBlanks(spe)

    # (4) add molecules
    spe <- .addMolecules(spe, recs)
    return(spe)
}

.centerCoords <- function(cdat)
{
    cdat$z <- factor(cdat$z, levels = unique(cdat$z))
    sids <- unique(cdat$sample_id)
    for(i in sids)
    {
        ind <- cdat$sample_id == i
        ccdat <- cdat[ind,]
        spl <- split(ccdat, ccdat$z)
        
        for(j in seq_along(spl)) 
        {   
            spl[[j]]$x <- scale(spl[[j]]$x, scale = FALSE)
            spl[[j]]$y <- scale(spl[[j]]$y, scale = FALSE)
        }
    
        cdat[ind,] <- do.call(rbind, spl)
    }
    cdat$z <- as.numeric(as.character(cdat$z))
    return(cdat)
}

# blanks (= negative controls)
.separateBlanks <- function(spe)
{
    blanks <- paste("Blank", 1:5, sep = "_")
    blank.exprs <- assay(spe)[blanks,]
    blank.spe <- SpatialExperiment::SpatialExperiment(
                    assays = list(exprs = blank.exprs)) 
    ind <- setdiff(rownames(spe), blanks)
    spe <- spe[ind,]
    SingleCellExperiment::altExps(spe)$Blank <- blank.spe   
    return(spe)
}

.separateAnalogs <- function(spe, mol.dat)
{
    sf.genes <- setdiff(rownames(spe), mol.dat$gene)
    sf.exprs <- assay(spe)[sf.genes,]
    sf.spe <- SpatialExperiment::SpatialExperiment(
                    assays = list(exprs = sf.exprs)) 
    ind <- setdiff(rownames(spe), sf.genes)
    spe <- spe[ind,]
    SingleCellExperiment::altExps(spe)$seqFISH <- sf.spe   
    SingleCellExperiment::mainExpName(spe) <- "smFISH"
    return(spe)
}

.addMolecules <- function(spe, recs)
{
    # extract relevant columns
    mol.dat <- .getResource(recs, "_molecules")
    rel.cols <- c("Gene_name", "Cell_name", "Centroid_X", "Centroid_Y", "Centroid_Z")
    mol.dat <- mol.dat[,rel.cols]
    colnames(mol.dat) <- c("gene", "cell", "x", "y", "z")
    
    # solve overlap (some cells in molecule data not
    # in segmentation table, and vice versa) 
    mol.dat <- subset(mol.dat, cell %in% spe$cell_id)
    ucell <- unique(mol.dat$cell)
    spe <- subset(spe, , spe$cell_id %in% ucell)
    
    # replace long cell barcodes by integer index
    # to speed up bumpy matrix creation
    ind <- match(mol.dat$cell, spe$cell_id)
    mol.dat$cell <- ind
    mol.dat$cell <- factor(mol.dat$cell, levels = seq_len(ncol(spe)))
    
    # there are three types of features: 
    # (1) digital (with molecules), combinatorial smFISH imaging
    # (2) analog (no molecules), non-combinatorial sequential FISH
    # (3) blanks (negative controls)

    # store molecules for blanks in the blank altExp
    blanks <- paste("Blank", 1:5, sep = "-")
    ind <- mol.dat$gene %in% blanks
    bmol.dat <- mol.dat[ind,]
    mol.dat <- mol.dat[!ind,] 
    bmol <- BumpyMatrix::splitAsBumpyMatrix(bmol.dat[, c("x", "y", "z")], 
                                            row = bmol.dat$gene,
                                            col = bmol.dat$cell)   
    rownames(bmol) <- paste("Blank", 1:5, sep = "_")
    colnames(bmol) <- NULL 
    SummarizedExperiment::assay(
        SingleCellExperiment::altExps(spe)$Blank, "molecules") <- bmol   
    
    # divide between 
    # (1) analog (no molecules), non-combinatorial sequential FISH 
    # (2) digital (with molecules), combinatorial smFISH imaging
    spe <- .separateAnalogs(spe, mol.dat)
    mol <- BumpyMatrix::splitAsBumpyMatrix(mol.dat[, c("x", "y", "z")], 
                                           row = mol.dat$gene,
                                           col = mol.dat$cell)
    colnames(mol) <- NULL 
    SummarizedExperiment::assay(spe, "molecules") <- mol    
    return(spe)
}
