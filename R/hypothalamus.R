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
    
    # (A) RAW
    # (1) molecule data
    #mol.dat <- .getResource(recs, "_molecules")
    
    # (2) image data
    #if(use.images) img.dat <- .getImageData(recs)

    # (B) PROCESSED
    seg.dat <- .getResource(recs, "_segmentation")

    # (3) expression matrix (assay)
    ind <- grep("Neuron_cluster_ID", colnames(seg.dat))
    genes <- colnames(seg.dat)[(ind + 1):ncol(seg.dat)]
    exprs <- t(as.matrix(seg.dat[,genes]))
    colnames(exprs) <- NULL 
   
    # (4) cell metadata (colData)
    cdat <- seg.dat[,seq_len(ind)]
    colnames(cdat)[c(2:3,5:7)] <- c("sample_id", "sex", "z", "x", "y")
    colnames(cdat) <- tolower(colnames(cdat))
    if(center.coords) cdat <- .centerCoords(cdat)

    spe <- SpatialExperiment::SpatialExperiment(assays = list(exprs = exprs),
                                                colData = cdat,
                                                spatialCoordsNames = c("x", "y", "z"))

    return(spe)
}

.centerCoords <- function(cdat)
{
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
    return(cdat)
}

