#' MERFISH mouse ileum dataset from Petukhov et al., 2021
#' @description Obtain the MERFISH mouse ileum dataset from Petukhov et al., 2021
#' @details Spatial transcriptomics protocols based on in situ sequencing or 
#' multiplexed RNA fluorescent hybridization can reveal detailed tissue organization.
#' Distinguishing the boundaries of individual cells in such data is challenging.
#' Current segmentation methods typically approximate cells positions using nuclei
#' stains. 
#'
#' Petukhov et al., 2021, describe Baysor, a segmentation method,
#' which optimizes 2D or 3D cell boundaries considering joint likelihood of
#' transcriptional composition and cell morphology.
#' Baysor can also perform segmentation based on the detected transcripts alone.
#'
#' Petukhov et al., 2021, compare the results of Baysor segmentation (mRNA-only)
#' to the results of a deep learning-based segmentation method called Cellpose from
#' Stringer et al., 2021. Cellpose applies a machine learning framework for the
#' segmentation of cell bodies, membranes and nuclei from microscopy images.
#' 
#' The function allows to obtain segmented MERFISH mouse ileum data for both
#' segmentation methods.
#'
#' A note on storing images within a \code{\linkS4class{SpatialExperiment}}:
#' The default \code{use.images = TRUE} reduces the 9-frame z-stack images
#' for DAPI stain and Membrane Na+/K+ - ATPase fluorecense to single-frame
#' images (taking the first frame). For working with the 9-frame z-stack
#' images it is recommended to load the images individually from ExperimentHub.
#' @param segmentation character. Should be either \code{"baysor"} or
#' \code{"cellpose"}. Defaults to \code{"baysor"}. See details.
#' @param use.images logical. Should DAPI and Membrane Na+/K+ - ATPase images
#' be loaded into memory and annotated to the \code{\link{imgData}} slot of
#' the returned \code{\linkS4class{SpatialExperiment}}? Defaults to \code{TRUE}.
#' See details.
#' @param use.polygons logical. Should polygon cell boundaries be annotated
#' to the \code{\link{metadata}} of the returned \code{\linkS4class{SpatialExperiment}}? 
#' Defaults to \code{TRUE}. Only available for Baysor segmentation.
#' @return An object of class \code{\linkS4class{SpatialExperiment}}.
#' @references Petukhov et al. (2021) Cell segmentation in imaging-based
#' spatial transcriptomics. Nat Biotechnol, 40(3), 345-54. 
#'
#' Stringer et al. (2021) Cellpose: a generalist algorithm for cellular
#' segmentation.  Nat Methods, 18(1), 100-6.  
#' @source \url{https://doi.org/10.5061/dryad.jm63xsjb2}
#' @examples spe <- MouseIleumPetukhov2021()
#' @export
MouseIleumPetukhov2021 <- function(segmentation = c("baysor", "cellpose"),
                                   use.images = TRUE,
                                   use.polygons = TRUE)
{
    segmentation <- match.arg(segmentation)

    eh <- ExperimentHub::ExperimentHub()
    recs <- AnnotationHub::query(eh, c("MERFISH", "Petukhov2021"))
    
    # (A) RAW
    # (1) molecule data
    mol.dat <- .getResource(recs, "_molecules")
    
    # (2) image data
    if(use.images) img.dat <- .getImageData(recs)

    # (B) PROCESSED
    # (3) counts
    suffix <- paste(segmentation, "counts", sep = "_")
    counts <- .getResource(recs, suffix)
    ass.dat <- list(counts = counts) 

    # (4) colData        
    suffix <- paste(segmentation, "coldata", sep = "_")
    col.dat <- .getResource(recs, suffix)
    scn <- c("x", "y", "z")

    # for baysor: add segmentation and polygons
    if(segmentation == "baysor")
    {
        # (5) segmentation
        mol <- .getSegmentation(recs, mol.dat)
        ass.dat <- c(ass.dat, molecules = mol)
        scn <- scn[1:2]    
    }
    
    # create SpatialExperiment object
    spe <- SpatialExperiment::SpatialExperiment(assays = ass.dat,
                                                colData = col.dat,
                                                sample_id = "ileum",
                                                spatialCoordsNames = scn)
    if(use.images) SpatialExperiment::imgData(spe) <- img.dat
    # (6) polygons
    if(use.polygons && segmentation == "baysor")
        metadata(spe)$polygons <- .annotatePolygons(recs)
    return(spe)
}

.getImageData <- function(recs)
{
    # (a) DAPI images
    dapi.img <- .getResource(recs, "_dapi")
    dapi.img <- as.raster(dapi.img)
    dapi.img <- SpatialExperiment::SpatialImage(dapi.img)    
    dapi.img <- SpatialExperiment::mirrorImg(dapi.img, "h")

    # (b) Membrane marker images
    mem.img <- .getResource(recs, "_membrane")
    mem.img <- as.raster(mem.img)
    mem.img <- SpatialExperiment::SpatialImage(mem.img)    
    mem.img <- SpatialExperiment::mirrorImg(mem.img, "h")

    img.list <- list(dapi = dapi.img, mem = mem.img)
    img.data <- S4Vectors::DataFrame(
                            sample_id = "ileum",
                            image_id = c("dapi", "membrane"),
                            data = I(img.list),
                            scaleFactor = NA_real_)
    return(img.data)
}

.getSegmentation <- function(recs, mol.dat)
{
        seg <- .getResource(recs, "_segmentation")
        ind <- match(seg$molecule_id, mol.dat$molecule_id)
        rel.cols <- c("gene", "x_pixel", "y_pixel", "z_pixel")
        seg.add <- mol.dat[ind, rel.cols]
        colnames(seg.add)[2:4] <- c("x", "y", "z")
        seg <- cbind(seg, seg.add)
        mol <- BumpyMatrix::splitAsBumpyMatrix(seg[, c("x", "y", "z")], 
                                               row = seg$gene, col = seg$cell)
        return(mol)
}

.annotatePolygons <- function(recs)
{
    poly <- .getResource(recs, "_polygons")
    return(poly)
}
