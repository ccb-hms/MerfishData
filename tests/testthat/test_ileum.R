context("Ileum")

checkSPE <- function(spe, 
                     seg = c("baysor", "cellpose"),
                     images = TRUE,
                     polygons = TRUE)
{
    # check overall object
    expect_is(spe, "SpatialExperiment")

    # check dimensions
    seg <- match.arg(seg)
    expect_equal(nrow(spe), 241) 
    expect_equal(ncol(spe), ifelse(seg == "baysor", 5800, 8439))
        
    # check colData
    cdat.cols <- c("sample_id", "leiden_final")
    expect_true(all(cdat.cols %in% colnames(colData(spe))))

    # check spatialCoords
    expect_equal(ncol(spatialCoords(spe)), ifelse(seg == "baysor", 2, 3)) 
    expect_true(all(colnames(spatialCoords(spe))[1:2] == c("x", "y")))
    expect_true(is.numeric(spatialCoords(spe)))

    # check images
    if(images)
    {
        expect_equal(ncol(imgData(spe)), 4)
        expect_equal(nrow(imgData(spe)), 2)
        expect_true(all(imgData(spe)$image_id == c("dapi", "membrane")))
        expect_is(imgData(spe)$data$dapi, "LoadedSpatialImage")
    }

    # check polygons
    if(seg == "baysor" && polygons)
    {
        expect_true("polygons" %in% names(metadata(spe)))    
        expect_is(metadata(spe)$polygons, "matrix")
        expect_equal(mode(metadata(spe)$polygons), "numeric")
        expect_equal(ncol(metadata(spe)$polygons), 4)
        pcols <- c("z", "cell", "x", "y")
        expect_true(all(colnames(metadata(spe)$polygons) == pcols))
    }
}

test_that("standard usage", {
    spe <- MouseIleumPetukhov2021()
    checkSPE(spe)
})

test_that("segmentation = cellpose", {
    spe <- MouseIleumPetukhov2021(segmentation = "cellpose", use.images = FALSE)
    checkSPE(spe, seg = "cellpose", images = FALSE)
})

test_that("use.images = FALSE", {
    spe <- MouseIleumPetukhov2021(use.images = FALSE)
    checkSPE(spe, images = FALSE)
})

test_that("use.polygons = FALSE", {
    spe <- MouseIleumPetukhov2021(use.images = FALSE, use.polygons = FALSE)
    checkSPE(spe, images = FALSE, polygons = FALSE)
})
