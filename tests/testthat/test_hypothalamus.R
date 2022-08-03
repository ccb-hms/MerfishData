context("Hypothalamus")

checkSPE <- function(spe)
{
    # check overall object
    expect_is(spe, "SpatialExperiment")
    expect_equal(nrow(spe), 161)
    expect_equal(ncol(spe), 73655)

    # check colData
    cdat.cols <- c("cell_id", "sample_id", "sex", "behavior", "cell_class")
    expect_true(all(cdat.cols %in% colnames(colData(spe))))

    # check spatialCoords
    expect_equal(ncol(spatialCoords(spe)), 3) 
    expect_true(all(colnames(spatialCoords(spe)) == c("x", "y", "z")))
    expect_true(is.numeric(spatialCoords(spe)))
}


test_that("standard usage", {
    spe <- MouseHypothalamusMoffitt2018()
    checkSPE(spe)
    spl <- split(spatialCoords(spe)[,"x"], spatialCoords(spe)[,"z"])
    expect_equal(round(mean(spl[[1]])), 0)
    spl <- split(spatialCoords(spe)[,"y"], spatialCoords(spe)[,"z"])
    expect_equal(round(mean(spl[[1]])), 0)
})

test_that("center.coords = FALSE", {
    spe <- MouseHypothalamusMoffitt2018(center.coords = FALSE)
    checkSPE(spe)
    spl <- split(spatialCoords(spe)[,"x"], spatialCoords(spe)[,"z"])
    expect_false(round(mean(spl[[1]])) == 0)
    spl <- split(spatialCoords(spe)[,"y"], spatialCoords(spe)[,"z"])
    expect_false(round(mean(spl[[1]])) == 0)
})
