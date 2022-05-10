context("metadata validity")

test_that("metadata is valid",
{
    if(!requireNamespace("ExperimentHubData", quietly = TRUE))
        BiocManager::install("ExperimentHubData")

    path <- find.package("")
    metadata <- system.file("extdata", "metadata.csv", package = "MerfishData")
    expect_true(ExperimentHubData::makeExperimentHubMetadata(path, metadata))
})
