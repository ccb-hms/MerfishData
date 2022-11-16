# MerfishData
Collection of public MERFISH datasets

## Installation

### Installation from Bioconductor

The package is available from [Bioconductor](https://bioconductor.org/).
Please follow the installation instructions on the
[package landing page](https://bioconductor.org/packages/MerfishData)
to install the package.

### Installation from GitHub

Make sure to have the latest release version of 
[R](https://cran.r-project.org/) installed.

Then proceed from within R via:

```
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "devel")
BiocManager::valid()
```

to install the 
[BiocManager](https://cran.r-project.org/web/packages/BiocManager/index.html) package
and to use the
[development version of Bioconductor](https://contributions.bioconductor.org/use-devel.html).

This is necessary as the `MerfishData` package uses resources from ExperimentHub
that are currently only available in the development version of Bioconductor.

The package itself can then be installed via:

```
BiocManager::install("ccb-hms/MerfishData")
```

NOTE: you will need the `remotes` package to install from github. 

To build the package vignettes upon installation use:

```
BiocManager::install("ccb-hms/MerfishData",
                     build_vignettes = TRUE,
                     dependencies = TRUE)
```

Once you have the package installed, you can inspect the vignettes from within
R via:

```
browseVignettes("MerfishData")
```

