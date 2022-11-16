#' Plot raster image
#' @description Small helper function to plot a raster image.
#' @param img a raster object representing a bitmap image.
#' @return A ggplot object.
#' @examples
#'  hgrid <- hcl(0, 80, seq(50, 80, 10))
#'  img <- as.raster(matrix(hgrid, nrow = 4, ncol = 5))
#'  plotRasterImage(img)
#' @export
plotRasterImage <- function(img)
{
    if(!requireNamespace("grid"))
        stop("Please install the 'grid' package to use 'plotRasterImage'")
    if(!requireNamespace("ggplot2"))
        stop("Please install the 'ggplot2' package to use 'plotRasterImage'")

    grb <- grid::rasterGrob(img,
                      interpolate = FALSE,
                      width = grid::unit(1, "npc"),
                      height = grid::unit(1, "npc"))
    p <- ggplot2::ggplot() + 
         ggplot2::annotation_custom(grob = grb,
                           xmin = 0,
                           xmax = ncol(grb$raster),
                           ymin = 0,
                           ymax = nrow(grb$raster)) + 
        ggplot2::coord_fixed(xlim = c(0, ncol(grb$raster)),
                             ylim = c(0, nrow(grb$raster)))
    return(p)
}

#' Plot spatial image with data overlay
#' @description A helper function to overlay data onto a spatial image.
#' @param df A \code{data.frame} storing the data to plot.
#' @param col character. A column of \code{df} to use for overlay onto the image.
#' @param img a raster object representing a bitmap image.
#' @return A ggplot.
#' @examples 
#'  gene <- rep(c("Cd44", "Cd8b1", "Cd79b"), each = 2)
#'  x <- c(1693, 1701, 1820, 3188, 1631, 1881)
#'  y <- c(1831, 1666, 1855, 6612, 1533, 942)
#'  df <- data.frame(gene = gene, x = x, y = y)
#'  plotXY(df, "gene")
#' @export
plotXY <- function(df, col, img = NULL)
{
    if(!requireNamespace("ggpubr"))
        stop("Please install the 'ggpubr' package to use 'plotXY'")

    .data <- x <- y <- NULL
    if(!is.null(img)) p <- plotRasterImage(img)
    else p <- ggplot2::ggplot()
    p <- p + ggplot2::geom_point(data = df,
                        shape = 16,
                        size = 0.5,
                        ggplot2::aes(x, y, col = .data[[col]])) +
                        ggplot2::theme_void()
    if(is.numeric(df[[col]])) p + ggplot2::scale_color_viridis_c()
    else
    {
        ucsc <- function(n) ggpubr::get_palette("ucscgb", n)
        cols <- ucsc(length(levels(factor(df[[col]]))))
        p + ggplot2::theme(legend.key.size = grid::unit(0.5, "lines")) +
            ggplot2::scale_color_manual(values = cols) +
            ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 2)))
    }
}

#' Plot a tabset 
#' @description Plot a tabset of colData annotations of one or more
#' SpatialExperiment objects over an image.
#' @param spe.list A named \code{list} of \code{\linkS4class{SpatialExperiment}}
#' objects.
#' @param img a raster object representing a bitmap image.
#' @return None. Produces a tabset for rendering with \code{rmarkdown}.
#' @examples
#'     # create simulated data as described in the SpatialExperiment man page
#'     example("SpatialExperiment", package = "SpatialExperiment", echo = FALSE)
#'     spe <- spe_mol   
#' 
#'     # add simulated cell centroids
#'     s <- cbind(x = runif(20), y = runif(20))    
#'     spatialCoords(spe) <- s 
#'
#'     # add simulated cell type and cell cycle annotation
#'     ct <- c("ct1", "ct2", "ct3")
#'     cc <- c("G1", "G2", "S", "M") 
#'     spe$type <- sample(ct, ncol(spe), replace = TRUE)
#'     spe$cycle <- sample(cc, ncol(spe), replace = TRUE)
#'
#'     # create an example image
#'     hgrid <- hcl(0, 80, seq(50, 80, 10))
#'     img <- as.raster(matrix(hgrid, nrow = 4, ncol = 5))
#'
#'     # plotTabset
#'     spe.list <- list(myseg = spe)
#'     plotTabset(spe.list, img)  
#' @export
plotTabset <- function(spe.list, img)
{
    for (n in names(spe.list))
    {
        cat("### ", n, " {.tabset} \n\n")
        spe <- spe.list[[n]]
        df <- data.frame(SummarizedExperiment::colData(spe),
                         SpatialExperiment::spatialCoords(spe))
        rel.cols <- c("sample_id", "cell", "seg", "x", "y")
        rel.cols <- setdiff(names(df), rel.cols)
        # ...plot each variable
        for (col in rel.cols)
        {
            p <- plotXY(df, col, img)
            cat("#### ", col, '{-}', "\n")
            print(p)
            cat("\n\n")
        }
    }
}

#' Add holes to polygons
#' @description Add holes to a `data.frame` of cell polygon coordinates 
#' @param poly A \code{data.frame} storing cell polygon coordinates.
#' Expected columns include \itemize{
#' \item \code{"cell"} storing the cell ID,
#' \item \code{"x"} storing x-coordinates of the polygon corners,
#' \item \code{"y"} storing the y-coordinates of the polygon corners.}
#' @return A \code{data.frame}
#' @export
#' @examples
#'      x <- c(2053, 2053, 2053, 2056, 2059, 2059)
#'      y <- c(51, 54, 57, 57, 57, 54)
#'      poly <- data.frame(cell = 1, x = x, y = y)
#'      poly <- addHolesToPolygons(poly)
addHolesToPolygons <- function(poly)
{
    .f <- function(df) 
    {
        df$x <- df$x + 0.5 * (mean(df$x) - df$x)
        df$y <- df$y + 0.5 * (mean(df$y) - df$y)
        return(df)
    }
    spl <- split(poly, poly$cell)
    dl <- lapply(spl, .f)
    holes <- do.call(rbind, dl)
    poly$subid <- 1L
    holes$subid <- 2L
    poly <- rbind(poly, holes)
    return(poly)
}

.getResource <- function(recs, suffix)
{
    suffix <- paste0(suffix, "$")
    ind <- grep(suffix, recs$title)
    id <- recs$ah_id[ind]
    suppressMessages(rec <- recs[[id]])
    return(rec)
}
