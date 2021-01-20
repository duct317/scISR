#' @title plotscISR: Visualization of scRNA-seq data.
#' @description Plotting scRNA-seq data to visualize the cell transciptome landscape.
#'
#' @param data Input matrix or data frame. Rows represent genes while columns represent samples.
#' @param label A vector that contains the cell type labels.
#' @param perplexity Number of closed neighbors for each data point.
#
#' @details
#' This function utilizes the dimension reduction feature of t-SNE package to viusualize high throughput data over two dimensional space.
#'
#' @return
#' \code{plotscISR} Return the visualization of the cell populations in the scRNA-seq dataset.
#'
#' @examples
#'
#' # Load the sample dataset scISRExample
#' data(scISRExample)
#'
#' # Perform the imputation
#' imputed <- scISR(data = scISRExample$dropout)
#'
#' # Plot the complete data
#' plot_raw <- plotscISR(scISRExample$raw, label = scISRExample$celltype)
#'
#' # Plot the dropout data
#' plot_dropout <- plotscISR(scISRExample$dropout, label = scISRExample$celltype)
#'
#' # Plot the imputed data
#' plot_imputed <- plotscISR(imputed, label = scISRExample$celltype)
#'
#' @import Rtsne ggplot2
#' @export

plotscISR <- function(data, label, perplexity=30){
    set.seed(1)
    pca <- Rtsne(t(data),dims = 2, perplexity= perplexity)$Y
    pca <- data.frame(pca)
    pca$group <- label
    colnames(pca) <- c("PC1","PC2","Group")
    colors <- c("#BC3C28", "#E18726", "#0072B5", "#1F854E", "#7776B1", "#00B0C6", "#857C19", "#D2B48C", "#1B1F2A", "#FFED4F", "#8A6A5F", "#CEC277", "#19227E", "#AE7F39", "#800000", "#aaffc3", "#66BA69", "#ffd8b1", "#000075", "#a9a9a9")
    # plot for raw data
    plot <- ggplot(pca, aes(x = PC1, y = PC2, colour = Group)) + labs(x = "t-SNE1", y = "t-SNE2") +
        theme(plot.title = element_text(hjust = 0.5)) +
        geom_point(size=2, alpha = 0.9) +
        scale_color_manual(values = colors)+
        theme_classic()+
        theme(axis.title=element_text(size=14),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              legend.position="bottom",
              legend.box = "horizontal",
              legend.title=element_text(size=18),
              legend.margin = margin(0),
              legend.spacing.x = unit(0.1, "cm"),
              legend.text=element_text(size=18))
    return(plot)
}
