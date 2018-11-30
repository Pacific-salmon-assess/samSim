#' Minimalist ggplot theme
#'
#' This function is a copy of S. Anderson's \code{ggsidekick} which is not
#' supported in the newest R version. Adds an argument for top, bottom, middle
#' for multipanel functionality.
#'
#' @import ggplot2
#'
#' @param base_size Defaults to 11 (don't change).
#' @param base_family Default font type (don't change).
#' @param position A character vector that can take the values \code{standard}
#' (default), \code{bottom}, \code{middle}, \code{top}, or \code{topWithX}.
#' Major differences between them are the presence of x-axis and facet labels.
#' @param axisSize A numeric representing font size for the axis labels;
#' defaults to 10.
#' @param legendSize A numeric scalar for the size of the legend; defaults to
#' 1.
#' @param facetSize A numeric scalar for the relative size of facet labels.
#' @return Returns a minimalist ggplot (i.e. no grid lines, more white and less
#' grey).
#'
#' @export

theme_sleekX <- function(base_size = 11, base_family = "", position = "standard",
                         axisSize = 10, legendSize = 1, facetSize = 1.1) {
  half_line <- base_size/2
  q <- theme_light(base_size = 11, base_family = "") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks.length = unit(half_line / 2.2, "pt"),
      strip.background = element_rect(fill = NA, colour = NA),
      strip.text.x = element_text(colour = "grey30", size = rel(facetSize)),
      strip.text.y = element_text(colour = "grey30"),
      axis.text = element_text(colour = "grey30", size = 0.9*axisSize),
      axis.title = element_text(colour = "grey30", size = axisSize),
      legend.title = element_text(colour = "grey30", size = rel(1.1 * legendSize)),
      panel.border = element_rect(fill = NA, colour = "grey70", size = 1),
      legend.key.size = unit(1, "lines"),
      legend.text = element_text(size = rel(legendSize), colour = "grey30"),
      legend.key = element_rect(colour = NA, fill = NA),
      legend.background = element_rect(colour = NA, fill = NA),
      plot.title = element_text(colour = "grey30", size = rel(1)),
      plot.subtitle = element_text(colour = "grey30", size = rel(.85))
    )
  if (position == "bottom") {
    q <- q + theme(strip.background = element_blank(),
                   strip.text.x = element_blank(),
                   strip.text.y = element_blank(),
                   axis.text.y = element_text(size = 0.9*axisSize),
                   axis.text.x = element_text(size = 0.9*axisSize),
                   axis.title = element_text(size = axisSize)
    )
  }
  if (position == "top") {
    q <- q + theme(strip.text = element_text(size = axisSize),
                   axis.text.y = element_text(size = 0.9*axisSize),
                   axis.text.x = element_blank(),
                   axis.title.y = element_text(size = axisSize),
                   axis.title.x = element_blank()
    )
  }
  if (position == "topWithX") {
    q <- q + theme(strip.text = element_text(size = axisSize),
                   axis.text.y = element_text(size = 0.9*axisSize),
                   axis.text.x = element_text(size = 0.9*axisSize),
                   axis.title.y = element_text(size = axisSize),
                   axis.title.x = element_blank()
    )
  }
  if (position == "mid") {
    q <- q + theme(strip.background = element_blank(),
                   strip.text.x = element_blank(),
                   strip.text.y = element_blank(),
                   axis.text.y = element_text(size = 0.9*axisSize),
                   axis.text.x = element_blank(),
                   axis.title.y = element_text(size = axisSize),
                   axis.title.x = element_blank()
    )
  }
  return(q)
}
