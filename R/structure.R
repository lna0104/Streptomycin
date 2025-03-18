
#' Visualising the protein structure of rpoB
#'
#' @param file_pdb pdb file containing the protein structure
#' @param pos positions to label
#' @param mut_colour_variable values used for colour of labels
#' @param n_cols number of colours in colour gradient
#' @param file_html name of output html file
#'
#' @return returns the visualisation in html format
#' 
visualise_RpoB_structure <- function(file_pdb, 
                                     pos,
                                     mut_colour_variable,
                                     n_cols = 100,
                                     file_html) {
  
  max_value <- max(mut_colour_variable[pos], na.rm = TRUE)
  colourScale <- hsv(0.66, seq(0, 1, length.out = n_cols), 1)
  
  ngl <- NGLVieweR("./data/5uac_chainC.pdb", width = "800px", height = "800px")
  
  # Color the entire chain by position
  ngl <- ngl |>
    stageParameters(backgroundColor = "white", zoomSpeed = 1) |>
    addRepresentation("cartoon",
                      param = list(name = "cartoon", colorScheme = "residueindex"))
  
  # add labels at known mutation sites:
  for(p in pos) {
    col_i <- round(mut_colour_variable[p]/max_value * n_cols)
    ngl <- ngl |> 
      addRepresentation("label",
                        param = list(
                          sele = paste(p, "and .CA"),
                          #labelType = "format",
                          #labelFormat = "%(resname)s%(resno)s",
                          labelGrouping = "residue",
                          showBackground = TRUE,
                          showBorder = TRUE,
                          borderColor = "black",
                          colorValue = "black",
                          backgroundColor = colourScale[col_i],
                          backgroundOpacity = 1,
                          radius = 3,
                          fixedSize = FALSE
                        )) #|>
#      addRepresentation("line",
#                        param = list(sele = as.character(p), colorValue = "grey30"))
  }
  # export as html:
  htmlwidgets::saveWidget(ngl, file_html, selfcontained = TRUE)
  
  # legends:
  png("plots/distance_legend.png", width = 200, height = 100)
  par(mar = c(4,1,1,1))
  plot(NA, xlim = c(0, max_value), ylim = c(0,1),
       xlab = "Mean amino acid distance", ylab = NA,
       axes = FALSE, cex.lab = 1.2)
  axis(1, at = c(0, 10, 20, 30, 40, round(max_value)), 
       pos = 0)
  rect(xleft = 0:99 * max_value/100, xright = 1:100 * max_value/100, 
       ytop = 1, ybottom = 0,
       col = colourScale, border = NA)
  dev.off()
  return(ngl)
}

