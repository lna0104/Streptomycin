
#' Visualising the protein structure of rpoB
#'
#' @param file_pdb pdb file containing the protein structure
#' @param chain_id Character. Protein chain identifier to display (default = "L").
#' @param ligand_resid Character. Residue name of the ligand to include (default = "5I0").
#' @param pos positions to label
#' @param mut_colour_variable values used for colour of labels
#' @param n_cols number of colours in colour gradient
#' @param file_html name of output html file
#'
#' @return returns the visualisation in html format
#' 
visualise_rpsL_structure <- function(file_pdb, 
                                     chain_ids = c("L","A"),
                                     ligand_resid = "5I0",
                                     pos,
                                     mut_colour_variable,
                                     n_cols = 100,
                                     file_html) {
  
  max_value <- max(mut_colour_variable[pos], na.rm = TRUE)
  colourScale <- hsv(0.66, seq(0, 1, length.out = n_cols), 1)
  # colourScale <- colorRampPalette(brewer.pal(9, "Spectral"))(n_cols)
 
  ngl <- NGLVieweR(file_pdb, width = "800px", height = "800px") |>
    stageParameters(backgroundColor = "white", zoomSpeed = 1)
  
  # Color the entire chain and ligand by position
  # ngl <- ngl |>
  #   stageParameters(backgroundColor = "white", zoomSpeed = 1) |>
  #   addRepresentation("cartoon",
  #                     param = list(name = "cartoon", colorScheme = "residueindex"))

  # Add the two chains with different colors
  chain_colors <- brewer.pal(8, name = "Set2")

  for (i in seq_along(chain_ids)) {
    ngl <- ngl |>
      addRepresentation("cartoon",
                        param = list(
                          name = paste0("chain_", chain_ids[i]),
                          sele = paste0(":", chain_ids[i]),
                          colorScheme = "uniform", 
                          colorValue = chain_colors[i]
                        ))
  }

  # Add the ligand as ball+stick
  ligand_sele <- paste0(ligand_resid, " and :", chain_ids[2]) 

  ngl <- ngl |>
    addRepresentation("ball+stick",
                      param = list(
                        name = "ligand",
                        sele = ligand_sele,
                        colorScheme = "element"
                      ))

  # add labels at known mutation sites:
  for(p in pos) {
    col_i <- round(mut_colour_variable[p]/max_value * n_cols)
    ngl <- ngl |> 
      addRepresentation("label",
                        param = list(
                          sele = paste(p, "and .CA and", paste0(":", chain_ids[1])),
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
  
  # legend for mutation color scale 

  png("plots/distance_legend.pdf", width = 200, height = 100)
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

  # legend for chain colors 
  png("plots/chain_legend.pdf", width = 300, height = 150)
  par(mar = c(4, 4, 2, 2))
  plot.new()
  # Draw legend with chain names and the corresponding pastel colors
  legend("center", legend = chain_ids, fill = chain_colors,
         title = "Chains", cex = 1.2, bty = "n")
  dev.off()

  return(ngl)
}

