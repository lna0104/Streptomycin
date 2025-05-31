# library(igraph)
# library(tidyverse)
# 
# 
plot_full_codon_network <- function() {
  ns <- c("A", "C", "T", "G")
  nodes <- expand.grid(ns, ns, ns) |>
    apply(1, paste0, collapse = "")

  adj <- matrix(NA,
                nrow = length(nodes),
                ncol = length(nodes),
                dimnames = list(nodes, nodes))
  for(i in 1:length(nodes)) {
    for(j in 1:length(nodes)) {
      adj[i,j] <- as.integer(sum((substr(nodes[i], 1, 1) == substr(nodes[j], 1, 1)) ||
                               (substr(nodes[i], 2, 2) == substr(nodes[j], 2, 2)) ||
                               (substr(nodes[i], 3, 3) == substr(nodes[j], 3, 3))) == 1L)
    }
  }
  diag(adj) <- 0L
  g <- graph_from_adjacency_matrix(adj, mode = "undirected")
  plot(g)
}


get_codon_network10_coords <- function(alpha = 0,      # angle between the central upper node and a vertical line
                                       beta = pi/6,    # angle between nodes in a group of three
                                       c = 10,         # distance between side peripheral and central node
                                       d = 12)         # distance between middle peripheral and central node
{
  x_coords <- array(NA, dim = c(4, 4, 4))
  y_coords <- x_coords
  
  # do first 4 nodes "by hand":
  x_coords[1, 1, 1] <- 0
  y_coords[1, 1, 1] <- 0
  x_coords[2, 1, 1] <- sin(alpha - beta) * c
  y_coords[2, 1, 1] <- cos(alpha - beta) * c
  x_coords[3, 1, 1] <- sin(alpha) * d
  y_coords[3, 1, 1] <- cos(alpha) * d
  x_coords[4, 1, 1] <- sin(alpha + beta) * c
  y_coords[4, 1, 1] <- cos(alpha + beta) * c

  # other nodes by rotation by +-120 degrees:
  x_coords[1, 2, 1] <- sin(alpha - beta + 2 * pi /3) * c
  y_coords[1, 2, 1] <- cos(alpha - beta + 2 * pi /3) * c
  x_coords[1, 3, 1] <- sin(alpha + 2 * pi /3) * d
  y_coords[1, 3, 1] <- cos(alpha + 2 * pi /3) * d
  x_coords[1, 4, 1] <- sin(alpha + beta + 2 * pi /3) * c
  y_coords[1, 4, 1] <- cos(alpha + beta + 2 * pi /3) * c
  
  x_coords[1, 1, 2] <- sin(alpha - beta - 2 * pi /3) * c
  y_coords[1, 1, 2] <- cos(alpha - beta - 2 * pi /3) * c
  x_coords[1, 1, 3] <- sin(alpha - 2 * pi /3) * d
  y_coords[1, 1, 3] <- cos(alpha - 2 * pi /3) * d
  x_coords[1, 1, 4] <- sin(alpha + beta - 2 * pi /3) * c
  y_coords[1, 1, 4] <- cos(alpha + beta - 2 * pi /3) * c
  
  # turn into data frame:
  nodes <- expand.grid(1:4, 1:4, 1:4) |>
    mutate(coord_x = NA, coord_y = NA)
  for(i in 1:nrow(nodes)) {
    nodes$coord_x[i] <- x_coords[nodes$Var1[i], nodes$Var2[i], nodes$Var3[i]]
    nodes$coord_y[i] <- y_coords[nodes$Var1[i], nodes$Var2[i], nodes$Var3[i]]
  }
  return(nodes)
}


get_codon_network28_coords <- function(alpha = -pi/24,   # angle between central and peripheral node
                                       beta = pi/9,    # angle between peripheral and intermediate node
                                       c = 6,          # distance of central nodes from origin
                                       d = 15,          # distance between central and peripheral nodes
                                       e = 10) {        # distance between central and intermediate nodes

  x_coords <- array(NA, dim = c(4, 4, 4))
  y_coords <- x_coords
  
  # do first 7 nodes "by hand":
  x_coords[1, 1, 1] <- c
  y_coords[1, 1, 1] <- c
  x_coords[1, 2, 1] <- sin(alpha) * d + c
  y_coords[1, 2, 1] <- cos(alpha) * d + c
  x_coords[1, 3, 1] <- sin(alpha + beta) * e + c
  y_coords[1, 3, 1] <- cos(alpha + beta) * e + c
  x_coords[1, 4, 1] <- sin(alpha + 2*beta) * d + c
  y_coords[1, 4, 1] <- cos(alpha + 2*beta) * d + c
  x_coords[1, 1, 4] <- sin(pi/2 - alpha) * d + c
  y_coords[1, 1, 4] <- cos(pi/2 - alpha) * d + c
  x_coords[1, 1, 3] <- sin(pi/2 - alpha - beta) * e + c
  y_coords[1, 1, 3] <- cos(pi/2 - alpha - beta) * e + c
  x_coords[1, 1, 2] <- sin(pi/2 - alpha - 2*beta) * d + c
  y_coords[1, 1, 2] <- cos(pi/2 - alpha - 2*beta) * d + c

  # other nodes by rotation:
  x_coords[2, ,] <- y_coords[1, ,]
  y_coords[2, ,] <- - x_coords[1, ,]
  x_coords[3, ,] <- -x_coords[1, ,]
  y_coords[3, ,] <- -y_coords[1, ,]
  x_coords[4, ,] <- - y_coords[1, ,]
  y_coords[4, ,] <- x_coords[1, ,]
  
  # turn into data frame:
  nodes <- expand.grid(1:4, 1:4, 1:4) |>
    mutate(coord_x = NA, coord_y = NA)
  for(i in 1:nrow(nodes)) {
    nodes$coord_x[i] <- x_coords[nodes$Var1[i], nodes$Var2[i], nodes$Var3[i]]
    nodes$coord_y[i] <- y_coords[nodes$Var1[i], nodes$Var2[i], nodes$Var3[i]]
  }
  return(nodes)
}


get_codon_network64_coords <- function(alpha = -pi/24,   # angle between central and peripheral node
                                       beta = pi/9,    # angle between peripheral and intermediate node
                                       c = 6,          # distance of central nodes from origin
                                       d = 15,          # distance between central and peripheral nodes
                                       e = 10,         # distance between central and intermediate nodes
                                       f = 25) {       # distance between central and outer nodes
  
  x_coords <- array(NA, dim = c(4, 4, 4))
  y_coords <- x_coords
  
  # do first nodes "by hand":
  x_coords[1, 1, 1] <- c
  y_coords[1, 1, 1] <- c
  
  x_coords[1, 2, 1] <- sin(alpha) * d + c
  y_coords[1, 2, 1] <- cos(alpha) * d + c
  x_coords[1, 3, 1] <- sin(alpha + beta) * e + c
  y_coords[1, 3, 1] <- cos(alpha + beta) * e + c
  x_coords[1, 4, 1] <- sin(alpha + 2*beta) * d + c
  y_coords[1, 4, 1] <- cos(alpha + 2*beta) * d + c
  x_coords[1, 1, 4] <- sin(pi/2 - alpha) * d + c
  y_coords[1, 1, 4] <- cos(pi/2 - alpha) * d + c
  x_coords[1, 1, 3] <- sin(pi/2 - alpha - beta) * e + c
  y_coords[1, 1, 3] <- cos(pi/2 - alpha - beta) * e + c
  x_coords[1, 1, 2] <- sin(pi/2 - alpha - 2*beta) * d + c
  y_coords[1, 1, 2] <- cos(pi/2 - alpha - 2*beta) * d + c
  
  gamma <- (pi/2 - 2 * alpha) / 8
  x_coords[1, 2, 2] <- sin(alpha) * f + c
  y_coords[1, 2, 2] <- cos(alpha) * f + c
  x_coords[1, 2, 3] <- sin(alpha + 1 * gamma) * f + c
  y_coords[1, 2, 3] <- cos(alpha + 1 * gamma) * f + c
  x_coords[1, 3, 2] <- sin(alpha + 2 * gamma) * f + c
  y_coords[1, 3, 2] <- cos(alpha + 2 * gamma) * f + c
  x_coords[1, 2, 4] <- sin(alpha + 3 * gamma) * f + c
  y_coords[1, 2, 4] <- cos(alpha + 3 * gamma) * f + c
  x_coords[1, 3, 3] <- sin(alpha + 4 * gamma) * f + c
  y_coords[1, 3, 3] <- cos(alpha + 4 * gamma) * f + c
  x_coords[1, 4, 2] <- sin(alpha + 5 * gamma) * f + c
  y_coords[1, 4, 2] <- cos(alpha + 5 * gamma) * f + c
  x_coords[1, 3, 4] <- sin(alpha + 6 * gamma) * f + c
  y_coords[1, 3, 4] <- cos(alpha + 6 * gamma) * f + c
  x_coords[1, 4, 3] <- sin(alpha + 7 * gamma) * f + c
  y_coords[1, 4, 3] <- cos(alpha + 7 * gamma) * f + c
  x_coords[1, 4, 4] <- sin(alpha + 8 * gamma) * f + c
  y_coords[1, 4, 4] <- cos(alpha + 8 * gamma) * f + c
  
  # other nodes by rotation:
  x_coords[2, ,] <- y_coords[1, ,]
  y_coords[2, ,] <- - x_coords[1, ,]
  x_coords[3, ,] <- -x_coords[1, ,]
  y_coords[3, ,] <- -y_coords[1, ,]
  x_coords[4, ,] <- - y_coords[1, ,]
  y_coords[4, ,] <- x_coords[1, ,]
  
  # turn into data frame:
  nodes <- expand.grid(1:4, 1:4, 1:4) |>
    mutate(coord_x = NA, coord_y = NA)
  for(i in 1:nrow(nodes)) {
    nodes$coord_x[i] <- x_coords[nodes$Var1[i], nodes$Var2[i], nodes$Var3[i]]
    nodes$coord_y[i] <- y_coords[nodes$Var1[i], nodes$Var2[i], nodes$Var3[i]]
  }
  return(nodes)
}



get_codon_freqs <- function(output, pos) {
  return(output |>
           filter(AA_pos_Ecoli == pos) |>
           distinct(species, .keep_all = TRUE) |>
           group_by(Codon_target) |>
           summarise(n = n()) |>
           mutate(freq = n / sum(n)))
}


# Function to add a horizontal gradient legend to a base R plot (chatGPT solution)
add_gradient_legend <- function(xmin, xmax, ymin, ymax, max_value) {
  gradient <- colorRampPalette(c("white", "yellow"))(100)
  n <- length(gradient)
  
  # Draw the gradient rectangles horizontally
  for (i in 1:n) {
    rect(xmin + (xmax - xmin) * (i - 1) / n, ymin, xmin + (xmax - xmin) * i / n, ymax,
         col = gradient[i], border = NA)
  }
  
  # Add minimum and maximum labels using text()
  text(x = xmin, y = ymin - (ymax-ymin)/200, labels = 0, pos = 1)
  text(x = xmax, y = ymin - (ymax-ymin)/200, labels = max_value, pos = 1)
}


plot_codon_network <- function(type = "type28",
                               site_order = c(1, 2, 3),
                               pos,
                               focal_codon,
                               mutations_list,
                               output,
                               node_radius = 2.2,
                               lines_occupied_codons = TRUE,
                               freq_threshold = 0.001,
                               file_path = NULL)
{
  par(mar=c(0,0,0,0) + 0.1)
  if (type == "type64") {
    nodes <- get_codon_network64_coords()
  } else if (type == "type28") {
    nodes <- get_codon_network28_coords()
  } else if (type == "type10") {
    nodes <- get_codon_network10_coords()
  } else {
    stop("Unknown type argument, must be either \"type10\", \"type28\" or \"type64\".")
  }
  
  nts <- c("A", "T", "G", "C")
  nts_pos1 <- c(substr(focal_codon, 1, 1), nts[nts != substr(focal_codon, 1, 1)])
  nts_pos2 <- c(substr(focal_codon, 2, 2), nts[nts != substr(focal_codon, 2, 2)])
  nts_pos3 <- c(substr(focal_codon, 3, 3), nts[nts != substr(focal_codon, 3, 3)])
  
  nodes$Codon_target <- paste0(nts_pos1[nodes[, paste0("Var", site_order[1])]], 
                               nts_pos2[nodes[, paste0("Var", site_order[2])]], 
                               nts_pos3[nodes[, paste0("Var", site_order[3])]])
  
  nodes$AA1 <- sapply(nodes$Codon_target, 
                      function(x) {as.character(translate(DNAString(x), no.init.codon = TRUE))} )
  nodes$AA3 <- c(AMINO_ACID_CODE, "*"="*")[nodes$AA1]
  nodes <- nodes |>
    left_join(get_codon_freqs(output, pos), by = join_by(Codon_target)) |>
    mutate(freq = ifelse(is.na(freq), 0, freq))
  
  nodes$resistant <- paste(pos, nodes$AA1, sep = "_") %in% 
    (paste(mutations_list$AA_pos_Ecoli, mutations_list$AA_mutation, sep = "_"))
  
  if (!is.null(file_path)) {
    file_name <- paste0(file_path, "codon_network_", type, "_AApos", pos, "_focal", focal_codon, ".pdf")
    if (type == "type10") pdf(file_name, width = 4, height = 4)
    if (type == "type28") pdf(file_name, width = 6, height = 6)
    if (type == "type64") pdf(file_name, width = 9, height = 9)
  }
  par(mar = c(0,0,0,0) + 0.1)
  minx <- min(nodes$coord_x - node_radius, na.rm = TRUE)
  maxx <- max(nodes$coord_x + node_radius, na.rm = TRUE)
  miny <- min(nodes$coord_y - node_radius, na.rm = TRUE)
  maxy <- max(nodes$coord_y + node_radius, na.rm = TRUE)
  plot(NA, type="n", axes=FALSE, ann=FALSE, 
       xlim = c(minx, maxx), ylim = c(miny, maxy))
  for(i in 1:nrow(nodes)) {
    if (!is.na(nodes$coord_x[i])) {
      for(j in 1:nrow(nodes)) {
        if ((hamming(nodes[i, 1:3], nodes[j, 1:3]) == 1L) &&
            ((nodes$freq[i] > freq_threshold) || (!lines_occupied_codons)))
          lines(x = c(nodes$coord_x[i], nodes$coord_x[j]),
                y = c(nodes$coord_y[i], nodes$coord_y[j]))
      }
    }
  }
  
  for(i in 1:nrow(nodes)) {
    if (!is.na(nodes$coord_x[i])) {
      draw.circle(nodes$coord_x[i], 
                  nodes$coord_y[i], 
                  node_radius, 
                  col = ifelse(nodes$freq[i] < freq_threshold, "lightgrey",
                               hsv(0.16, nodes$freq[i]/max(nodes$freq), 1)))

      # draw.circle(nodes$coord_x[i], 
      #             nodes$coord_y[i], 
      #             node_radius, 
      #             col = ifelse(nodes$freq[i] < freq_threshold, "lightgrey",
      #                          hsv(0.16, nodes$freq[i]/0.7, 1)))

      if (nodes$resistant[i]) {
        draw.circle(nodes$coord_x[i], 
                    nodes$coord_y[i], 
                    node_radius, 
                    border = "#C93312",
                    lwd = 3)
      }
      text(nodes$coord_x[i], nodes$coord_y[i] + node_radius/3, nodes$Codon_target[i])
      text(nodes$coord_x[i], nodes$coord_y[i] - node_radius/3, nodes$AA3[i])
    }
  }
  text(max(nodes$coord_x + node_radius, na.rm = TRUE),
       max(nodes$coord_y + node_radius, na.rm = TRUE),
      paste("pos =", pos),
      adj = c(1, 1))
  
  add_gradient_legend(minx, 
                      minx + (maxx-minx)/6, 
                      miny, 
                      miny + (maxy-miny)/30, 
                      round(max(nodes$freq, na.rm = TRUE), digits = 2))
  if (!is.null(file_path)) {
    dev.off()
  }
}
