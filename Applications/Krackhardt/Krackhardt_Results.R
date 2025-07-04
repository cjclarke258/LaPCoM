####################################################################################################
# KRACKHARDT DATA APPLICATION (RESULTS)
####################################################################################################

rm(list = ls())
setwd("/Volumes/My Passport/PhD/PhD Project 1/LaPCoM_Code/Applications/Krackhardt/")

# ================================================================================================
# LOAD THE DATA
# ================================================================================================
load("Krackhardt_Multiplex_Directed.Rdata")
# load("pp_res_krack.Rdata")
load("pp_res_krack_PPMethodB.Rdata")

# ================================================================================================
# OPTIMAL CLUSTERING
# ================================================================================================
G_plus_hat = krack_pp$krack_G
optimal_clustering_networks = krack_pp$krack_C
table(optimal_clustering_networks)

# ------------------------------------------------------------------------------------------------
# INVESTIGATE A LITTLE
# ------------------------------------------------------------------------------------------------
for (g in 1:G_plus_hat) {
  g_inds = which(optimal_clustering_networks == g)
  dens = rep(NA, length(g_inds))
  for (m in 1:length(g_inds)) {
    dens[m] = network::network.density(network::as.network(krack_multi[, , g_inds[m]]))
  } # end m for loop
  cat("The average density of a network in cluster", g, "is: ", round(mean(dens), 2), "\n")
} # end g for loop 

# ================================================================================================
# LOOK AT RESULTS
# ================================================================================================

# ------------------------------------------------------------------------------------------------
# COMPARE TO SIGNORELLI OUT OF CURIOSITY
# ------------------------------------------------------------------------------------------------
clust_SigWit = c(1, 2, 1, 1, 1, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1)
table(optimal_clustering_networks, clust_SigWit)
mclust::adjustedRandIndex(optimal_clustering_networks, clust_SigWit)

# ================================================================================================
# SELECT A COLOURBLIND FRIENDLY PALETTE
# ================================================================================================
# use Paul Tol's bright colour scheme
cols = khroma::color("muted")(9)
M = dim(krack_multi)[3]

# ------------------------------------------------------------------------------------------------
# PLOT THE NETWORKS WITH COLOUR CORRESPONDING TO CLUSTER
# ------------------------------------------------------------------------------------------------
# jpeg("Krackhardt_Data_Clustered_Arbitrary_Coords.png", width = 480 * 2, height = 480)
# par(mfrow = c(3, 7))
# par(mar = c(1.5, 1, 1.5, 1)) # bltr
# for (m in 1:M) {
#   plot(network::as.network(krack_multi[, , m], directed = T), main = LETTERS[m], cex.main = 1.8,
#        edge.col = "grey", vertex.col = cols[optimal_clustering_networks[m]], vertex.cex = 3, arrowhead.cex = 3)
# } # end m for loop (plotting)
# dev.off()

# ------------------------------------------------------------------------------------------------
# PLOT THE NETWORKS WITH COLOUR CORRESPONDING TO CLUSTER 
# ------------------------------------------------------------------------------------------------
ind = which.max(apply(krack_multi, 3, sum))
library(statnet); coords = gplot(network::as.network(krack_multi[, , ind], directed = T))
pdf("Results_PDF_Plots/Krackhardt_Data_Clustered_Same_Coords.pdf", width = 16, height = 8)
par(mfrow = c(3, 7))
par(mar = c(1.5, 1, 1.5, 1)) # bltr
for (m in 1:M) {
  plot(network::as.network(krack_multi[, , m], directed = T), main = LETTERS[m], cex.main = 1.8,
       edge.col = "grey", vertex.cex = 3, arrowhead.cex = 3,
       vertex.col = cols[optimal_clustering_networks[m]], 
       # vertex.border = cols[optimal_clustering_networks[m]], 
       coord = coords)
} # end m for loop (plotting)
dev.off()

# ------------------------------------------------------------------------------------------------
# PLOT THE NETWORKS WITH COLOUR CORRESPONDING TO CLUSTER AND NODE SIZE CORRESPONDING TO DEGREE
# ------------------------------------------------------------------------------------------------
# jpeg("Krackhardt_Data_Clustered_LS_Coords_Degree.png", width = 480 * 2, height = 480)
# par(mfrow = c(3, 7))
# par(mar = c(1.5, 1, 1.5, 1)) # bltr
# for (m in 1:M) {
#   plot(network::as.network(krack_multi[, , m], directed = T), main = LETTERS[m], cex.main = 1.8,
#        edge.col = "grey", arrowhead.cex = 3,
#        vertex.col = cols[optimal_clustering_networks[m]], 
#        vertex.cex = igraph::degree(igraph::graph_from_adjacency_matrix(krack_multi[, , m]))/6, 
#        coord = krack_pp$krack_Z[[optimal_clustering_networks[m]]])
# } # end m for loop (plotting)
# dev.off()

# ------------------------------------------------------------------------------------------------
# OVERLAY THE NETWORKS WITHIN THEIR LATENT SPACE  
# ------------------------------------------------------------------------------------------------
# jpeg("Krackhardt_Data_Clustered_LS_Coords_Overlay.png", width = 480 * 2, height = 480)
# par(mfrow = c(3, 7))
# # par(mar = c(1.5, 1, 1.5, 1)) # bltr
# par(mar = c(1.75, 1.75, 1.5, 0.5) + 0.4) # bltr
# for (m in 1:M) {
#   
#   z1 = krack_pp$krack_Z[[optimal_clustering_networks[m]]][, 1]
#   z2 = krack_pp$krack_Z[[optimal_clustering_networks[m]]][, 2]
# 
#   ones = which(krack_multi[, , m] != 0)
#   
#   plot(krack_pp$krack_Z[[optimal_clustering_networks[m]]], type = "n", las = 1, 
#        ylim = c(min(krack_pp$krack_Z[[optimal_clustering_networks[m]]][, 2]) - 0.2, max(krack_pp$krack_Z[[optimal_clustering_networks[m]]][, 2]) + 0.2),
#        xlim = c(min(krack_pp$krack_Z[[optimal_clustering_networks[m]]][, 1]) - 0.2, max(krack_pp$krack_Z[[optimal_clustering_networks[m]]][, 1]) + 0.2),
#        main = LETTERS[m], cex.main = 1.8)
#   grid()
#   
#   segments(z1[row(krack_multi[, , m])[ones]], z2[row(krack_multi[, , m])[ones]], 
#            z1[col(krack_multi[, , m])[ones]], z2[col(krack_multi[, , m])[ones]], 
#            lwd = 1, col = adjustcolor("grey10", 0.4))
#   
#   points(krack_pp$krack_Z[[optimal_clustering_networks[m]]], 
#          bg = cols[optimal_clustering_networks[m]], 
#          pch = 21, cex = 1.75)
#   
# } # end m for loop (plotting)
# dev.off()

# ------------------------------------------------------------------------------------------------
# OVERLAY THE NETWORKS WITHIN THEIR LATENT SPACE, REORDERED BY CLUSTER
# ------------------------------------------------------------------------------------------------
# jpeg("Krackhardt_Data_Clustered_LS_Coords_Overlay_Ordered.png", width = 480 * 2, height = 480)
# par(mfrow = c(3, 7))
# # par(mar = c(1.5, 1, 1.5, 1)) # bltr
# par(mar = c(1.75, 1.75, 1.5, 0.5) + 0.4) # bltr
# for (g in 1:G_plus_hat) {
#   
#   for (m in 1:M) {
#     
#     if (optimal_clustering_networks[m] == g) {
#       
#       z1 = krack_pp$krack_Z[[optimal_clustering_networks[m]]][, 1]
#       z2 = krack_pp$krack_Z[[optimal_clustering_networks[m]]][, 2]
#       
#       ones = which(krack_multi[, , m] != 0)
#       
#       plot(krack_pp$krack_Z[[optimal_clustering_networks[m]]], type = "n", las = 1, 
#            ylim = c(min(krack_pp$krack_Z[[optimal_clustering_networks[m]]][, 2]) - 0.2, max(krack_pp$krack_Z[[optimal_clustering_networks[m]]][, 2]) + 0.2),
#            xlim = c(min(krack_pp$krack_Z[[optimal_clustering_networks[m]]][, 1]) - 0.2, max(krack_pp$krack_Z[[optimal_clustering_networks[m]]][, 1]) + 0.2),
#            main = LETTERS[m], cex.main = 1.8)
#       grid()
#       
#       segments(z1[row(krack_multi[, , m])[ones]], z2[row(krack_multi[, , m])[ones]], 
#                z1[col(krack_multi[, , m])[ones]], z2[col(krack_multi[, , m])[ones]], 
#                lwd = 1, col = adjustcolor("grey10", 0.4))
#       
#       points(krack_pp$krack_Z[[optimal_clustering_networks[m]]], 
#              bg = cols[optimal_clustering_networks[m]], 
#              pch = 21, cex = 1.75)
#       
#     } # end if statement
#     
#   } # end m for loop 
#   
# } # end g for loop
# 
# dev.off()

# 

library(ggplot2)
library(dplyr)

# Create an empty list to store plots for each network (m)
plot_list = list()

model_cols = unname(cols[1:3])
  
for (m in 1:M) {
  
  # Extract latent space coordinates (Z1, Z2) for current network
  z1 = krack_pp$krack_Z[[optimal_clustering_networks[m]]][, 1]
  z2 = krack_pp$krack_Z[[optimal_clustering_networks[m]]][, 2]
  
  # Create a data frame for nodes
  nodes_df = data.frame(
    ID = 1:length(z1),
    Z1 = z1,
    Z2 = z2,
    Network = rep("A", length(z1)),
    Cluster = rep(optimal_clustering_networks[m], length(z1))
  )
  
  # Identify non-zero edges from the adjacency matrix (krack_multi)
  ones = which(krack_multi[, , m] != 0, arr.ind = TRUE)
  
  # Create a data frame for edges
  edges_df = data.frame(
    From = ones[, 1],
    To = ones[, 2],
    Z1_from = z1[ones[, 1]],
    Z2_from = z2[ones[, 1]],
    Z1_to = z1[ones[, 2]],
    Z2_to = z2[ones[, 2]],
    Network = rep("A", length(ones[, 1]))
    # Cluster = rep(optimal_clustering_networks[m], length(ones[, 1]))
  )
  
  # Create the plot for the current network
  # Create the plot with centered titles, same-sized boxes, and enhanced readability
  p = ggplot() +
    # Plot the edges (connections)
    geom_segment(data = edges_df, aes(x = Z1_from, y = Z2_from, xend = Z1_to, yend = Z2_to),
                 color = "grey10", linewidth = 0.5, alpha = 0.4) +
    # Plot the nodes (latent space points) with fully colored circles
    # geom_point(data = nodes_df, aes(x = Z1, y = Z2, color = as.factor(Cluster)),
    geom_point(data = nodes_df, aes(x = Z1, y = Z2),
                          size = 3, shape = 21, fill = model_cols[optimal_clustering_networks[m]]) +  # 'fill' will be colored based on clusters
    # scale_y_continuous(limits = c(cluster_limits[optimal_clustering_networks[m], "y_min"], cluster_limits[optimal_clustering_networks[m], "y_max"])) +
    # scale_x_continuous(limits = c(cluster_limits[optimal_clustering_networks[m], "x_min"], cluster_limits[optimal_clustering_networks[m], "x_max"])) +
    scale_y_continuous(breaks = seq(c(-1.5, -1, -1)[optimal_clustering_networks[m]], 
                                    c(1, 1, 1)[optimal_clustering_networks[m]], 
                                    c(0.5, 0.1, 0.5)[optimal_clustering_networks[m]])) +
    scale_x_continuous(breaks = seq(c(-1, -1, -1.5)[optimal_clustering_networks[m]], 
                                    c(2, 1, 0.5)[optimal_clustering_networks[m]],
                                    c(1, 0.1, 1)[optimal_clustering_networks[m]])) +
    # Set plot title, centered
    ggtitle(paste("Network", LETTERS[m])) +
    # Set theme with improved readability
    theme_minimal(base_size = 16) +
    theme(
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16),
      strip.text = element_text(size = 18, face = "bold", hjust = 0.5),  # Center titles
      panel.grid = element_blank(),  # Remove default grid
      panel.background = element_rect(fill = "white"),
      legend.position = "none",  # No legend
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),  # Center plot title
      panel.spacing = unit(1, "lines"),  # Space between panels
      strip.background = element_blank(),  # No background for facet labels
      plot.margin = margin(10, 10, 10, 10)  # Add margin around plot
    ) +
    labs(x = NULL, y = NULL) +  # Remove axis labels
    # scale_color_manual(values = model_cols) +  # Adjust colors based on clusters
    theme(panel.grid.major = element_line(color = "gray90", size = 0.5),  # Add grid lines
          panel.grid.minor = element_line(color = "gray95", size = 0.25))  # Add minor grid lines
  
  # Ensure fixed panel sizes
  p + facet_wrap(~ Scenario, nrow = 1, scales = "free_x") 
  
  # Store the plot in the list
  plot_list[[m]] = p
    
} # end m for loop
rm(m)

inds = c(which(optimal_clustering_networks == 1), which(optimal_clustering_networks == 2), which(optimal_clustering_networks == 3))
plot_list_reordered = plot_list[inds]

library(gridExtra)
# Combine the plots into a single multi-panel plot and save as PDF
pdf("Results_PDF_Plots/Krackhardt_Data_Clustered_LS_Coords_Overlay_Ordered_ggplot.pdf", width = 16, height = 8)
grid.arrange(grobs = plot_list_reordered, ncol = 7)  # Adjust the number of columns to 7
dev.off()

