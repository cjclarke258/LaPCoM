####################################################################################################
# PRIMARY SCHOOL (COUNT) DATA APPLICATION (RESULTS)
####################################################################################################
rm(list = ls())
setwd("/Volumes/My Passport/PhD/PhD Project 1/LaPCoM_Code/Applications/Primary_School_Count/")

# ================================================================================================
# LOAD THE DATA
# ================================================================================================
primary_school_multi = readRDS("primary_school_adj_matrices_count.rds")

load("pp_res_psc.Rdata")

# ================================================================================================
# OPTIMAL CLUSTERING OF NETWORKS
# ================================================================================================
G_plus_hat = primary_school_pp$primary_school_G
optimal_clustering_networks = primary_school_pp$primary_school_C
table(optimal_clustering_networks)
cbind(optimal_clustering_networks, dimnames(primary_school_multi)[[3]])

# ------------------------------------------------------------------------------------------------
# INVESTIGATE A LITTLE
# ------------------------------------------------------------------------------------------------
for (g in 1:G_plus_hat) {
  g_inds = which(optimal_clustering_networks == g)
  dens = rep(NA, length(g_inds))
  for (m in 1:length(g_inds)) {
    dens[m] = network::network.density(network::as.network(primary_school_multi[, , g_inds[m]]))
  } # end m for loop
  cat("The average density of a network in cluster", g, "is: ", round(mean(dens), 2), "\n")
} # end g for loop 

# ================================================================================================
# OPTIMAL CLUSTERING OF NODES IN EACH LATENT SPACE
# ================================================================================================
Kg_plus_hat = primary_school_pp$primary_school_Kg_plus
Kg_plus_hat
optimal_clustering_nodes = vector("list", length = G_plus_hat)
for (g in 1:G_plus_hat) {
  optimal_clustering_nodes[[g]] = primary_school_pp$primary_school_S[[g]]
} # end g for loop
optimal_clustering_nodes

lapply(optimal_clustering_nodes, table)

# ================================================================================================
# SELECT A COLOURBLIND FRIENDLY PALETTE
# ================================================================================================
# use Paul Tol's bright colour scheme
cols = khroma::color("muted")(9)[3:9]
cols_boxes = khroma::color("muted")(9)[1:2]
M = dim(primary_school_multi)[3]
set.seed(123)

# ------------------------------------------------------------------------------------------------
# PLOT THE NETWORKS WITH COLOUR CORRESPONDING TO CLUSTER
# ------------------------------------------------------------------------------------------------
pdf("PDF_Plots/PSC_Clustered_Arbitrary_Coords.pdf", width = 8, height = 8)
par(mfrow = c(4, 4))
par(mar = c(1, 1, 3, 1)) # bltr
for (m in 1:M) {
  plot(network::as.network(primary_school_multi[, , m], directed = F), 
       main = dimnames(primary_school_multi)[[3]][m], cex.main = 1,
       edge.col = "grey", vertex.col = cols[optimal_clustering_networks[m]], vertex.cex = 3)
} # end m for loop (plotting)
dev.off()

# ------------------------------------------------------------------------------------------------
# PLOT THE NETWORKS WITH COLOUR CORRESPONDING TO CLUSTER (SAME COORDS)
# ------------------------------------------------------------------------------------------------
ind = which.max(apply(primary_school_multi, 3, sum))
library(statnet); coords = gplot(network::as.network(primary_school_multi[, , ind], directed = T))
pdf("PDF_Plots/PSC_Clustered_Same_Coords.pdf", width = 8, height = 8)
par(mfrow = c(4, 4))
par(mar = c(1, 1, 3, 1)) # bltr
for (m in 1:M) {
  plot(network::as.network(primary_school_multi[, , m], directed = F), 
       main = dimnames(primary_school_multi)[[3]][m], cex.main = 1,
       edge.col = "grey", vertex.cex = 3,
       vertex.col = cols[optimal_clustering_networks[m]], 
       coord = coords)
} # end m for loop
dev.off()

# ------------------------------------------------------------------------------------------------
# OVERLAY THE NETWORKS WITHIN THEIR LATENT SPACE  
# ------------------------------------------------------------------------------------------------
pdf("PDF_Plots/PSC_Clustered_LS_Coords_Overlay.pdf", width = 8, height = 8)
par(mfrow = c(4, 4))
# par(mar = c(1.5, 1, 1.5, 1)) # bltr
par(mar = c(1.75, 1.75, 1.5, 0.5) + 0.4) # bltr
for (m in 1:M) {
  
  z1 = primary_school_pp$primary_school_Z[[optimal_clustering_networks[m]]][, 1]
  z2 = primary_school_pp$primary_school_Z[[optimal_clustering_networks[m]]][, 2]
  
  ones = which(primary_school_multi[, , m] != 0)
  
  plot(primary_school_pp$primary_school_Z[[optimal_clustering_networks[m]]], type = "n", las = 1, 
       ylim = c(min(primary_school_pp$primary_school_Z[[optimal_clustering_networks[m]]][, 2]) - 0.2, max(primary_school_pp$primary_school_Z[[optimal_clustering_networks[m]]][, 2]) + 0.2),
       xlim = c(min(primary_school_pp$primary_school_Z[[optimal_clustering_networks[m]]][, 1]) - 0.2, max(primary_school_pp$primary_school_Z[[optimal_clustering_networks[m]]][, 1]) + 0.2),
       main = dimnames(primary_school_multi)[[3]][m], cex.main = 1)
  grid()
  
  segments(z1[row(primary_school_multi[, , m])[ones]], z2[row(primary_school_multi[, , m])[ones]], 
           z1[col(primary_school_multi[, , m])[ones]], z2[col(primary_school_multi[, , m])[ones]], 
           lwd = 1, col = adjustcolor("grey60", 0.4))
  
  points(primary_school_pp$primary_school_Z[[optimal_clustering_networks[m]]], 
         bg = cols[optimal_clustering_networks[m]], 
         pch = 21, cex = 1)
  
} # end m for loop (plotting)
dev.off()

# ------------------------------------------------------------------------------------------------
# OVERLAY THE NETWORKS WITHIN THEIR LATENT SPACE, REORDERED BY CLUSTER
# ------------------------------------------------------------------------------------------------
pdf("PDF_Plots/PSC_Clustered_LS_Coords_Overlay_Ordered_and_Coloured.pdf", width = 10, height = 12)
# par(mfrow = c(1, 5))
# par(mar = c(1.5, 1, 1.5, 1)) # bltr
cols = rep(khroma::color("light")(8), 2)

layout(matrix(1:20, nrow = 5, ncol = 4, byrow = TRUE), widths = rep(1, 20), heights = 1)
par(mar = c(1.75, 1.75, 1.5, 0.5) + 0.4) # bltr
for (g in 1:G_plus_hat) {
  
  for (m in 1:M) {
    
    if (optimal_clustering_networks[m] == g) {
      
      z1 = primary_school_pp$primary_school_Z[[optimal_clustering_networks[m]]][, 1]
      z2 = primary_school_pp$primary_school_Z[[optimal_clustering_networks[m]]][, 2]
      
      ones = which(primary_school_multi[, , m] != 0)
      
      plot(primary_school_pp$primary_school_Z[[optimal_clustering_networks[m]]], type = "n", las = 1, 
           ylim = c(min(primary_school_pp$primary_school_Z[[optimal_clustering_networks[m]]][, 2]) - 0.2, max(primary_school_pp$primary_school_Z[[optimal_clustering_networks[m]]][, 2]) + 0.2),
           xlim = c(min(primary_school_pp$primary_school_Z[[optimal_clustering_networks[m]]][, 1]) - 0.2, max(primary_school_pp$primary_school_Z[[optimal_clustering_networks[m]]][, 1]) + 0.2),
           main = dimnames(primary_school_multi)[[3]][m], cex.main = 1.5, cex.axis = 1)
      grid()
      
      segments(z1[row(primary_school_multi[, , m])[ones]], z2[row(primary_school_multi[, , m])[ones]], 
               z1[col(primary_school_multi[, , m])[ones]], z2[col(primary_school_multi[, , m])[ones]], 
               lwd = 1, col = adjustcolor("grey60", 0.4))
      
      points(primary_school_pp$primary_school_Z[[optimal_clustering_networks[m]]], 
             bg = cols[optimal_clustering_nodes[[g]]], 
             pch = 21, cex = 1.5)
      
      # box(col = cols[c(1:2)[g]], lwd = 3)
      box(col = cols_boxes[c(1:2)[g]], lwd = 3)
      
    } # end if statement
    
  } # end m for loop 
  
} # end g for loop

plot.new()  # Create an empty plot in the 6th space

# Add the legend for box colors (to the right of the plots)
legend("center", legend = c("Network-level Cluster 1", "Network-level Cluster 2"), lwd = 3, cex = 1.5,
       title = expression(atop(underline("Box Colour"))), col = cols_boxes[1:2], bty = "n", xpd = TRUE, inset = c(0, 0))

plot.new()  # Create an empty plot in the 6th space

# Add the legend for point colors (roles)
legend("center", legend = c(paste0("Node-level Cluster ", 1:6)), pch = 19, cex = 1.5,
       col = cols[1:6], title = expression(atop(underline("Point Colour"))), bty = "n", xpd = TRUE, inset = c(0, 0))

plot.new()  # Create an empty plot in the 6th space

# Add the legend for point colors (roles)
legend("center", legend = c(paste0("Node-level Cluster ", 7:max(Kg_plus_hat))), pch = 19, cex = 1.5,
       col = cols[7:13], title = "", bty = "n", xpd = TRUE, inset = c(0, 0))

dev.off()

#

# ------------------------------------------------------------------------------------------------
# IS THERE ANY CORRESPONDENCE WITH THE ROLE/CLASS AND THE NODE-LEVEL CLUSTERING?
# ------------------------------------------------------------------------------------------------
meta = read.table("metadata_primaryschool.txt", col.names = c("ID", "Class", "Gender"))
meta$Class = as.factor(meta$Class)
meta$Gender = as.factor(meta$Gender)

mclust::adjustedRandIndex(meta$Class, primary_school_pp$primary_school_S[[1]])
mclust::adjustedRandIndex(meta$Class, primary_school_pp$primary_school_S[[2]])

mclust::adjustedRandIndex(meta$Gender, primary_school_pp$primary_school_S[[1]])
mclust::adjustedRandIndex(meta$Gender, primary_school_pp$primary_school_S[[2]])

nrow(meta)

N = dim(primary_school_multi)[1]

pdf("PDF_Plots/PSC_Clustered_LS_Coords_Overlay_Ordered_and_Coloured_by_Class.pdf", width = 10, height = 12)
layout(matrix(1:20, nrow = 5, ncol = 4, byrow = TRUE), widths = rep(1, 20), heights = 1)
par(mar = c(1.75, 1.75, 1.5, 0.5) + 0.4) # bltr
for (g in 1:G_plus_hat) {
  
  for (m in 1:M) {
    
    if (optimal_clustering_networks[m] == g) {
      
      z1 = primary_school_pp$primary_school_Z[[optimal_clustering_networks[m]]][, 1]
      z2 = primary_school_pp$primary_school_Z[[optimal_clustering_networks[m]]][, 2]
      
      ones = which(primary_school_multi[, , m] != 0)
      
      plot(primary_school_pp$primary_school_Z[[optimal_clustering_networks[m]]], type = "n", las = 1, 
           ylim = c(min(primary_school_pp$primary_school_Z[[optimal_clustering_networks[m]]][, 2]) - 0.2, max(primary_school_pp$primary_school_Z[[optimal_clustering_networks[m]]][, 2]) + 0.2),
           xlim = c(min(primary_school_pp$primary_school_Z[[optimal_clustering_networks[m]]][, 1]) - 0.2, max(primary_school_pp$primary_school_Z[[optimal_clustering_networks[m]]][, 1]) + 0.2),
           main = dimnames(primary_school_multi)[[3]][m], cex.main = 1.5, cex.axis = 1)
      grid()
      
      segments(z1[row(primary_school_multi[, , m])[ones]], z2[row(primary_school_multi[, , m])[ones]], 
               z1[col(primary_school_multi[, , m])[ones]], z2[col(primary_school_multi[, , m])[ones]], 
               lwd = 1, col = adjustcolor("grey60", 0.4))
      
      points(primary_school_pp$primary_school_Z[[optimal_clustering_networks[m]]], 
             bg = cols[meta$Class], 
             # pch = c(21, 24, 18)[meta$Gender], 
             pch = 21,
             cex = 1.5)
      
      # box(col = cols[c(1:2)[g]], lwd = 3)
      box(col = cols_boxes[c(1:2)[g]], lwd = 3)
      
    } # end if statement
    
  } # end m for loop 
  
} # end g for loop

plot.new()  # Create an empty plot in the 6th space

# Add the legend for box colors (to the right of the plots)
legend("center", legend = c("Network-level Cluster 1", "Network-level Cluster 2"), lwd = 3, cex = 1.5,
       title = expression(atop(underline("Box Colour"))), col = cols_boxes[1:2], bty = "n", xpd = TRUE, inset = c(0, 0))

plot.new()  # Create an empty plot in the 6th space

# Add the legend for point colors (class)
# legend("center", legend = levels(meta$Class), pch = 19, cex = 1.2,
#        col = cols[1:length(levels(meta$Class))], title = expression(atop(underline("Point Colour"))), bty = "n", xpd = TRUE, inset = c(0, 0))

# Add the legend for point colors (roles)
legend("center", legend = levels(meta$Class)[1:6], pch = 19, cex = 1.5,
       col = cols[1:6], title = expression(atop(underline("Point Colour"))), bty = "n", xpd = TRUE, inset = c(0, 0))

plot.new()  # Create an empty plot in the 6th space

# Add the legend for point colors (roles)
legend("center", legend = levels(meta$Class)[7:11], pch = 19, cex = 1.5,
       col = cols[7:11], title = "", bty = "n", xpd = TRUE, inset = c(0, 0))

plot.new()  # Create an empty plot in the 6th space

# Add the legend for point shapes (gender)
# legend("center", legend = levels(meta$Gender), pch = c(21, 24, 18), cex = 1.2,
#        col = c(1, 1, 1), title = expression(atop(underline("Point Shape"))), bty = "n", xpd = TRUE, inset = c(0, 0))

dev.off()

# 

# ------------------------------------------------------------------------------------------------
# CLUSTER 1 ZOOMED IN
# ------------------------------------------------------------------------------------------------
pdf("PDF_Plots/PSC_Clustered_LS_Coords_Overlay_Ordered_and_Coloured_by_Class_C1_Zoomed.pdf", width = 15, height = 3)
layout(matrix(1:5, nrow = 1, ncol = 5, byrow = TRUE), widths = rep(1, 5), heights = 1)
par(mar = c(1.75, 1.75, 1.5, 0.5) + 0.4) # bltr
g = 1  
for (m in 1:M) {
  
  if (optimal_clustering_networks[m] == g) {
    
    z1 = primary_school_pp$primary_school_Z[[optimal_clustering_networks[m]]][, 1]
    z2 = primary_school_pp$primary_school_Z[[optimal_clustering_networks[m]]][, 2]
    
    ones = which(primary_school_multi[, , m] != 0)
    
    plot(primary_school_pp$primary_school_Z[[optimal_clustering_networks[m]]], type = "n", las = 1, 
         ylim = c(-1.5, 1.5), xlim = c(-2.5, 0),
         main = dimnames(primary_school_multi)[[3]][m], cex.main = 1.5, cex.axis = 1)
    grid()
    
    segments(z1[row(primary_school_multi[, , m])[ones]], z2[row(primary_school_multi[, , m])[ones]], 
             z1[col(primary_school_multi[, , m])[ones]], z2[col(primary_school_multi[, , m])[ones]], 
             lwd = 1, col = adjustcolor("grey60", 0.4))
    
    points(primary_school_pp$primary_school_Z[[optimal_clustering_networks[m]]], 
           bg = cols[meta$Class], 
           # pch = c(21, 24, 18)[meta$Gender], 
           pch = 21,
           cex = 1)
    
    # box(col = cols[c(1:2)[g]], lwd = 3)
    box(col = cols_boxes[c(1:2)[g]], lwd = 3)
    
  } # end if statement
  
} # end m for loop 
  
plot.new()  # Create an empty plot in the 6th space

# Add the legend for point colors (class)
legend("center", legend = levels(meta$Class), pch = 19, cex = 1.2,
       col = cols[1:length(levels(meta$Class))], title = expression(atop(underline("Point Colour"))), bty = "n", xpd = TRUE, inset = c(0, 0))

dev.off()