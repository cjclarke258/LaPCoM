####################################################################################################
# AARHUS DATA APPLICATION (RESULTS)
####################################################################################################
rm(list = ls())
setwd("/Volumes/My Passport/PhD/PhD Project 1/LaPCoM_Code/Applications/Aarhus/")

# ================================================================================================
# LOAD THE DATA
# ================================================================================================
aarhus_multi = readRDS("aarhus_adjacency_matrices_updated.rds")

# load("pp_res_aarhus.Rdata")
load("pp_res_aarhus_PPMethodB.Rdata")

# ================================================================================================
# OPTIMAL CLUSTERING OF NETWORKS
# ================================================================================================
G_plus_hat = aarhus_pp$aarhus_G
optimal_clustering_networks = aarhus_pp$aarhus_C
table(optimal_clustering_networks)
cbind(optimal_clustering_networks, dimnames(aarhus_multi)[[3]])

# ------------------------------------------------------------------------------------------------
# INVESTIGATE A LITTLE
# ------------------------------------------------------------------------------------------------
for (g in 1:G_plus_hat) {
  g_inds = which(optimal_clustering_networks == g)
  dens = rep(NA, length(g_inds))
  for (m in 1:length(g_inds)) {
    dens[m] = network::network.density(network::as.network(aarhus_multi[, , g_inds[m]]))
  } # end m for loop
  cat("The average density of a network in cluster", g, "is: ", round(mean(dens), 2), "\n")
} # end g for loop 

# ================================================================================================
# OPTIMAL CLUSTERING OF NODES IN EACH LATENT SPACE
# ================================================================================================
Kg_plus_hat = aarhus_pp$aarhus_Kg_plus
Kg_plus_hat
optimal_clustering_nodes = vector("list", length = G_plus_hat)
for (g in 1:G_plus_hat) {
  optimal_clustering_nodes[[g]] = aarhus_pp$aarhus_S[[g]]
} # end g for loop
optimal_clustering_nodes

# ================================================================================================
# SELECT A COLOURBLIND FRIENDLY PALETTE
# ================================================================================================
# use Paul Tol's bright colour scheme
cols = khroma::color("muted")(9)
M = dim(aarhus_multi)[3]
set.seed(123)

# ------------------------------------------------------------------------------------------------
# PLOT THE NETWORKS WITH COLOUR CORRESPONDING TO CLUSTER
# ------------------------------------------------------------------------------------------------
# pdf("PDF_Plots/PPMethodC_Aarhus_Data_Clustered_Arbitrary_Coords.pdf", width = 10, height = 2)
# par(mfrow = c(1, 5))
# par(mar = c(1, 1, 3, 1)) # bltr
# for (m in 1:M) {
#   plot(network::as.network(aarhus_multi[, , m], directed = F), main = dimnames(aarhus_multi)[[3]][m], cex.main = 1.8,
#        edge.col = "grey", vertex.col = cols[optimal_clustering_networks[m]], vertex.cex = 3)
# } # end m for loop (plotting)
# dev.off()

# ------------------------------------------------------------------------------------------------
# PLOT THE NETWORKS 
# ------------------------------------------------------------------------------------------------
ind = which.max(apply(aarhus_multi, 3, sum))
library(statnet); coords = gplot(network::as.network(aarhus_multi[, , ind], directed = T))
pdf("PDF_Plots/Aarhus_Data_Unclustered_Same_Coords.pdf", width = 10, height = 2)
par(mfrow = c(1, 5))
par(mar = c(1, 1, 3, 1)) # bltr
for (m in 1:M) {
  plot(network::as.network(aarhus_multi[, , m], directed = F), main = dimnames(aarhus_multi)[[3]][m], cex.main = 1.5,
       edge.col = "grey", vertex.cex = 3,
       vertex.col = 1, 
       coord = coords)
} # end m for loop
dev.off()

# ------------------------------------------------------------------------------------------------
# PLOT THE NETWORKS WITH COLOUR CORRESPONDING TO CLUSTER (SAME COORDS)
# ------------------------------------------------------------------------------------------------
ind = which.max(apply(aarhus_multi, 3, sum))
library(statnet); coords = gplot(network::as.network(aarhus_multi[, , ind], directed = T))
pdf("PDF_Plots/PPMethodC_Aarhus_Data_Clustered_Same_Coords.pdf", width = 10, height = 2)
par(mfrow = c(1, 5))
par(mar = c(1, 1, 3, 1)) # bltr
for (g in 1:G_plus_hat) {
  for (m in 1:M) {
    if (optimal_clustering_networks[m] == g) {
      plot(network::as.network(aarhus_multi[, , m], directed = F), main = dimnames(aarhus_multi)[[3]][m], cex.main = 1.5,
           edge.col = "grey", vertex.cex = 3,
           vertex.col = cols[optimal_clustering_networks[m]], 
           coord = coords)
    } # end if statement
  } # end m for loop
} # end g for loop
dev.off()

# ------------------------------------------------------------------------------------------------
# OVERLAY THE NETWORKS WITHIN THEIR LATENT SPACE  
# ------------------------------------------------------------------------------------------------
pdf("PDF_Plots/Aarhus_Data_Clustered_LS_Coords_Overlay.pdf", width = 10, height = 2)
par(mfrow = c(1, 5))
# par(mar = c(1.5, 1, 1.5, 1)) # bltr
par(mar = c(1.75, 1.75, 1.5, 0.5) + 0.4) # bltr
for (g in 1:G_plus_hat) {
  for (m in 1:M) {
    if (optimal_clustering_networks[m] == g) {
      z1 = aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][, 1]
      z2 = aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][, 2]
      
      ones = which(aarhus_multi[, , m] != 0)
      
      plot(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]], type = "n", las = 1,
           ylim = c(min(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][, 2]) - 0.2, max(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][, 2]) + 0.2),
           xlim = c(min(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][, 1]) - 0.2, max(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][, 1]) + 0.2),
           main = dimnames(aarhus_multi)[[3]][m], cex.main = 1.8)
      grid()
      
      segments(z1[row(aarhus_multi[, , m])[ones]], z2[row(aarhus_multi[, , m])[ones]],
               z1[col(aarhus_multi[, , m])[ones]], z2[col(aarhus_multi[, , m])[ones]],
               lwd = 1, col = adjustcolor("grey10", 0.4))
      
      points(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]],
             bg = cols[optimal_clustering_networks[m]],
             pch = 21, cex = 1.5)
    } # end if statement
  } # end m for loop (plotting)
} # end g for loop
dev.off()

# ------------------------------------------------------------------------------------------------
# OVERLAY THE NETWORKS WITHIN THEIR LATENT SPACE, REORDERED BY CLUSTER
# ------------------------------------------------------------------------------------------------
pdf("PDF_Plots/Aarhus_Data_Clustered_LS_Coords_Overlay_Ordered_and_Coloured.pdf", width = 10, height = 6.5)
# par(mfrow = c(1, 5))
# par(mar = c(1.5, 1, 1.5, 1)) # bltr
layout(matrix(1:6, nrow = 2, ncol = 3, byrow = TRUE), widths = rep(1, 5), heights = 1)
par(mar = c(1.75, 1.75, 1.5, 0.5) + 0.4) # bltr
for (g in 1:G_plus_hat) {
  
  for (m in 1:M) {
    
    if (optimal_clustering_networks[m] == g) {
      
      z1 = aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][, 1]
      z2 = aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][, 2]
      
      ones = which(aarhus_multi[, , m] != 0)
      
      plot(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]], type = "n", las = 1, 
           ylim = c(min(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][, 2]) - 0.2, max(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][, 2]) + 0.2),
           xlim = c(min(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][, 1]) - 0.2, max(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][, 1]) + 0.2),
           main = dimnames(aarhus_multi)[[3]][m], cex.main = 1.5, cex.axis = 0.8)
      grid()
      
      segments(z1[row(aarhus_multi[, , m])[ones]], z2[row(aarhus_multi[, , m])[ones]], 
               z1[col(aarhus_multi[, , m])[ones]], z2[col(aarhus_multi[, , m])[ones]], 
               lwd = 1, col = adjustcolor("grey10", 0.4))
      
      points(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]], 
             bg = cols[3:9][optimal_clustering_nodes[[g]]], 
             pch = 21, cex = 1.5)
      
      box(col = cols[c(1:2)[g]], lwd = 3)
      
    } # end if statement
    
  } # end m for loop 
  
} # end g for loop

plot.new()  # Create an empty plot in the 6th space

# Add the legend for box colors (to the right of the plots)
legend("top", legend = c("Network-level Cluster 1", "Network-level Cluster 2"), lwd = 3, cex = 1.2,
       title = expression(atop(underline("Box Colour"))), col = cols[1:2], bty = "n", xpd = TRUE, inset = c(0, 0))

# Add the legend for point colors (roles)
legend("bottom", legend = c(paste0("Node-level Cluster ", 1:max(Kg_plus_hat))), pch = 19, cex = 1.2,
       col = cols[3:9], title = expression(atop(underline("Point Colour"))), bty = "n", xpd = TRUE, inset = c(0, 0))

dev.off()

# ------------------------------------------------------------------------------------------------
# IS THERE ANY CORRESPONDENCE WITH THE EMPLOYEE ROLE AND THE NODE-LEVEL CLUSTERING?
# ------------------------------------------------------------------------------------------------
load("aarhus_cs_network_Michael.RData"); aarhus_full_Michael = out; rm(out)

inds_remove = which(aarhus_full_Michael$meta$id %in% setdiff(aarhus_full_Michael$meta$id, dimnames(aarhus_multi)[[1]]))
aarhus_full_Michael$meta = aarhus_full_Michael$meta[-inds_remove, ]

aarhus_full_Michael$meta$value = ifelse(aarhus_full_Michael$meta$value == "Phd (visiting)", "PhD", aarhus_full_Michael$meta$value)

employee_roles = aarhus_full_Michael$meta$value
employee_roles_factor = factor(aarhus_full_Michael$meta$value, levels = unique(employee_roles), labels = c(1:length(unique(employee_roles))))

N = dim(aarhus_multi)[1]

# 

# pdf("PDF_Plots/Aarhus_Data_Clustered_LS_Coords_Overlay_Ordered_and_Coloured_with_Shape_by_Role.pdf", width = 10, height = 6.5)
# # par(mfrow = c(1, 5))
# layout(matrix(1:6, nrow = 2, ncol = 3, byrow = TRUE), widths = rep(1, 5), heights = 1)
# # par(mar = c(1.5, 1, 1.5, 1)) # bltr
# par(mar = c(1.75, 1.75, 1.5, 0.5) + 0.4) # bltr
# for (g in 1:G_plus_hat) {
#   
#   for (m in 1:M) {
#     
#     if (optimal_clustering_networks[m] == g) {
#       
#       z1 = aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][, 1]
#       z2 = aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][, 2]
#       
#       ones = which(aarhus_multi[, , m] != 0)
#       
#       plot(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]], type = "n", las = 1, 
#            ylim = c(min(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][, 2]) - 0.2, max(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][, 2]) + 0.2),
#            xlim = c(min(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][, 1]) - 0.2, max(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][, 1]) + 0.2),
#            main = dimnames(aarhus_multi)[[3]][m], cex.main = 1.5, cex.axis = 0.8,
#            xlab = "Dimension 1", ylab = "Dimension 2")
#       grid()
#       
#       segments(z1[row(aarhus_multi[, , m])[ones]], z2[row(aarhus_multi[, , m])[ones]], 
#                z1[col(aarhus_multi[, , m])[ones]], z2[col(aarhus_multi[, , m])[ones]], 
#                lwd = 1, col = adjustcolor("grey10", 0.4))
#       
#       points(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]], 
#              bg = cols[3:9][optimal_clustering_nodes[[g]]], 
#              pch = c(21:25, 4)[employee_roles_factor], cex = 1.5)
#       
#       box(col = cols[c(1:2)[g]], lwd = 3)
#       
#     } # end if statement
#     
#   } # end m for loop 
#   
# } # end g for loop
# 
# plot.new()  # Create an empty plot in the 6th space
# 
# # Add the legend for box colors (to the right of the plots)
# legend("top", legend = c("Network-level Cluster 1", "Network-level Cluster 2"), lwd = 3, cex = 1.2,
#        title = expression(atop(underline("Box Colour"))), col = cols[1:2], bty = "n", xpd = TRUE, inset = c(0, 0))
# 
# # Add the legend for point colors (cluster)
# legend("bottomleft", legend = paste0("Node-level Cluster ", 1:max(Kg_plus_hat)), cex = 1.2,
#        fill = cols[3:9], title = expression(atop(underline("Point Colour"))), bty = "n", xpd = TRUE, inset = c(0, 0))
# 
# # Add the legend for point shapes (roles)
# legend("bottomright", legend = unique(employee_roles), pch = c(21:25, 4), cex = 1.2,
#        title = expression(atop(underline("Point Shape"))), bty = "n", xpd = TRUE, inset = c(0, 0))
# 
# dev.off()

# 

employee_roles_broader = aarhus_full_Michael$meta$value
employee_roles_broader = ifelse(employee_roles_broader %in% c("Professor", "Associate"), "Professor", employee_roles_broader)
employee_roles_broader = ifelse(employee_roles_broader %in% c("PhD", "Postdoc"), "PhD/Postdoc", employee_roles_broader)

employee_roles_factor_broader = factor(employee_roles_broader, levels = unique(employee_roles_broader), labels = c(1:length(unique(employee_roles_broader))))

pdf("PDF_Plots/Aarhus_Data_Clustered_LS_Coords_Overlay_Ordered_and_Coloured_with_Shape_by_Role_Broader.pdf", width = 10, height = 6.5)
# par(mfrow = c(1, 5))
layout(matrix(1:6, nrow = 2, ncol = 3, byrow = TRUE), widths = rep(1, 5), heights = 1)
# par(mar = c(1.5, 1, 1.5, 1)) # bltr
par(mar = c(1.75, 1.75, 1.5, 0.5) + 0.4) # bltr
for (g in 1:G_plus_hat) {
  
  for (m in 1:M) {
    
    if (optimal_clustering_networks[m] == g) {
      
      z1 = aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][, 1]
      z2 = aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][, 2]
      
      ones = which(aarhus_multi[, , m] != 0)
      
      plot(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]], type = "n", las = 1, 
           ylim = c(min(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][, 2]) - 0.2, max(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][, 2]) + 0.2),
           xlim = c(min(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][, 1]) - 0.2, max(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][, 1]) + 0.2),
           main = dimnames(aarhus_multi)[[3]][m], cex.main = 1.5, cex.axis = 0.8,
           xlab = "Dimension 1", ylab = "Dimension 2")
      grid()
      
      segments(z1[row(aarhus_multi[, , m])[ones]], z2[row(aarhus_multi[, , m])[ones]], 
               z1[col(aarhus_multi[, , m])[ones]], z2[col(aarhus_multi[, , m])[ones]], 
               lwd = 1, col = adjustcolor("grey10", 0.4))
      
      points(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]], 
             bg = cols[3:9][optimal_clustering_nodes[[g]]], 
             pch = c(21:23, 4)[employee_roles_factor_broader], 
             cex = 1.5)
      
      box(col = cols[c(1:2)[g]], lwd = 3)
      
    } # end if statement
    
  } # end m for loop 
  
} # end g for loop

plot.new()  # Create an empty plot in the 6th space

# Add the legend for box colors (to the right of the plots)
legend("top", legend = c("Network-level Cluster 1", "Network-level Cluster 2"), lwd = 3, cex = 1.2,
       title = expression(atop(underline("Box Colour"))), col = cols[1:2], bty = "n", xpd = TRUE, inset = c(0, 0))

# Add the legend for point colors (cluster)
legend("bottomleft", legend = paste0("Node-level Cluster ", 1:max(Kg_plus_hat)), cex = 1.2,
       fill = cols[3:9], title = expression(atop(underline("Point Colour"))), bty = "n", xpd = TRUE, inset = c(0, 0))

# Add the legend for point shapes (roles)
legend("bottomright", legend = unique(employee_roles_broader), pch = c(21:23, 4), cex = 1.2,
       title = expression(atop(underline("Point Shape"))), bty = "n", xpd = TRUE, inset = c(0, 0.1175))

dev.off()

# 

# USE WG OPTIAL S2 FORM CHAIN 31 TO SHOW UNCONNECTED VS CONNECTED CLUSTER
optimal_clustering_nodes[[2]] = readRDS(file = "Chain_31_WG_Optimal_S2.RDS")

pdf("PDF_Plots/Aarhus_Data_Clustered_LS_Coords_Overlay_Ordered_and_Coloured_with_Shape_by_Role_Broader_WG_S2.pdf", width = 10, height = 6.5)
# par(mfrow = c(1, 5))
layout(matrix(1:6, nrow = 2, ncol = 3, byrow = TRUE), widths = rep(1, 5), heights = 1)
# par(mar = c(1.5, 1, 1.5, 1)) # bltr
par(mar = c(1.75, 1.75, 1.5, 0.5) + 0.4) # bltr
for (g in 1:G_plus_hat) {
  
  for (m in 1:M) {
    
    if (optimal_clustering_networks[m] == g) {
      
      z1 = aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][, 1]
      z2 = aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][, 2]
      
      ones = which(aarhus_multi[, , m] != 0)
      
      plot(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]], type = "n", las = 1, 
           ylim = c(min(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][, 2]) - 0.025, max(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][, 2]) + 0.025),
           xlim = c(min(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][, 1]) - 0.025, max(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][, 1]) + 0.025),
           main = dimnames(aarhus_multi)[[3]][m], cex.main = 1.5, cex.axis = 0.8,
           xlab = "Dimension 1", ylab = "Dimension 2")
      grid()
      
      segments(z1[row(aarhus_multi[, , m])[ones]], z2[row(aarhus_multi[, , m])[ones]], 
               z1[col(aarhus_multi[, , m])[ones]], z2[col(aarhus_multi[, , m])[ones]], 
               lwd = 1, col = adjustcolor("grey10", 0.4))
      
      points(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]], 
             bg = cols[3:9][optimal_clustering_nodes[[g]]], 
             pch = c(21:23, 4)[employee_roles_factor_broader], 
             cex = 1.5)
      
      box(col = cols[c(1:2)[g]], lwd = 3)
      
    } # end if statement
    
  } # end m for loop 
  
} # end g for loop

plot.new()  # Create an empty plot in the 6th space

# Add the legend for box colors (to the right of the plots)
legend("top", legend = c("Network-level Cluster 1", "Network-level Cluster 2"), lwd = 3, cex = 1.2,
       title = expression(atop(underline("Box Colour"))), col = cols[1:2], bty = "n", xpd = TRUE, inset = c(0, 0))

# Add the legend for point colors (cluster)
legend("bottomleft", legend = paste0("Node-level Cluster ", 1:max(Kg_plus_hat)), cex = 1.2,
       fill = cols[3:9], title = expression(atop(underline("Point Colour"))), bty = "n", xpd = TRUE, inset = c(0, 0))

# Add the legend for point shapes (roles)
legend("bottomright", legend = unique(employee_roles_broader), pch = c(21:23, 4), cex = 1.2,
       title = expression(atop(underline("Point Shape"))), bty = "n", xpd = TRUE, inset = c(0, 0.1175))

dev.off()

# ------------------------------------------------------------------------------------------------
# LOCATE WHERE CERTAIN POINTS ARE IN BOTH CLUSTERS
# ------------------------------------------------------------------------------------------------
# pdf("PDF_Plots/PPMethodC_Aarhus_Data_Clustered_LS_Coords_Overlay_Ordered_Test_Unconnected_Point.pdf", width = 10, height = 2)
# par(mfrow = c(1, 5))
# par(mar = c(1.75, 1.75, 1.5, 0.5) + 0.4) # bltr
# g = 2
# for (m in 1:M) {
#   
#   if (optimal_clustering_networks[m] == g) {
#     
#     z1 = aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][, 1]
#     z2 = aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][, 2]
#     
#     ones = which(aarhus_multi[, , m] != 0)
#     
#     plot(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]], type = "n", las = 1, 
#          ylim = c(min(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][, 2]) - 0.2, max(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][, 2]) + 0.2),
#          xlim = c(min(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][, 1]) - 0.2, max(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][, 1]) + 0.2),
#          main = dimnames(aarhus_multi)[[3]][m], cex.main = 1.5, cex.axis = 0.8)
#     grid()
#     
#     segments(z1[row(aarhus_multi[, , m])[ones]], z2[row(aarhus_multi[, , m])[ones]], 
#              z1[col(aarhus_multi[, , m])[ones]], z2[col(aarhus_multi[, , m])[ones]], 
#              lwd = 1, col = adjustcolor("grey10", 0.4))
#     
#     points(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]], 
#            bg = cols[optimal_clustering_networks[m]], 
#            pch = 21, cex = 1.5)
#     
#     point = which.max(z1)
#     points(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][point, 1],
#            aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][point, 2],
#            bg = "green",
#            pch = 21, cex = 1.5)
#     
#   } # end if statement
#   
# } # end m for loop
# 
# g = 1
# for (m in 1:M) {
#   
#   if (optimal_clustering_networks[m] == g) {
#     
#     z1 = aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][, 1]
#     z2 = aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][, 2]
#     
#     ones = which(aarhus_multi[, , m] != 0)
#     
#     plot(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]], type = "n", las = 1, 
#          ylim = c(min(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][, 2]) - 0.2, max(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][, 2]) + 0.2),
#          xlim = c(min(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][, 1]) - 0.2, max(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][, 1]) + 0.2),
#          main = dimnames(aarhus_multi)[[3]][m], cex.main = 1.5, cex.axis = 0.8)
#     grid()
#     
#     segments(z1[row(aarhus_multi[, , m])[ones]], z2[row(aarhus_multi[, , m])[ones]], 
#              z1[col(aarhus_multi[, , m])[ones]], z2[col(aarhus_multi[, , m])[ones]], 
#              lwd = 1, col = adjustcolor("grey10", 0.4))
#     
#     points(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]], 
#            bg = cols[optimal_clustering_networks[m]], 
#            pch = 21, cex = 1.5)
#     
#     points(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][point, 1], 
#            aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][point, 2], 
#            bg = "green",
#            pch = 21, cex = 1.5)
#     
#   } # end if statement
#   
# } # end m for loop
# dev.off()
# 
# # 
# 
# pdf("PDF_Plots/PPMethodC_Aarhus_Data_Clustered_LS_Coords_Overlay_Ordered_Test_Connected_Blob_Facebook.pdf", width = 10, height = 2)
# par(mfrow = c(1, 5))
# par(mar = c(1.75, 1.75, 1.5, 0.5) + 0.4) # bltr
# g = 1
# for (m in 1:M) {
#   
#   if (optimal_clustering_networks[m] == g) {
#     
#     z1 = aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][, 1]
#     z2 = aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][, 2]
#     
#     ones = which(aarhus_multi[, , m] != 0)
#     
#     plot(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]], type = "n", las = 1, 
#          ylim = c(min(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][, 2]) - 0.2, max(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][, 2]) + 0.2),
#          xlim = c(min(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][, 1]) - 0.2, max(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][, 1]) + 0.2),
#          main = dimnames(aarhus_multi)[[3]][m], cex.main = 1.5, cex.axis = 0.8)
#     grid()
#     
#     segments(z1[row(aarhus_multi[, , m])[ones]], z2[row(aarhus_multi[, , m])[ones]], 
#              z1[col(aarhus_multi[, , m])[ones]], z2[col(aarhus_multi[, , m])[ones]], 
#              lwd = 1, col = adjustcolor("grey10", 0.4))
#     
#     points(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]], 
#            bg = cols[optimal_clustering_networks[m]], 
#            pch = 21, cex = 1.5)
#     
#     blob = which((z1 < -0.4) & (z1 > -1) & (z2 < -0.5) & (z2 > -1) )
#     points(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][blob, ], 
#            bg = "green",
#            pch = 21, cex = 1.5)
#     
#   } # end if statement
#   
# } # end m for loop
# 
# g = 2
# for (m in 1:M) {
#   
#   if (optimal_clustering_networks[m] == g) {
#     
#     z1 = aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][, 1]
#     z2 = aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][, 2]
#     
#     ones = which(aarhus_multi[, , m] != 0)
#     
#     plot(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]], type = "n", las = 1, 
#          ylim = c(min(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][, 2]) - 0.2, max(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][, 2]) + 0.2),
#          xlim = c(min(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][, 1]) - 0.2, max(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][, 1]) + 0.2),
#          main = dimnames(aarhus_multi)[[3]][m], cex.main = 1.5, cex.axis = 0.8)
#     grid()
#     
#     segments(z1[row(aarhus_multi[, , m])[ones]], z2[row(aarhus_multi[, , m])[ones]], 
#              z1[col(aarhus_multi[, , m])[ones]], z2[col(aarhus_multi[, , m])[ones]], 
#              lwd = 1, col = adjustcolor("grey10", 0.4))
#     
#     points(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]], 
#            bg = cols[optimal_clustering_networks[m]], 
#            pch = 21, cex = 1.5)
#     
#     points(aarhus_pp$aarhus_Z[[optimal_clustering_networks[m]]][blob, ], 
#            bg = "green",
#            pch = 21, cex = 1.5)
#     
#   } # end if statement
#   
# } # end m for loop
# dev.off()