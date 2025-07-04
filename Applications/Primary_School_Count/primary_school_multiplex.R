# ######################################################################################################################
# DATASET: Primary school temporal network data
# 
# Release data: Sep 30, 2015
# 
# This data set contains the temporal network of contacts between the children and teachers used in the study published 
# in BMC Infectious Diseases 2014, 14:695. The file contains a tab-separated list representing the active contacts 
# during 20-second intervals of the data collection. Each line has the form “t i j Ci Cj”, where i and j are the 
# anonymous IDs of the persons in contact, Ci and Cj are their classes, and the interval during which this contact was 
# active is [ t – 20s, t ]. If multiple contacts are active in a given interval, you will see multiple lines starting 
# with the same value of t. Time is measured in seconds.
# 
# Terms and conditions
# The data are distributed to the public under a Creative Commons Attribution-NonCommercial-ShareAlike license. When
# this data is used in published research or for visualization purposes, please cite the following papers:
#   
# Mitigation of infectious disease at school: targeted class closure vs school closure, 
# BMC Infectious Diseases 14:695 (2014)
# 
# primary-Resolution Measurements of Face-to-Face Contact Patterns in a Primary School, PLOS ONE 6(8): e23176 (2011)
# 
# Please also acknowledge the SocioPatterns collaboration and provide a link to http://www.sociopatterns.org.
# 
# 
#  
# The time-resolved contact network considered here, analyzed in Ref. [38], describes the contacts among 232 children 
# and 10 teachers in a primary school in Lyon, France, and covers two days of school activity (Thursday, October 1st and 
# Friday, October 2nd 2009). The school is composed by 5 grades, each of them comprising two classes, for a total of 10 
# classes. Contacts events are individually resolved, and their starting and ending times are known up to the 20-second 
# resolution of the measurement system.
# ######################################################################################################################

rm(list = ls())
setwd("/Volumes/My Passport/PhD/PhD Project 1/LaPCoM_Code/Applications/Primary_School_Count/")

# ======================================================================================================================
# LOAD THE DATA AND FORMAT IT 
# ======================================================================================================================
contacts = read.table("primaryschool.csv", col.names = c("Time", "ID_i", "ID_j", "Class_i", "Class_j"))
metadata = read.table("metadata_primaryschool.txt", col.names = c("ID", "Class", "Gender"))

# ----------------------------------------------------------------------------------------------------------------------
# FORMAT THE TIME COLUMN (interval during which this contact was active is [ t – 20s, t ])
# ----------------------------------------------------------------------------------------------------------------------

# define the custom origin
custom_origin = as.POSIXct("2009-10-01 00:00:00", tz = "UTC")

# convert the Time to POSIXct directly (without using origin)
contacts$Time = as.POSIXct(contacts$Time, tz = "UTC")

# adjust Time by adding the difference from the custom origin
contacts$Time = custom_origin + as.numeric(difftime(contacts$Time, as.POSIXct("1970-01-01", tz = "UTC"), units = "secs"))

max(contacts$Time) - min(contacts$Time)

unique_days = as.numeric(unique(format(contacts$Time, format = "%d")))

start_times = end_times = rep(NA, length(unique_days))

for (day_ind in 1:length(unique_days)) {
  
  day = unique_days[day_ind]
  day_df = contacts[as.numeric(unique(format(contacts$Time, format = "%d"))) == day, ]
  start_times[day_ind] = min(format(day_df$Time, format = "%H:%M"))
  end_times[day_ind] = max(format(day_df$Time, format = "%H:%M"))
  
} # end day_ind for loop
rm(day_ind, day, day_df)

start_times
end_times

rm(start_times, end_times)
# ----------------------------------------------------------------------------------------------------------------------
# WE SHOULD ONLY USE HOURS WHICH ARE FULL, I.E. FROM 9AM TO 5PM EVERY DAY
# I.E. WE HAVE DATA FOR 8 HOURS FOR EACH OF 2 DAYS -> 16 NETWORKS
# ----------------------------------------------------------------------------------------------------------------------
from_9am = as.numeric(format(contacts$Time, format = "%H")) >= 09
until_5pm = as.numeric(format(contacts$Time, format = "%H")) < 17

contacts = contacts[(from_9am & until_5pm), ]

rm(from_9am, until_5pm)

# ----------------------------------------------------------------------------------------------------------------------
# ALTERNATIVE TIME COLUMN
# ----------------------------------------------------------------------------------------------------------------------
contacts$Time = paste0("Oct", as.numeric(format(contacts$Time, format = "%d")), "_H.", as.numeric(format(contacts$Time, format = "%H")))

# ----------------------------------------------------------------------------------------------------------------------
# BUILD BINARY AND WEIGHTED ADJACENCY MATRICES WHERE A NETWORK COVERS 1 HOUR -> 16 NETWORKS
# ----------------------------------------------------------------------------------------------------------------------
M = 16
M_names = unique(contacts$Time)

N = length(unique(metadata$ID))
N_names = unique(metadata$ID)

# 

primary_school_adj_matrices_binary = array(0, c(N, N, M), dimnames = list(N_names, N_names, M_names))
primary_school_adj_matrices_count = array(0, c(N, N, M), dimnames = list(N_names, N_names, M_names))

for (m in 1:M) {
  
  network_m = M_names[m]
  
  net_df = contacts[contacts$Time == network_m, ]
  
  for (i in 1:N) {
    
    node_i = N_names[i]
    
    net_node_df = net_df[net_df$ID_i == node_i, ]
    
    for (j in 1:N) {
      
      if (i < j) { # fill upper triangle
        
        node_j = N_names[j]
        
        net_node_node_df = net_node_df[net_node_df$ID_j == node_j, ]
        
        primary_school_adj_matrices_binary[i, j, m] = ifelse(nrow(net_node_node_df > 0), 1, 0)
        primary_school_adj_matrices_count[i, j, m] = nrow(net_node_node_df > 0)
        
        rm(node_j, net_node_node_df)
        
      } # end if statement
      
    } # end j for loop
    rm(j)
    
  } # end i for loop
  rm(i, node_i, net_node_df)
  
  primary_school_adj_matrices_binary[, , m] = primary_school_adj_matrices_binary[, , m] + t(primary_school_adj_matrices_binary[, , m])
  primary_school_adj_matrices_count[, , m] = primary_school_adj_matrices_count[, , m] + t(primary_school_adj_matrices_count[, , m])
  
  cat("Network", m, "done!\n")
} # end m for loop
rm(m, network_m, net_df)

# save the adjacency arrays
saveRDS(primary_school_adj_matrices_binary, "primary_school_adj_matrices_binary.rds")
saveRDS(primary_school_adj_matrices_count, "primary_school_adj_matrices_count.rds")

# ----------------------------------------------------------------------------------------------------------------------
# ARE THERE ANY COMPLETELY UNCONNECTED N ODES THROUGHOUT THE ENTIRE MULTIPLEX?
# ----------------------------------------------------------------------------------------------------------------------
hist(apply(primary_school_adj_matrices_binary, 3, sum)) # number of connections in each network

num_interactions_binary = apply(primary_school_adj_matrices_binary, 1, sum) # degree of each node across the entire multiplex
hist(num_interactions_binary)

# KEEP ONLY NODES THAT HAVE AT LEAST 20 CONNECTIONS THROUGHOUT THE ENTIRE MULTIPLEX
nodes_remove = as.vector(which(num_interactions_binary < 20))

# ----------------------------------------------------------------------------------------------------------------------
# ARE THERE ANY EMPTY NETWORKS?
# ----------------------------------------------------------------------------------------------------------------------
any(apply(primary_school_adj_matrices_binary, 3, sum) == 0)
nets_remove = as.vector(which(apply(primary_school_adj_matrices_binary, 3, sum) < 10)) 

save.image("Priamry_School_Workspace.Rdata")

# ================================================================================================
# PLOTS
# ================================================================================================
set.seed(123)

rm(list = ls())

setwd("/Volumes/My Passport/PhD/PhD Project 1/LaPCoM_Code/Applications/Primary_School_Count/")

primary_school_adj_matrices_count = readRDS("primary_school_adj_matrices_count.rds")
M = dim(primary_school_adj_matrices_count)[3]
M_names = dimnames(primary_school_adj_matrices_count)[[3]]

cols = khroma::color("muted")(9)

count_palette = grDevices::colorRampPalette(c("white", "black"))(171)

# ------------------------------------------------------------------------------------------------
# EXTRA PLOT
# ------------------------------------------------------------------------------------------------
ind = which.max(apply(primary_school_adj_matrices_count > 0, 3, sum))
ind = 9
library(statnet); coords = gplot(network::as.network(primary_school_adj_matrices_count[, , ind], directed = T))

pdf("Data_PDF_Plots/PSC_Data_Unclustered.pdf", width = 12, height = 8)
par(mfrow = c(3, 6))
par(mar = c(1.5, 1, 1.5, 1)) # bltr
# par(mar = c(2.5, 2, 3.5, 2)) # bltr
for (m in 1:M) { 
  plot(network::as.network(primary_school_adj_matrices_count[, , m], directed = F), 
       main = M_names[m], cex.main = 1.5, edge.col = "grey", vertex.col = 1, vertex.cex = 1.5, 
       coord = coords, edge.lwd = primary_school_adj_matrices_count[, , m] / 10)
} # end m for loop
dev.off()

# 

# ------------------------------------------------------------------------------------------------
# VIEW ARBITRARILY AS HEATMAPS
# ------------------------------------------------------------------------------------------------
pdf("Data_PDF_Plots/Primary_School_Undirected_Binary_Unclustered_Heatmaps.pdf", width = 10, height = 5)
par(mfrow = c(2, 8))
par(mar = c(2.5, 2, 3.5, 2)) # bltr
for (m in 1:M) { 
  image(primary_school_adj_matrices_binary[, , m], main = M_names[m], cex.main = 1.8, xaxt = "n", yaxt = "n") } 
dev.off()

pdf("Data_PDF_Plots/Primary_School_Undirected_Count_Unclustered_Heatmaps.pdf", width = 10, height = 5)
par(mfrow = c(2, 8))
par(mar = c(2.5, 2, 3.5, 2)) # bltr
for (m in 1:M) { 
  image(primary_school_adj_matrices_count[, , m], main = M_names[m], cex.main = 1.8, xaxt = "n", yaxt = "n", col = count_palette) } 
dev.off()

# ------------------------------------------------------------------------------------------------
# VIEW AS NETWORKS WITH ARBITRARY COORDIANTES
# ------------------------------------------------------------------------------------------------
pdf("Data_PDF_Plots/Primary_School_Undirected_Binary_Unclustered_Arbitrary_Coords.pdf", width = 10, height = 5)
par(mfrow = c(2, 8))
par(mar = c(1.5, 1, 1.5, 1)) # bltr
for (m in 1:M) {
  plot(network::as.network(primary_school_adj_matrices_binary[, , m], directed = T), main = M_names[m], cex.main = 1.8,
       edge.col = "grey", vertex.col = cols[1], vertex.cex = 3)
} # end m for loop (plotting)
dev.off()

pdf("Data_PDF_Plots/Primary_School_Undirected_Count_Unclustered_Arbitrary_Coords.pdf", width = 10, height = 5)
par(mfrow = c(2, 8))
par(mar = c(1.5, 1, 1.5, 1)) # bltr
for (m in 1:M) {
  plot(network::as.network(primary_school_adj_matrices_count[, , m], directed = T), main = M_names[m], cex.main = 1.8,
       edge.col = "grey", vertex.col = cols[1], vertex.cex = 3, edge.lwd = primary_school_adj_matrices_count[, , m])
} # end m for loop (plotting)
dev.off()

# ------------------------------------------------------------------------------------------------
# VIEW AS NETWORKS WITH THE SAME COORDINATES
# ------------------------------------------------------------------------------------------------
ind = which.max(apply(primary_school_adj_matrices_binary, 3, sum))
library(statnet); coords = gplot(network::as.network(primary_school_adj_matrices_binary[, , ind], directed = T))

pdf("Data_PDF_Plots/Primary_School_Undirected_Binary_Unclustered_Same_Coords.pdf", width = 10, height = 5)
par(mfrow = c(2, 8))
par(mar = c(1.5, 1, 2.5, 1)) # bltr
for (m in 1:M) {
  plot(network::as.network(primary_school_adj_matrices_binary[, , m], directed = T), main = M_names[m], cex.main = 1.8,
       edge.col = "grey", vertex.col = cols[1], vertex.cex = 3, coord = coords)
} # end m for loop (plotting)
dev.off()

pdf("Data_PDF_Plots/Primary_School_Undirected_Count_Unclustered_Same_Coords.pdf", width = 10, height = 5)
par(mfrow = c(2, 8))
par(mar = c(1.5, 1, 2.5, 1)) # bltr
for (m in 1:M) {
  plot(network::as.network(primary_school_adj_matrices_count[, , m], directed = T), main = M_names[m], cex.main = 1.8,
       edge.col = "grey", vertex.col = cols[1], vertex.cex = 3, coord = coords, edge.lwd = primary_school_adj_matrices_count[, , m])
} # end m for loop (plotting)
dev.off()

# ------------------------------------------------------------------------------------------------
# COMBINED PLOT
# ------------------------------------------------------------------------------------------------
pdf("Data_PDF_Plots/Primary_School_Undirected_Binary_Unclustered_Combined_Arbitrary_Coords.pdf", width = 10, height = 6)
par(mfrow = c(4, 8))
par(mar = c(2, 1, 2, 1) + 0.1)
for (m in 1:M) { 
  image(primary_school_adj_matrices_binary[, , m], xaxt = "n", yaxt = "n", main = M_names[m], cex.main = 1.5)
  plot(network::as.network(primary_school_adj_matrices_binary[, , m], directed = F), main = M_names[m], cex.main = 1.5,
       edge.col = "grey", vertex.col = cols[1], vertex.cex = 2)
} # end m for loop
dev.off()

pdf("Data_PDF_Plots/Primary_School_Undirected_Count_Unclustered_Combined_Arbitrary_Coords.pdf", width = 10, height = 6)
par(mfrow = c(4, 8))
par(mar = c(2, 1, 2, 1) + 0.1)
for (m in 1:M) { 
  image(primary_school_adj_matrices_count[, , m], xaxt = "n", yaxt = "n", main = M_names[m], cex.main = 1.5, col = count_palette)
  plot(network::as.network(primary_school_adj_matrices_binary[, , m], directed = F), main = M_names[m], cex.main = 1.5,
       edge.col = "grey", vertex.col = cols[1], vertex.cex = 2, edge.lwd = primary_school_adj_matrices_count[, , m])
} # end m for loop
dev.off()

# ------------------------------------------------------------------------------------------------
# COMBINED PLOT - SAME COORDS
# ------------------------------------------------------------------------------------------------
pdf("Data_PDF_Plots/Primary_School_Undirected_Binary_Unclustered_Combined_Same_Coords.pdf", width = 10, height = 6)
par(mfrow = c(4, 8))
par(mar = c(2, 1, 2, 1) + 0.1)
for (m in 1:M) { 
  image(primary_school_adj_matrices_binary[, , m], xaxt = "n", yaxt = "n", main = M_names[m], cex.main = 1.5)
  plot(network::as.network(primary_school_adj_matrices_binary[, , m], directed = F), main = M_names[m], cex.main = 1.5,
       edge.col = "grey", vertex.col = cols[1], vertex.cex = 2, coord = coords)
} # end m for loop
dev.off()

pdf("Data_PDF_Plots/Primary_School_Undirected_Count_Unclustered_Combined_Same_Coords.pdf", width = 10, height = 6)
par(mfrow = c(4, 8))
par(mar = c(2, 1, 2, 1) + 0.1)
for (m in 1:M) { 
  image(primary_school_adj_matrices_count[, , m], xaxt = "n", yaxt = "n", main = M_names[m], cex.main = 1.5, col = count_palette)
  plot(network::as.network(primary_school_adj_matrices_count[, , m], directed = F), main = M_names[m], cex.main = 1.5,
       edge.col = "grey", vertex.col = cols[1], vertex.cex = 2, coord = coords, edge.lwd = primary_school_adj_matrices_count[, , m])
} # end m for loop
dev.off()
