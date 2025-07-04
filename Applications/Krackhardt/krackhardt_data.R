####################################################################################################
# KRACKHARDT MANAGEMENT DATA
####################################################################################################

rm(list = ls())
setwd("/Volumes/My Passport/PhD/PhD Project 1/LaPCoM_Code/Applications/Krackhardt/")

set.seed(1230)

# ================================================================================================
# LOAD AND FORMAT THE DATA
# ================================================================================================
krack_data = XML::xmlToList(XML::xmlParse("krackhardt.xml"))

krack_data = krack_data$MetaNetwork$networks

N = M = 21 
krack_multi = array(NA, c(N, N, M))

for (m in 1:M) {
  network_m_el = matrix(NA, nrow = length(krack_data[[m]]) - 2, ncol = 2)
  for (i in 1:(length(krack_data[[m]]) - 2)) { network_m_el[i, ] = as.numeric(krack_data[[m]][[i]][c(1, 2)]) }
  graph_m = igraph::graph_from_edgelist(network_m_el)
  krack_multi[, , m] = igraph::as_adjacency_matrix(graph_m, sparse = F)
} # end m for loop
rm(m, network_m_el, graph_m)

all(apply(krack_multi, 3, isSymmetric)) # not symmetric

save(krack_multi, file = "Krackhardt_Multiplex_Directed.Rdata")

# ------------------------------------------------------------------------------------------------
# IGNORE DIRECTIONALITY BECAUSE LaPCoM ONLY WORKS WITH UNDIRECTED NETWORKS
# ------------------------------------------------------------------------------------------------
# for (m in 1:M) {
#   graph_m = igraph::as.undirected(igraph::graph_from_adjacency_matrix(krack_multi[, , m]), mode = "collapse")
#   krack_multi[, , m] = igraph::as_adjacency_matrix(graph_m, sparse = F)
# } # end m for loop 
# rm(m, graph_m)
# 
# save(krack_multi, file = "Krackhardt_Multiplex_Undirected.Rdata")

# ================================================================================================
# PLOTS
# ================================================================================================
rm(list = ls())
setwd("/Volumes/My Passport/PhD/PhD Project 1/LaPCoM_Code/Applications/Krackhardt/")

set.seed(1230)

load("Krackhardt_Multiplex_Directed.Rdata")

# cols = khroma::color("muted")(9) # no colours because no clustering so "no need"
cols = 1:9
N = M = 21 

# ------------------------------------------------------------------------------------------------
# VIEW ARBITRARILY AS HEATMAPS
# ------------------------------------------------------------------------------------------------
pdf("Data_PDF_Plots/Krackhardt_Directed_Unclustered_Heatmaps.pdf", width = 10, height = 5)
par(mfrow = c(3, 7))
par(mar = c(2.5, 2, 3.5, 2)) # bltr
for (m in 1:M) { image(krack_multi[, , m], main = LETTERS[m], cex.main = 1.8,
                       xaxt = "n", yaxt = "n") } 
dev.off()

library(ggplot2)
library(dplyr)

n = dim(krack_multi)[1]
M = dim(krack_multi)[3]

heat_df = do.call(rbind, lapply(1:M, function(m) {
  mat = krack_multi[,,m]
  df = as.data.frame(as.table(mat))
  names(df) = c("Row", "Col", "Value")
  df$Network = LETTERS[m]
  df$Row = as.numeric(df$Row)
  df$Col = as.numeric(df$Col)
  df$Value = as.factor(df$Value)  # = treat as factor for discrete scale
  df
}))

# Plot with discrete black/white fill
p = ggplot(heat_df, aes(x = Col, y = Row, fill = Value)) +
  geom_tile(color = "black", linewidth = 0.2) +
  scale_fill_manual(values = c("0" = "white", "1" = "black")) +
  scale_y_continuous(trans = "reverse") +
  coord_fixed() +
  facet_wrap(~ Network, nrow = 3) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    strip.text = element_text(size = 18, face = "bold"),
    legend.position = "bottom",
    plot.margin = margin(2, 2, 2, 2)
  )

# Save to PDF
ggsave("Data_PDF_Plots/Krackhardt_Directed_Unclustered_Heatmaps.pdf", plot = p, width = 8, height = 5)

# ------------------------------------------------------------------------------------------------
# VIEW AS NETWORKS WITH ARBITRARY COORDIANTES
# ------------------------------------------------------------------------------------------------
pdf("Data_PDF_Plots/Krackhardt_Directed_Unclustered_Arbitrary_Coords.pdf", width = 10, height = 5)
par(mfrow = c(3, 7))
par(mar = c(1.5, 1, 1.5, 1)) # bltr
for (m in 1:M) {
  plot(network::as.network(krack_multi[, , m], directed = T), main = LETTERS[m], cex.main = 1.8,
       edge.col = "grey", vertex.col = cols[1], vertex.cex = 3)
} # end m for loop (plotting)
dev.off()

# ------------------------------------------------------------------------------------------------
# VIEW AS NETWORKS WITH THE SAME COORDINATES
# ------------------------------------------------------------------------------------------------
ind = which.max(apply(krack_multi, 3, sum))
library(statnet); coords = gplot(network::as.network(krack_multi[, , ind], directed = T))

pdf("Data_PDF_Plots/Krackhardt_Directed_Unclustered_Same_Coords.pdf", width = 10, height = 5)
par(mfrow = c(3, 7))
par(mar = c(1.5, 1, 2.5, 1)) # bltr
for (m in 1:M) {
  plot(network::as.network(krack_multi[, , m], directed = T), main = LETTERS[m], cex.main = 1.8,
       edge.col = "grey", vertex.col = cols[1], vertex.cex = 3, coord = coords)
} # end m for loop (plotting)
dev.off()