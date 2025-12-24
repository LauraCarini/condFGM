library(igraph)
library(ggraph)
library(extrafont)
library(patchwork)
library (pracma)
folder = "EEG_dataset/"
folder_res = "EEG_dataset/seed_1/results/"


plot_single_graph <- function(adjm, title){
  
  graph <- graph.adjacency(adjm, weighted = TRUE, mode = "undirected")
  
  node_names <- read.csv(
    paste0(folder, "position_list.txt"),
    header = FALSE,
    stringsAsFactors = FALSE
  )
  V(graph)$name <- unlist(node_names, use.names = FALSE)
  
  gp <- ggraph(graph, layout = "linear", circular = TRUE) +
    
    geom_edge_arc(
      aes(edge_colour = ifelse(weight == 0, NA, weight)),
      show.legend = TRUE
    ) +
    
    geom_node_point(shape = 21, size = 3, fill = "darkgrey") +
    
    theme_graph(base_family = "Calibri") +
    
    scale_edge_colour_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      na.value = "transparent",
      oob = scales::squish_infinite,
      breaks = scales::pretty_breaks(n = 5),
      name = "Edge weight"
    ) +
    
    guides(
      fill = guide_colorbar(
        title.position = "top",
        barheight = unit(4, "cm"),
        barwidth  = unit(0.4, "cm")
      )
    ) +
    
    coord_fixed(xlim = c(-1.4, 1.4), ylim = c(-1.4, 1.4)) +
    
    geom_node_text(
      aes(
        label = name,
        x = x * 1.2,
        y = y * 1.2,
        hjust = ifelse(
          atan2(y, x) > -pi/2 & atan2(y, x) < pi/2, 0, 1
        )
      ),
      size = 2
    ) +
    
    ggtitle(title)
  
  return(gp)
}



load(paste0(folder, "Seed_1/results/Test1_Adj_estimation.rda"))
adj_group = G.our.symm.weighted$group
adj_group = log(adj_group)
adj_group[is.na(adj_group)] = 0

# remove diagonal if needed
diag(adj_group) <- 0

# vector of existing edges
vals <- adj_group[adj_group != 0 & !is.na(adj_group)]
lower_thr <- quantile(vals, probs = 0.05, na.rm = TRUE)
upper_thr <- quantile(vals, probs = 0.95, na.rm = TRUE)
adj_extreme <- adj_group

adj_extreme[
  adj_group > lower_thr & adj_group < upper_thr
] <- 0


plot = plot_single_graph(adj_extreme, "Differential network")
ggsave(plot, filename = paste0(folder_res, "Differential_network.png"), width = 8, height = 8, dpi = 300)


