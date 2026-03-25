# SIOSSAR
Spatial Transcriptomics Analysis for Pathway Activity and Hotspot Detection

## Introduction
**SIOSSAR** is an R package designed for analyzing spatial transcriptomics data to identify pathway activity and functional hotspots based on cell‑cell interactions.  
It integrates spatial network construction, ligand‑receptor expression, random walk with restart, and permutation testing to provide a comprehensive view of pathway activation in tissue architecture.
<div align=center>
<img width="430" height="452" alt="e5365d3f8702092611734aa5e1ae1aac" src="https://github.com/user-attachments/assets/dab46e84-3645-4e40-b8a6-7e18c4aeba0c" />
</div>

## Installation

```
#install from github
remotes::install_github("https://github.com/qing-HMU/SIOSSAR")
#devtools::install_github("qing-HMU/SIOSSAR")

### Development version (recommended)
install.packages('SIOSSAR')

```

## Usage
```shell
#RUN examples
library(SIOSSAR)

# Load example data
file_path <- system.file("extdata", "siossar_example.rds", package = "SIOSSAR")
data <- readRDS(file_path)

# Define pathways to analyse
pathway <- c("antigen_processing_and_presentation")

# Run SIOSSAR with default parameters
# This may take several minutes depending on dataset size
data <- compute_siossar(
  object = data,
  pathway = pathway,
  verbose = F
)

# Plot activity scores for the pathway
p1 <- score_plot(data, "antigen_processing_and_presentation")
p1

# Plot activation regions
p2 <- hot_region_plot(data, "antigen_processing_and_presentation")
p2

# Plot hierarchical layers (show layers 1-10)
p3 <- hplot(data, "antigen_processing_and_presentation", show_layers = 1:10)
p3

# Create Sankey diagram for pathway
p4 <- plot_PPI_sankey(data, pathway[1])
p4


# Plot heatmap for top pathways
p5 <- plot_top_pathway_heatmap(data, pathway, top_n = 1)
p5
```





