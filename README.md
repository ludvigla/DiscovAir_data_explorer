# DiscovAir spatial transcriptomics data explorer
With this app you can explore spatial transcriptomics data of sections from lung tissue.

# Installation
Clone the repo from the terminal by running:

`$git clone https://github.com/ludvigla/DiscovAir_lung_explorer`

To run the app, you first need to install the following R packages:

CRAN:
- ggplot2
- zeallot
- shiny
- shinydashboard
- RColorBrewer
- viridis
- scales
- shinyWidgets
- leaflet
- leaflet.extras
- htmlwidgets
- shinyBS

GitHub:
- d3heatmap

You can install the packages directly with the install-packages.R script:

`$Rscript install-packages.R`

# Run the app
From RStudio, navigate to the cloned repository:

`setwd("~/DiscovAir_data_explorer")`

You can then activate the app by running:

`library(shiny)`

`runApp()`

Or alternatively you can open the app.R file File->Open file->.../app.R and then click on the Run App icon at the top of the script.

# How to use
When you open the app, you will see an H&E image of the first lung tissue section (`section 1`) overlaid by spots. The color of the spots correspond either to factor activity values, normalized gene expression or cell type proportions (will be added soon) which you can select from in the left panel. You can also switch tabs to visualize the other lung tissue section. 

Formatting options:
  * opacity : sets the (maximum) spot opacity to control transparanecy
  * maxcutoff : trim values based on a specifiec quantile, eg. a value of 0.99 will clip the values of the displayed value vector to the 99th quantile
  * reverse palette : should the color palette be reversed?
  * opacity gradient : should the spot opacity be scaled with the displayed value?
  * Factor gene weights : show a bar plot of the top contributing genes for a selected factor
  * Value histogram : show a histogram of the selected value vector
  * palette : select a color palette
