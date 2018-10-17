# https://satijalab.org/loomR/loomR_tutorial.html
# Install devtools from CRAN
install.packages("devtools")
# Use devtools to install hdf5r and loomR from GitHub
devtools::install_github(repo = "hhoeflin/hdf5r")
devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")
# Load loomR
library(loomR)
# Download an example file from http://loom.linnarssonlab.org/
download.file(url = "http://loom.linnarssonlab.org/clone/osmFISH/osmFISH_SScortex_mouse_all_cells.loom", destfile = "osmFISH_SScortex_mouse_all_cells.loom")
# Connect to the loom file in read/write mode
lfile <- connect(filename = "osmFISH_SScortex_mouse_all_cells.loom", mode = "r+")
lfile
lfile[["matrix"]][1:5, 1:5]



