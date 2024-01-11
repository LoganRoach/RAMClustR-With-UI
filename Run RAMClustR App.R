#### Load Packages ####
# loads the packages used for the script to load packages
packages_to_install <- c("shiny")
for (package_name in packages_to_install) {
  if(!requireNamespace(package_name, quietly = TRUE)) {
    install.packages(package_name)
  }
  library(package_name, character.only = TRUE)
}
# run App
runApp("RamClustR App.R")
