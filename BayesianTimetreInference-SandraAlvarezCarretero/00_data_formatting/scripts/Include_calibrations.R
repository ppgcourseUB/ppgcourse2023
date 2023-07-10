#-------------------#
# CLEAN ENVIRONMENT #
#-------------------#
rm( list = ls( ) )

#-----------------------#
# SET WORKING DIRECTORY #
#-----------------------#
# Load package to help find the path to this source file 
library(rstudioapi) 
# Get the path of current open file
path_to_file <- getActiveDocumentContext()$path 
# Get working directory path
wd      <- paste( dirname( path_to_file ), "/", sep = "" )
wd.name <- dirname( path_to_file )
# Set wd. Note that the wd will be the same directory where this 
# script is saved.
setwd( wd )

#----------------------------#
# DEFINE GLOBAL VARS BY USER #
#----------------------------#
# Name of output calibrated tree file ready to be used by `MCMCtree`.
# Note that the file name will have the following format
# "<your_selected_out_name>_calib_MCMCtree.tree".
#
# NOTE: Make always sure that there is at least one blank line at the 
# end of the this text file! Otherwise, you will get an error telling you that 
# there is an incomplete final line in these files. This file needs to be
# already in PHYLIP format. Please follow the same format as used in the 
# example tree file provided.
out_name1 <- c( "tree_fosscal" )
out_name2 <- c( "tree_seccal" )

# Path to your input tree with calibrations. If in the same directory,
# you only need to write the name. If in a different directory, please
# type the absolute or relative path to this file. You need to include the
# flags within square brackets (e.g., [Mammalia]) and write them on the node
# that is to be calibrated.
#
# NOTE: Make always sure that there is at least one blank line at the 
# end of the this text file! Otherwise, you will get an error telling you that 
# there is an incomplete final line in these files.
path_tree <- c( "../calibs/cytb_calibnames.tree" )

# Path to your input text file that allows you to match the flags you have 
# used to label the nodes that are to be calibrated with the calibration you 
# want to use in `MCMCtree` format. If in the same directory,
# you only need to write the name. If in a different directory, please
# type the absolute or relative path to this file.
# The format of this file meets the following criteria:
#
#   - No header.
#   - One row per calibration.
#   - No spaces at all.
#   - The pipe character, "|", is used to separate the flag you are using on
#     the node you want to calibrate (left) from the `MCMCtree` calibration
#     that you will use to calibrate the node (right).
#   - The `MCMCtree` calibrations are written in `MCMCtree` format and with 
#     single quotation marks. No spaces.
#   - No spaces when writing the name of your flag either.
# 
# E.g.: row in this text file to calibrate node "Mammalia". The flag used in the
# tree to locate the node "Mammalia" is "Mammalia" and the `MCMCtree`
# calibration is a soft-bound calibration:
#
# ```
# Mammalia|'B(1.649,2.51254)'
# ```
#`
# If in doubt, please follow the same format as used in the example text 
# file provided.
path_textconv1 <- c( "../calibs/Calib_converter_1.txt" )
path_textconv2 <- c( "../calibs/Calib_converter_2.txt" )

#---------------------------------#
# READ TREE AND CALIBRATIONS FILE #
#---------------------------------#
# Read tree and get phylip header
# NOTE: Make always sure that there is at least one blank line at the 
# end of the tree file! Otherwise, you will get an error telling you that 
# there is an incomplete final line in these files.
tt_name       <- path_tree
tt            <- readLines( tt_name )
phylip.header <- tt[1]
tt1           <- tt2 <- tt[2]

# Read file converter that will be used to rename calibrations
calibrations1 <- read.table( file = path_textconv1,
                             stringsAsFactors = F, sep = "|",
                             blank.lines.skip = T )
calibrations2 <- read.table( file = path_textconv2,
                             stringsAsFactors = F, sep = "|",
                             blank.lines.skip = T )
colnames( calibrations1 ) <- colnames( calibrations1 ) <- 
  c( "name", "MCMCtree calib" )

#--------------------------------#
# REPLACE TAGS WITH CALIBRATIONS #
#--------------------------------#
# Replace calibration names with corresponding calibration
for( j in 1:dim( calibrations1 )[1] ){
  tt1 <- gsub( pattern = paste0("\\[",calibrations1[j,1],"\\]"),
               x = tt1,
               replacement = paste( "'", calibrations1[j,2], "'", sep = "" ) )
  tt2 <- gsub( pattern = paste0("\\[",calibrations2[j,1],"\\]"),
               x = tt2,
               replacement = paste( "'", calibrations2[j,2], "'", sep = "" ) )
}


#-------------------------------#
# WRITE CALIBRATED TREE IN FILE #
#-------------------------------#
out_dir <- "../01_inp_data/"
write( x = phylip.header, file = paste( out_dir, out_name1, ".tree", sep = "" ) )
write( x = tt1, file = paste( out_dir, out_name1, ".tree", sep = "" ),
       append = TRUE )
write( x = phylip.header, file = paste( out_dir, out_name2, ".tree", sep = "" ) )
write( x = tt2, file = paste( out_dir, out_name2, ".tree", sep = "" ),
       append = TRUE )
  




