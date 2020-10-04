library(roxygen2)
# where shall we assemble the package?
PKG_ROOT <- "C:/My/tmp"
# where are the source files?
PKG_SOURCE <- "C:/My/code/robsurvey"

#-------------------------------------------------------------------------------
# package skeleton
setwd(paste0(PKG_SOURCE, "/R"))
utils::package.skeleton(name = "robsurvey", 
   path = PKG_ROOT, 
   force = TRUE, # i.e. overwrite existing directory
   code_files = list.files(paste0(PKG_SOURCE, "/R"), pattern = "\\.R$"))

# delete readme
file.remove(paste0(PKG_ROOT, "/robsurvey/Read-and-delete-me"))

# delete all *.Rd-files in /man
rd_files <- list.files(paste0(PKG_ROOT, "/robsurvey/man"), full.names = TRUE)
file.remove(rd_files)

# copy existing *.Rd from the source to the package
rd_files <- list.files(paste0(PKG_SOURCE, "/man"), pattern = "\\.Rd$", 
   full.names = TRUE)
file.copy(rd_files, paste0(PKG_ROOT, "/robsurvey/man"), overwrite = TRUE)

#-------------------------------------------------------------------------------
# generate file 'DESCRIPTION'
description <- 'Package: robsurvey
Type: Package
Title: Robust Survey Statistics Estimation
Version: 0.2
Authors@R: 
   c(person(
      given = "Beat", 
      family = "Hulliger", 
      email = "beat.hulliger@fhnw.ch", 
      role = c("aut", "cph")),
    person(
      given ="Tobias", 
      family = "Schoch", 
      email = "tobias.schoch@fhnw.ch", 
      role = c("aut","cre", "cph")),
    person(
      given = "Martin", 
      family = "Sterchi", 
      email = "martin.sterchi@fhnw.ch", 
      role = c("aut")))
Author: Beat Hulliger [aut, cph],
  Tobias Schoch [aut, cre, cph],
  Martin Sterchi [aut]
Maintainer: Tobias Schoch <tobias.schoch@fhnw.ch>
Description: Functions to compute robust (outlier-resistant) estimates of finite 
    population characteristics. The package supports the computations of robust 
    means, totals, ratios, etc. Available methods are Huber M-estimators, 
    trimming, and winsorization. The package robsurvey complements the survey.
License: MIT + file LICENSE
URL: https://github.com/tobiasschoch/robsurvey
BugReports: https://github.com/tobiasschoch/robsurvey/issues
Encoding: UTF-8
LazyData: true
Depends: 
    R (>= 3.5.0)
Imports:
    grDevices,
    stats,
    survey (>= 3.35-1),
    KernSmooth
'

cat(description, file = paste0(PKG_ROOT, "/robsurvey/DESCRIPTION"))

#-------------------------------------------------------------------------------
# generate file 'LICENSE' 
lic <- paste0("YEAR: ", format(Sys.time(), "%Y"), 
   "\nCOPYRIGHT HOLDER: Beat Hulliger and Tobias Schoch")
cat(lic, file = paste0(PKG_ROOT, "/robsurvey/LICENSE"))

#-------------------------------------------------------------------------------
# copy /src
src_folder <- paste0(PKG_ROOT, "/robsurvey/src")
if (!dir.exists(src_folder)) dir.create(src_folder)
src_files_c <- list.files(paste0(PKG_SOURCE, "/src"), pattern = "\\.c$", 
   full.names = TRUE)
src_files_h <- list.files(paste0(PKG_SOURCE, "/src"), pattern = "\\.h$", 
   full.names = TRUE)
# remove C files that are not needed
# src_files_c <- src_files_c[-grep("sctutils.c", src_files_c)]
file.copy(c(src_files_c, src_files_h), src_folder, overwrite = TRUE)
# copy file 'Makevars'
file.copy(paste0(PKG_SOURCE, "/src/Makevars"), src_folder, overwrite = TRUE)

#-------------------------------------------------------------------------------
# copy /data
dir.create(paste0(PKG_ROOT, "/robsurvey/data"))
data_files <- list.files(paste0(PKG_SOURCE, "/data"), pattern = "\\.RData$", 
   full.names = TRUE)
file.copy(data_files, paste0(PKG_ROOT, "/robsurvey/data"), overwrite = TRUE)

#-------------------------------------------------------------------------------
# roxygenize
setwd(PKG_ROOT)
roxygenize(package.dir = "robsurvey")

# delete file 'NAMESPACE' (we let roxygen2 do it for us)
file.remove(paste0(PKG_ROOT, "/robsurvey/NAMESPACE"))
roxygenize(package.dir = "robsurvey", namespace_roclet())

# clean up after roxygenize...
object_files <- list.files(paste0(PKG_ROOT, "/robsurvey/src"), pattern = "\\.o$", 
   full.names = TRUE)
file.remove(object_files)

dll_files <- list.files(paste0(PKG_ROOT, "/robsurvey/src"), pattern = "\\.dll$", 
   full.names = TRUE)
file.remove(dll_files)

#-------------------------------------------------------------------------------
# generate the manual 'by hand' (for inspection) 
# system("R CMD Rd2pdf robsurvey")

#-------------------------------------------------------------------------------
# check 
# system("R CMD check robsurvey")

#-------------------------------------------------------------------------------
# install and build
system("R CMD INSTALL --build robsurvey")

#-------------------------------------------------------------------------------
# only build 
# system("R CMD build robsurvey")


