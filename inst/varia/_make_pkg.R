#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
    mode <- "help"
} else {
    mode <- switch(args,
        fast    = "fast",
        check   = "check",
        full    = "full",
        tar     = "tar",
        h       = "help",
        help    = "help",
        pdf     = "pdf")
}

if (mode == "help") {
    cat("\nOptions:\n")
    cat("  fast (install only x64, w/o html)\n")
    cat("  check\n")
    cat("  tar (generate tar ball)\n")
    cat("  full (standard installation)\n\n")
    q(save = "no")
}

#-------------------------------------------------------------------------------
PKG <- "robsurvey"
if (.Platform$OS.type == "unix") {
	PKG_SOURCE <- "/mnt/c/my/code"
	PKG_ROOT <- "/mnt/c/my/tmp"
} else {
	PKG_SOURCE <- "C:/My/code"
	PKG_ROOT <- "C:/My/tmp"
}
setwd(PKG_ROOT)

# delete package folder if it already exists
if (dir.exists(PKG))
    unlink(PKG, recursive = TRUE)

# copy entire directory (excl. files/folders with a leading dot, e.g. '.git')
dir.create(PKG)
pkg_files <- list.files(paste0(PKG_SOURCE, "/", PKG), full.names = TRUE)
file.copy(pkg_files, paste0(PKG_ROOT, "/", PKG), recursive = TRUE)

# copy .Rbuildignore and .Rinstignore
f_R_build_ignore <- paste0(PKG_SOURCE, "/", PKG, "/.Rbuildignore")
if (file.exists(f_R_build_ignore))
    file.copy(f_R_build_ignore, paste0(PKG_ROOT, "/", PKG))
f_R_inst_ignore <- paste0(PKG_SOURCE, "/", PKG, "/.Rinstignore")
if (file.exists(f_R_inst_ignore))
    file.copy(f_R_inst_ignore, paste0(PKG_ROOT, "/", PKG))

# clean src folder (remove binary files)
binary_files <- list.files(paste0(PKG_ROOT, "/", PKG, "/src"),
    pattern = "\\.o$|\\.dll$|\\.so$")
if (length(binary_files) > 0)
    file.remove(paste0(PKG_ROOT, "/", PKG, "/src/", binary_files))

#-------------------------------------------------------------------------------
setwd(PKG_ROOT)

# fast install (only x64 arch; without html help files and vignette)
if (mode == "fast")
    system(paste0("R CMD INSTALL ", PKG, " --no-html --no-multiarch"))

# build the tar ball
if (mode == "tar") {
    system(paste0("R CMD build ", PKG))
}

# build the tar ball and check
if (mode == "check") {
    vers <- trimws(strsplit(readLines(paste0(PKG,
        "/DESCRIPTION"))[4], ":")[[1]][2])
    pkg_tar <- paste0(PKG, "_", vers, ".tar.gz")
    if (file.exists(pkg_tar))
        file.remove(pkg_tar)
    system(paste0("R CMD build ", PKG))
    system(paste0("R CMD check --as-cran ", pkg_tar))
}

# full build and install (incl. html, vignette, etc)
if (mode == "full") {
    vers <- trimws(strsplit(readLines(paste0(PKG,
        "/DESCRIPTION"))[4], ":")[[1]][2])
    pkg_tar <- paste0(PKG, "_", vers, ".tar.gz")
    if (file.exists(pkg_tar))
        file.remove(pkg_tar)
    system(paste0("R CMD build ", PKG))
    system(paste0("R CMD INSTALL ", pkg_tar))
}

# render pdf manual
#system("R CMD Rd2pdf robsurvey")
