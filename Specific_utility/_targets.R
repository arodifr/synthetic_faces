## this could be a stand alone script
## It could be run after global utility and safety have been evaluated
### this code is for selecting the subsample that will best replace
## the original sample. 


library(targets)
source("./R/functions.R")


# Packages needed
tar_option_set(packages = c("geomorph", "tidyverse", "webmorphR", "ggplot2"))

path_data = "./Data/"


list(
  #1. reading original data
  tar_target(
    #target name, function
    dataO, imag_tem(i_path = path_data)
  ),
  ##2 creating synthetic faces using 2, 3, and 4 parent faces
  tar_target(
    synt3, synt_face(nparent = 3, 
                     facesO = dataO[[1]], 
                     o_path = "./Data/synt_tem3/", 
                     nm = "S3_")
  ),
  #3 useful samples
  tar_target(
    useful, util_subsamp(facesO = dataO[[1]],
                         facesidO = as.character(unlist(dataO[[2]])),
                         facesS = synt3,
                         facesidS =  attributes(synt3)$names,
                         o_path = "./Data/synt_tem3/",
                         lm_n = 70, seed = 1, nreplica = 30, nparent = 3)
  )
)

#to read the output information
#tar_read(useful) 
#formato(tar_read(useful)[[1]][[1]])
#formato(tar_read(useful)[[3]][[1]])
#formato(tar_read(useful)[[2]][[1]])
