### this code compares the utility and safety of synthetic databases 
## created with 2, 3, and 4 parent faces 

library(targets)
source("./R/functions.R")

# Packages needed
tar_option_set(packages = c("geomorph", "tidyverse", "webmorphR", "ggplot2"))

path_data = "./Data/"

paste(path_data, "synt_tem2/", sep = "")

# List of target objects.
list(
  
  #1 reading original data
  tar_target(
    dataO, imag_tem(i_path = path_data)
  ),
  
  #2 creating synthetic faces using 2, 3, and 4 parent faces
  tar_target(
    synt2, synt_face(nparent = 2, 
                     facesO = dataO[[1]], 
                     o_path = "./Data/synt_tem2/", 
                     nm = "S2_")
  ),
  tar_target(
    synt3, synt_face(nparent = 3, 
                     facesO = dataO[[1]], 
                     o_path = "./Data/synt_tem3/", 
                     nm = "S3_")
  ),
  tar_target(
    synt4, synt_face(nparent = 4, 
                     facesO = dataO[[1]], 
                     o_path = "./Data/synt_tem4/", 
                     nm = "S4_")
  ),
  
  #3 assessing utility 
  tar_target(
    util2, utility(facesO = dataO[[1]],
                   facesidO = as.character(unlist(dataO[[2]])),
                   facesS = synt2,
                   facesidS = attributes(synt2)$names,
                   o_path = "./Data/synt_tem2/",
                   lm_n = 70, seed = 1, nreplica = 30, nparent = 2, nm = "S2_")
  ),
  tar_target(
    util3, utility(facesO = dataO[[1]],
                   facesidO = as.character(unlist(dataO[[2]])),
                   facesS = synt3,
                   facesidS = attributes(synt2)$names,
                   o_path = "./Data/synt_tem3/",
                   lm_n = 70, seed = 1, nreplica = 30, nparent = 3, nm = "S3_")
  ),
  tar_target(
    util4, utility(facesO = dataO[[1]],
                   facesidO = as.character(unlist(dataO[[2]])),
                   facesS = synt4,
                   facesidS = attributes(synt2)$names,
                   o_path = "./Data/synt_tem4/",
                   lm_n = 70, seed = 1, nreplica = 30, nparent = 4, nm = "S4_")
  ),
  
  #4 assesing safety 
  tar_target(
    safe2, safety(facesO = dataO[[1]],
                  facesidO = as.character(unlist(dataO[[2]])),
                  facesS = synt2,
                  facesidS = attributes(synt2)$names,
                  o_path = "./Data/synt_tem2/",
                  lm_n = 70, seed = 1, nreplica = 30, nparent = 2, nsampSize = 12, nm = "S2_")
    ),
  
  tar_target(
    safe3, safety(facesO = dataO[[1]],
                  facesidO = as.character(unlist(dataO[[2]])),
                  facesS = synt3,
                  facesidS = attributes(synt3)$names,
                  o_path = "./Data/synt_tem3/",
                  lm_n = 70, seed = 1, nreplica = 30, nparent = 3, nsampSize = 12, nm = "S3_")
  ),
  
  tar_target(
    safe4, safety(facesO = dataO[[1]],
                  facesidO = as.character(unlist(dataO[[2]])),
                  facesS = synt4,
                  facesidS = attributes(synt4)$names,
                  o_path = "./Data/synt_tem4/",
                  lm_n = 70, seed = 1, nreplica = 30, nparent = 4, nsampSize = 12, nm = "S4_")
  ),
  
  #5 list of the utility results for each parent type
  tar_target(
    utility_lists, list(util2[util2$Type == "S2",], util3[util3$Type == "S3",], util4[util4$Type == "S4",])
  ),
  
  #6 list of the safety results for each parent type
  tar_target(
    safety_lists, list(safe2[safe2$Type == "S2",], safe3[safe3$Type == "S3",], safe4[safe4$Type == "S4",])
  ),
  
  #7 plotting the results (see Figure 3 within the paper)
  tar_target(
    allplots, all_plots(Outil = util2[util2$Type == "Original",], 
                        list_Sutil = utility_lists, 
                        list_Ssafe = safety_lists, 
                        o_path = "./Data/")
  )
)

