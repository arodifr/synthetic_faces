### FUNCTIONS FOR PAPER TAL

# 1. Data setup
imag_tem <- function(i_path){
  # This function reads your stimuli (images and templates)
  
  ## Requires: the directory (input path) where the stimuli is located ##
  
  faces <- read_stim(i_path, pattern = NULL)
  facesid <- data.frame(Folio = as.character(names(faces)))
  
  # Returns a list containing the stimuli and their ID
  return(list(faces, facesid))
}


synt_face <- function(nparent, facesO, o_path, nm, noriginal){
  # This function creates synthetic faces from parent images
  
  ## Requires:
  #          the number of parents used (nparent)
  #          the stim_list of images with their landmarks (facesO)
  #          the directory where the synthetic images will be saved (o_path)
  #          the name that will be included as name (i.e. "s2_" ) (nm)
  #                    use "_" in this name to separate nm from the iterator (S2_i, S2_i+1, S2_i+2...)
  
  
  # DF with al possible unique combinations between n parent faces
  synt_all <- as.data.frame(t(combn(sort(names(facesO)),nparent)))
  # Creates and saves the average faces depending on how many parent images are selected 
  synt_temp <- c()
  nm1 <- c()
  if(nparent == 2){
    for (i in 1:length(synt_all$V1)) {
      synt_temp[i] <- as_stimlist(avg(as_stimlist(c(
        facesO[synt_all$V1[i]],
        facesO[synt_all$V2[i]]))
      ))
      nm1[i] <- paste0(nm, i)
      write_stim(synt_temp[i], dir = o_path, names = nm1[i])
    }
    
  } else if(nparent == 3){
    for (i in 1:length(synt_all$V1)) {
      synt_temp[i] <- as_stimlist(avg(as_stimlist(c(
        facesO[synt_all$V1[i]],
        facesO[synt_all$V2[i]],
        facesO[synt_all$V3[i]]))
      ))
      nm1[i] <- paste0(nm, i)
      write_stim(synt_temp[i], dir = o_path, names = nm1[i])
    }
  } else if(nparent == 4){
    for (i in 1:length(synt_all$V1)) {
      synt_temp[i] <- as_stimlist(avg(as_stimlist(c(
        facesO[synt_all$V1[i]],
        facesO[synt_all$V2[i]],
        facesO[synt_all$V3[i]],
        facesO[synt_all$V4[i]]))
      ))
      nm1[i] <- paste0(nm, i)
      write_stim(synt_temp[i], dir = o_path, names = nm1[i])
    }
  } else{
    print("nparent must be 2, 3, or 4")
  }
  
  # Returns a list containing the synthetic faces, templates and their ID
  synth <- as_stimlist(synt_temp)
  attributes(synth)$names <- nm1
  
  
  return(synth)
}


formato <- function(data){
  format(data, digit = 2)
}

# euclidean <- function(a, b) sqrt(sum((a - b)^2))


util_subsamp <- function(facesO, facesidO, facesS, facesidS, o_path, lm_n,  seed, nreplica, nparent){
    # This function computes inter-individual distances (INDEPENDENTLY) for:
    #             the original database
    #             the synthetic database
    
    
    ## Requires:
    #          the stim_list of original images with their landmarks (facesO)
    #          the ID of the original images (facesidO)  
    #          the stim_list of synthetic images with their landmarks (facesS)
    #          the ID of the synthetic images (facesidS)
    #          the directory where the results will be saved (o_path)
    #          the number of landmarks in the templates (lm_n)
    #          the seed to replicate the results (seed)
    #          the number of random synthetic databases to be generated (nreplica)
    #          the number of parents used (nparent)

    
  #read synt images
  #  facesS <- read_stim(i_path, pattern = NULL)
  fs <- facesS
  
  # Get landmarks
  lmO <- get_point(facesO, 0:(lm_n-1))
  lmO <- as.matrix(lmO[,c(3,4)]) #col 1 y 2: image y point
  lmO <- arrayspecs(lmO, lm_n, 2)
  
  lmS <- get_point(fs, 0:(lm_n-1))
  lmS <- as.matrix(lmS[,c(3,4)]) #col 1 y 2: image y point
  lmS <- arrayspecs(lmS, lm_n, 2)
  
  # changing array format
  lmO <- two.d.array(lmO) 
  attributes(lmO)$dimnames[[1]] <- facesidO
  
  lmS <- two.d.array(lmS)
  attributes(lmS)$dimnames[[1]] <- facesidS
  
  ## THE MOST DIVERSE SYNTHETIC DATABASE RELATIVE TO THE ORIGINAL DIVERSITY ## 
  set.seed(seed)
  
  ## ORIGINAL DATA ###
  # Generalized Procrustes analysis for the original data and all the synthetic replicas
  original_gpa <- gpagen(arrayspecs(lmO, lm_n, 2), print.progress = F, verbose = T, ProcD = F)
  original_alig <- rotate.coords(original_gpa$coords, type = "rotateC")
  original_pca <- gm.prcomp(original_alig)
  original_dist <- unlist(as.vector(dist(two.d.array(original_alig))))
  
  minOO <- min(original_dist)
  maxOO <- max(original_dist)
  
  # Sampling the synthetic database "nreplica" times
  nsample = dim(lmO)[1]
  sample_list <- replicate(nreplica, sample(attributes(lmS)$dimnames[[1]], nsample), simplify=FALSE)

  ## SYNTHETIC DATA ###
  # Generalized Procrustes analysis for all the synthetic replicas
  sample_list <- lapply(sample_list, function(x) lmS[x,])
  sample_list <- lapply(sample_list, function(x) arrayspecs(x, lm_n, 2))
  sample_gpa <- lapply(sample_list, function(x) gpagen(x, print.progress = FALSE,
                                                       verbose = TRUE, ProcD = F))
  sample_aligned <- lapply(sample_gpa, function(x) rotate.coords(x$coords, type = "rotateC"))
  sample_pca <- lapply(sample_aligned, function(x) gm.prcomp(x))
  # Using PC scores to compute the Procrustes distances among synthetic images
  sample_distS <- lapply(sample_pca, function(x, i) {unlist(as.vector(dist(x$x)))})
  
  
  ## SIMILARITY IN TYPE OF SHAPE VARIATION (PC LOADINGS) ##
  PCs <- as.data.frame(matrix(NA, nrow = nreplica, ncol = 3))
  colnames(PCs) <- c("pc1", "pc2", "pc3") 
  for (i in 1:nreplica) {
    PCs$pc1[i] <- coef(lm(sample_pca[[i]]$rotation[,1]~original_pca$rotation[,1]))[[2]]
    PCs$pc2[i] <- coef(lm(sample_pca[[i]]$rotation[,2]~original_pca$rotation[,2]))[[2]]
    PCs$pc3[i] <- coef(lm(sample_pca[[i]]$rotation[,3]~original_pca$rotation[,3]))[[2]]
  }
  
  sim_loading <- PCs[order( -PCs[,1], -PCs[,2]),]
  # Subsample with similar type of shape diversity (i.e. PC1 and PC2 loadings) to the original type of diversity
  sim_loadings <- sim_loading[1,]
  
  ## SIMILARITY IN SHAPE DIVERSITY (INTER-INDIVIDUAL PROCRUSTES DISTANCES) ##
  # Descriptive statistics
  sample_descr <- lapply(sample_distS, function(x) c(min(x), max(x)))
  sample_df <- do.call(rbind.data.frame, sample_descr)
  colnames(sample_df) <- c("min", "max")
  sample_df <- sample_df[order( -sample_df[,1], -sample_df[,2]),]
  # Subsample with the highest shape diversity (inter-individual distance variation)
  sim_shape <- sample_df[1,] 
  
  
  ## SIMILARITY IN IN TYPE OF SHAPE VARIATION (PC LOADINGS)  & SHAPE DIVERSITY (INTER-INDIVIDUAL PROCRUSTES DISTANCES) ##
  
  # Subsample with similar type of shape diversity (i.e. PC1 and PC2 loadings) to the original type of diversity
  ## but also, with the highest shape diversity (inter-individual distance variation)
  load_shape <- cbind(PCs, sample_df)
  load_shapes <- load_shape[order( -load_shape$min, -load_shape$max, load_shape$pc1, load_shape$pc2), ]
  # selecting the first row with positive PC1 values, after considering the ordering
  sim_loadshape <- load_shapes[which(load_shapes$pc1>0), ] 
  sim_loadshape <- sim_loadshape[1,]
  
  # names of the synthetic images used to create a useful subsample
  sim_shape_n <- attributes(sample_list[[as.numeric(rownames(sample_df[1,]))]])$dimnames[[3]]
  sim_loading_n <- as.numeric(rownames(sim_loadings))
  sim_loading_n <- attributes(sample_list[[sim_loading_n]])$dimnames[[3]]
  sim_loadshape_n<- as.numeric(rownames(sim_loadshape))
  sim_loadshape_n <- attributes(sample_list[[sim_loadshape_n]])$dimnames[[3]]
  
  
  #saving the images of the useful samples
  #sample with the highest distance diversity
  f <- fs[unlist(sim_shape_n)]
  f <- f |> remove_tem() |> as_ggplot()
  
  ggsave(filename = "subsample_shapeDiversity.png", plot = print(f), 
         path = o_path, device = "png")
  
  #sample with most similar shape diversity to the original
  g <- fs[unlist(sim_loading_n)]
  g <- g |> remove_tem() |> as_ggplot()
  
  ggsave(filename = "subsample_shapeType.png", plot = print(g),
         path = o_path, device = "png")
  
  #sample with highest distance diversity and most similar shape
  h <- fs[unlist(sim_loadshape_n)]
  h <- h |> remove_tem() |> as_ggplot()
  
  ggsave(filename = "subsample_DiversityType.png", plot = print(h), 
         path = o_path, device = "png")
  
  
  sim_shape_list <- list(shapes = sim_shape[1,], shape_imag = sim_shape_n)
  sim_loadings_list <- list(types = sim_loading[1,], type_imag = sim_loading_n)
  sim_loadshape_list <- list(shapes_types = sim_loadshape, shapes_types_imag = sim_loadshape_n)
  
  return(list(sim_shape_list, sim_loadings_list, sim_loadshape_list))
  
}



formato <- function(data){
  format(data, digit = 3)
}

