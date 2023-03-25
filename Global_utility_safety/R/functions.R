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

euclidean <- function(a, b) sqrt(sum((a - b)^2))



utility <- function(facesO, facesidO, facesS, facesidS, o_path, 
                    lm_n, seed, nreplica, nparent, nm){
  
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
  #          the name included as name (e.g. "s2_" ) (nm)
  #                    use "_" in this name to separate nm from the iterator (S2_i, S2_i+1, S2_i+2...)
  
  # getting synthetic images
  fs <- facesS
  
  ## LANDMARKS ##
  # getting landmarks for original and synthetic data
  lmO <- get_point(facesO, 0:(lm_n-1))
  lmO <- as.matrix(lmO[,c(3,4)]) 
  lmO <- arrayspecs(lmO, lm_n, 2)
  
  lmS <- get_point(fs, 0:(lm_n-1))
  lmS <- as.matrix(lmS[,c(3,4)]) 
  lmS <- arrayspecs(lmS, lm_n, 2)
  
  # changing array format 
  lmO <- two.d.array(lmO) 
  attributes(lmO)$dimnames[[1]] <- facesidO
  
  lmS <- two.d.array(lmS)
  attributes(lmS)$dimnames[[1]] <- names(fs)
  
  ## Choosing the most diverse synthetic sample 
  ## the one that matches closely the original diversity
  set.seed(seed)
  
  ## ORIGINAL DATA ##
  # Generalized Procrustes analysis
  original_gpa <- gpagen(arrayspecs(lmO, lm_n, 2), print.progress = F, verbose = T, ProcD = F)
  original_alig <- rotate.coords(original_gpa$coords, type = "rotateC")
  original_pca <- gm.prcomp(original_alig)
  
  # Using PC scores to compute the Procrustes distances among original images
  original_dist <- unlist(as.vector(dist(two.d.array(original_alig))))
  # Minimum and Maximum of the Procrustes distances
  minOO <- min(original_dist)
  maxOO <- max(original_dist)
  
  ## SYNTHETIC DATA ##
  # Sampling the synthetic database "nreplica" times
  nsample = dim(lmO)[1]
  sample_list0 <- replicate(nreplica, sample(attributes(lmS)$dimnames[[1]], nsample), simplify=FALSE)

  # Generalized Procrustes analysis in the nreplicas generated
  sample_list1 <- lapply(sample_list0, function(x) lmS[x,])
  sample_list2 <- lapply(sample_list1, function(x) arrayspecs(x, lm_n, 2))
  sample_gpa <- lapply(sample_list2, function(x) gpagen(x, print.progress = FALSE,
                                                        verbose = TRUE, ProcD = F))
  sample_aligned <- lapply(sample_gpa, function(x) rotate.coords(x$coords, type = "rotateC"))
  sample_pca <- lapply(sample_aligned, function(x) gm.prcomp(x))
  
  # Using PC scores to compute the Procrustes distances among synthetic images in the nreplicas generated
  sample_distS <- lapply(sample_pca, function(x, i) {unlist(as.vector(dist(x$x)))})
  # Minimum, Maximum and mean of those distances
  sample_descr <- lapply(sample_distS, function(x) c(min(x), mean(x), max(x)))
  sample_df <- do.call(rbind.data.frame, sample_descr)
  colnames(sample_df) <- c("min", "mean", "max")
  
  # choosing the replica with the highest distance values
  higher <- rownames(sample_df[order( -sample_df[,1], -sample_df[,3] ),][1,]) 
  # Minimum and Maximum of the Procrustes distances in the selected replica
  synt_dist <- sample_distS[[as.numeric(higher)]]
  minSS <- min(synt_dist)
  maxSS <- max(synt_dist)
  
  ## INTER-INDIVIDUAL DISTANCE DATABASE ##
  disDF <- data.frame(distance = c(unlist(original_dist),
                                   unlist(synt_dist)),
                      Type = c(rep("Original", length(original_dist)),
                               #using "_" in the name (nm) is very useful in the next function
                               rep(gsub('_','',nm), length(synt_dist))))
  

  ## VISUAL COMPARISON OF THE INTER-INDIVIDUAL DISTANCES ##
  title <- paste("Original (min = ",formato(minOO),", max = ", formato(maxOO),
                 ") and \n synth (min = ", formato(minSS),",max = ", 
                 formato(maxSS), ") distances.", "\n Using ", nparent, 
                 " parent faces.")

  # Generates Utility plot (FIGURE & IN THE PAPER)
  paleta <- c("#F7B267", "#8DB580")
  alls <- ggplot(disDF, aes(x=distance, y = Type, fill = type)) +
    geom_boxplot(color=paleta, fill=paleta, alpha=0.2) + theme_bw()+
    xlim(0, 0.3) + coord_flip()+
    ggtitle(title) + theme(text = element_text(size = 10)) 
  

  ## SAVES DE PLOT ##
  f_name <- paste("utilityPlot_", nparent, "faces.png", sep = "")
  ggsave(filename = f_name, plot = print(alls), path = o_path, device = "png", width = 7, height = 7)
  
  return(disDF)
}



safety <- function(facesO, facesidO, facesS, facesidS, o_path, 
                   lm_n, seed, nreplica, nparent, nsampSize, nm){
  
  
  # This function computes inter-individual distances (DEPENDENTLY):
  #             for the original database (as in the previous function)
  #             between parents and synthetic children
  
  
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
  #          the sample size of the original database (nsampleSize)
  #          the name included as name (e.g. "s2_" ) (nm)
  #                    use "_" in this name to separate nm from the iterator (S2_i, S2_i+1, S2_i+2...)

  
  ### COMPARING PARENTS WITH DESCENDENTS
  #read synthetic images
  fs <- facesS
  
  # Get landmarks
  lmO <- get_point(facesO, 0:(lm_n-1))
  lmO <- as.matrix(lmO[,c(3,4)])
  lmO <- arrayspecs(lmO, lm_n, 2)
  
  lmS <- get_point(fs, 0:(lm_n-1))
  lmS <- as.matrix(lmS[,c(3,4)])
  lmS <- arrayspecs(lmS, lm_n, 2)
  
  # changing array format
  lmO <- two.d.array(lmO) 
  attributes(lmO)$dimnames[[1]] <- facesidO
  
  lmS <- two.d.array(lmS)
  attributes(lmS)$dimnames[[1]] <- names(fs)
  
  set.seed(seed)
  
  # Sampling the synthetic database "nreplica" times
  nsample = dim(lmO)[1]
  sample_list0 <- replicate(nreplica, sample(attributes(lmS)$dimnames[[1]], nsample), simplify=FALSE)
  
  ## ORIGINAL AND SYNTHETIC DATA ###
  # Generalized Procrustes analysis for the original data and all the synthetic replicas
  or_syn0 <- lapply(sample_list0, function(x) rbind(lmO, lmS[x,]))
  or_syn1 <- lapply(or_syn0, function(x) arrayspecs(x, lm_n, 2))
  or_syn_gpa <- lapply(or_syn1, function(x) gpagen(x, print.progress = FALSE,
                                                   verbose = TRUE, ProcD = F))
  or_syn_align <- lapply(or_syn_gpa, function(x) rotate.coords(x$coords, type = "rotateC"))
  or_syn_pca <- lapply(or_syn_align, function(x) gm.prcomp(x))
  
  # gets and orders the Procrustes (inter-individual) distance generated
  or_syn_dist <- lapply(or_syn_pca, function(x, i) {dist(x$x)})
  or_syn_dist2 <- lapply(or_syn_pca, function(x, i) {unlist(as.vector(dist(x$x)))})
  or_syn_nom <- c(attributes(lmO)$dimnames[[1]],unlist(sample_list0))
  or_syn_order <- lapply(or_syn_dist, function(x) which(lower.tri((x), diag=TRUE), arr.ind=TRUE))
  or_syn_order2 <- or_syn_order[[1]][or_syn_order[[1]][,1] != or_syn_order[[1]][,2], ]
  or_syn_DF <- lapply(or_syn_dist2, function(x) data.frame(distance = x))
  or_syn_DF <- lapply(or_syn_DF, function(x) cbind(x, or_syn_order2))
  

  synt_all <- as.data.frame(t(combn(sort(names(facesO)[1:nsampSize]),nparent)))

  #using "_" in the name (nm) is very useful in the next function
  sample_list0.n <- lapply(sample_list0, function(x)  as.numeric(gsub(nm,'',x)))
  parent <- lapply(sample_list0.n, function(x) synt_all[x,])
 
  if(nparent == 2){
    
    for (i in 1:length(or_syn_pca)) {
      for (j in 1:length(sample_list0[[1]])) {
        #first parent
        parent[[i]]$distP1[j] <- euclidean(or_syn_pca[[i]]$x[rownames(or_syn_pca[[i]]$x) == parent[[i]]$V1[j],],
                                           or_syn_pca[[i]]$x[rownames(or_syn_pca[[i]]$x) == sample_list0[[i]][j],])
        #second parent
        parent[[i]]$distP2[j] <- euclidean(or_syn_pca[[i]]$x[rownames(or_syn_pca[[i]]$x) == parent[[i]]$V2[j],],
                                           or_syn_pca[[i]]$x[rownames(or_syn_pca[[i]]$x) == sample_list0[[i]][j],])
      }
    }
    
  } else if (nparent == 3){
    for (i in 1:length(or_syn_pca)) {
      for (j in 1:length(sample_list0[[1]])) {
        #first parent
        parent[[i]]$distP1[j] <- euclidean(or_syn_pca[[i]]$x[rownames(or_syn_pca[[i]]$x) == parent[[i]]$V1[j],],
                                           or_syn_pca[[i]]$x[rownames(or_syn_pca[[i]]$x) == sample_list0[[i]][j],])
        #second parent
        parent[[i]]$distP2[j] <- euclidean(or_syn_pca[[i]]$x[rownames(or_syn_pca[[i]]$x) == parent[[i]]$V2[j],],
                                           or_syn_pca[[i]]$x[rownames(or_syn_pca[[i]]$x) == sample_list0[[i]][j],])
        #third parent
        parent[[i]]$distP3[j] <- euclidean(or_syn_pca[[i]]$x[rownames(or_syn_pca[[i]]$x) == parent[[i]]$V3[j],],
                                           or_syn_pca[[i]]$x[rownames(or_syn_pca[[i]]$x) == sample_list0[[i]][j],])
      }
    }
    
    
    
  } else if (nparent == 4){
    
    for (i in 1:length(or_syn_pca)) {
      for (j in 1:length(sample_list0[[1]])) {
        #first parent
        parent[[i]]$distP1[j] <- euclidean(or_syn_pca[[i]]$x[rownames(or_syn_pca[[i]]$x) == parent[[i]]$V1[j],],
                                           or_syn_pca[[i]]$x[rownames(or_syn_pca[[i]]$x) == sample_list0[[i]][j],])
        #second parent
        parent[[i]]$distP2[j] <- euclidean(or_syn_pca[[i]]$x[rownames(or_syn_pca[[i]]$x) == parent[[i]]$V2[j],],
                                           or_syn_pca[[i]]$x[rownames(or_syn_pca[[i]]$x) == sample_list0[[i]][j],])
        #third parent
        parent[[i]]$distP3[j] <- euclidean(or_syn_pca[[i]]$x[rownames(or_syn_pca[[i]]$x) == parent[[i]]$V3[j],],
                                           or_syn_pca[[i]]$x[rownames(or_syn_pca[[i]]$x) == sample_list0[[i]][j],])
        
        #fourth parent
        parent[[i]]$distP4[j] <- euclidean(or_syn_pca[[i]]$x[rownames(or_syn_pca[[i]]$x) == parent[[i]]$V4[j],],
                                           or_syn_pca[[i]]$x[rownames(or_syn_pca[[i]]$x) == sample_list0[[i]][j],])
      }
    }
    
  } else {print("nparent must be 2, 3, or 4")}
  
  
  
  parent <- (do.call(rbind.data.frame, parent))
  descent <- (unlist(sample_list0))
  pd <- parent
  pd$Synt <- descent

  
  if(nparent == 2){
    minPD <- min(c(pd$distP1, pd$distP2))
    maxPD <- max(c(pd$distP1, pd$distP2))
    pd_l <- gather(pd, "value", "distance", -c("V1", "V2", "Synt"))
    
  } else if (nparent == 3){
    minPD <- min(c(pd$distP1, pd$distP2, pd$distP3))
    maxPD <- max(c(pd$distP1, pd$distP2, pd$distP3))
    pd_l <- gather(pd, "value", "distance", -c("V1", "V2", "V3", "Synt"))
    
  } else if (nparent == 4){
    minPD <- min(c(pd$distP1, pd$distP2, pd$distP3, pd$distP4))
    maxPD <- max(c(pd$distP1, pd$distP2, pd$distP3, pd$distP4))
    pd_l <- gather(pd, "value", "distance", -c("V1", "V2", "V3", "V4", "Synt"))
    
  } else {print("nparent must be 2, 3, or 4")}
  
  
  
  ## Distances between originals
  or_syn_DF_OO <- lapply(or_syn_DF, function(x) x[x$col <= dim(lmO)[1] & 
                                                    x$row <= dim(lmO)[1],])
  OO <- do.call(rbind.data.frame, or_syn_DF_OO)
  #Removes duplicate rows/columns ignoring column order
  OO <- OO[!duplicated(apply(OO[,2:3], 1, function(row) paste(sort(row), collapse=""))),]
  minOO <- min(OO$distance)
  maxOO <- max(OO$distance)
 
   ### COMPARING SAFETY
  title <- paste("Similarity between original images (min = ",formato(minOO),", max = ", formato(maxOO),
                 ") and \n between parent and synthetic images (min = ", formato(minPD),",max = ", 
                 formato(maxPD), ").", "\n Using ", nparent, 
                 " parent faces.")
  
  

  OO_DF <- data.frame(distance = OO[,"distance"])
  OO_DF$Type <- "Original"
  PD_DF <- data.frame(distance = pd_l[,"distance"])
  #using "_" in the name (nm) is very useful in the next function
  PD_DF$Type <- gsub('_','',nm)
  
  priv_risk <- rbind(OO_DF[, c("distance", "Type")], PD_DF[, c("distance", "Type")])
  paleta <- c("#F7B267", "#8DB580")
  p <- ggplot(priv_risk, aes(x=distance, y = Type, fill = Type)) +
    geom_boxplot(color=paleta, fill=paleta, alpha=0.2) + theme_bw()+
    xlim(0, 0.3) + coord_flip()+
    ggtitle(title) + theme(text = element_text(size = 10))
  
  f_name <- paste("safetyPlot_", nparent, "faces.png", sep = "")
  ggsave(filename = f_name, plot = print(p), path = o_path, device = "png", width = 7, height = 7)
  
  return(priv_risk)
  
}



all_plots <- function(Outil, list_Sutil, list_Ssafe, o_path){
  utilDF <- rbind(Outil, bind_rows(list_Sutil))
  safeDF <- rbind(Outil, bind_rows(list_Ssafe))
  
  paleta <- c("#F7B267", "#8DB580","#8DB580","#8DB580")
  
  alls <- ggplot(utilDF, aes(x=distance, y = Type, fill = Type)) +
    geom_boxplot(color=paleta, fill=paleta, alpha=0.2) + theme_bw()+
    xlim(0, 0.3) + coord_flip()
  
  allr <- ggplot(safeDF, aes(x=distance, y = Type, fill = Type)) +
    geom_boxplot(color=paleta, fill=paleta, alpha=0.2) + theme_bw()+
    xlim(0, 0.3) + coord_flip()
  
  util_name <- "all_utilityPlot_faces.png"
  ggsave(filename = util_name, plot = print(alls), path = o_path, device = "png", width = 7, height = 7)
  
  safe_name <- "all_safetyPlot_faces.png"
  ggsave(filename = safe_name, plot = print(allr), path = o_path, device = "png", width = 7, height = 7)
  
}
