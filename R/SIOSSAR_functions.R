#' @title RandomWalk function
#'
#' @param igraphM  matrix of adjacency,the rows and columns are spot
#' @param VertexWeight  a vector v with a length equal to the total number of spots in the slice
#' @param EdgeWeight  logical value, default is TRUE
#' @param gamma  restart probability,default is 0.7
#'
#' @return returns a vector ,reflecting the probability each spot is labeled during the propagation
#' @export
#' @examples
#' \dontrun{
#' # This example requires a valid adjacency matrix and seed vector
#' # visProbs <- RandomWalk2igraph(matrix1, Vector1, EdgeWeight = TRUE, gamma = 0.7)
#' }
#'
RandomWalk2igraph<-function(igraphM,VertexWeight,EdgeWeight=TRUE,gamma=0.7){
  W<-igraphM
  p0 <- VertexWeight

  p0<-t(p0)
  p0 <- p0/sum(p0)
  PT <- p0
  k <- 0
  delta <- 1
  Ng <- dim(W)[2]
  for (i in 1:Ng) {
    sumr<-sum(W[,i])
    if(sumr==0){
      W[,i] <-numeric(length=length(W[,i]))
    }
    if(sumr>0){
      W[,i] <- W[,i]/sum(W[,i])
    }
  }
  W<-as.matrix(W)
  while(delta>1e-10) {
    PT1 <- (1-gamma)*W
    PT2 <- PT1 %*% t(PT)
    PT3 <- (gamma*p0)
    PT4 <- t(PT2) + PT3
    delta <- sum(abs(PT4 - PT))
    PT <- PT4
    k <- k + 1
  }
  PT<-t(PT)
  rownames(PT)<-NULL
  return(drop(PT))
}

#' @title The distance matrix undirected network was constructed to calculate the distance between any two spots
#'
#' @param M  A data frame containing the coordinates of the idle data and the barcode
#' @return A data frame containing the Euclidean distance between any two spots
#' @export
#' @importFrom stats dist
#' @importFrom reshape2 melt
#' @importFrom dplyr group_by filter
#' @importFrom tidyr spread
#'
#' @examples
#' \dontrun{
#' dist_long <- construct_adj_long(meta)
#' }
construct_adj_long <- function(M){
  st_dist <- as.matrix(stats::dist(x = cbind(M$x, M$y)))
  rownames(st_dist) <- M$barcode
  colnames(st_dist) <- M$barcode
  st_dist <- as.data.frame(st_dist)

  dist_long <- st_dist
  dist_long <- as.matrix(dist_long)
  dist_long <- reshape2::melt(dist_long)
  colnames(dist_long) <- c("from","to","distance")
  dist_long <- dist_long[-which(dist_long$from==dist_long$to),]
  return(dist_long)
}


#' @title Construct the weighted adjacency matrix
#'
#' @param dist_long_sub  A data frame including the distance between any two spots
#' @export
#' @importFrom tidyr spread
#' @examples
#' \dontrun{
#' Network <- construct_adj_wide(dist_long)
#' }
#'
construct_adj_wide <-function(dist_long_sub){
  dist_long_sub <- as.data.frame(dist_long_sub)
  dist_long_sub$from <- as.character(dist_long_sub$from)
  dist_long_sub$to <- as.character(dist_long_sub$to)
  dist_long_sub$value <- 1/dist_long_sub$distance
  dist_long_sub$distance <- NULL
  wide_data <- spread(dist_long_sub, key = from, value = value,fill = 0)
  rownames(wide_data) <- wide_data$to
  wide_data$to <- NULL
  row_no <- setdiff(colnames(wide_data),rownames(wide_data))
  if(length(row_no)!=0){
    fill_data <- as.data.frame(matrix(0,length(row_no),ncol(wide_data)))
    rownames(fill_data) <- row_no
    colnames(fill_data) <- colnames(wide_data)
    network <- rbind(wide_data,fill_data)
  }else{
    network <-wide_data
  }

  return(network)
}


#' @title Normalization of expression profiles
#'
#' @param object  a seurat object
#' @param norm  normalization method, NormalizeData or SCTransform
#' @return returns a dataframe with genes as rows and genomic loci as columns
#' @export
#' @importFrom Seurat GetAssayData
#' @examples
#' \dontrun{
#' exp <- normdata(object)
#' }
#'
normdata <- function(object, norm = "NormalizeData") {
  ver <- packageVersion("Seurat") >= "5.0.0"

  if (ver) {
    sp_counts <- GetAssayData(
      object = object, assay = "Spatial", layer = "counts"
    )
  } else {
    sp_counts <- GetAssayData(
      object = object, assay = "Spatial", slot = "counts"
    )
  }

  rawdata <- CreateSeuratObject(counts = sp_counts)

  if (norm == "NormalizeData") {
    nordata <- NormalizeData(rawdata, verbose = FALSE)
    if (ver) {
      exp <- as.data.frame(GetAssayData(
        nordata, assay = "RNA", layer = "data"
      ))
    } else {
      exp <- as.data.frame(GetAssayData(
        nordata, assay = "RNA", slot = "data"
      ))
    }
  } else if (norm == "SCTransform") {
    nordata <- SCTransform(rawdata, vst.flavor = "v1", verbose = FALSE)
    if (ver) {
      exp <- as.data.frame(GetAssayData(
        nordata, assay = "SCT", layer = "data"
      ))
    } else {
      exp <- as.data.frame(GetAssayData(
        nordata, assay = "SCT", slot = "data"
      ))
    }
  } else {
    stop("norm must be either 'NormalizeData' or 'SCTransform'")
  }

  return(exp)
}


#' @title Expression profile of secreted proteins
#'
#' @param exp  Expression profiling data with genes as rows and genomic loci as columns
#' @param gene  A vector, a secreted protein gene
#' @export
#' @examples
#' \dontrun{
#' sm <- secret_mean(exp,gene=c("CCL21","ITGB1"))
#' }
#'
secret_mean<- function(exp, gene){
  secretion_exp <- exp[which(rownames(exp) %in% gene),]
  sm <- as.data.frame(apply(secretion_exp,2,mean))
  colnames(sm) <- "Smean"
  sm$barcode <- rownames(sm)
  sm <- sm[order(-sm$Smean),]
  return(sm)
}

#' @title Membrane protein expression profile
#'
#' @param exp  Expression profiling data with genes as rows and genomic loci as columns
#' @param gene  A vector, a membrane protein gene
#' @export
#' @examples
#' \dontrun{
#' mm <- membrane_mean(exp, gene = c("CCR7", "SRC"))
#' }
#'
membrane_mean<- function(exp, gene){
  membrane_exp <- exp[which(rownames(exp) %in% gene),]
  mm <- as.data.frame(apply(membrane_exp,2,mean))
  colnames(mm) <- "Mmean"
  return(mm)
}

#' @title Filter the number of seed nodes
#'
#' @param SM_exp  A data frame containing information on the average secreted protein expression for each spot and whether the expression content is zero
#' @export
#' @examples
#' \dontrun{
#' topn_sm <- top_secret(SM_exp)
#' }
#'
top_secret <- function(SM_exp){
  seed_node <- ceiling(nrow(SM_exp) * 0.01)
  topn_sm <- SM_exp[c(1:seed_node), ]


  zero_count <- sum(topn_sm$SType == 0)

  if (zero_count == seed_node) {

    return(list(status = "break_loop", message = "No expressive seed detected"))

  } else if (zero_count == 0) {

    topn_sm$Sper <- topn_sm$Smean / sum(topn_sm$Smean)
    return(list(status = "success", data = topn_sm))

  } else {

    topn_sm <- topn_sm[topn_sm$SType != 0, ]
    topn_sm$Sper <- topn_sm$Smean / sum(topn_sm$Smean)
    return(list(status = "success", data = topn_sm))
  }
}


#' @title Run restart random walk
#'
#' @param Network  A data frame, adjacency matrix with weights, rows and columns with spots
#' @param all_spot  All spots in the expression profile
#' @param topn_sm  A data frame containing information on the mean expression of secreted proteins used as seed spots
#' @export
#' @examples
#' \dontrun{
#'  aa <- run_rwr(Network,all_spot,topn_sm)
#' }
#'
run_rwr<- function(Network,all_spot,topn_sm){
  matrix1 <- Network[all_spot,all_spot]
  Vector1<-rep(0,length=length(all_spot))
  Vector1[match(topn_sm$barcode,all_spot)]<-topn_sm[topn_sm$barcode,]$Sper
  visProbs<-RandomWalk2igraph(matrix1,Vector1,EdgeWeight=TRUE,gamma=0.7)
  aa<-cbind(all_spot,visProbs)
  colnames(aa)<-c("spotId","globalScore")
  aa <- as.data.frame(aa)
  rownames(aa) <- aa$spotId
  aa$globalScore <- as.numeric(aa$globalScore)
  return(aa)
}


#' @title Integrated scores
#'
#' @param aa  A data frame, run_rwr function run results
#' @param SM_exp  A data frame containing information on the average secreted protein expression for each spot and whether the expression content is zero
#' @param topn_sm  A data frame containing information on the mean expression of secreted proteins used as seed spots
#' @export
#'
#' @examples
#' \dontrun{
#' # These objects are typically generated by the SIOSSAR workflow.
#' # Example :
#' # aa <- run_rwr(Network, all_spot, topn_sm)
#' # topn_sm <- top.secret(SM_exp)
#' # SM_exp <- secret.mean(exp, genes)
#' IS <- ISs(aa, topn_sm, SM_exp)
#' }
#'
#'
ISs<- function(aa,topn_sm,SM_exp){
  aa1 <- aa[-which(aa$spotId %in% topn_sm$barcode),]
  IS <- SM_exp[-which(SM_exp$barcode %in% topn_sm$barcode),]
  IS$globalScore <- aa1[IS$barcode,]$globalScore
  IS$Score_rank <- rank(IS$globalScore)
  IS$Mexp_rank <- rank(IS$Mmean)
  IS$mean_rank <- (IS$Score_rank+IS$Mexp_rank)/2
  if(length(which(IS$globalScore==0|IS$Mmean==0))!=0){
    IS[which(IS$globalScore==0|IS$Mmean==0),]$mean_rank <- 1
  }
  return(IS)
}




#' @title Define hierarchy
#'
#' @param plot_all Results for the scores as well as p-values
#' @return A data frame containing hierarchical information
#' @importFrom stats quantile
#' @export
#'
#' @examples
#' \dontrun{
#' # This function expects a data frame with specific columns from compute_siossar.
#' # See the vignette for a complete example.
#' plot_all <- hierarchy(plot_all)
#' }
#'
hierarchy <- function(plot_all){

  plot_all$level <- NA

  plot_all$level[plot_all$group != "others"] <- "layer 1"


  plot_all_other <- plot_all[plot_all$group == "others", ]


  plot_all_other <- plot_all_other[order(-plot_all_other$scale_mean_rank), ]


  q <- quantile(plot_all_other$scale_mean_rank,
                probs = seq(0, 1, length.out = 10),
                na.rm = TRUE)

  for(i in 1:9){

    idx <- which(
      plot_all$group == "others" &
        plot_all$scale_mean_rank >= q[10-i] &
        plot_all$scale_mean_rank < q[11-i]
    )

    if(length(idx) > 0){
      plot_all$level[idx] <- paste0("layer ", i + 1)
    }
  }

  return(plot_all)
}
#1. When dividing into uniform intervals (seq) based on scale_mean_rank, the score distribution may be uneven, resulting in some intervals where no distinct tier is identified.
#2. When using quantiles, the score differences may be small or negligible, and the distribution may be relatively flat; however, distinct tiers can still be identified. Currently, method 2 is being used.


#' @title Pathway activation maps were calculated as well as scores for each spot
#'
#' @param object  a seurat object
#' @param LR_PPI  Threshold for background interaction pairs (confidence level),e.g. low-, medium-, and high-confidence
#' @param pathway  A vector that needs to run the pathway
#' @param norm  normalization method,NormalizeData or SCTransform
#' @param verbose  Logical value, whether to print detailed progress information
#' @return returns a dataframe with cells as rows and mca coordinates as columns
#' @export
#' @import Seurat
#' @importFrom dplyr group_by filter
#' @importFrom dplyr sample_n
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom magrittr %>%
#'
#' @usage compute_siossar(object, LR_PPI = "high_ppi_list", pathway,
#'                         norm="NormalizeData" ,verbose = TRUE)
#'
#'
compute_siossar <- function(object, LR_PPI = "high_ppi_list", pathway, norm="NormalizeData" ,verbose = TRUE){

  data(list = LR_PPI, package = "SIOSSAR", envir = environment())
  LR_PPI <- get(LR_PPI, envir = environment())


  PPI_list <- list()   # Save the PPI data for all pathways

  total_start_time <- Sys.time()
  if(verbose) cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  if(verbose) cat("Starting SIOSSAR analysis for", length(pathway), "pathways\n")
  if(verbose) cat(paste(rep("=", 60), collapse = ""), "\n\n")

  #Location information
  if(verbose) cat("Step 1/3  Preparing spatial coordinates and network...\n")
  if(packageVersion("Seurat")<"5.0.0"){
    meta <- object@images[["image"]]@coordinates
    meta$barcode <- rownames(meta)
    colnames(meta) <-c("tissue","array_row","array_col","x","y","barcode")
  }else{
    meta <- GetTissueCoordinates(object)
    meta$barcode <- rownames(meta)
  }
  dist_long <- construct_adj_long(meta)
  dist_long_sub <- dist_long %>% group_by(from) %>% dplyr::filter(distance %in% sort(distance)[1:10])
  Network <- construct_adj_wide(dist_long_sub)
  if(verbose) cat("  [OK] Network constructed with", nrow(Network), "spots\n")

  if(verbose) cat("Step 2/3  Normalizing expression data...\n")
  exp <- normdata(object,norm)
  if(verbose) cat("  [OK] Expression data normalized with", nrow(exp), "genes\n")

  if(verbose) cat("Step 3/3  Processing pathways...\n\n")

  #Cyclic running pathway
  for (n in 1:length(pathway)) {
    pathway_start_time <- Sys.time()

    if(verbose){
      cat("  ", paste(rep("-", 50), collapse = ""), "\n")
      cat("  [", n, "/", length(pathway), "] Processing pathway  ", pathway[n], "\n", sep = "")
      cat("  ", paste(rep("-", 50), collapse = ""), "\n")
    }

    # Step 1  Extract PPI for current pathway
    if(verbose) cat("    Step 1/6  Extracting PPI interactions... ")
    PPI <- LR_PPI[[pathway[n]]]
    #PPI <- LR_PPI[which(LR_PPI$Description ==pathway[n]),]
    PPI_sub <-PPI[which(PPI$source %in% rownames(exp) & PPI$target %in% rownames(exp)),]
    PPI_sub <- unique(PPI_sub[,c(1,2)])
    if(verbose) cat("Found", nrow(PPI_sub), "interactions\n")


    PPI_list[[pathway[n]]] <- PPI_sub

    if(nrow(PPI_sub)!=0){
      # Step 2  Calculate expression means
      if(verbose) cat("    Step 2/6  Calculating secreted and membrane protein expression... ")
      sm <- secret_mean(exp,PPI_sub$source)
      mm <- membrane_mean(exp,PPI_sub$target)
      SM_exp <- sm
      SM_exp$Mmean <- mm[SM_exp$barcode,1]
      #The expression of secreted protein was judged to be 0
      SM_exp$SType <- ifelse((SM_exp$Smean>0),1,0)
      if(verbose) cat("Done\n")

      # Step 3  Select seed nodes
      if(verbose) cat("    Step 3/6  Selecting seed nodes (top 1%)... ")
      topn_sm_result <- top_secret(SM_exp)

      if (topn_sm_result$status == "break_loop") {
        cat("\n    [WARNING]", topn_sm_result$message, "\n")
        cat("    Skipping pathway  ", pathway[n], " due to no expressive seed\n", sep = "")
        next
      } else if (topn_sm_result$status == "success") {
        topn_sm <- topn_sm_result$data
        if(verbose) cat("Found", nrow(topn_sm), "seed nodes\n")
      }

      all_spot <- SM_exp$barcode

      # Step 4  Run random walk
      if(verbose) cat("    Step 4/6  Running random walk with restart... ")
      aa <- run_rwr(Network,all_spot,topn_sm)
      if(verbose) cat("Done\n")

      # Step 5  Calculate integrated scores
      if(verbose) cat("    Step 5/6  Calculating integrated scores... ")
      IS <- ISs(aa,topn_sm,SM_exp)

      ####Random 100 times to generate a null distribution
      if(verbose) cat("Running 100x randomization...\n")
      random <- as.data.frame(matrix(NA,nrow(IS),100))
      rownames(random) <-IS$barcode
      set.seed(123)

      if(verbose) {
        pb <- txtProgressBar(min = 0, max = 100, style = 3, char = ".")
      }

      for (h in 1:100) {
        Rdist_long_sub <- dist_long %>% group_by(from) %>% sample_n(10, replace = FALSE)
        RNetwork <- construct_adj_wide(Rdist_long_sub)
        Random_aa <- run_rwr(RNetwork,all_spot,topn_sm)

        Random_IS <-SM_exp
        Random_IS$random_Mexp <- sample(Random_IS$Mmean)

        Random_IS <- ISs(Random_aa,topn_sm,Random_IS)
        random[,h] <-Random_IS[rownames(random),]$mean_rank

        if(verbose && h %% 10 == 0) setTxtProgressBar(pb, h)
      }

      if(verbose) {
        close(pb)
        cat("\n")
      }

      random$real<-IS$mean_rank
      random$p<-apply(random,1,function(x){length(which(x[1:100] > x[101]))})
      IS$p_value<-random$p/100
      if(length(which(IS$mean_rank==1))!=0){
        IS[which(IS$mean_rank==1),]$p_value <- 1
      }
      IS_signif <- IS[which(IS$p_value<0.1),]

      # Step 6  Summarize scores and define hierarchy
      if(verbose) cat("    Step 6/6  Summarizing scores and defining hierarchy... ")
      plot_all <- meta
      rownames(plot_all) <- plot_all$barcode
      plot_all$mean_rank<- IS[plot_all$barcode,]$mean_rank
      plot_all$p_value <- IS[plot_all$barcode,]$p_value

      plot_all[which(plot_all$barcode %in% topn_sm$barcode),(ncol(plot_all)+1)] <- "seed node"
      plot_all[which(plot_all$barcode %in% IS_signif$barcode),ncol(plot_all)] <- "activated spots"
      plot_all[which(is.na(plot_all[,ncol(plot_all)])),ncol(plot_all)] <- "others"
      colnames(plot_all)[ncol(plot_all)] <- "group"

      plot_all$modif_mean_rank <- ifelse(plot_all$group=="seed node",nrow(plot_all),plot_all$mean_rank)
      plot_all$scale_mean_rank <- plot_all$modif_mean_rank/nrow(plot_all)

      plot_all <- hierarchy(plot_all)

      object@meta.data[[pathway[n]]] <- plot_all[colnames(object),]$scale_mean_rank
      object@meta.data[[paste0(pathway[n],"_region")]] <- plot_all[colnames(object),]$group
      object@meta.data[[paste0(pathway[n],"_level")]] <- plot_all[colnames(object),]$level

      if(verbose) cat("Done\n")

      pathway_end_time <- Sys.time()
      pathway_elapsed <- difftime(pathway_end_time, pathway_start_time, units = "mins")

      if(verbose){
        cat("    [OK] Pathway", pathway[n], "completed in",
            round(pathway_elapsed, 2), "minutes\n")


        n_seed <- nrow(topn_sm)
        n_activated <- nrow(IS_signif)
        cat("      Summary ", n_seed, "seed nodes,", n_activated, "activated spots\n")
      }

    } else {
      if(verbose) cat("    [WARNING] Pathway not found in background interaction data\n")
    }

    if(verbose && n < length(pathway)){
      elapsed_total <- difftime(Sys.time(), total_start_time, units = "mins")
      avg_time <- elapsed_total / n
      remaining <- avg_time * (length(pathway) - n)
      cat("    Progress ", n, "/", length(pathway),
          "pathways processed (", round(n/length(pathway)*100, 1), "%)\n", sep = "")
      cat("    Estimated remaining time ", round(remaining, 2), "minutes\n\n")
    }
  }

  total_end_time <- Sys.time()
  total_elapsed <- difftime(total_end_time, total_start_time, units = "mins")

  if(verbose){
    cat("\n", paste(rep("=", 60), collapse = ""), "\n")
    cat("SIOSSAR analysis completed!\n")
    cat("Total time ", round(total_elapsed, 2), "minutes\n")
    cat("Processed", length(pathway), "pathways\n")
    cat(paste(rep("=", 60), collapse = ""), "\n")
  }

  object@misc$SIOSSAR_PPI <- PPI_list
  return(object)
}





#' @title Output hierarchical results
#'
#' @param out_data output of 'compute_siossar'
#' @param pathway pathway
#' @param show_layers numeric vector of layers to display,Default is 10
#' @param layer_colors optional named color vector
#' @return Output graph
#' @export
#'
#' @importFrom Seurat SpatialDimPlot
#' @importFrom ggplot2 scale_fill_manual theme element_text element_rect
#' @examples
#' \dontrun{
#' p <- hplot(object, "glycolytic_process",1:3)
#' }
hplot <- function(out_data, pathway, show_layers = 1:10, layer_colors = NULL){

  m <- out_data@meta.data
  layer_col <- paste0(pathway, "_level")


  if(!layer_col %in% colnames(m)){
    stop(paste("Layer column", layer_col, "not found in metadata"))
  }


  m$layer_num <- as.numeric(gsub("layer ", "", m[, layer_col]))


  acti <- m[m$layer_num %in% show_layers, ]


  if(nrow(acti) == 0){
    warning("No spots found in the specified layers")
    return(NULL)
  }

  # subset Seurat object
  subset_obj <- out_data[, rownames(acti)]


  if(is.null(layer_colors)){
    layer_colors <- c(
      "#583C88FF", "#764395FF", "#924CA1FF",
      "#AB57A9FF", "#BE68ACFF", "#CF7BAFFF",
      "#DD8EB3FF", "#EAA2B8FF", "#F4B7C0FF", "#FBCCCCFF"
    )
  }


  layer_colors <- layer_colors[1:length(show_layers)]
  names(layer_colors) <- paste0("layer ", show_layers)


  subset_barcodes <- rownames(subset_obj@meta.data)


  subset_layer_nums <- m[subset_barcodes, "layer_num"]


  if(any(is.na(subset_layer_nums))){
    warning("Some spots have NA layer numbers")

    valid_spots <- !is.na(subset_layer_nums)
    subset_obj <- subset_obj[, valid_spots]
    subset_layer_nums <- subset_layer_nums[valid_spots]
  }


  subset_obj@meta.data$plot_layer <- factor(
    paste0("layer ", subset_layer_nums),
    levels = paste0("layer ", show_layers)
  )


  p <- SpatialDimPlot(
    subset_obj,
    group.by = "plot_layer",
    pt.size.factor = 2,
    crop = TRUE,
    image.alpha = 0.7,
    stroke = NA
  ) +
    scale_fill_manual(values = layer_colors, na.value = "grey50") +
    theme(
      #legend.position = "top",
      #legend.title.position = "top",
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 10),
      plot.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white")
    )

  return(p)
}



#' @title Output score results
#'
#' @param out_data output of 'compute_siossar'
#' @param pathway pathway
#' @return Output graph
#' @export
#' @importFrom Seurat SpatialDimPlot
#' @importFrom ggplot2 scale_fill_gradientn theme element_text element_rect
#' @importFrom scales rescale
#' @importFrom Seurat SpatialFeaturePlot
#' @examples
#' \dontrun{
#' p <- score_plot(object, "glycolytic_process")
#' }
score_plot <- function(out_data,pathway){

  p <- SpatialFeaturePlot(out_data, features =pathway, pt.size.factor =2, alpha = c(0.3, 1),
                          crop =T,
                          image.alpha =0.7,stroke =NA)+
    theme(#legend.position = "top",
          #legend.title.position="top",
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 10),
          plot.background = element_rect(fill = "white"),
          legend.background = element_rect(fill = "white")) +
    scale_fill_gradientn(colors = c("#FCF8BAFF", "#FECA8DFF", "#C93E73FF", "#331067FF"),
                         values = rescale(c(0, 0.5, 0.75, 1)))
}




#' @title Output region results
#'
#' @param out_data output of 'compute_siossar'
#' @param pathway pathway
#' @return Output graph
#' @export
#' @importFrom Seurat SpatialDimPlot
#' @importFrom ggplot2 scale_fill_manual theme element_text element_rect
#' @examples
#' \dontrun{
#' p <- hot_region_plot(object, "glycolytic_process")
#' }
hot_region_plot <- function(out_data,pathway){

  p <- SpatialDimPlot(out_data, label = F, label.size = 2, pt.size.factor =1.5,alpha = 1,
                      group.by = paste0(pathway,"_region"), image.alpha = 1,stroke =NA)+
    theme(#legend.position = "top",
          #legend.title.position="top",
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          plot.background = element_rect(fill = "white"),
          legend.background = element_rect(fill = "white"))+
    scale_fill_manual(values = c("seed node"="#FFB900",
                                 "activated spots"="#5773CC",
                                 "others"="grey90"))

  return(p)
}



# Visualization
#

#' @title Output Contribution of SP-MP
#'
#' @param object output of 'compute_siossar'
#' @param pathway pathway
#' @return Output graph
#' @export
#' @importFrom ggalluvial geom_alluvium geom_stratum
#' @importFrom dplyr group_by summarise mutate
#' @importFrom ggplot2 ggplot aes geom_text after_stat scale_x_discrete theme_bw theme element_blank
#' @importFrom tidyr gather spread
#' @importFrom magrittr %>%
#' @importFrom dplyr n
#' @examples
#' \dontrun{
#' p <- plot_PPI_sankey(object,"glycolytic_process")
#' }
#'
plot_PPI_sankey <- function(object, pathway){

  PPI_list <- object@misc$SIOSSAR_PPI

  all_df <- data.frame()

  for(p in pathway){

    df <- PPI_list[[p]]

    if(!is.null(df) && nrow(df)>0){

      df$pathway <- p
      colnames(df)[1:2] <- c("ligand","receptor")

      all_df <- rbind(all_df, df)

    }
  }

  sankey_df <- all_df %>%
    group_by(pathway, ligand, receptor) %>%
    summarise(freq = n(), .groups="drop")

  sankey_df <- sankey_df %>%
    mutate(id = paste(pathway, ligand, receptor, sep="_"))

  p <- ggplot(
    sankey_df,
    aes(axis1 = pathway,
        axis2 = ligand,
        axis3 = receptor,
        y = freq)
  ) +
    geom_alluvium(aes(fill = pathway), width = 1/12, alpha = 0.8) +
    geom_stratum(width = 1/12, fill = "grey90", color = "black") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
    scale_x_discrete(
      limits = c("Pathway","Ligand","Receptor"),
      expand = c(.05, .05)
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      axis.title = element_blank()
    )

  return(p)

}




#' @title Output top 20 pathway
#'
#' @param object output of 'compute_siossar'
#' @param pathway pathway
#' @param top_n Select the number of top pathways to display based on the count of activated spots
#' @return Output graph
#' @export
#' @importFrom pheatmap pheatmap
#' @examples
#' \dontrun{
#' p <- plot_top_pathway_heatmap(object, pathway, top_n = 1)
#' }
plot_top_pathway_heatmap <- function(object, pathway, top_n = 20){


  stat <- data.frame()

  for(p in pathway){

    region <- object@meta.data[[paste0(p,"_region")]]

    n_active <- sum(region != "others")

    ratio <- n_active / length(region)

    stat <- rbind(
      stat,
      data.frame(
        pathway=p,
        active_spot=n_active,
        ratio=ratio
      )
    )
  }

  stat <- stat[order(-stat$active_spot),]

  top_pathway <- stat$pathway[1:top_n]

  mat <- matrix(
    NA,
    nrow=length(top_pathway),
    ncol=ncol(object)
  )

  rownames(mat) <- top_pathway
  colnames(mat) <- colnames(object)

  for(i in 1:length(top_pathway)){

    p <- top_pathway[i]

    mat[i,] <- object@meta.data[[p]]

  }

  p <- pheatmap(
    mat,
    cluster_rows=F,
    cluster_cols=TRUE,
    show_colnames=FALSE
  )

  return(p)

}
