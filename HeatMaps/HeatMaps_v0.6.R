# HeatMaps v0.6
    # Thomas Ashhurst
    # Sydney Cytometry Facility

###################################################### 1. INSTALL AND-OR LOAD PACKAGES ###################################################### 

### 1.1 - Install packages (if not already installed)
    if(!require('flowCore')) {install.packages('flowCore')}
    if(!require('plyr')) {install.packages('plyr')}
    if(!require('data.table')) {install.packages('data.table')}
    if(!require('rstudioapi')) {install.packages('rstudioapi')}
    if(!require('R.utils')) {install.packages('R.utils')}
    if(!require('ggplot2')) {install.packages('ggplot2')}
    if(!require('gplots')) {install.packages('gplots')}
    if(!require('RColorBrewer')) {install.packages('RColorBrewer')}
    if(!require('tidyr')) {install.packages('tidyr')}
    if(!require('Biobase')) {install.packages('Biobase')}

### 1.2 Load packages       
    library('flowCore')
    library('plyr')
    library('data.table')
    library('rstudioapi')
    library('R.utils')
    library('ggplot2')
    library('gplots')
    library('RColorBrewer')
    library('tidyr')
    library('Biobase')

### 1.3 - Set working directory and assign as 'PrimaryDirectory'
    
    ## In order for this to work, a) rstudioapi must be installed and b) the location of this .r script must be in your desired working directory
    dirname(rstudioapi::getActiveDocumentContext()$path)            # Finds the directory where this script is located
    setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Sets the working directory to where the script is located
    getwd()
    PrimaryDirectory <- getwd()
    PrimaryDirectory
    
    ## Create an output directory
    dir.create("Output_HeatMap")
    
    ## Create colour schemes
    colour.palette <- (colorRampPalette(brewer.pal(9, "YlGnBu"))(31)) # 256
    group.col <- (colorRampPalette(brewer.pal(11, "Spectral")))  
    fold.palette <- colorRampPalette(rev(c("#ffeda0","#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#bd0026","#800026","black","#023858","#045a8d","#0570b0","#3690c0","#74a9cf","#a6bddb","#d0d1e6","#ece7f2")))

    
### 1.4 - READ DATA
   
    ## Read files
    list.files(path=PrimaryDirectory, pattern = ".csv")
    heatmap.input <- read.csv("demo_cellspertissue_fold_log2.csv")
    
    ## What is the name of the column that denotes the ROW names (usually sample names)
    as.matrix(names(heatmap.input))
    samp.col <- "Sample"
    
    ## Select any columns that don't reflect cellular parameter (sample name, group name etc)
    names(heatmap.input)
    col.rmv <- c(1:2) # Choose which columns to REMOVE (remove columns that are Sample names, Group names, etc, leaving only actual measured parameters)
    
    heatmap.input[c(col.rmv)] # check columns to be REMOVED
    heatmap.input[-c(col.rmv)] # check columns to be KEPT
    

### 1.5 - Preferences

    ## Basics
        plot.title        <- "Cells per cluster (fold-change)"
        do.transpose      <- 0        # Default = 0. Turns rows into columns and columns into rows BEFORE the normalisation or clustering is performed. By default, rows are samples, and columns are clusters.
        
        plot.width        <- 11.69    # Default = 11.69 inches, target for A4 landscape.
        plot.height       <- 8.27     # Default = 8.26 inches, target for A4 landscape.
        
    ## Fold-change options
        is.foldchange     <- 1        # Default = 1. Is the data setup as a fold-change data in Log2? IF THIS IS SELECTED, THEN do.normalise WILL NOT BE PERFORMED.

            # IF data IS fold-change
            max.range         <- 3        # Maximum value of symmetrical range. Default = 3
            min.range         <- -3       # Minimum value of symmetrical range. Default = -3
            (max.range + min.range)       # For fold-change, must return 0.
              
            # IF data IS NOT fold-change
            do.normalise      <- 0        # default = 0. This will normalise data PER COLUMN. If any NaN's are present, normalise function will prevent the full script from working.

    ## Clustering of columns and rows
        do.dendrogram     <- 1        # This will create row and/or column dendrograms. Default = 1.
        dendro.set        <- "both"   # Do you want to create dendrograms for rows and/or columns? "none", "both", "row", "column"
            
        row.clustering    <- TRUE     # TRUE or FALSE
        n.row.groups      <- 3        # choose number of row groupings (by default, rows are samples, unless do.transpose = 1). Must be >1.
            
        col.clustering    <- TRUE     # TRUE or FALSE
        n.col.groups      <- 3        # choose number of column groupings (by default, columns are populations, unless do.transpose = 1). Must be >1.
        
    
        
#########################################################################################
################################### END USER INPUT ######################################
#########################################################################################

        
###################################################### 2. DATA PREPARATION ###################################################### 
    
### 2.1 - Embed cluster or population name as row name
    heatmap.data <- heatmap.input
    rownames(heatmap.data) <- t(heatmap.input[samp.col])
    heatmap.data
    
    heatmap.data <- heatmap.data[-c(col.rmv)] # remove columns
    heatmap.data

### 2.2 - Transpose (ONLY IF REQUIRED) -- the longest set (clusters or parameters) on x-axis -- by default MARKERS are columns, CLUSTERS are rows -- transpose to flip these defaults
    if(do.transpose == 1){
      heatmap.data.t <- as.data.frame(t(heatmap.data))
      heatmap.data <- heatmap.data.t
    }
    
### 2.3 - NORMALISE BY COLUMN (i.e. each column/parameter has a max of 1 and a minimum of 0) # This is optional, but allows for better comparison between markers
    if(do.normalise == 1){
      if(is.foldchange == 0){
        row.nam <- row.names(heatmap.data)
        
        norm.fun <- function(x) {(x - min(x, na.rm=TRUE))/(max(x,na.rm=TRUE) -min(x, na.rm=TRUE))}
        heatmap.data.norm <- as.data.frame(lapply(heatmap.data, norm.fun)) # by default, removes the names of each row
        max(heatmap.data.norm)
        max(heatmap.data.norm[,5])
        heatmap.data.norm <- as.matrix(heatmap.data.norm)
        heatmap.data <- heatmap.data.norm
        rownames(heatmap.data) <- row.nam # add row names back
      }
    }
    
    # convert to matrix
    heatmap.data <- as.matrix(heatmap.data)


### 2.4 - Set up clustering
    
    if(do.dendrogram == 1){
      
      # set the custom distance and clustering functions, per your example
      hclustfunc <- function(x) hclust(x, method="complete")
      distfunc <- function(x) dist(x, method="euclidean")
      
      # perform clustering on rows and columns
      cl.row <- hclustfunc(distfunc(heatmap.data))
      cl.col <- hclustfunc(distfunc(t(heatmap.data)))
      
      # work out no cols and rows
      nrow(heatmap.data)
      ncol(heatmap.data)
      
      # extract cluster assignments; i.e. k=8 (rows) k=5 (columns)
      gr.row <- cutree(cl.row, n.row.groups)
      gr.col <- cutree(cl.col, n.col.groups)
      
      col1 <- group.col(n.row.groups)
      col2 <- group.col(n.col.groups)
    }


### 2.6 - Set up fold-change or normal preferences

    if(is.foldchange == 1){
      map.colour <- fold.palette(31)
      sym.key <- FALSE # TRUE need if NOT creating own breaks
      sym.breaks <- TRUE
      my.breaks <- seq(min.range, max.range, length.out = 32)
    }

    if(is.foldchange == 0){
      map.colour <- colour.palette
      sym.key <- FALSE
      sym.breaks <- FALSE
      my.breaks <- seq(min(heatmap.data), max(heatmap.data), length.out = 32)
    }

    scale.set <- "none" # Can be: "none", "both", "row", "column"

    
### 2.7 - Plot heatmap and dendrograms
    
    dev.off()    
    par(cex.main=1.3)
    par(xpd = TRUE)
    
    pdf("Output_HeatMap/HeatMap.pdf",width=plot.width,height=plot.height)
    
    HeatMap <- heatmap.2(as.matrix(heatmap.data),  #### http://stanstrup.github.io/heatmaps.html # https://stackoverflow.com/questions/22278508/how-to-add-colsidecolors-on-heatmap-2-after-performing-bi-clustering-row-and-co
                         main=plot.title,
                         notecol= "black", # font colour of all cells to black
                         key=TRUE,
                         #keysize=0.75,#key.title = NA, 
                         #key.par = list(cex=1),
                         
                         hclust=hclustfunc, distfun=distfunc,
                         
                         dendrogram =dendro.set, # "both", "row", "column" etc
                         Rowv = row.clustering, # turns on row clustering # = as.dendrogram(cluster) 
                         Colv = col.clustering, # turns on column clustering # = reorder(as.dendrogram(cluster.row), 10:1) to flip around, or maybe dictage
                         RowSideColors = col1[gr.row],
                         ColSideColors = col2[gr.col], # as.character(clusters)                     
                         
                         revC = FALSE, # default FALSE
                         symm=FALSE, # default FALSE
                         symkey= sym.key, # default FALSE, TRUE FOR SYMMETRY IN FOLD CHANGE SITUATIONS
                         symbreaks=sym.breaks, # default FALSE
                         
                         breaks=my.breaks,
                         
                         scale=scale.set, # "none" or "row" or "column". default is to scale by row, I think? # scale argument only for heat data, not dendrogram
                         trace="none", # trace lines inside the heatmap
                         cexRow=1.1,
                         cexCol=1.1, 
                         col=map.colour,
                         #breaks=palette.breaks, # pairs.breaks or "none"
                         
                         #margins = c(9,18), # extra Y and X margins respectively # 9 and 18
                         #lhei = c(1,4), # lhei = c(1,7), # relative height of the legend and plot respectively # default is c(1.5,4,1)
                         #lwid = c(1,4), # lwid = c(1,4), # relativewidth of the legend and plot respectively # default is c(1.5,4)
                         
                         sepwidth = c(0.1, 0.1),
                         sepcolor = "white",
                         #rowsep = c(6, 8, 10, 14, 17, 20, 21, 23, 24),
                         #colsep = c(1, 9, 10, 14, 20, 22, 24, 28, 31),
                         
                         density.info="histogram") # "histogram" "density" "none"
    dev.off() 
    
  
  