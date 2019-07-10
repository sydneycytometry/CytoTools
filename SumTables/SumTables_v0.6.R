### SumTables v0.5
    # Thomas Ashhurst
    # 2018-09-25
    # thomas.ashhurst@sydney.edu.au
    # www.github.com/SydneyCytometry

###################################################### 1. INSTALL AND-OR LOAD PACKAGES ###################################################### 
    
    ### 1.1 - Install packages (if not already installed)
        if(!require('flowCore')) {install.packages('flowCore')}
        if(!require('plyr')) {install.packages('plyr')}
        if(!require('data.table')) {install.packages('data.table')}
        if(!require('rstudioapi')) {install.packages('rstudioapi')}
        if(!require('R.utils')) {install.packages('R.utils')}
        if(!require('tidyr')) {install.packages('tidyr')}
        if(!require('Biobase')) {install.packages('Biobase')}
    
    ### 1.2 Load packages       
        library('flowCore')
        library('plyr')
        library('data.table')
        library('rstudioapi')
        library('R.utils')
        library('tidyr')
        library('Biobase')

    ### 1.3 - Set working directory and assign as 'PrimaryDirectory'

        ## In order for this to work, a) rstudioapi must be installed and b) the location of this .r script must be in your desired working directory
        dirname(rstudioapi::getActiveDocumentContext()$path)            # Finds the directory where this script is located
        setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Sets the working directory to where the script is located
        getwd()
        PrimaryDirectory <- getwd()
        PrimaryDirectory
        
        ## Use this to manually set the working directory
        #setwd("/Users/Tom/Desktop/Experiment")                          # Set your working directory here (e.g. "/Users/Tom/Desktop/") -- press tab when selected after the '/' to see options
        #getwd()                                                         # Check your working directory has changed correctly
        #PrimaryDirectory <- getwd()                                     # Assign the working directory as 'PrimaryDirectory'
        #PrimaryDirectory

    ### 1.4 - Preferences
    
        # Basis
        use.groups            <- 1        # Add group keyword to output files (currently only works for cell number summaries)
        
        # Cell numbers
        do.CELLNUM            <- 1        # Calculate cell numbers per cluster/sample/file
        do.CELLPROPORTIONS    <- 1        # Calculate cell proprotions per cluster/sample/file
        do.CELLSPERTISSUE     <- 1        # Calculate cell numbers per cluster/sample/tissue        # Turn OFF if you don't have any cell counts
        do.FOLDCHANGE         <- 1        # Calculate cell number/proportion changes as fold-change         # Turn OFF if you don't have group keywords

        # Expression levels
        fun.type              <- "mean"    # Choose summary function for expression levels. Default = "mean". Can be "mean" or "median".
        
        do.MFI.per.sample     <- 1        # Creates a table of the MFI of each marker on each cluster, one file per sample.
        do.MFI.per.marker     <- 1        # Creates a table of the MFI of each marker on each samples, one file per marker. WARNING - runs slowly, only works when clusters are numbers (not population names etc)



###################################################### 2. READ FILES ###################################################### 
    
    # DO ON MAIN FlowSOM FILES, preferably not tSNE files -- numbers too small, differential downsampling
    
    ## Use to list the .fcs or .csv files in the working directory
        as.matrix(list.files(path=PrimaryDirectory, pattern = ".fcs"))     # see a list of FCS files
        as.matrix(list.files(path=PrimaryDirectory, pattern = ".csv"))     # see a list of CSV files
    
    ## Specify general options
        File_Type             <- ".csv"              # Specify the type of input file. Use ".fcs" OR ".csv", in lower case. Fcs files not current supported
        FileNos               <- "single"            # Can be "single" for one file, or "multi" for multiple files. If using 'multi', then only the relevant files should be in the directory. Each row in will be embedded with the name of the name of the file it came from. 

        file.single           <- c(1)               # ONLY REQUIRED IF READING FROM SINGLE FILE. Enter the number of the single file between the brackets
        
    ## Loop for single files

        if(FileNos == "single"){
          if(File_Type == ".csv"){
            files <- list.files(path=PrimaryDirectory, pattern = ".csv")
            data <- fread(files[file.single])
          }
          
          ## FCS files functionality not currently compatible - do not use
          #if(File_Type == ".fcs"){
          #  files <- list.files(path=PrimaryDirectory, pattern = ".fcs")
          #  data <- exprs(read.FCS(files[file.single], transformation = FALSE)) 
          #  data <- parameters(read.FCS(File)) 
          #}
        }
    
    ## Loop for multiple files
        if(FileNos == "multi"){
          DataList=list() # Creates and empty list to start 
          if (File_Type == ".csv"){
            FileNames <- list.files(path=PrimaryDirectory, pattern = ".csv")
            for (File in FileNames) { # Loop to read files into the list
              tempdata <- fread(File)
              File <- gsub(".csv", "", File)
              DataList[[File]] <- tempdata
            }
          }
          
          ## FCS files functionality not currently compatible - do not use
          #if (File_Type == ".fcs"){
          #  FileNames <- list.files(path=PrimaryDirectory, pattern = ".fcs")
          #  ParaList=list()
          #  for (File in FileNames) { # Loop to read files into the list
          #    tempdata <- exprs(read.FCS(FileNames[1], transformation = FALSE)) 
          #    File <- gsub(".fcs", "", File)
          #    DataList[[File]] <- tempdata
          #    ParaList[[File]] <- parameters(read.FCS(File)) 
          #  }
          #}
          data <- rbindlist(DataList) ## Concatenate into one large data frame -- what if they have a column conflict??
        }
        
    ## Review data
        data <- as.data.frame(data)
        dim(data)                             # Review data dimensionality (number of cells, number of parameters)
        head(data)                            # Review first 6 rows of the data

    ## Choose the column that defines the sample names/numbers 
        as.matrix(names(data))
        samp.col <- "SampleName"              # Enter name of the column that indicates SAMPLE names here
        as.matrix(names(data[samp.col]))
    
    ## Choose the column that defines the GROUPS (not required if not using groups and not using fold-change)
        as.matrix(names(data))
        grp.col <- "GroupName"                # Enter name of the column that indicates GROUP names here. Not required if not using groups (i.e. use.groups = 0)
        as.matrix(names(data[grp.col]))    
        
    ## Choose the column that defines the cluster/population names/numbers 
        as.matrix(names(data))
        clust.col <- "FlowSOM_meta_cluster"   # enter the NAME of the column that represents clusters or populations
        clust.col.num <- c(38)                # enter the NUMBER of the column that represents clusters or populations
        
        as.matrix(names(data[clust.col]))

    ## Choose useful annotation columns (i.e. choose which columns show CLUSTER, POPULATION, SAMPLE, and GROUP names or numbers, etc)
        
        as.matrix(names(data))
        annot.col <- c(33:37) # DO NOT INCLUDE THE 'CLUSTER' COLUMN USED ABOVE!
        as.matrix(names(data[annot.col]))
        
    ## Choose control group. ONLY required if doing fold-change comparisons (i.e. do.FOLDCHANGE = 1)
        unique(data[grp.col])
        ctrls <- "Mock_D7" # indicate which group represents the control group
       
    ## Add cell counts. ONLY required if running calculating cells per cluster per tissue (i.e. do.CELLSPERTISSUE = 1)
        
        as.matrix(unique(data$SampleName)) # Review sample names -- counts MUST be in the same order as the sample names
        
        cell.counts           <- c(2e+07,
                                   2e+07,
                                   2e+07,
                                   2e+07,
                                   2e+07,
                                   2e+07,
                                   1.8e+07,
                                   1.8e+07,
                                   1.8e+07,
                                   1.8e+07,
                                   1.8e+07,
                                   1.8e+07) 
        
        length(cell.counts) == length(unique(data[[samp.col]])) # Check to see if the number of cell counts equals the number of sampless. Should return TRUE.
        
  
############################################################################################################################          
###################################################### END USER INPUT ######################################################           
############################################################################################################################        
        
        
###################################################### 3. DATA PREPARATION  ######################################################      
    
    AllSampleNames <- as.matrix(unique(data[samp.col]))
    if(use.groups == 1){AllGroupNames <- as.matrix(unique(data[grp.col]))}

    
###################################################### 4. GENERATE CELL NUMBER SUMMARY TABLES  ######################################################          

    if(do.CELLNUM == 1){
    
      dir.create("Output_CellNum")
      
    ### 4.1 - Work out cells per cluster for each sample    
            
        # Specify which column denotes the sample names, and check sample names -- 'table' function will show how many rows (cells) exist per sample
        data.cellnum <- data
        data.cellnum <- data.cellnum[order(data.cellnum[samp.col]),] # re-order by sample name (alphabetically)
        
        # subset data to cluster number/name, sample number/name, and other parameters that are helpful (e.g. group name, sample name, etc) -- can remove all marker columns (e.g. 141Nd-Ly6G)
        data.pops <- data.cellnum[c(annot.col, clust.col.num)]
        samp.name <- data.pops[samp.col]
        
            #Samples <- unique(data.pops[samp.col])
            #as.matrix(Samples)
            
        samp.list = list()
        
        # Split data into separate 'samples'  - specific sample name, all clusters
        for(i in AllSampleNames) {
          temp <- subset(data.pops, samp.name == i)
          samp.list[[i]] <- temp
        }
        
        # some checks
            samp.list[[1]]
            nrow(samp.list[[1]])
            
            all.samps <- rbindlist(samp.list)
            all.samps
            
        ClustName <- all.samps[[clust.col]]
        
        # number of cells per cluster in each sample and order
        length(as.matrix(unique(ClustName))) # find number of clusters
        
        num.clusters <- as.matrix(unique(ClustName)) # list of cluster numbers or names
        num.clusters <- sort(num.clusters, decreasing = FALSE)

        new.sample.list = list()
        
        for(a in AllSampleNames){
          samp.res = data.frame(Cluster = numeric(0), NumCells= numeric(0), Group=character(0)) # empty dataframe with names 

          for(i in num.clusters){ 
            num.cells <- sum(samp.list[[a]][clust.col] == i) # number of cells per cluster
            #group.label <- samp.list[[a]]$GroupName[1] # picks the group name from the first group of that sample
            
                if(use.groups == 1){
                  group.label <- samp.list[[a]][grp.col][1,1] # picks the group name from the first group of that sample
                  res <- data.frame(Cluster = i, NumCells = num.cells, Group = group.label) # added group label
                }
                
                if(use.groups == 0){
                  res <- data.frame(Cluster = i, NumCells = num.cells) # added group label
                  }

            samp.res <- rbind(samp.res, res)
          }
          
          new.sample.list[[a]] <- samp.res
          new.sample.list[[a]][samp.col] <- a
        }
        
        merged.df <- rbindlist(new.sample.list) ## Concatenate into one large data frame -- what if they have a column conflict??
        
        ## Arrange table (long to wide)
        merged.df <- spread(merged.df, Cluster, NumCells)
        merged.df <- as.data.frame(merged.df)
        
        rownames(merged.df)
        colnames(merged.df)
        merged.df[samp.col] # determine which column has the desired 'row names' -- might be first or second column, depending on whether group ID was added
        
        ## Save table
        write.csv(merged.df, "Output_CellNum/SumTable_CellNum.csv")
        
        # Resulting table represents samples (rows) and clusters (columns). Numbers in each cell represent the number of cells from each cluster in each sample
        ## Data can now be examined
    
    ### 4.2 - Calculate cell proportions
        
        merged.df
 
        # create ROW totals (i.e. total cells per sample)
        x <- merged.df
        x[samp.col] <- NULL
        
        if(use.groups == 1){
          x["Group"] <- NULL
        }
        
        row.totals <- rowSums(x)
        as.matrix(row.totals) # sum(row.totals)
            
            if(use.groups == 1){
              cells.per.file <- cbind(merged.df[c(1:2)], row.totals)
            }
            
            if(use.groups == 0){
              cells.per.file <- cbind(merged.df[c(1)], row.totals)
            }

        cells.per.file
        
        if(use.groups == 1){
          per <- apply(merged.df[-c(1:2)], 1, function(i) i/sum(i)) # https://stackoverflow.com/questions/35691205/divide-each-each-cell-of-large-matrix-by-sum-of-its-row
          per <- t(per)                                                       # 1 = rows, 2 = columns
          sumtables.percentages <- cbind(merged.df[c(1:2)], per)
        }
        
        if(use.groups == 0){
          per <- apply(merged.df[-c(1)], 1, function(i) i/sum(i)) # https://stackoverflow.com/questions/35691205/divide-each-each-cell-of-large-matrix-by-sum-of-its-row
          per <- t(per)                                                       # 1 = rows, 2 = columns
          sumtables.percentages <- cbind(merged.df[c(1)], per)
        }
        
        write.csv(sumtables.percentages, "Output_CellNum/SumTable_CellProportions.csv")
        
        
    ### 4.3 - Calculate CELLS PER TISSUE        
        
        if(do.CELLSPERTISSUE == 1){
        
            # times by cell count # input string of 'cell counts'
            as.matrix(cell.counts)
          
            cellsper <- sweep(per,MARGIN=1,cell.counts,`*`)
            
                if(use.groups == 1){
                  sumtables.cellspertissue <- cbind(merged.df[c(1:2)], cellsper)
                }
                if(use.groups == 0){
                  sumtables.cellspertissue <- cbind(merged.df[c(1)], cellsper)
                }

            write.csv(sumtables.cellspertissue, "Output_CellNum/SumTable_CellsPerTissue.csv")
            }
        
    ### 4.4 - Calculate FOLD CHANGE (PROPORTIONS)
        
        if(do.FOLDCHANGE == 1){
          if(use.groups == 1){
            sumtables.percentages
          
            ctrl.grp <- subset(sumtables.percentages, Group == ctrls)
            ctrl.grp
            
            ctrl.grp.means <- colMeans(ctrl.grp[-c(1:2)])
            as.matrix(ctrl.grp.means)
            
            # As above, divide each 'cell' to the vector (by column...)
            fold.raw <- t(t(per) / ctrl.grp.means)
            fold.raw
            
            sumtables.fold.raw <- cbind(merged.df[c(1:2)], fold.raw)
            write.csv(sumtables.fold.raw, "Output_CellNum/SumTable_Proportions_FoldChangeRaw.csv")
            
            # convert table to log2
            
            fold.log2 <- log(x = fold.raw, 2)
            
            sumtables.fold.log2 <- cbind(merged.df[c(1:2)], fold.log2)
            write.csv(sumtables.fold.log2, "Output_CellNum/SumTable_Proportions_FoldChangeLog2.csv")
            }
          }
        

    ### 4.5 - Calculate FOLD CHANGE (CELLS PER TISSUE)
  
        if(do.FOLDCHANGE == 1){
          if(use.groups == 1){
            if(do.CELLSPERTISSUE == 1){
              
              # average of the columns for rows with specific group name
              sumtables.cellspertissue
              
              ctrl.grp <- subset(sumtables.cellspertissue, Group == ctrls)
              ctrl.grp
              
              ctrl.grp.means <- colMeans(ctrl.grp[-c(1:2)])
              as.matrix(ctrl.grp.means)
          
              # As above, divide each 'cell' to the vector (by column...)
              fold.raw <- t(t(cellsper) / ctrl.grp.means)
              fold.raw
              
              sumtables.fold.raw <- cbind(merged.df[c(1:2)], fold.raw)
              write.csv(sumtables.fold.raw, "Output_CellNum/SumTable_CellsPerTissue_FoldChangeRaw.csv")
          
              # convert table to log2
                
              fold.log2 <- log(x = fold.raw, 2)
            
              sumtables.fold.log2 <- cbind(merged.df[c(1:2)], fold.log2)
              write.csv(sumtables.fold.log2, "Output_CellNum/SumTable_CellsPerTissue_FoldChangeLog2.csv")
            }
          }
        }
        
    } 


        
###################################################### 5. GENERATE MFI SUMMARY TABLES  ######################################################   
        
        if(do.MFI.per.sample == 1){   
          
          setwd(PrimaryDirectory)
          dir.create("Output_MFI_per_sample")
          
          ### 5.1 - MFI SUMMARY FOR WHOLE DATASET # Use aggregate to determine MFI of each column, by CLUSTER
          temp <- data
          colnames(temp)[which(names(temp) == clust.col)] <- "CLUSTER"
          fun.type <- noquote(fun.type)
          
          temp <- temp[-c(annot.col)]
          
          data.MFI <- aggregate(. ~CLUSTER,
                                data = temp, 
                                FUN=fun.type)   # Choose the function (mean, median)
          
          write.csv(data.MFI, "Output_MFI_per_sample/SumTable_MFI.csv") 
          
          ## 5.2 - LOOP TO REPEAT FOR EACH SAMPLE
          for(a in AllSampleNames){
            temp <- data
            colnames(temp)[which(names(temp) == clust.col)] <- "CLUSTER"
            data.subset <- subset(temp, SampleName == a)
            
            data.subset <- data.subset[-c(annot.col)]
            
            data.subset.MFI <- aggregate(. ~ CLUSTER, data = data.subset, FUN=fun.type) 
            write.csv(data.subset.MFI, paste("Output_MFI_per_sample/SumTable_MFI_", a, ".csv", sep="")) 
          }
          
          ## 5.3 - LOOP TO REPEAT FOR EACH GROUP
          if(use.groups == 1){
            for(a in AllGroupNames){
              temp <- data
              colnames(temp)[which(names(temp) == clust.col)] <- "CLUSTER"
              data.subset <- subset(temp, GroupName == a)
              
              data.subset <- data.subset[-c(annot.col)]
              
              data.subset.MFI <- aggregate(. ~ CLUSTER, data = data.subset, FUN=fun.type) 
              write.csv(data.subset.MFI, paste("Output_MFI_per_sample/SumTable_MFI_", a, ".csv", sep="")) 
            }
          }
        }
        
        
        
###################################################### 6. GENERATE MFI PER MARKER SUM TABLES  ######################################################          
        
        
        if(do.MFI.per.marker == 1){
          
          setwd(PrimaryDirectory)
          dir.create("Output_MFI_per_marker")
          
          ### 6.1 - MFI SUMMARY FOR WHOLE DATASET # Use aggregate to determine MFI of each column, by CLUSTER
          
          temp <- data
          
          markers <- names(temp)[-c(annot.col, clust.col.num)]
          colnames(temp)[which(names(temp) == clust.col)] <- "CLUSTER"
          
          num.clusters <- as.matrix(unique(temp["CLUSTER"])) # find number of clusters
          num.clusters <- sort(num.clusters, decreasing = FALSE)
          
          
          ### 6.2 - Cycle through markers
          
          for (i in c(1:length(unique(markers)))){
            #i <- 1
            z <- markers[i]
            temp[[z]]
            list.per.marker <- list()
            
            ### 6.3 - Cycle through samples
            
            for (a in c(AllSampleNames)){
              #a <- "TA147_01_Air_12_PBS"
              data.subset <- subset(temp, SampleName == a)
              data.subset <- data.subset[-c(annot.col)]
              
              #data.subset <- aggregate(. ~ CLUSTER, data = data.subset, FUN=fun.type)  # works without this line, but doesn't work with this line if some samples are missing clusters
              ## add something here to deal with missing cluster value. c(1:length_of_clusters) -- 
              ## if the number of row doesnt not equal the expected number of clusters, add a row for which ever
              
              ## Replacement for aggregate:
              median.per.cluster <- list()
              
              for(c in num.clusters){ 
                r <- data.subset[data.subset["CLUSTER"] == c,] # selects rows that belong to cluster 'x'
                r
                #r <- r[-c(1:11), ] 
                if(fun.type == "mean"){r <- colMeans(r)} # for mean
                if(fun.type == "median"){r <- apply(r, 2, FUN = median)} # for median
                median.per.cluster[[c]] <- r
              }
              
              median.per.cluster
              
              #TEST <- rbindlist(list(median.per.cluster), use.names=TRUE)    # doesn't work, don't use
              TEST <- do.call(rbind, unname(median.per.cluster))
              rownames(TEST) <- TEST[,"CLUSTER"]
              TEST <- TEST[,-c(which( colnames(TEST)=="CLUSTER" ))]
              data.subset <- as.data.frame(TEST)
              
              data.subset[[i]]
              
              ## Have subseted a marker from one sample
              list.per.marker[[a]] <- data.subset[[i]]
            }
            
            lst <- rbindlist(list(list.per.marker)) ## rbind list won't work with samples that are missing some clusters -- tried to use rbind in this case, didn't help
            #lst <- rbind(list(list.per.marker))
            
            lst <- as.data.frame(t(lst))
            lst
            
            write.csv(lst, paste0("Output_MFI_per_marker/SumTable_MFI_per_marker", "_", z, ".csv")) 
          }
        }
        
        
        ### CURRENTLY DOESN'T INCLUDE THE GROUP
        
