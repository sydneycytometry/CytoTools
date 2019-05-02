# SampleSplit_v0.3
# Thomas Ashhurst
# Sydney Cytometry Facility

###################################################### 1. INSTALL AND-OR LOAD PACKAGES ###################################################### 
    
    ### 1.1 - Install packages (if not already installed) from CRAN
            
            ## From CRAN (required)
                if(!require('plyr')) {install.packages('plyr')}
                if(!require('data.table')) {install.packages('data.table')}
                if(!require('rstudioapi')) {install.packages('rstudioapi')}
                if(!require('devtools')){install.packages("devtools")}
            
            ## From Bioconductor (required)
                if(!require('flowViz')) {source("https://bioconductor.org/biocLite.R") 
                  biocLite('flowViz')}
                if(!require('flowCore')) {source("https://bioconductor.org/biocLite.R") 
                  biocLite('flowCore')}
                if(!require('Biobase')) {source("https://bioconductor.org/biocLite.R") 
                  biocLite('Biobase')}
            
    ### 1.2 Load packages       
                library('plyr')
                library('data.table')
                library('rstudioapi')
                library('devtools')
                library('flowCore') 
                library('Biobase')
                library('flowViz')

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
            
            write.merged.file     <- 1    # 1 = yes, 0 = no   
            write.sep.files       <- 1    # 1 = yes, 0 = no
            write.group.files     <- 1    # 1 = yes, 0 = no
            
            writeCSV              <- 1    # 1 = yes, 0 = no
            writeFCS              <- 1    # 1 = yes, 0 = no
                
###################################################### 2. DATA INPUT and PREPARATION ###################################################### 
         
    ### 2.1 - Read data into workspace
            
            ## Use to list the .fcs or .csv files in the working directory
            ## IMPORTANT! The only FCS/CSV file in the directory should be the one desired for analysis. If more than one are found, only the first file will be used
                list.files(path=PrimaryDirectory, pattern = ".fcs")     # see a list of FCS files
                list.files(path=PrimaryDirectory, pattern = ".csv")     # see a list of CSV files
                
            ## File details
                file.type               <- ".csv"         # Specfiy file type (".csv" or ".fcs") -- readings .fcs files not currently functional
                data.name               <- "DEMO.DATASET"    # a new name for the data - suggest name is sampletype_date_time (e.g. liver_20180203_1400)
                
            ## Check the list of files
                FileNames <- list.files(path=PrimaryDirectory, pattern = file.type) # Generate list of files
                as.matrix(FileNames) # See file names in a list

            ## Read data from Files into list of data frames
                data <- fread(FileNames[1], check.names = FALSE)
                data <- as.data.frame(data)
                dim(data)   # Review dimensionality of the data (cells and parameters)
                head(data)  # Review the first 6 rows of the data
                
    ### 2.2 - Identify columns you would like to use to separate samples and groups
                
            ## Choose the column that defines the sample NAMES
                as.matrix(names(data))
                samp.col <- "SampleName"
                head(data[samp.col]) # check that this is actually the sample names
                
            ## Choose the column that defines the GROUPS NAMES
                as.matrix(names(data))
                grp.col <- "GroupName"
                head(data[grp.col]) # check that this is actually the group names
                
                
                
###################################################### 3. USER INPUT - CLUSTERING AND DIMENSIONALITY REDUCTION PREFERENCES ###################################################### 
            

############################################################################################################################
###################################################### END USER INPUT ###################################################### 
############################################################################################################################ 

    ## Create output directory
        setwd(PrimaryDirectory)
        dir.create("Output_splitfiles", showWarnings = FALSE)
        setwd("Output_splitfiles")
        OutputDirectory <- getwd()
        OutputDirectory

        AllSampleNames <- unique(data$SampleName)
        group.names <- unique(data$GroupName)
        
        
    ### 4.1 - Export data in single large file -- both .csv and .fcs format
       
        if(write.merged.file == 1){
        
          ## write .csv
          if(writeCSV == 1){
            csv.filename <- paste(paste0(data.name), "_ALL.csv", sep = "")
            fwrite(x = data, file = csv.filename, row.names=FALSE)
          }
  
          ## write .fcs
          if(writeFCS == 1){
            metadata <- data.frame(name=dimnames(data)[[2]],desc=paste('column',dimnames(data)[[2]],'from dataset'))
          
            ## Create FCS file metadata - ranges, min, and max settings
            #metadata$range <- apply(apply(data,2,range),2,diff)
            metadata$minRange <- apply(data,2,min)
            metadata$maxRange <- apply(data,2,max)
          
            data.ff <- new("flowFrame",exprs=as.matrix(data), parameters=AnnotatedDataFrame(metadata)) # in order to create a flow frame, data needs to be read as matrix by exprs
            head(data.ff)
            write.FCS(data.ff, paste0(data.name, "_ALL.fcs"))
          }
        }
        
        
    ### 4.2 - Export data in individual files -- both .csv and .fcs format

        if(write.sep.files == 1){
          
          for (a in AllSampleNames) {
            
            data_subset <- subset(data, SampleName == a)
            dim(data_subset)
            
            ## write .csv
            if(writeCSV == 1){
              fwrite(data_subset, file = paste(data.name, "_", a,".csv", sep = ""), row.names=FALSE)
            }
            
            ## write .fcs
            if(writeFCS == 1){
              metadata <- data.frame(name=dimnames(data_subset)[[2]],desc=paste('column',dimnames(data_subset)[[2]],'from dataset'))
              
              ## Create FCS file metadata - ranges, min, and max settings
              #metadata$range <- apply(apply(data_subset,2,range),2,diff)
              metadata$minRange <- apply(data_subset,2,min)
              metadata$maxRange <- apply(data_subset,2,max)
              
              data_subset.ff <- new("flowFrame",exprs=as.matrix(data_subset), parameters=AnnotatedDataFrame(metadata)) # in order to create a flow frame, data needs to be read as matrix by exprs
              head(data_subset.ff)
              write.FCS(data_subset.ff, paste0(data.name, "_", a, ".fcs"))
            }
          }
        }
         
    ### 4.3 - Export data as grouped files -- both .csv and .fcs format
                  
        if(write.group.files == 1){
          
          for(a in group.names){
            data_subset <- subset(data, GroupName == a)
            dim(data_subset)
            
            ## write .csv
            if(writeCSV == 1){
              fwrite(data_subset, file = paste(data.name, "_", a,"_GROUP.csv", sep = ""), row.names=FALSE)
            }
            
            ## write .fcs
            if(writeFCS == 1){
              metadata <- data.frame(name=dimnames(data_subset)[[2]],desc=paste('column',dimnames(data_subset)[[2]],'from dataset'))
              metadata
              
              ## Create FCS file metadata - ranges, min, and max settings
              #metadata$range <- apply(apply(data_subset,2,range),2,diff)
              metadata$minRange <- apply(data_subset,2,min)
              metadata$maxRange <- apply(data_subset,2,max)
              
              data_subset.ff <- new("flowFrame",exprs=as.matrix(data_subset), parameters=AnnotatedDataFrame(metadata)) # in order to create a flow frame, data needs to be read as matrix by exprs
              head(data_subset.ff)
              write.FCS(data_subset.ff, paste0(data.name,"_",a,"_GROUP.fcs"))
            }
          }
        }

        