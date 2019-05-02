# SampleMerge_v0.2
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
            
            writeCSV <- 1       # 1 = yes, 0 = no
            writeFCS <- 1       # 1 = yes, 0 = no
              
                
###################################################### 2. DATA INPUT and PREPARATION ###################################################### 
                
            ### 2.1 - Read data into workspace
                
                ## Use to list the .fcs or .csv files in the working directory
                list.files(path=PrimaryDirectory, pattern = ".fcs")     # see a list of FCS files
                list.files(path=PrimaryDirectory, pattern = ".csv")     # see a list of CSV files
                
                ## File details
                file.type               <- ".csv"         # Specfiy file type (".csv" or ".fcs") -- readings .fcs files not currently functional
                data.name               <- "DEMO.DATASET"    # a new name for the data - suggest name is sampletype_date_time (e.g. liver_20180203_1400)
                
                ## Check the list of files
                FileNames <- list.files(path=PrimaryDirectory, pattern = file.type) # Generate list of files
                as.matrix(FileNames) # See file names in a list
                
                ## Read data from Files into list of data frames
                DataList=list() # Creates and empty list to start 
                Length_check = list() # creates an empty list
                ColName_check = list() 
                nrow_check = list()
                
                if (file.type == ".csv"){
                  for (File in FileNames) { # Loop to read files into the list
                    tempdata <- fread(File, check.names = FALSE)
                    File <- gsub(".csv", "", File)
                    DataList[[File]] <- tempdata
                  }
                  for(i in c(1:(length(DataList)))){Length_check[[i]] <- length(names(DataList[[i]]))} # creates a list of the number of columns in each sample
                  for(i in c(1:(length(DataList)))){ColName_check[[i]] <- names(DataList[[i]])}
                  name.table <- data.frame(matrix(unlist(ColName_check), nrow = length(DataList), byrow = T))
                  for(i in c(1:(length(DataList)))){nrow_check[[i]] <- nrow(DataList[[i]])}
                }
                
                rm(tempdata)

                
            ### 2.2 - Data review
                
                ## Review the data from the first sample
                head(DataList[[1]])
                
                ## Review number of columns (parameters, or features), and the number of rows (cells) in each sample
                as.matrix(Length_check)    # number of columns in each sample
                as.matrix(nrow_check)      # number of rows in each sample
                
            ### 2.3 - Review column names, and remove troubleshom columns (if required)
                
                ## Review number of columns, and then all column names
                as.matrix(Length_check)
                name.table
                
                ## Remove troublesom columns (if required)
                    ############### ONLY IF REQUIRED ###############
                    ## Remove any troublesome columns (if required)
                    #for (i in c(1:(length(DataList)))) {
                    #  DataList[[i]]$SampleID <- NULL # after the '$', put column name here
                    #}
                    ################################################ 
                
                ## Final check -- ensure the number of columns in each file is consistent
                Length_check = list() # creates an empty list
                for(i in c(1:(length(DataList)))){Length_check[[i]] <- length(names(DataList[[i]]))} # creates a list of the number of columns in each sample
                as.matrix(Length_check) # check that the number of columns in each sample is the same length
                
                ## Final check -- check all column names
                for(i in c(1:(length(DataList)))){ColName_check[[i]] <- names(DataList[[i]])}
                name.table <- data.frame(matrix(unlist(ColName_check), nrow = length(DataList), byrow = T))
                name.table
                
            ### 2.4 - Remove any columns that represent empty channels -- leave in all cellular markers, live/dead, etc
                
                ## Column selection
                  #as.matrix(names(DataList[[1]])) # review the list of columns
                  #col.rmv <- c() # select the columns to remove
                
                ## Column review    
                  #as.matrix(names(DataList[[1]][-c(col.rmv)])) # Check columns to KEEP
                  #as.matrix(names(DataList[[1]][c(col.rmv)])) # Check columns to REMOVE
                
                ## Remove columns and check data
                  #for (i in c(1:(length(DataList)))) {
                  #  DataList[[i]] <- DataList[[i]][-c(col.rmv)]
                  #}
                
                ## Review columns remaining
                  #as.matrix(names(DataList[[1]]))
                
            ### 2.5 - Add sample identifiers (necessary to merge sample)
                
                ## Create a list of 'SampleName' and SampleNo' entries -- 1:n of samples, these will be matched in order
                AllSampleNames <- c(names(DataList))
                AllSampleNames # Character (words)
                
                AllSampleNos <- c(1:(length(DataList)))       ### <-- changed to 4+ to match with sample
                AllSampleNos # Intiger/numeric (numbers)
                
                ## Add 'SampleNo' and SampleName' to each 
                for (i in AllSampleNos) {DataList[[i]]$SampleNo <- i}
                for (a in AllSampleNames) {DataList[[a]]$SampleName <- a} # Inserting names doesn't stop the rest working, just messes with some auto visual stuff
                
                head(DataList[[1]]) # check to see that the sample name and number keywords have been entered on each row
                
                as.matrix(AllSampleNos)
                as.matrix(AllSampleNames)
                
            ### 2.6 - Add 'GroupNo' and 'GroupName' to each (required if generating 'group' files) (can skip, if no groups present in sample)
                
                ## REQUIRED IF YOU NEED TO SEPARATE DIFFERENT TREATMENT GROUPS -- OTHERWISE, CAN IGNORE
                
                ## Create empty lists
                group.names = list()
                group.nums = list()
                
                ## Setup group names for each sample [[1]] for sample 1, [[2]] for sample 2, etc
                group.names[[1]] <- "Mock_D7"                    # name of group 1
                group.names[[2]] <- "WNV_D7"     # name of group 2 (repeat line with [[3]] for a third group, continue for other groups)
                
                group.names # check group names
                
                ## Specify which samples belong to each group
                as.matrix(AllSampleNames)
                
                group.nums[[1]] <- c(1:6)       # samples that belong to group 1
                group.nums[[2]] <- c(7:12)      # samples that belong to group 2 (repeat line with [[3]] for a third group, continue for other groups)
                
                group.nums # check group names
                
                ## Add 'GroupName' and 'GroupNo' keywords to dataset
                num.of.groups <- c(1:length(group.names))
                for(a in num.of.groups){
                  for(i in c(group.nums[[a]])){
                    DataList[[i]]$GroupNo <- a   
                    DataList[[i]]$GroupName <- group.names[[a]]
                  }
                }
                
            ### 2.7 - Merge samples
                
                ## Check column names
                head(DataList[[1]]) # check the entries for one of the samples in group 1
                head(DataList[[7]]) # check the entries for one of the samples in group 2
                
                ## Concatenate into one large data frame
                data <- rbindlist(DataList) 
                data
                
                
############################################################################################################################
###################################################### END USER INPUT ###################################################### 
############################################################################################################################ 
                

###################################################### 3. WRITE FILES ###################################################### 
            
    ### 3.1 - Identify columns you would like to use to separate samples and groups
                
        ## Choose the column that defines the sample NAMES
            as.matrix(names(data))
            samp.col <- "SampleName"
            head(data[samp.col]) # check that this is actually the sample names
            
        ## Choose the column that defines the GROUPS NAMES
            as.matrix(names(data))
            grp.col <- "GroupName"
            head(data[grp.col]) # check that this is actually the group names
         
        ## Create output directory
            setwd(PrimaryDirectory)
            dir.create("Output_SampleMerge", showWarnings = FALSE)
            setwd("Output_SampleMerge")
            OutputDirectory <- getwd()
            OutputDirectory
    
            AllSampleNames <- unique(data$SampleName)
            group.names <- unique(data$GroupName)
        
    ### 3.2 - Export data in single large file -- both .csv and .fcs format

          ## write .csv
            if(writeCSV == 1){
              csv.filename <- paste(paste0(data.name), ".csv", sep = "")
              fwrite(x = data, file = csv.filename, row.names=FALSE)
            }
            
  
          ## write .fcs
            if(writeFCS == 1){
              metadata <- data.frame(name=dimnames(data)[[2]],desc=paste('column',dimnames(data)[[2]],'from dataset'))
              
              # Create FCS file metadata - ranges, min, and max settings
              #metadata$range <- apply(apply(data,2,range),2,diff)
              metadata$minRange <- apply(data,2,min)
              metadata$maxRange <- apply(data,2,max)
              
              data.ff <- new("flowFrame",exprs=as.matrix(data), parameters=AnnotatedDataFrame(metadata)) # in order to create a flow frame, data needs to be read as matrix by exprs
              head(data.ff)
              write.FCS(data.ff, paste0(data.name, ".fcs"))
            }
            
