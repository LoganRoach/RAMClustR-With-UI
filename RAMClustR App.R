#### Load Packages ####
# loads the packages used for the script to load packages
packages_to_install <- c("tidyverse", "ggrepel", "shiny",  "readxl", "RAMClustR")
for (package_name in packages_to_install) {
  if(!requireNamespace(package_name, quietly = TRUE)) {
    install.packages(package_name)
  }
  library(package_name, character.only = TRUE)
}
ui <- fluidPage(
  
  #App title
  headerPanel("Select CSV Files and Parameter Sheet For RAMClustR"),
  h4("**To Work Properly: Be sure the CSV file you are using, the ExpDes.csv file and the associated R files are present in the same folder**"),
  #File Inputs
      h4("Browse for Files"),
      fileInput("file1", "Select a CSV File:", multiple = FALSE, accept = ".csv"),
      fileInput("file2", "Upload a Template Excel File:"),
      h4("Selecting the Checkbox Includes Internal Standards"),
      checkboxInput("checkbox1", "Include Internal Standards"),
      checkboxInput("checkbox2", "Return Summaries of the PsuedoClusters (Top 2 Features in Each Psuedocluster and Filtered PseudoClusters for Dominate Feature)"),
      # action button that submits the user inputs into the Shiny App function
      actionButton("Submit_Button", "Submit")
      
      
    ) 

shinyServer <- function(input, output){
  observeEvent(input$Submit_Button, {
    if (is.null(input$file1)) {
      print("Need a CSV File")
      break
    } 
    else {
      csvfile_name <- input$file1$name
      file_name <- sub("\\.csv$","", csvfile_name)
      
    }
    
    if (is.null(input$file2)) {
      print("Need a Parameter Excel Sheet To Identify Sample/ Blank Names")
      break
    } 
    else {
      parameters_template <- read_excel(input$file2$datapath)
      
    }
    
    ##### FUNCTIONS ####
    #Renaming RT columns to Retention Time
    RAMCLUST_NO_IS <- function(csv) {
      print("RUNNING RAMCLUSTR WITH NO INTERNAL STANDARDS INCLUDED")
      csvfile_name <- input$file1$name
      file_name <- sub("\\.csv$","", csvfile_name)
      input_file <- read_csv(input$file1$datapath)
      #Renaming RT columns to Retention Time
      if ("RT" %in% colnames (input_file)) {
        input_file <- input_file %>% rename(`Retention_Time` = RT)
      }
      if (!"row.ID" %in% colnames (input_file)) {
        input_file <- input_file %>%  
          mutate(`row.ID` = row_number())
      }
      samples <- as.character(parameters_template[1,2])
      blanks <- as.character(parameters_template[2,2])
      # creating a max area column
      maxareas <- input_file %>%
        pivot_longer(matches(blanks), names_to = "blanks", values_to ="blank_vol") %>%
        pivot_longer(matches(samples),names_to = "samples", values_to = "sample_vol") %>%
        group_by(Mass, Retention_Time) %>%
        summarize(blankvol = mean(blank_vol, na.rm = TRUE),
                  maxArea = max(sample_vol, na.rm = TRUE))
      input_file1 <- left_join(input_file, maxareas)
      missed_by_software <- input_file1 %>% 
        pivot_longer(matches(blanks),names_to = "Blanks", values_to ="blanks_vol") %>%
        pivot_longer(matches(samples), 
                     names_to = "Sample", values_to = "sample_features") %>%
        group_by(Mass, Retention_Time, Sample, Blanks, maxArea,) %>% 
        mutate(sample_features = mean(sample_features, na.rm = TRUE),
               blanks_vol = mean(blanks_vol, na.rm = TRUE)) %>% 
        distinct(row.ID, Mass, Retention_Time, Sample, sample_features, Blanks, blanks_vol)
      wide_missed <- missed_by_software %>% 
        pivot_wider(names_from = "Sample", values_from = "sample_features") %>% 
        pivot_wider(names_from = "Blanks", values_from = "blanks_vol") %>% group_by(Mass, Retention_Time,)
      input_file2 <- wide_missed %>% distinct( Mass, Retention_Time, .keep_all = TRUE)
      #Isolate feature annotation information and number rows for merging later
      feature_info <- input_file2 %>%
        select(-matches(samples),`maxArea`, "row.ID")
      #Isolate matrix of abundances and number rows for merging later
      abundance_info <- input_file2 %>%
        select(matches(samples),`maxArea`, "row.ID")
      #Subset abundance matrix to just blank files
      blank_abundance <- abundance_info %>%
        pivot_longer(matches(samples), names_to = "Sample", values_to = "vol") %>%   #switch to long form for efficient filtering
        filter(grepl(blanks, Sample)) %>%                                      #filter to select all blank types
        group_by(row.ID) %>%
        summarize(blank_vol = median(vol))
      #write_csv(blank_abundance, "NO IS BA.csv")
      #Subset CD input to only internal standards for IS correction calculation
      #IS_only <- input_file1 
      
      #IS_abundance <- IS_only %>% 
      # pivot_longer(matches("VP1468"), names_to = "Sample", values_to = "vol") %>%        #pivot to long form for efficient filter/summarize
      # group_by(Sample,Mass, Retention_Time) %>%  summarize(IS_medabun = median(vol)) %>%        #Determine median abundance for each IS per sample
      #ungroup() %>% group_by(Mass) %>% left_join(summarize(.,mean(IS_medabun))) %>% #Determine average of all samples and merge back to original file inline
      #mutate(scalar = IS_medabun / `mean(IS_medabun)`) %>% ungroup() %>%            #Calculate deviation of each IS from grand average per-sample
      #group_by(Sample,Mass, Retention_Time) %>% summarize(adjust = mean(scalar))                 #Calculate adjustment factor per sample based on average deviation of all IS compounds
      
      #Adjust features for IS levels, filter for 5x blank abundance
      filter_frame <- blank_abundance %>% left_join(abundance_info) %>%               #Combine measurement data matrix and blank info
        rename(max_vol = `maxArea`) %>%
        rename(b_v = `blank_vol`) %>% 
        select(-matches(as.character(blanks))) %>% 
        pivot_longer(matches(as.character(samples)), names_to = "Sample", values_to = "vol") %>%
        #left_join(IS_abundance) %>%                                                   #Combine IS measurement data as well
        #mutate(vol = vol / adjust) %>%                                                #Adjust intensity information using IS scalar value
        left_join(feature_info) %>%                                                   #Merge back feature annotation information
        rename(blank_vol = `b_v`) %>% 
        filter(max_vol > 5*blank_vol)
      
      #Try to merge isomers to check variability and filter by feature CV
      #Would need to be rewritten for non-PFAS due to the type of isomers expected
      CV_filter <- filter_frame %>% filter(!grepl(as.character(blanks), Sample)) %>%      #select only samples
        mutate(ifelse(vol < 1e5, 0, vol)) %>%                                         #zero out minor compounds
        mutate(Sample = gsub(perl = TRUE, as.character(samples),"",Sample)) %>%                 #convert filenames to sample names by stripping replicate annotations. This needs to be custom (re)-written until I develop a fixed naming scheme
        group_by(Sample, Mass, Retention_Time) %>% mutate(vol_sum = sum(vol)) %>%                 #add abundance from all major isomers of the same formula
        summarize(sample_vol = median(vol_sum),
                  sample_CV = sd(vol_sum)/sample_vol) %>%                             #calculate median feature abundance and RSD
        mutate(keep = sample_vol > 1000000 & sample_CV < 1) %>%                       #filter if less than abundance threshold or greater than CV threshold
        group_by(Mass, Retention_Time) %>% summarize(drop = sum(keep)) %>% filter(drop == 0)       #keep compounds that pass filter
      #Prepare midpoint data matrix for later subfiltering by RAMClust, enviohomolog etc.
      final_frame <- filter_frame %>% 
        mutate(Sample = gsub(perl = TRUE, samples,"",Sample)) %>%                 #convert column labels to sample names, retaining replicate tags
        filter(!Mass %in% CV_filter$Mass) %>% 
        filter(!Retention_Time %in% CV_filter$Retention_Time) %>%                                             #drop based on CV_filtering logic
        select(Mass, Retention_Time, row.ID, Sample, vol) %>%                                        #select only abundance data and minimal feature info
        group_by(Mass, Retention_Time, Sample) %>%
        summarize(vol = mean(vol),
                  row.ID = min(row.ID)) %>%
        left_join(feature_info) %>%
        mutate(Sample = case_when(grepl("Samlpe",Sample) ~ gsub("Samlpe","Sample",Sample),
                                  TRUE ~ Sample)) %>%
        filter(`Retention_Time` > 0.5) %>%                                                        #remove early eluters ( > 0.5)
        arrange(desc(`maxArea`)) %>% unique() %>% ungroup()                           # arrange by maxArea
      #Prepare "traditional" formatted table of sample x abundance
      wide_frame <- final_frame %>%
        pivot_wider(names_from = Sample, values_from = vol)
      MIDPOINT_ANALYSIS2 <- paste0("Analysis_midpoint_",file_name,".csv")
      MIDPOINT <- write_csv(wide_frame, MIDPOINT_ANALYSIS2, na = "")
      #Convert to csv file for RAMClust, have not rewritten to take direct object inputs
      ramclustR_file <- input_file2 %>%
        select(matches(samples)) %>%
        t() %>%
        as_tibble(.name_repair = make.names)
      #Get column and row names to prepare csv file for RAMClustR input
      ramClustR_colnames <- input_file2 %>%
        ungroup() %>% 
        mutate(colnames = paste0(`Mass`,"_",`Retention_Time`*60)) %>%                     
        select(colnames) %>% as_vector() 
      ramClustR_rownames <- colnames(select(input_file2,matches(samples))) %>%
        gsub(pattern = "Area: ", replacement ="")
      colnames(ramclustR_file) <- ramClustR_colnames
      ramclustR_file["sample"] <- ramClustR_rownames
      write_csv(select(ramclustR_file, sample, contains("_")), "RAMClustR_file.csv")
      ## Run RamClustR
      experiment <- defineExperiment(csv = "ExpDes.csv")                              #experiment info is required input, but is placeholder, not modified
      RC_cluster <- ramclustR(ms = "RAMClustR_file.csv",                                    #using csv input because I am too lazy to generate xcmsObj, features as columns + samples as rows
                              featdelim = "_", 
                              st = 4,                                                        #half-value of chromatographic peak width
                              maxt = 45,                                                     #maximum  shifted for cluster
                              deepSplit = FALSE,
                              ExpDes=experiment, 
                              sampNameCol = 1,
                              mspout = FALSE)
      #Expected Adducts for CAMERA annotation
      ads_list <- c("[M-H]-", "[M+Na-2H]-", "[M+K-2H]-", "[M+CH2O2-H]-", 
                    "[2M-H]-", "[2M+Na-2H]-", "[2M+K-2H]-", "[2M+CH2O2- H]-", 
                    "[3M-H]-", "[3M+Na-2H]-", "[3M+K-2H]-", "[3M+CH2O2- H]-")
      #Expected Neutral Losses for CAMERA annotation
      nls_list <- c("[M-H-NH3]-", "[M-H-H2O]-", "[M-H-COCH2]-", 
                    "[M-H-CO2]-", "[M-H-NH3-CO2]-", "[M-H-NH3-HCOOH]-", 
                    "[M-H-NH3-H2O]-", "[M-H-NH3-COCH2]-", "[M-H-S]-", 
                    "[M-H-S-NH3-HCOOH]-", "[M-H-H4O2]-", "[M-H-CH2]-", 
                    "[M-H-O]-", "[M-H-C2H2]-", "[M-H-C2H4]-", "[M-H-CO]-",
                    "[M-H-CO2CF2]-", "[M-H-CF2]-", "[M-H-HF]-")
      #Tried to assign M-H in pseduoclusters
      RC <- do.findmain(RC_cluster, mode = "negative", mzabs.error = 0.02, ppm.error = 5,
                        ads = ads_list,
                        nls = nls_list,
                        writeMat = FALSE,
                        writeMS = TRUE,)
      #Bind cluster information from RAMClust output back to input data frame
      clustering_order <- tibble(cluster = RC$featclus, order = RC$xcmsOrd) %>%
        arrange(order)
      working_matrix <- cbind(input_file2, clustering_order) 
      working_matrix_Cluster_Order_Long <- working_matrix %>% 
        pivot_longer(matches(samples),names_to = "samples", values_to = "sample_vol") %>% 
        group_by( Mass, ) %>% 
        summarise(MeanArea = mean(sample_vol, na.rm = TRUE))
      working_matrix2 <- left_join(working_matrix, working_matrix_Cluster_Order_Long)
      working_matrix_Clustering_Order <- working_matrix2 %>% 
        group_by(cluster) %>% arrange(desc(MeanArea), .by_group = TRUE) %>% 
        mutate(ClusterOrder = row_number())
      RAMCLUST_FILENAME <- paste0("RamClustR_Output_",file_name,".csv")
      RAMCLUST_OUTPUT <- write_csv(working_matrix_Clustering_Order, RAMCLUST_FILENAME)
      #return(RAMCLUST_OUTPUT)
      #break
      #stopApp()
      print("RAMClustR has finished")
    }
    #
    RAMCLUST_WITH_IS <- function(csv) {
      print("RUNNING RAMCLUSTR WITH INTERNAL STANDARDS INCLUDED")
      csvfile_name <- input$file1$name
      file_name <- sub("\\.csv$","", csvfile_name)
      input_file <- read_csv(input$file1$datapath)
      #Renaming RT columns to Retention Time and adding extra necessary cols
      if ("RT" %in% colnames (input_file)) {
        input_file <- input_file %>% rename(`Retention_Time` = RT)
      }
      if (!"row.ID" %in% colnames (input_file)) {
        input_file <- input_file %>%  
          mutate(`row.ID` = row_number())
      }
      samples <- as.character(parameters_template[1,2])
      blanks <- as.character(parameters_template[2,2])
      # creating a max area column
      maxareas <- input_file %>%
        pivot_longer(matches(blanks), names_to = "blanks", values_to ="blank_vol") %>%
        pivot_longer(matches(samples),names_to = "samples", values_to = "sample_vol") %>%
        group_by(Mass, Retention_Time) %>%
        summarize(blankvol = mean(blank_vol, na.rm = TRUE),
                  maxArea = max(sample_vol, na.rm = TRUE))
      input_file1 <- left_join(input_file, maxareas)
      missed_by_software <- input_file1 %>% 
        pivot_longer(matches(blanks),names_to = "Blanks", values_to ="blanks_vol") %>%
        pivot_longer(matches(samples), 
                     names_to = "Sample", values_to = "sample_features") %>%
        group_by(Mass, Retention_Time, Sample, Blanks, maxArea,) %>% 
        mutate(sample_features = mean(sample_features, na.rm = TRUE),
               blanks_vol = mean(blanks_vol, na.rm = TRUE)) %>% 
        distinct(row.ID, Mass, Retention_Time, Sample, sample_features, Blanks, blanks_vol)
      wide_missed <- missed_by_software %>% 
        pivot_wider(names_from = "Sample", values_from = "sample_features") %>% 
        pivot_wider(names_from = "Blanks", values_from = "blanks_vol") %>% group_by(Mass, Retention_Time,)
      input_file2 <- wide_missed %>% distinct( Mass, Retention_Time, .keep_all = TRUE)
      #Isolate feature annotation information and number rows for merging later
      feature_info <- input_file2 %>%
        select(-matches(samples),`maxArea`, "row.ID")
      #Isolate matrix of abundances and number rows for merging later
      abundance_info <- input_file2 %>%
        select(matches(samples),`maxArea`, "row.ID")
      #Subset abundance matrix to just blank files
      blank_abundance <- abundance_info %>%
        pivot_longer(matches(samples), names_to = "Sample", values_to = "vol") %>%   #switch to long form for efficient filtering
        filter(grepl(blanks, Sample)) %>%                                      #filter to select all blank types
        group_by(row.ID) %>%
        summarize(blank_vol = median(vol))
      #write_csv(blank_abundance, "blank_abundance.csv")
      #Subset CD input to only internal standards for IS correction calculation
      IS_only <- input_file1 
      IS_abundance <- IS_only %>% 
      pivot_longer(matches(samples), names_to = "Sample", values_to = "vol") %>%        #pivot to long form for efficient filter/summarize
      group_by(Sample,Mass, Retention_Time) %>%  summarize(IS_medabun = median(vol)) %>%        #Determine median abundance for each IS per sample
      ungroup() %>% group_by(Mass) %>% left_join(summarize(.,mean(IS_medabun))) %>% #Determine average of all samples and merge back to original file inline
      mutate(scalar = IS_medabun / `mean(IS_medabun)`) %>% ungroup() %>%            #Calculate deviation of each IS from grand average per-sample
      group_by(Sample,Mass, Retention_Time) %>% summarize(adjust = mean(scalar))                 #Calculate adjustment factor per sample based on average deviation of all IS compounds
      #Adjust features for IS levels, filter for 5x blank abundance
      filter_frame <- blank_abundance %>% left_join(abundance_info) %>%               #Combine measurement data matrix and blank info
        rename(max_vol = `maxArea`) %>%
        rename(b_v = `blank_vol`) %>% 
        select(-matches(as.character(blanks))) %>% 
        pivot_longer(matches(as.character(samples)), names_to = "Sample", values_to = "vol") %>%
        #left_join(IS_abundance) %>%                                                   #Combine IS measurement data as well
        #mutate(vol = vol / adjust) %>%                                                #Adjust intensity information using IS scalar value
        left_join(feature_info) %>%                                                   #Merge back feature annotation information
        rename(blank_vol = `b_v`) %>% 
        filter(max_vol > 5*blank_vol)
      #Try to merge isomers to check variability and filter by feature CV
      #Would need to be rewritten for non-PFAS due to the type of isomers expected
      CV_filter <- filter_frame %>% filter(!grepl(blanks, Sample)) %>%      #select only samples
        mutate(ifelse(vol < 1e5, 0, vol)) %>%                                         #zero out minor compounds
        mutate(Sample = gsub(perl = TRUE, samples,"",Sample)) %>%                 #convert filenames to sample names by stripping replicate annotations. This needs to be custom (re)-written until I develop a fixed naming scheme
        group_by(Sample, Mass, Retention_Time) %>% mutate(vol_sum = sum(vol)) %>%                 #add abundance from all major isomers of the same formula
        summarize(sample_vol = median(vol_sum),
                  sample_CV = sd(vol_sum)/sample_vol) %>%                             #calculate median feature abundance and RSD
        mutate(keep = sample_vol > 1000000 & sample_CV < 1) %>%                       #filter if less than abundance threshold or greater than CV threshold
        group_by(Mass, Retention_Time) %>% summarize(drop = sum(keep)) %>% filter(drop == 0)       #keep compounds that pass filter
      #Prepare midpoint data matrix for later subfiltering by RAMClust, enviohomolog etc.
      final_frame <- filter_frame %>% 
        mutate(Sample = gsub(perl = TRUE, samples,"",Sample)) %>%                 #convert column labels to sample names, retaining replicate tags
        filter(!Mass %in% CV_filter$Mass) %>% 
        filter(!Retention_Time %in% CV_filter$Retention_Time) %>%                                             #drop based on CV_filtering logic
        select(Mass, Retention_Time, row.ID, Sample, vol) %>%                                        #select only abundance data and minimal feature info
        group_by(Mass, Retention_Time, Sample) %>%
        summarize(vol = mean(vol),
                  row.ID = min(row.ID)) %>%
        left_join(feature_info) %>%
        mutate(Sample = case_when(grepl("Samlpe",Sample) ~ gsub("Samlpe","Sample",Sample),
                                  TRUE ~ Sample)) %>%
        filter(`Retention_Time` > 0.5) %>%                                                        #remove early eluters ( > 0.5)
        arrange(desc(`maxArea`)) %>% unique() %>% ungroup()                           # arrange by maxArea
      #Prepare "traditional" formatted table of sample x abundance
      wide_frame <- final_frame %>%
        pivot_wider(names_from = Sample, values_from = vol)
      MIDPOINT_ANALYSIS2 <- paste0("Analysis_midpoint_",file_name,".csv")
      MIDPOINT <- write_csv(wide_frame, MIDPOINT_ANALYSIS2, na = "")
      #Convert to csv file for RAMClust, have not rewritten to take direct object inputs
      ramclustR_file <- input_file2 %>%
        select(matches(samples)) %>%
        t() %>%
        as_tibble(.name_repair = make.names)
      #Get column and row names to prepare csv file for RAMClustR input
      ramClustR_colnames <- input_file2 %>%
        ungroup() %>% 
        mutate(colnames = paste0(`Mass`,"_",`Retention_Time`*60)) %>%                     
        select(colnames) %>% as_vector() 
      ramClustR_rownames <- colnames(select(input_file2,matches(samples))) %>%
        gsub(pattern = "Area: ", replacement ="")
      colnames(ramclustR_file) <- ramClustR_colnames
      ramclustR_file["sample"] <- ramClustR_rownames
      write_csv(select(ramclustR_file, sample, contains("_")), "RAMClustR_file.csv")
      ## Run RamClustR
      experiment <- defineExperiment(csv = "ExpDes.csv")                              #experiment info is required input, but is placeholder, not modified
      RC_cluster <- ramclustR(ms = "RAMClustR_file.csv",                                    #using csv input because I am too lazy to generate xcmsObj, features as columns + samples as rows
                              featdelim = "_", 
                              st = 4,                                                        #half-value of chromatographic peak width
                              maxt = 45,                                                     #maximum  shifted for cluster
                              deepSplit = FALSE,
                              ExpDes=experiment, 
                              sampNameCol = 1,
                              mspout = FALSE)
      #Expected Adducts for CAMERA annotation
      ads_list <- c("[M-H]-", "[M+Na-2H]-", "[M+K-2H]-", "[M+CH2O2-H]-", 
                    "[2M-H]-", "[2M+Na-2H]-", "[2M+K-2H]-", "[2M+CH2O2- H]-", 
                    "[3M-H]-", "[3M+Na-2H]-", "[3M+K-2H]-", "[3M+CH2O2- H]-")
      #Expected Neutral Losses for CAMERA annotation
      nls_list <- c("[M-H-NH3]-", "[M-H-H2O]-", "[M-H-COCH2]-", 
                    "[M-H-CO2]-", "[M-H-NH3-CO2]-", "[M-H-NH3-HCOOH]-", 
                    "[M-H-NH3-H2O]-", "[M-H-NH3-COCH2]-", "[M-H-S]-", 
                    "[M-H-S-NH3-HCOOH]-", "[M-H-H4O2]-", "[M-H-CH2]-", 
                    "[M-H-O]-", "[M-H-C2H2]-", "[M-H-C2H4]-", "[M-H-CO]-",
                    "[M-H-CO2CF2]-", "[M-H-CF2]-", "[M-H-HF]-")
      #Tried to assign M-H in pseduoclusters
      RC <- do.findmain(RC_cluster, mode = "negative", mzabs.error = 0.02, ppm.error = 5,
                        ads = ads_list,
                        nls = nls_list,
                        writeMat = FALSE,
                        writeMS = TRUE,)
      #Bind cluster information from RAMClust output back to input data frame
      clustering_order <- tibble(cluster = RC$featclus, order = RC$xcmsOrd) %>%
        arrange(order)
      working_matrix <- cbind(input_file2, clustering_order) 
      working_matrix_Cluster_Order_Long <- working_matrix %>% 
        pivot_longer(matches(samples),names_to = "samples", values_to = "sample_vol") %>% 
        group_by( Mass, ) %>% 
        summarise(MeanArea = mean(sample_vol, na.rm = TRUE))
      working_matrix2 <- left_join(working_matrix, working_matrix_Cluster_Order_Long)
      working_matrix_Clustering_Order <- working_matrix2 %>% 
        group_by(cluster) %>% arrange(desc(MeanArea), .by_group = TRUE) %>% 
        mutate(ClusterOrder = row_number())
      RAMCLUST_FILENAME <- paste0("RamClustR_Output_",file_name,".csv")
      RAMCLUST_OUTPUT <- write_csv(working_matrix_Clustering_Order, RAMCLUST_FILENAME)
      #return(RAMCLUST_OUTPUT)
      #break
      #stopApp()
      print("RAMClustR has finished")
    }
    #
    RAMCLUST_Check2 <- function(csv) {
      print("RUNNING RAMCLUSTR WITH ADDITIONAL SUMMARIES")
      csvfile_name <- input$file1$name
      file_name <- sub("\\.csv$","", csvfile_name)
      input_file <- read_csv(input$file1$datapath)
      #Renaming RT columns to Retention Time and adding extra necessary cols
      if ("RT" %in% colnames (input_file)) {
        input_file <- input_file %>% rename(`Retention_Time` = RT)
      }
      if (!"row.ID" %in% colnames (input_file)) {
        input_file <- input_file %>%  
          mutate(`row.ID` = row_number())
      }
      samples <- as.character(parameters_template[1,2])
      blanks <- as.character(parameters_template[2,2])

      # creating a max area column
      maxareas <- input_file %>%
        pivot_longer(matches(blanks), names_to = "blanks", values_to ="blank_vol") %>%
        pivot_longer(matches(samples),names_to = "samples", values_to = "sample_vol") %>%
        group_by(Mass, Retention_Time) %>%
        summarize(blankvol = mean(blank_vol, na.rm = TRUE),
                  maxArea = max(sample_vol, na.rm = TRUE))
      input_file1 <- left_join(input_file, maxareas)
      missed_by_software <- input_file1 %>% 
        pivot_longer(matches(blanks),names_to = "Blanks", values_to ="blanks_vol") %>%
        pivot_longer(matches(samples), 
                     names_to = "Sample", values_to = "sample_features") %>%
        group_by(Mass, Retention_Time, Sample, Blanks, maxArea,) %>% 
        mutate(sample_features = mean(sample_features, na.rm = TRUE),
               blanks_vol = mean(blanks_vol, na.rm = TRUE)) %>% 
        distinct(row.ID, Mass, Retention_Time, Sample, sample_features, Blanks, blanks_vol)
      wide_missed <- missed_by_software %>% 
        pivot_wider(names_from = "Sample", values_from = "sample_features") %>% 
        pivot_wider(names_from = "Blanks", values_from = "blanks_vol") %>% group_by(Mass, Retention_Time,)
      input_file2 <- wide_missed %>% distinct( Mass, Retention_Time, .keep_all = TRUE)
      #Isolate feature annotation information and number rows for merging later
      feature_info <- input_file2 %>%
        select(-matches(samples),`maxArea`, "row.ID")
      #Isolate matrix of abundances and number rows for merging later
      abundance_info <- input_file2 %>%
        select(matches(samples),`maxArea`, "row.ID")
      #Subset abundance matrix to just blank files
      #print("starting blank_abundance")
      blank_abundance <- abundance_info %>%
        pivot_longer(matches(samples), names_to = "Sample", values_to = "vol") %>%   #switch to long form for efficient filtering
        filter(grepl(blanks, Sample)) %>%                                      #filter to select all blank types
        group_by(row.ID) %>%
        summarize(blank_vol = median(vol))
      #write_csv(blank_abundance,"blank_abundance.csv")
      #Subset CD input to only internal standards for IS correction calculation
      #IS_only <- input_file1 
      
      #IS_abundance <- IS_only %>% 
      # pivot_longer(matches("VP1468"), names_to = "Sample", values_to = "vol") %>%        #pivot to long form for efficient filter/summarize
      # group_by(Sample,Mass, Retention_Time) %>%  summarize(IS_medabun = median(vol)) %>%        #Determine median abundance for each IS per sample
      #ungroup() %>% group_by(Mass) %>% left_join(summarize(.,mean(IS_medabun))) %>% #Determine average of all samples and merge back to original file inline
      #mutate(scalar = IS_medabun / `mean(IS_medabun)`) %>% ungroup() %>%            #Calculate deviation of each IS from grand average per-sample
      #group_by(Sample,Mass, Retention_Time) %>% summarize(adjust = mean(scalar))                 #Calculate adjustment factor per sample based on average deviation of all IS compounds
      
      #Adjust features for IS levels, filter for 5x blank abundance
      filter_frame <- blank_abundance %>% left_join(abundance_info) %>%               #Combine measurement data matrix and blank info
        rename(max_vol = `maxArea`) %>%
        rename(b_v = `blank_vol`) %>% 
        select(-matches(as.character(blanks))) %>% 
        pivot_longer(matches(as.character(samples)), names_to = "Sample", values_to = "vol") %>%
        #left_join(IS_abundance) %>%                                                   #Combine IS measurement data as well
        #mutate(vol = vol / adjust) %>%                                                #Adjust intensity information using IS scalar value
        left_join(feature_info) %>%                                                   #Merge back feature annotation information
        rename(blank_vol = `b_v`) %>% 
        filter(max_vol > 5*blank_vol)
      
      #Try to merge isomers to check variability and filter by feature CV
      #Would need to be rewritten for non-PFAS due to the type of isomers expected
      CV_filter <- filter_frame %>% filter(!grepl(blanks, Sample)) %>%      #select only samples
        mutate(ifelse(vol < 1e5, 0, vol)) %>%                                         #zero out minor compounds
        mutate(Sample = gsub(perl = TRUE, "VP","",Sample)) %>%                 #convert filenames to sample names by stripping replicate annotations. This needs to be custom (re)-written until I develop a fixed naming scheme
        group_by(Sample, Mass, Retention_Time) %>% mutate(vol_sum = sum(vol)) %>%                 #add abundance from all major isomers of the same formula
        summarize(sample_vol = median(vol_sum),
                  sample_CV = sd(vol_sum)/sample_vol) %>%                             #calculate median feature abundance and RSD
        mutate(keep = sample_vol > 1000000 & sample_CV < 1) %>%                       #filter if less than abundance threshold or greater than CV threshold
        group_by(Mass, Retention_Time) %>% summarize(drop = sum(keep)) %>% filter(drop == 0)       #keep compounds that pass filter
      #Prepare midpoint data matrix for later subfiltering by RAMClust, enviohomolog etc.
      final_frame <- filter_frame %>% 
        mutate(Sample = gsub(perl = TRUE, samples,"",Sample)) %>%                 #convert column labels to sample names, retaining replicate tags
        filter(!Mass %in% CV_filter$Mass) %>% 
        filter(!Retention_Time %in% CV_filter$Retention_Time) %>%                                             #drop based on CV_filtering logic
        select(Mass, Retention_Time, row.ID, Sample, vol) %>%                                        #select only abundance data and minimal feature info
        group_by(Mass, Retention_Time, Sample) %>%
        summarize(vol = mean(vol),
                  row.ID = min(row.ID)) %>%
        left_join(feature_info) %>%
        mutate(Sample = case_when(grepl("Samlpe",Sample) ~ gsub("Samlpe","Sample",Sample),
                                  TRUE ~ Sample)) %>%
        filter(`Retention_Time` > 0.5) %>%                                                        #remove early eluters ( > 0.5)
        arrange(desc(`maxArea`)) %>% unique() %>% ungroup()                           # arrange by maxArea
      #Prepare "traditional" formatted table of sample x abundance
      wide_frame <- final_frame %>%
        pivot_wider(names_from = Sample, values_from = vol)
      MIDPOINT_ANALYSIS2 <- paste0("Analysis_midpoint_",file_name,".csv")
      MIDPOINT <- write_csv(wide_frame, MIDPOINT_ANALYSIS2, na = "")
      #Convert to csv file for RAMClust, have not rewritten to take direct object inputs
      ramclustR_file <- input_file2 %>%
        select(matches(samples)) %>%
        t() %>%
        as_tibble(.name_repair = make.names)
      #Get column and row names to prepare csv file for RAMClustR input
      ramClustR_colnames <- input_file2 %>%
        ungroup() %>% 
        mutate(colnames = paste0(`Mass`,"_",`Retention_Time`*60)) %>%                     
        select(colnames) %>% as_vector() 
      ramClustR_rownames <- colnames(select(input_file2,matches(samples))) %>%
        gsub(pattern = "Area: ", replacement ="")
      colnames(ramclustR_file) <- ramClustR_colnames
      ramclustR_file["sample"] <- ramClustR_rownames
      write_csv(select(ramclustR_file, sample, contains("_")), "RAMClustR_file.csv")
      ## Run RamClustR
      experiment <- defineExperiment(csv = "ExpDes.csv")                              #experiment info is required input, but is placeholder, not modified
      RC_cluster <- ramclustR(ms = "RAMClustR_file.csv",                                    #using csv input because I am too lazy to generate xcmsObj, features as columns + samples as rows
                              featdelim = "_", 
                              st = 4,                                                        #half-value of chromatographic peak width
                              maxt = 45,                                                     #maximum  shifted for cluster
                              deepSplit = FALSE,
                              ExpDes=experiment, 
                              sampNameCol = 1,
                              mspout = FALSE)
      #Expected Adducts for CAMERA annotation
      ads_list <- c("[M-H]-", "[M+Na-2H]-", "[M+K-2H]-", "[M+CH2O2-H]-", 
                    "[2M-H]-", "[2M+Na-2H]-", "[2M+K-2H]-", "[2M+CH2O2- H]-", 
                    "[3M-H]-", "[3M+Na-2H]-", "[3M+K-2H]-", "[3M+CH2O2- H]-")
      #Expected Neutral Losses for CAMERA annotation
      nls_list <- c("[M-H-NH3]-", "[M-H-H2O]-", "[M-H-COCH2]-", 
                    "[M-H-CO2]-", "[M-H-NH3-CO2]-", "[M-H-NH3-HCOOH]-", 
                    "[M-H-NH3-H2O]-", "[M-H-NH3-COCH2]-", "[M-H-S]-", 
                    "[M-H-S-NH3-HCOOH]-", "[M-H-H4O2]-", "[M-H-CH2]-", 
                    "[M-H-O]-", "[M-H-C2H2]-", "[M-H-C2H4]-", "[M-H-CO]-",
                    "[M-H-CO2CF2]-", "[M-H-CF2]-", "[M-H-HF]-")
      #Tried to assign M-H in pseduoclusters
      RC <- do.findmain(RC_cluster, mode = "negative", mzabs.error = 0.02, ppm.error = 5,
                        ads = ads_list,
                        nls = nls_list,
                        writeMat = FALSE,
                        writeMS = TRUE,)
      #Bind cluster information from RAMClust output back to input data frame
      clustering_order <- tibble(cluster = RC$featclus, order = RC$xcmsOrd) %>%
        arrange(order)
      working_matrix <- cbind(input_file2, clustering_order) 
      working_matrix_Cluster_Order_Long <- working_matrix %>% 
        pivot_longer(matches(samples),names_to = "samples", values_to = "sample_vol") %>% 
        group_by( Mass, ) %>% 
        summarise(MeanArea = mean(sample_vol, na.rm = TRUE))
      working_matrix2 <- left_join(working_matrix, working_matrix_Cluster_Order_Long)
      working_matrix_Clustering_Order <- working_matrix2 %>% 
        group_by(cluster) %>% arrange(desc(MeanArea), .by_group = TRUE) %>% 
        mutate(ClusterOrder = row_number())
      RAMCLUST_FILENAME <- paste0("RamClustR_Output_",file_name,".csv")
      RAMCLUST_OUTPUT <- write_csv(working_matrix_Clustering_Order, RAMCLUST_FILENAME)
      #For all features not already manually assigned, take top 2 most abundance features from the pseudocluster
      decimate <- working_matrix_Clustering_Order %>%
        filter(cluster != 0) %>%
        group_by(cluster) %>%
        top_n(2, wt = maxArea)
      Top2Summary <- paste("Top_2_Features_",file_name,".csv")
      TOP2 <- write_csv(decimate, Top2Summary)
      #For features where top feature for pseudocluster is more than 30% higher than second highest, use top feature as exemplar for pseudocluster, otherwise keep both
      dominate <- decimate %>%
        summarize(top = max(maxArea)) %>%
        full_join(decimate) %>%
        filter(maxArea/top > 0.7) %>%
        select(-top) %>%
        full_join(filter(working_matrix, cluster == 0|order == 1)) %>% ungroup() %>%
        select(one_of(colnames(decimate)), -order) %>%
        group_by(cluster) %>%
        arrange(desc(maxArea), .by_group = TRUE) %>% ungroup()
      #Store as csv
      TopFeaturesSummary <- paste("Dominate_Feature_Summary_",file_name,".csv")
      TOPFEATS <- write_csv(dominate, TopFeaturesSummary)
      print("RAMClustR has finished")
      #return(RAMCLUST_OUTPUT)
      #break
      #stopApp()
      
    }
    #
    RAM_ALL <- function(csv) {
      print("RUNNING RAMCLUSTR WITH INTERNAL STANDARDS AND ADDTIONAL SUMMARIES")
      csvfile_name <- input$file1$name
      file_name <- sub("\\.csv$","", csvfile_name)
      input_file <- read_csv(input$file1$datapath)
      #Renaming RT columns to Retention Time and adding extra necessary cols
      if ("RT" %in% colnames (input_file)) {
        input_file <- input_file %>% rename(`Retention_Time` = RT)
      }
      if (!"row.ID" %in% colnames (input_file)) {
        input_file <- input_file %>%  
          mutate(`row.ID` = row_number())
      }
      samples <- as.character(parameters_template[1,2])
      blanks <- as.character(parameters_template[2,2])
      # creating a max area column
      maxareas <- input_file %>%
        pivot_longer(matches(blanks), names_to = "blanks", values_to ="blank_vol") %>%
        pivot_longer(matches(samples),names_to = "samples", values_to = "sample_vol") %>%
        group_by(Mass, Retention_Time) %>%
        summarize(blankvol = mean(blank_vol, na.rm = TRUE),
                  maxArea = max(sample_vol, na.rm = TRUE))
      input_file1 <- left_join(input_file, maxareas)
      missed_by_software <- input_file1 %>% 
        pivot_longer(matches(blanks),names_to = "Blanks", values_to ="blanks_vol") %>%
        pivot_longer(matches(samples), 
                     names_to = "Sample", values_to = "sample_features") %>%
        group_by(Mass, Retention_Time, Sample, Blanks, maxArea,) %>% 
        mutate(sample_features = mean(sample_features, na.rm = TRUE),
               blanks_vol = mean(blanks_vol, na.rm = TRUE)) %>% 
        distinct(row.ID, Mass, Retention_Time, Sample, sample_features, Blanks, blanks_vol)
      wide_missed <- missed_by_software %>% 
        pivot_wider(names_from = "Sample", values_from = "sample_features") %>% 
        pivot_wider(names_from = "Blanks", values_from = "blanks_vol") %>% group_by(Mass, Retention_Time,)
      input_file2 <- wide_missed %>% distinct( Mass, Retention_Time, .keep_all = TRUE)
      #Isolate feature annotation information and number rows for merging later
      feature_info <- input_file2 %>%
        select(-matches(samples),`maxArea`, "row.ID")
      #Isolate matrix of abundances and number rows for merging later
      abundance_info <- input_file2 %>%
        select(matches(samples),`maxArea`, "row.ID")
      #Subset abundance matrix to just blank files
      blank_abundance <- abundance_info %>%
        pivot_longer(matches(samples), names_to = "Sample", values_to = "vol") %>%   #switch to long form for efficient filtering
        filter(grepl(blanks, Sample)) %>%                                      #filter to select all blank types
        group_by(row.ID) %>%
        summarize(blank_vol = median(vol))
      #Subset CD input to only internal standards for IS correction calculation
      IS_only <- input_file1 
      IS_abundance <- IS_only %>% 
        pivot_longer(matches(samples), names_to = "Sample", values_to = "vol") %>%        #pivot to long form for efficient filter/summarize
        group_by(Sample,Mass, Retention_Time) %>%  summarize(IS_medabun = median(vol)) %>%        #Determine median abundance for each IS per sample
        ungroup() %>% group_by(Mass) %>% left_join(summarize(.,mean(IS_medabun))) %>% #Determine average of all samples and merge back to original file inline
        mutate(scalar = IS_medabun / `mean(IS_medabun)`) %>% ungroup() %>%            #Calculate deviation of each IS from grand average per-sample
        group_by(Sample,Mass, Retention_Time) %>% summarize(adjust = mean(scalar))                 #Calculate adjustment factor per sample based on average deviation of all IS compounds
      #Adjust features for IS levels, filter for 5x blank abundance
      filter_frame <- blank_abundance %>% left_join(abundance_info) %>%               #Combine measurement data matrix and blank info
        rename(max_vol = `maxArea`) %>%
        rename(b_v = `blank_vol`) %>% 
        select(-matches(as.character(blanks))) %>% 
        pivot_longer(matches(as.character(samples)), names_to = "Sample", values_to = "vol") %>%
        #left_join(IS_abundance) %>%                                                   #Combine IS measurement data as well
        #mutate(vol = vol / adjust) %>%                                                #Adjust intensity information using IS scalar value
        left_join(feature_info) %>%                                                   #Merge back feature annotation information
        rename(blank_vol = `b_v`) %>% 
        filter(max_vol > 5*blank_vol)
      #Try to merge isomers to check variability and filter by feature CV
      #Would need to be rewritten for non-PFAS due to the type of isomers expected
      CV_filter <- filter_frame %>% filter(!grepl(blanks, Sample)) %>%      #select only samples
        mutate(ifelse(vol < 1e5, 0, vol)) %>%                                         #zero out minor compounds
        mutate(Sample = gsub(perl = TRUE, samples,"",Sample)) %>%                 #convert filenames to sample names by stripping replicate annotations. This needs to be custom (re)-written until I develop a fixed naming scheme
        group_by(Sample, Mass, Retention_Time) %>% mutate(vol_sum = sum(vol)) %>%                 #add abundance from all major isomers of the same formula
        summarize(sample_vol = median(vol_sum),
                  sample_CV = sd(vol_sum)/sample_vol) %>%                             #calculate median feature abundance and RSD
        mutate(keep = sample_vol > 1000000 & sample_CV < 1) %>%                       #filter if less than abundance threshold or greater than CV threshold
        group_by(Mass, Retention_Time) %>% summarize(drop = sum(keep)) %>% filter(drop == 0)       #keep compounds that pass filter
      #Prepare midpoint data matrix for later subfiltering by RAMClust, enviohomolog etc.
      final_frame <- filter_frame %>% 
        mutate(Sample = gsub(perl = TRUE, samples,"",Sample)) %>%                 #convert column labels to sample names, retaining replicate tags
        filter(!Mass %in% CV_filter$Mass) %>% 
        filter(!Retention_Time %in% CV_filter$Retention_Time) %>%                                             #drop based on CV_filtering logic
        select(Mass, Retention_Time, row.ID, Sample, vol) %>%                                        #select only abundance data and minimal feature info
        group_by(Mass, Retention_Time, Sample) %>%
        summarize(vol = mean(vol),
                  row.ID = min(row.ID)) %>%
        left_join(feature_info) %>%
        mutate(Sample = case_when(grepl("Samlpe",Sample) ~ gsub("Samlpe","Sample",Sample),
                                  TRUE ~ Sample)) %>%
        filter(`Retention_Time` > 0.5) %>%                                                        #remove early eluters ( > 0.5)
        arrange(desc(`maxArea`)) %>% unique() %>% ungroup()                           # arrange by maxArea
      #Prepare "traditional" formatted table of sample x abundance
      wide_frame <- final_frame %>%
        pivot_wider(names_from = Sample, values_from = vol)
      MIDPOINT_ANALYSIS2 <- paste0("Analysis_midpoint_",file_name,".csv")
      MIDPOINT <- write_csv(wide_frame, MIDPOINT_ANALYSIS2, na = "")
      #Convert to csv file for RAMClust, have not rewritten to take direct object inputs
      ramclustR_file <- input_file2 %>%
        select(matches(samples)) %>%
        t() %>%
        as_tibble(.name_repair = make.names)
      #Get column and row names to prepare csv file for RAMClustR input
      ramClustR_colnames <- input_file2 %>%
        ungroup() %>% 
        mutate(colnames = paste0(`Mass`,"_",`Retention_Time`*60)) %>%                     
        select(colnames) %>% as_vector() 
      ramClustR_rownames <- colnames(select(input_file2,matches(samples))) %>%
        gsub(pattern = "Area: ", replacement ="")
      colnames(ramclustR_file) <- ramClustR_colnames
      ramclustR_file["sample"] <- ramClustR_rownames
      write_csv(select(ramclustR_file, sample, contains("_")), "RAMClustR_file.csv")
      ## Run RamClustR
      experiment <- defineExperiment(csv = "ExpDes.csv")                              #experiment info is required input, but is placeholder, not modified
      RC_cluster <- ramclustR(ms = "RAMClustR_file.csv",                                    #using csv input because I am too lazy to generate xcmsObj, features as columns + samples as rows
                              featdelim = "_", 
                              st = 4,                                                        #half-value of chromatographic peak width
                              maxt = 45,                                                     #maximum  shifted for cluster
                              deepSplit = FALSE,
                              ExpDes=experiment, 
                              sampNameCol = 1,
                              mspout = FALSE)
      #Expected Adducts for CAMERA annotation
      ads_list <- c("[M-H]-", "[M+Na-2H]-", "[M+K-2H]-", "[M+CH2O2-H]-", 
                    "[2M-H]-", "[2M+Na-2H]-", "[2M+K-2H]-", "[2M+CH2O2- H]-", 
                    "[3M-H]-", "[3M+Na-2H]-", "[3M+K-2H]-", "[3M+CH2O2- H]-")
      #Expected Neutral Losses for CAMERA annotation
      nls_list <- c("[M-H-NH3]-", "[M-H-H2O]-", "[M-H-COCH2]-", 
                    "[M-H-CO2]-", "[M-H-NH3-CO2]-", "[M-H-NH3-HCOOH]-", 
                    "[M-H-NH3-H2O]-", "[M-H-NH3-COCH2]-", "[M-H-S]-", 
                    "[M-H-S-NH3-HCOOH]-", "[M-H-H4O2]-", "[M-H-CH2]-", 
                    "[M-H-O]-", "[M-H-C2H2]-", "[M-H-C2H4]-", "[M-H-CO]-",
                    "[M-H-CO2CF2]-", "[M-H-CF2]-", "[M-H-HF]-")
      #Tried to assign M-H in pseduoclusters
      RC <- do.findmain(RC_cluster, mode = "negative", mzabs.error = 0.02, ppm.error = 5,
                        ads = ads_list,
                        nls = nls_list,
                        writeMat = FALSE,
                        writeMS = TRUE,)
      #Bind cluster information from RAMClust output back to input data frame
      clustering_order <- tibble(cluster = RC$featclus, order = RC$xcmsOrd) %>%
        arrange(order)
      working_matrix <- cbind(input_file2, clustering_order) 
      working_matrix_Cluster_Order_Long <- working_matrix %>% 
        pivot_longer(matches(samples),names_to = "samples", values_to = "sample_vol") %>% 
        group_by( Mass, ) %>% 
        summarise(MeanArea = mean(sample_vol, na.rm = TRUE))
      working_matrix2 <- left_join(working_matrix, working_matrix_Cluster_Order_Long)
      working_matrix_Clustering_Order <- working_matrix2 %>% 
        group_by(cluster) %>% arrange(desc(MeanArea), .by_group = TRUE) %>% 
        mutate(ClusterOrder = row_number())
      RAMCLUST_FILENAME <- paste0("RamClustR_Output_",file_name,".csv")
      RAMCLUST_OUTPUT <- write_csv(working_matrix_Clustering_Order, RAMCLUST_FILENAME)
      #For all features not already manually assigned, take top 2 most abundance features from the pseudocluster
      decimate <- working_matrix_Clustering_Order %>%
        filter(cluster != 0) %>%
        group_by(cluster) %>%
        top_n(2, wt = maxArea)
      Top2Summary <- paste("Top_2_Features_",file_name,".csv")
      TOP2 <- write_csv(decimate, Top2Summary)
      #For features where top feature for pseudocluster is more than 30% higher than second highest, use top feature as exemplar for pseudocluster, otherwise keep both
      dominate <- decimate %>%
        summarize(top = max(maxArea)) %>%
        full_join(decimate) %>%
        filter(maxArea/top > 0.7) %>%
        select(-top) %>%
        full_join(filter(working_matrix, cluster == 0|order == 1)) %>% ungroup() %>%
        select(one_of(colnames(decimate)), -order) %>%
        group_by(cluster) %>%
        arrange(desc(maxArea), .by_group = TRUE) %>% ungroup()
      #Store as csv
      TopFeaturesSummary <- paste("Dominate_Feature_Summary_",file_name,".csv")
      TOPFEATS <- write_csv(dominate, TopFeaturesSummary)
      #return(RAMCLUST_OUTPUT)
      #break
      #stopApp()
      print("RAMClustR has finished")
    }
    #### REST OF THE SCRIPT ####
    
    ##### if checkbox 1 is checked, turn on Include Internal Standards, else no IS included #####
    if (input$checkbox1) {
      csv <- input$file1
      suppressWarnings({
        RAMCLUST_WITH_IS(csv)
        })
      stopApp()
    } else if (input$checkbox2) {
      csv <- input$file1
      suppressWarnings({
        RAMCLUST_Check2(csv)
      })
      stopApp()
    }  else if (input$checkbox1 && input$checkbox2) {
      csv <- input$file1
      suppressWarnings({
        RAM_ALL(csv)
      })
      stopApp()
    } else  {
      csv <- input$file1
      suppressWarnings({
        RAMCLUST_NO_IS(csv)
      })
      stopApp()
    }
    
  })
 
}

shinyApp(ui,shinyServer)

#stopApp()

