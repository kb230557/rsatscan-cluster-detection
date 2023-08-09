library(dplyr)
library(readxl)
library(rgdal)
library(rsatscan)
library(readr)
library(lubridate)
library(rJava)
library(mailR)
library(keyring)
library(rmarkdown)
library(leaflet)
library(kableExtra)
library(purrr)
library(stringr)


####IMPORTANT NOTES:
#If the case history file is ever manipulated directly in Excel, be sure to save the RunDate column as a Cusom Format type 'm/d/yyyy' before closing
#Same as above for the cluster history file and the StartDate and EndDate columns


#Set email settings - only needs to be done once per computer
#key_set("satscan_sender")  ##Email address to appear in 'From' line of emails
#key_set("email_host")      ##Email smtp host name
#key_set("email")           ##Email smtp email address 
#key_set("email_pw")        ##Eamil smtp password


####QUESTIONS TO ANSWER:
#How to handle Salmonella serogroup?
#Keep or chuck addresses only verified to zip code centroid - currently chucking (bias?)


#=================SCRIPT SET-UP================#

#Setting working directory
setwd('S:/SatScan Project/Exact Locations Program')

#Setting options for Markdown reports
Sys.setenv(RSTUDIO_PANDOC="C:/Program Files/RStudio/bin/pandoc")



#=================PREPARING FILES================#

#TEST DATA SET - Business Objects Report
# events <- read_csv("S:/SatScan Project/Fake Data Set_Salmonella.csv", col_types = "ccccccccciccccccddccdcccc") %>%
#   mutate(`Onset Date` = mdy(str_remove(`Onset Date`, ' 0:00')),
#          `Event Date` = mdy(str_remove(`Event Date`, ' 0:00')))

#Importing Business Objects Report
events <- read_csv("C:/Users/rishi.kowalski/Data/satscan-data/SatScanData080723.csv", col_types = "ccccccccciccccccddccdcccc") %>%
              mutate(
                #`Onset Date` = ymd(str_remove(`Onset Date`, '[[:space:]].*')),
                #`Event Date` = ymd(str_remove(`Event Date`, '[[:space:]].*')))
                `Onset Date` = mdy(`Onset Date`),
                `Event Date` = mdy(`Event Date`)) #modify as BO changed date format 8/4/23


#Replacing spaces in column names
names(events) <- gsub(" |/","_",names(events))

#Combining disease types for legionella and GAS and cleaning disease names
events <- mutate(events, 
            Disease = case_when(
              grepl('Legionellosis', events$Disease) ~ 'Legionellosis',
              grepl('Streptococcal', events$Disease) ~ 'Group A Strep',
              grepl('^Shiga.*not cultured', events$Disease) ~ "STEC Not Cultured",
              grepl('^Shiga.*non-O157', events$Disease) ~ "STEC Non-O157",
              grepl('^Shiga.*O157:H7', events$Disease) ~ "STEC O157:H7",
              grepl("Varicella", events$Disease) ~ "Varicella",
              TRUE ~ Disease
              ))

#Keep only events with verified addresses
eventsVer <- filter(events, `Address_Verified_(Home)` == "YES")

#Adding day of week indicator
eventsVer$DOW <- wday(eventsVer$Event_Date)

#Importing disease parameter file
parameters <- read_tsv('SatScan Disease Parameters.csv', col_types = "cci")


###################COMMENT OUT WHEN DONE TESTING#################################
# eventsFilter <- filter(eventsVer, Disease == "Pertussis")
# paramFilter <- filter(parameters, Disease == "Pertussis")
# disease_run_name <- "Pertussis"
###################COMMENT OUT WHEN DONE TESTING#################################

#=================STARTING FUNCTION TO ITERATE THROUGH DISEASES================#

run_satscan <- function(disease, serogroup) {
  
  #Creating a disease name for incorporation in history files
  disease_run_name <- ifelse(!missing(serogroup), paste(disease, serogroup, sep = "-"), disease)
  
  #Filtering data set for disease specified
  eventsFilter <- filter(eventsVer, Disease == disease)
  
  if(!missing(serogroup)) {                             #====Added code to filter on serogroup but will start/end date need to be adjusted?====#
    eventsFilter <- filter(eventsFilter, Serogroup == serogroup)
  }
  
  #Exiting function early if no cases reported in last year
  if(nrow(eventsFilter) == 0) {
    sink(file = 'satscan_function_log.txt', append = T)
    cat('Zero cases reported for', disease_run_name, '\n')
    sink()
    
    #stop('Zero cases reported')
    #break
    return()
    
  }
  
  #Filtering paramaters for disease specified
  paramFilter <- filter(parameters, Disease == disease)
  
  
#=================RUNNING SATSCAN================#
  
#Creating cases file for SatScan
casfile <- eventsFilter %>%
            mutate(Count = 1) %>%
            select(EventID = State_Case_Number,
                   Count,
                   Event_Date,
                   DOW)
  
#Creating coordinates file for SatScan
geofile <- eventsFilter %>%
            select(EventID = State_Case_Number,
                   Latitude = `Latitude_(Home)`,
                   Longitude = `Longitude_(Home)`)

#Reset the options for the parameter file
invisible(ss.options(reset=TRUE))

#Set options related to the input tab 
ss.options(list(CaseFile = "SessionR.cas", 
                PrecisionCaseTimes = 3,
                StartDate = format(Sys.Date()-366, "%Y/%m/%d"), #366 initially
                EndDate = format(Sys.Date()-1, "%Y/%m/%d"),
                CoordinatesFile = "SessionR.geo",
                CoordinatesType = 1
                ))

#Set options related to the analysis tab 
ss.options(list(AnalysisType = 4,
                ModelType = 2,
                ScanAreas = 1,
                TimeAggregationUnits = 3,
                TimeAggregationLength = 1
))

#Set options related to the output tab 
ss.options(list(ResultsFile = "S:/SatScan Project/ScanResult.txt",  #doesn't appear to be writing
                OutputGoogleEarthKML = "n",
                OutputShapefiles = "y",
                MostLikelyClusterEachCentroidASCII = "y",
                MostLikelyClusterEachCentroidDBase = "y",
                MostLikelyClusterCaseInfoEachCentroidASCII = "n",
                MostLikelyClusterCaseInfoEachCentroidDBase = "n",
                CensusAreasReportedClustersASCII = "y",
                CensusAreasReportedClustersDBase = "y",
                IncludeRelativeRisksCensusAreasASCII = "n",
                IncludeRelativeRisksCensusAreasDBase = "n",
                SaveSimLLRsASCII = "n",
                SaveSimLLRsDBase = "n"
))

#Set advanced options 
ss.options(list(StudyPeriodCheckType=1, #relax geography bounds per cop
                GeographicalCoordinatesCheckType = 1,
                UseDistanceFromCenterOption = "n",
                MaxSpatialSizeInDistanceFromCenter = 1,
                NonCompactnessPenalty = 1,
                MaxTemporalSize = as.integer(paramFilter$Max_Size),
                MaxTemporalSizeInterpretation = 1,  
                ProspectiveStartDate = format(Sys.Date()-366, "%Y/%m/%d"),  #Not sure if this parameter is correct?
                CriteriaForReportingSecondaryClusters = 1,
                LogRunToHistoryFile = "y",
                Version = "9.6.1"  #actually running version 9.6 but when correct and enabled, script fails? Works when not specified. - need to specify last digit (rk 4/18/23)
))

#Write files
# workdir <- 'S:\\SatScan Project\\Aggregated Program\\WorkDir'
# write.ss.prm(workdir, "SessionR")
# write.cas(as.data.frame(casfile), workdir, "SessionR")
# write.geo(geofile, workdir, "SessionR")
# 
# #Run SatScan and store results
# session <- satscan(workdir, "SessionR", sslocation = "C:/Program Files (x86)/SaTScan", cleanup = FALSE, verbose = TRUE)

#Write files
td <- tempdir()
write.ss.prm(td, "SessionR")
write.cas(as.data.frame(casfile), td, "SessionR")
write.geo(as.data.frame(geofile), td, "SessionR")

print("prm written")

#Run SatScan and store results
session <- try(satscan(td, "SessionR", sslocation = "C:/SaTScan", cleanup = FALSE, verbose = TRUE), silent = T)

if (class(session) == "try-error"){
  sink(file = 'satscan_function_log.txt', append = T)
  print(paste0('No results returned for ', disease_run_name, ', run: ', Sys.time(), '\n'))
  sink()
  
  write(paste0(disease_run_name, "_", Sys.time(), " - ", session, "\n"), file = 'try_error_log.txt', append = T)
  
  return()
}

#=================PACKAGING SATSCAN RESULTS================#

#If/Else determines if a significant cluster has been detected and executes different steps depending on answer
if (max(session$col$RECURR_INT) < 100) {             ####REDUCE NUMBER IF TESTING OUTPUT/LOOPS
  
  #Add result to log file
  sink(file = 'satscan_function_log.txt', append = T)
  print(paste0('No significant cluster for ', disease_run_name, ', run: ', Sys.time(), '\n'))
  sink()
  
  #Creating log of most likely cluster to add to history file       
  cluster_temp <- as.data.frame(session$col) %>% 
    filter(CLUSTER == 1) %>%  
    mutate(ClusterID = paste(Sys.Date(), disease_run_name, "1", sep = "-"),
           StartDate = ymd(START_DATE),
           EndDate = ymd(END_DATE),
           Disease = disease_run_name) %>%   
    select(ClusterID, Disease, LATITUDE, LONGITUDE, RADIUS, StartDate, EndDate, P_VALUE, RECURR_INT, OBSERVED, EXPECTED, ODE)
  
  #Identifying case list of most likely cluster (current and past - filtered to just current involved in next step)
  all_involved <- session$gis %>% filter(CLUSTER == 1) %>% select(LOC_ID)
  
  #Saving cases from involved tracts to cluster log
  cluster_temp$CaseList <- eventsFilter %>% 
    filter(State_Case_Number %in% all_involved$LOC_ID & Event_Date >= cluster_temp$StartDate & Event_Date <= cluster_temp$EndDate) %>% 
    select(State_Case_Number) %>% 
    toString()
  
  #Adding cluster to history file
  write_csv(cluster_temp, 'Cluster History File.csv', append = TRUE)

  
}  else {
  
  
  #Determine number of significant clusters
  cluster_num <- sum(session$col$RECURR_INT > 99)   ####REDUCE NUMBER IF TESTING OUTPUT/LOOPS
  
  #Logging each cluster to history file
  for (i in 1:cluster_num) {
    
    #Creating log of cluster to add to history file      
    cluster_temp <- as.data.frame(session$col) %>% 
      filter(CLUSTER == i) %>%  
      mutate(ClusterID = paste(Sys.Date(), disease_run_name, i, sep = "-"),
             StartDate = ymd(START_DATE),
             EndDate = ymd(END_DATE),
             Disease = disease_run_name) %>%   
      select(ClusterID, Disease, LATITUDE, LONGITUDE, RADIUS, StartDate, EndDate, P_VALUE, RECURR_INT, OBSERVED, EXPECTED, ODE)
    
    #Identifying case list of most likely cluster (current and past - filtered to just current involved in next step)
    all_involved <- session$gis %>% filter(CLUSTER == 1) %>% select(LOC_ID)
    
    #Identifying just current cases involved in cluster
    cases_involved <- eventsFilter %>% 
      filter(State_Case_Number %in% all_involved$LOC_ID & Event_Date >= cluster_temp$StartDate & Event_Date <= cluster_temp$EndDate) %>% 
      select(State_Case_Number) 
    
    #Saving cases from involved tracts to cluster log
    cluster_temp$CaseList <- toString(cases_involved)
    
    #Adding cluster to history file
    write_csv(cluster_temp, 'Cluster History File.csv', append = TRUE)
    
    #Creating log for case history file
    cases_involved <- cases_involved %>% 
      mutate(ClusterID_FirstIdentified = paste(Sys.Date(), disease_run_name, i, sep = "-"),
             RunDate = Sys.Date()) %>% 
      select(RunDate, ClusterID_FirstIdentified, State_Case_Number)
    
    #Importing case history file to assist in determining if clusters are new or ongoing
    case_hist <- read_csv('Case History File.csv', col_types = "Dccc")
    
    #Determine if cluster is all new or new cases added to ongoing cluster
    if (any(cases_involved$State_Case_Number %in% case_hist$State_Case_Number) == FALSE) {
      
      #Add result to log file
      sink(file = 'satscan_function_log.txt', append = T)
      cat('Significant cluster detected for', disease_run_name, '- NEW \n')
      sink()
      
      #Add new cases to history file
      write_csv(cases_involved,'Case History File.csv', append = TRUE)
      
      #Select variables of interest for report
      line_list <- eventsFilter %>% 
        filter(State_Case_Number %in% cases_involved$State_Case_Number) %>%
        mutate(Attends_Resides_Name = ifelse(is.na(Day_Care_Facility_Name), Other_Facility, Day_Care_Facility_Name)) %>%
        select(State_Case_Number, 
               Serogroup:Event_Date, 
               Last_Name = Current_Last_Name,
               First_Name = Current_First_Name,
               Sex = Sex_at_Onset,
               Age = `Current_Age_(yrs)`,
               Race = Race_List,
               Ethnicity,
               Attends_Resides = Patient_Attends_Resides,
               Attends_Resides_Name,
               Epi_Comment,
               City = Patient_Home_City,
               Lat = `Latitude_(Home)`,
               Long = `Longitude_(Home)`)
      
      #Prepare line list
      line_list_table <- line_list %>%
        select(-Lat, -Long) %>%
        select_if(function(col) !(all(is.na(col)))) %>% #Will eliminate any columns that are empty for all cases, e.g. Serogroup
        arrange(Event_Date)
      
      #Prepare map
      cluster <- session$shapeclust[session$shapeclust$CLUSTER == i, ] 
      
      line_list_map <- SpatialPointsDataFrame(coords=line_list[, c("Long", "Lat")], line_list, 
                                              proj4string = CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"))
      
      cluster_map <- leaflet() %>% addProviderTiles(providers$CartoDB.Positron) %>%
      
        addPolygons(data = cluster,
                     popup = sprintf("%s<br/>%s<br/>%s", 
                                     paste("Date Range:", cluster$START_DATE, "to", cluster$END_DATE), 
                                     paste("Observed Cases:", cluster$OBSERVED), 
                                     paste("Expected Cases:", round(cluster$EXPECTED, 2))) %>% lapply(htmltools::HTML),
                    group = "SatScan Cluster")  %>%
        
        addCircleMarkers(data = line_list_map,
                   fillColor = 	'black',
                   fillOpacity = .8,
                   stroke = FALSE,
                   radius = 5,
                   label = sprintf("%s<br/>%s", line_list_map$State_Case_Number, 
                                   paste(line_list_map$Age, " year old", line_list_map$Sex)) %>% lapply(htmltools::HTML),
                   group = "Cases in Cluster") %>%
        
        addLayersControl(overlayGroups = c("SatScan Cluster", "Cases in Cluster"))
      
      #Packaging data for final products
      report_data <- list(line_list_table, cluster_map, cluster_temp)
      
      #Save data --- IF THIS ENDS UP TAKING UP TOO MUCH SPACE, IT'S PROBABLY NOT NECESSARY SINCE THE MARKDOWN FILE WILL BE SAVED
      save(report_data, file = paste0("Clusters/Cluster ", cluster_temp$ClusterID, " Files.Rdata"))
      
      #Render report on cluster in R Markdown
      rmarkdown::render(input = "New Cluster Report.Rmd", 
                        output_file = paste0("Clusters/Cluster Report ID ", cluster_temp$ClusterID, ".html"),
                        params = list(data = report_data))
      
      #Email report to epidemiologist
      send.mail(from = "noreply@cookcountyhhs.org",
                to = unlist(strsplit(paramFilter$Email, ",")),                       
                subject = paste("SECURELOCK: SatScan has detected a new cluster for", disease_run_name),
                body = "Please review the attached line list and determine if the cluster should be investigated further.",   
                smtp = list(host.name = key_get("cchhs_host"), #port = **,
                            user.name = key_get("rk_email"),
                            passwd = key_get("office365"), tls = TRUE),
                attach.files = paste0("Clusters/Cluster Report ID ", cluster_temp$ClusterID, ".html"),
                authenticate = TRUE,
                send = TRUE) 

      
    } else {      #new case(s) in ongoing cluster
      
      #Determine which cluster new cases are part of
      which_clust <- case_hist %>% 
        filter(case_hist$State_Case_Number %in% cases_involved$State_Case_Number) %>% 
        arrange(RunDate) %>%
        slice(1) %>%
        select(ClusterID_FirstIdentified) %>%
        as.character()
      
      #Determine which cases are new and assign ongoing cluster ID
      new_cases_involved <- cases_involved %>% 
        filter(!cases_involved$State_Case_Number %in% case_hist$State_Case_Number) %>%
        mutate(ClusterID_Assigned = which_clust)
      
      #Only generate output if new case was added to the cluster (as opposed to ongoing cluster, no new cases)
      if (nrow(new_cases_involved) >= 1) {
        
        #Add result to log file
        sink(file = 'satscan_function_log.txt', append = T)
        cat('Significant cluster detected for', disease_run_name, '- ONGOING - NEW CASE(S) \n')
        sink()
        
        #Add new cases to history file
        write_csv(new_cases_involved,'Case History File.csv', append = TRUE)
        
        #Prepare line list
        ongoing_case_list <- rbind(new_cases_involved, filter(case_hist, ClusterID_FirstIdentified == which_clust | ClusterID_Assigned == which_clust)) 
        
        line_list <- inner_join(eventsFilter, ongoing_case_list) %>% 
          mutate(Attends_Resides_Name = ifelse(is.na(Day_Care_Facility_Name), Other_Facility, Day_Care_Facility_Name)) %>%
          select(Date_Added = RunDate,
                 State_Case_Number,
                 Serogroup:Event_Date,
                 Last_Name = Current_Last_Name,
                 First_Name = Current_First_Name,
                 Sex = Sex_at_Onset,
                 Age = `Current_Age_(yrs)`,
                 Race = Race_List,
                 Ethnicity,
                 Attends_Resides = Patient_Attends_Resides,
                 Attends_Resides_Name,
                 Epi_Comment,
                 City = Patient_Home_City,
                 Original_Cluster = ClusterID_FirstIdentified,
                 Lat = `Latitude_(Home)`,
                 Long = `Longitude_(Home)`) 
        
        line_list_table <- line_list %>%
          select(-Lat, -Long) %>%
          select_if(function(col) !(all(is.nas(col)))) %>% #Will eliminate any columns that are empty for all cases, e.g. Serogroup
          arrange(Event_Date)
        
        #Prepare map
        cluster <- session$shapeclust[session$shapeclust$CLUSTER == i, ]   
        
        line_list_map <- SpatialPointsDataFrame(coords=line_list[, c("Long", "Lat")], line_list, 
                                                proj4string = CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"))
        
        line_list_map@data <- line_list_map@data %>%   #Flagging cases added today
          mutate(new_case_flag = ifelse(line_list_map$State_Case_Number %in% new_cases_involved$State_Case_Number, 1, 0))
        
        cluster_map <- leaflet() %>% addProviderTiles(providers$CartoDB.Positron) %>% 
          
          addPolygons(data = cluster,
                      popup = sprintf("%s<br/>%s<br/>%s", 
                                      paste("Date Range:", cluster$START_DATE, "to", cluster$END_DATE), 
                                      paste("Observed Cases:", cluster$OBSERVED), 
                                      paste("Expected Cases:", round(cluster$EXPECTED, 2))) %>% lapply(htmltools::HTML),
                      group = "New Area of Cluster")  %>%
          
          addCircleMarkers(data = line_list_map,
                           fillColor = 	~ifelse(line_list_map$new_case_flag == 1, "red", "black"),
                           fillOpacity = .8,
                           stroke = FALSE,
                           radius = 5,
                           label = sprintf("%s<br/>%s", line_list_map$State_Case_Number, 
                                           paste(line_list_map$Age, " year old", line_list_map$Sex, ", Onset:", format(line_list_map$Onset_Date, "%m/%d/%y"))) 
                           %>% lapply(htmltools::HTML),
                           group = "Cases in Cluster") %>%
          
          addLayersControl(overlayGroups = c("New Area of Cluster", "Cases in Cluster"))
        
        #Packaging data for final products
        report_data <- list(which_clust, line_list_table, cluster_map, cluster_temp, new_cases_involved[ , 3])
        
        #Saving data (currently replacing original file, might need to change later if all records wanted)
        save(report_data, file = paste0("Clusters/Cluster ", which_clust, " Files.Rdata"))
        
        #Render report on cluster in R Markdown
        rmarkdown::render(input = "Ongoing Cluster Report.Rmd", 
                          output_file = paste0("Clusters/Cluster Report ID ", which_clust, "_Updated.html"),
                          params = list(data = report_data))
        
        #Email report to epidemiologist
        send.mail(from = "noreply@cookcountyhhs.org",
                  to = unlist(strsplit(paramFilter$Email, ",")),                       
                  subject = paste("SECURELOCK: SatScan has detected new cases in an ongoing", disease_run_name, "cluster"),
                  body = paste0("New cases have been detected for cluster ", which_clust, ". Please review the updated line list and map attached."),   
                  smtp = list(host.name = key_get("cchhs_host"), #port = **,
                              user.name = key_get("rk_email"),
                              passwd = key_get("office365"), tls = TRUE),
                  attach.files = paste0("Clusters/Cluster Report ID ", which_clust, "_Updated.html"),
                  authenticate = TRUE,
                  send = TRUE) 
        
      } else {
        
        #Add result to log file
        sink(file = 'satscan_function_log.txt', append = T)
        cat('Significant cluster detected for', disease_run_name, '- ONGOING \n')
        sink()
        
        
      } #ongoing cluster, new cases if/else closure
      

      
    } #new or ongoing cluster if/else closure
    
    
    
  } #significant cluster for loop closure
  
  
  
} #is cluster significant if/else closure

} #run_satscan function closure


walk(parameters$Disease, run_satscan)  

# run_satscan("Campylobacteriosis")
run_satscan("Cryptosporidiosis")
# run_satscan("Cyclosporiasis") 
# run_satscan("Giardiasis")    
# run_satscan("Hepatitis A") 
# run_satscan("Legionellosis")
# run_satscan("Mumps") 
# run_satscan("Pertussis") 
# run_satscan("Salmonellosis") 
# run_satscan("STEC O157:H7")
# run_satscan("STEC Not Cultured")
# run_satscan("STEC Non-O157") 
# run_satscan("Shigellosis") 
# run_satscan("Group A Strep") 
# run_satscan("Varicella")




#=================NOTIFYING PROGRAMMER OF JOB RUN RESULTS================#

#Notify programmer via email of the scan results
if(as.Date(file.mtime("Cluster History File.csv")) == Sys.Date()){
  
  send.mail(from = "noreply@cookcountyhealth.org",
            to = key_get("rk_email"),
            subject = "Daily SatScan Report - Exact Locations",
            body = paste("Results of the run for", Sys.Date(), "are attached."),   
            smtp = list(host.name = key_get("cchhs_host"), #port = **,
                        user.name = key_get("rk_email"),
                        passwd = key_get("office365"), tls = TRUE),
            attach.files = c("satscan_function_log.txt", "try_error_log.txt"),
            authenticate = TRUE,
            send = TRUE) 
  
} else{ 
  
  send.mail(from = "noreply@cookcountyhealth.org",
            to = key_get("rk_email"),
            subject = "ATTENTION: Daily SatScan Cluster Detection Failed - Exact Locations",
            body = "Please check the script and troubleshoot.",
            smtp = list(host.name = key_get("cchhs_host"), #port = **,
                        user.name = key_get("rk_email"),
                        passwd = key_get("office365"), tls = TRUE),
            authenticate = TRUE,
            send = TRUE)
  
}

#Delete log file after emailed
unlink('satscan_function_log.txt')
unlink('try_error_log.txt')

#Delete SatScan Download
#unlink('C:/Users/kbemis/Downloads/SatScan_Data.csv')



