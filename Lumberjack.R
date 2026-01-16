
require(data.table)
require(devtools)
require(rstudioapi)
require(RCurl)
require(tools)

# Basics copied from Ben
today <- function(){
  format(Sys.time(), "%Y_%m_%d")
}

DateFileName <- function(x){
  name <-   paste0(today(), "_", x)
  print (name)
  return (name)
}

tdy <- function(){
  format(Sys.time(), "%y.%m.%d")
}

date_it <- function(x){
  name <- paste0(tdy(), "_", x)
  return (name)
}

moment <- function(){
  format(Sys.time(), "%y.%m.%d--%H.%M.%S")
}

moment_it <- function(x){
  name <- paste0(moment(), "_", x)
  return (name)
}

get_smpl <- function(x){
  return( tools::file_path_sans_ext( basename( x ) ) )
}

ScriptName <- function(scriptName = NULL){
  if(is.null(scriptName))
    scriptName <- rstudioapi::getActiveDocumentContext()$path
  if (is.null (scriptName) || scriptName == "")
    stop("No script name found -- you may need to save this file first")
  else 
    return(scriptName)
}

ScriptNamedDir <- function(scriptName = NULL){
  scriptName <- ScriptName(scriptName)
  outDir <- gsub(".[qR](md)?$", "_data", scriptName, ignore.case = TRUE)
  stopifnot( outDir != scriptName)
  if (!dir.exists(outDir) & grepl("[rR]md", tools::file_ext(scriptName)) ){
    message ("Creating directory associated with ", scriptName,", at ", outDir)
    dir.create(outDir)
  }
  return(outDir)
}

# Basics for logging 

get_logs_dir <- function(){
  # Check if there exists a _data folder with the same name as current script, 
  # return (and create if necessary) the logs folder inside <sript>_data or 
  # if there is no _data folder, return/create logs folder in parent directory
  scriptDir <- paste0("./",basename(ScriptNamedDir()))
  if (scriptDir %in% list.dirs(recursive = F)){
    logsFolder <- file.path(ScriptNamedDir(),"logs")
    if (  !file.exists(logsFolder) ){
      message ("Creating logs directory at ", logsFolder)
      dir.create(logsFolder)
    }
  } else{
    logsFolder <- file.path(dirname(ScriptName()), "logs")
    if (  !file.exists(logsFolder) ){
      message ("Creating logs directory at ", logsFolder)
      dir.create(logsFolder)
    }
  }
  return(logsFolder) 
}

get_git_commit_sha <- function(filepath = ".") {
  # Use 'git log' to get the last commit SHA for a specific file or the current directory
  # --pretty=format:"%H" ensures only the full SHA is returned
  # -n 1 limits the output to the most recent commit
  # --follow is useful if the file has been moved/renamed (works with 'git log' with a single file path)
  dir <- dirname(filepath)
  filename <- basename(filepath)
  command <- paste("cd ", dir,  "; git log -n 1 --pretty=format:\"%H\" --", filename)
  #message (command)
  
  # Execute the command and capture the output (intern=TRUE)
  sha <- system(command, intern = TRUE, ignore.stderr = TRUE)
  
  if (length(sha) == 0) {
    return("No commit found or not a git repository")
  } else {
    return(sha)
  }
}

get_git_repo_sha <- function(filepath = "."){
  dir <- dirname(filepath)
  command <- paste("cd ", dir,  "; git config --get remote.origin.url")
  #message (command)
  
  # Execute the command and capture the output (intern=TRUE)
  sha <- system(command, intern = TRUE, ignore.stderr = TRUE)
  
  if (length(sha) == 0) {
    return("No repo URL found or not a git repository")
  } else {
    return(sha)
  }
}

get_git_commit_hyperlink <- function(filepath = "."){
  repoUrl <- get_git_repo_sha(filepath)
  commitHash <- get_git_commit_sha(filepath)
  commitUrl <- paste0( tools::file_path_sans_ext( repoUrl ), "/commit/", commitHash  )
  #message(commitUrl)
  
  if ( RCurl::url.exists(commitUrl)){
    return( commitUrl )
  } else {
    message( sprintf("Unable to Find Git Commit URL associated with %s", filepath) )
    return( "NoCommitFound" )
  }
}

check_if_gitRemote_upToDate <- function(filepath = "."){
  return(F)
}

get_source_log_file <- function(logFileName = "00_sourced_Rcode_log.tsv"){
  logDir <- get_logs_dir()
  sourcelogFile <- file.path(logDir, logFileName)
  if (  !file.exists(sourcelogFile) ){
    message ("Creating sourced R code log file in folder ", logDir)
    headersOnly <- data.table::data.table("Date" = character(length = 0L), "Time" = character(length = 0L), 
                              "Source_Requested" = character(length = 0L), "Source_Used" = character(length = 0L), 
                              "Logged_As" = character(length = 0L), "Assc._Git_Commit" = character(length = 0L))
    data.table::fwrite(x = headersOnly, file = sourcelogFile, sep = "\t")
  }
  return(sourcelogFile)
}

log_source <- function(sourceRequested){
  message(sprintf("Sourcing and Logging R Script:   %s   ", basename(sourceRequested)))
  sourceFile <- copy(sourceRequested)
  recentSource <- NULL
  sourceOriginal <- T
  codeChanged <- T
  day <- moment() |> tstrsplit(split = "--", keep = 1) |> unlist()
  timeOfDay <- moment() |> tstrsplit(split = "--", keep = 2) |> unlist()
  
  # find log folder and set logged code filename
  logDir <- get_logs_dir()
  loggedCodeName <- paste0( "rSourced_", basename(sourceFile)) #get_smpl( ScriptName() ) )
  loggedCodePath <- file.path( logDir, date_it(loggedCodeName) )

  # find last sourced code in logs folder
  logFiles <- list.files( path = get_logs_dir(), full.names = T)
  
  if ( length(logFiles) > 0 )   logFiles <- logFiles[ grepl(loggedCodeName, logFiles) ] 
  if ( length(logFiles) > 0 )   recentSource <- logFiles[ file.info(logFiles)$ctime |> which.max() ]
  
  # if the target source file cannot be found, it will check the logs directory for last sourced file
  if ( !file.exists(sourceFile) & !is.null(recentSource) ){
    warning( sprintf("File Not Found.  %s \nSource Code File Found in logs Directory.  %s \n\nSourcing Previously Logged Code Instead...", sourceFile, recentSource))
    sourceFile <- recentSource
    sourceOriginal <- F
  } else if (!file.exists(sourceFile) & is.null(recentSource)) {
    stop( sprintf("File Not Found.  %s \nSource Code File Not Found in logs Directory.  %s. \n\nCheck path and try again.", sourceFile, logDir))
  }
  
  # make sure source code has no syntax errors before continuing -- parse() is the canary
  parse(sourceFile) |> capture.output() |> invisible()
    
  # check if source identical to last sourced version, if not, maintain both sources, if they are identical, do not write new code to avoid cluttering directory
  if (!is.null(recentSource))   codeChanged <- !identical( as.character(tools::md5sum(sourceFile)), as.character(tools::md5sum(recentSource)) )
  
  # write new sourced code
  if (codeChanged) {
    file.copy(from = sourceFile, to = loggedCodePath)
    warning("Sourced code content is changed from previous logged source")
    message(sprintf("New code source logged to %s", loggedCodePath))
  } else {
    loggedCodePath <- sprintf("Content unchanged, see previous logged %s", basename(sourceRequested) )
    message(sprintf("Code source content identical to previous log %s", recentSource))
  }
  
  # write entry in log
  sourceInfo <- data.table::data.table("Date" = c(day), "Time" = c(timeOfDay), 
                           "Source_Requested" = c(sourceRequested), "Source_Used" = sourceFile,
                           "Logged_As" = c(loggedCodePath), "Assc._Git_Commit" = get_git_commit_hyperlink(sourceFile) )
  data.table::fwrite(x = sourceInfo, file = get_source_log_file(), append = T, sep = "\t")

  # finalllyyyy source it!
  source(sourceFile)
}