getGeneScore_ArchR <- function(
  ArchRProj = projHeme1,
  name = markerGenes,
  imputeWeights = getImputeWeights(projHeme1),
  threads = getArchRThreads(),
  logFile = createLogFile("plotEmbedding"),
  ...
){
  colorBy = "GeneScoreMatrix"
  
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = name, name = "name", valid = c("character"))
  .validInput(input = imputeWeights, name = "imputeWeights", valid = c("list", "null"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))
  
  .requirePackage("ggplot2", source = "cran")
  
  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "Input-Parameters", logFile=logFile)
  
  suppressMessages(message(logFile))
  
  
  colorMat <- .getMatrixValues(
    ArchRProj = ArchRProj, 
    name = name, 
    matrixName = colorBy, 
    log2Norm = FALSE, 
    threads = threads,
    logFile = logFile
  )
  
  
  .logThis(colorMat, "colorMat-Before-Impute", logFile = logFile)
  
  if(!is.null(imputeWeights)){
    message("Imputing Matrix")
    colorMat <- imputeMatrix(mat = as.matrix(colorMat), imputeWeights = imputeWeights, logFile = logFile)
    if(!inherits(colorMat, "matrix")){
      colorMat <- matrix(colorMat, ncol = nrow(df))
      colnames(colorMat) <- rownames(df)
    }
  }
  
  .logThis(colorMat, "colorMat-After-Impute", logFile = logFile)
  
  
  gene_score <- as.data.frame(colorMat)
  new_col_names <- colnames(gene_score)
  new_col_names <- unlist(lapply(new_col_names, function(x) gsub(".*#","", x)))
  new_col_names <- unlist(lapply(new_col_names, function(x) gsub("-.*","", x)))
  colnames(gene_score) <- new_col_names
  gene_score <- log2(gene_score + 1)
    

  message("")
  
  .endLogging(logFile = logFile)
  
  return(gene_score)

}


# h5read implementation for optimal reading
.h5read <- function(
  file = NULL,
  name = NULL,
  method = "fast",
  index = NULL,
  start = NULL,
  block = NULL,
  count = NULL
){
  
  if(tolower(method) == "fast" & is.null(index) & is.null(start) & is.null(block) & is.null(count)){
    fid <- H5Fopen(file)
    dapl <- H5Pcreate("H5P_DATASET_ACCESS")
    did <- .Call("_H5Dopen", fid@ID, name, dapl@ID, PACKAGE='rhdf5')
    res <- .Call("_H5Dread", did, NULL, NULL, NULL, TRUE, 0L, FALSE, fid@native, PACKAGE='rhdf5')
    invisible(.Call("_H5Dclose", did, PACKAGE='rhdf5'))   
  }else{
    res <- h5read(file = file, name = name, index = index, start = start, block = block, count = count)
  }
  o <- h5closeAll()
  return(res)
}


.getMatrixValues <- function(
  ArchRProj = NULL, 
  name = NULL, 
  matrixName = NULL, 
  log2Norm = FALSE, 
  threads = getArchRThreads(),
  logFile = NULL
){
  
  o <- h5closeAll()
  
  .logMessage("Getting Matrix Values...", verbose = TRUE, logFile = logFile)
  
  featureDF <- .getFeatureDF(head(getArrowFiles(ArchRProj), 2), matrixName)
  .logThis(featureDF, "FeatureDF", logFile = logFile)
  
  matrixClass <- h5read(getArrowFiles(ArchRProj)[1], paste0(matrixName, "/Info/Class"))
  
  if(matrixClass == "Sparse.Assays.Matrix"){
    if(!all(unlist(lapply(name, function(x) grepl(":",x))))){
      .logMessage("When accessing features from a matrix of class Sparse.Assays.Matrix it requires seqnames\n(denoted by seqnames:name) specifying to which assay to pull the feature from.\nIf confused, try getFeatures(ArchRProj, useMatrix) to list out available formats for input!", logFile = logFile)
      stop("When accessing features from a matrix of class Sparse.Assays.Matrix it requires seqnames\n(denoted by seqnames:name) specifying to which assay to pull the feature from.\nIf confused, try getFeatures(ArchRProj, useMatrix) to list out available formats for input!")
    }
  }
  
  if(grepl(":",name[1])){
    
    sname <- stringr::str_split(name,pattern=":",simplify=TRUE)[,1]
    name <- stringr::str_split(name,pattern=":",simplify=TRUE)[,2]
    
    idx <- lapply(seq_along(name), function(x){
      ix <- intersect(which(tolower(name[x]) == tolower(featureDF$name)), BiocGenerics::which(tolower(sname[x]) == tolower(featureDF$seqnames)))
      if(length(ix)==0){
        .logStop(sprintf("FeatureName (%s) does not exist! See getFeatures", name[x]), logFile = logFile)
      }
      ix
    }) %>% unlist
    
  }else{
    
    idx <- lapply(seq_along(name), function(x){
      ix <- which(tolower(name[x]) == tolower(featureDF$name))[1]
      if(length(ix)==0){
        .logStop(sprintf("FeatureName (%s) does not exist! See getFeatures", name[x]), logFile = logFile)
      }
      ix
    }) %>% unlist
    
  }
  .logThis(idx, "idx", logFile = logFile)
  
  if(any(is.na(idx))){
    .logStop(sprintf("FeatureName (%s) does not exist! See getFeatures", paste0(name[which(is.na(idx))], collapse=",")), logFile = logFile)
  }
  
  featureDF <- featureDF[idx, ,drop=FALSE]
  .logThis(featureDF, "FeatureDF-Subset", logFile = logFile)
  
  #Get Values for FeatureName
  cellNamesList <- split(rownames(getCellColData(ArchRProj)), getCellColData(ArchRProj)$Sample)
  
  values <- .safelapply(seq_along(cellNamesList), function(x){
    message(x, " ", appendLF = FALSE)
    valuesx <- tryCatch({
      o <- h5closeAll()
      ArrowFile <- getSampleColData(ArchRProj)[names(cellNamesList)[x],"ArrowFiles"]
      valuesx <- .getMatFromArrow(
        ArrowFile = ArrowFile, 
        featureDF = featureDF,
        binarize = FALSE, 
        useMatrix = matrixName, 
        cellNames = cellNamesList[[x]],
        threads = 1
      )
      colnames(valuesx) <- cellNamesList[[x]]
      valuesx
    }, error = function(e){
      errorList <- list(
        x = x,
        ArrowFile = ArrowFile,
        ArchRProj = ArchRProj, 
        cellNames = ArchRProj$cellNames, 
        cellNamesList = cellNamesList, 
        featureDF = featureDF
      )
      .logError(e, fn = ".getMatFromArrow", info = "", errorList = errorList, logFile = logFile)  
    })
    valuesx
  }, threads = threads) %>% Reduce("cbind", .)
  values <- values[, ArchRProj$cellNames, drop = FALSE]
  message("")
  gc()
  .logThis(values, "Feature-Matrix", logFile = logFile)
  
  if(!inherits(values, "matrix")){
    values <- matrix(as.matrix(values), ncol = nCells(ArchRProj))
    colnames(values) <- ArchRProj$cellNames
  }
  
  #Values Summary
  if(!is.null(log2Norm)){
    if(log2Norm){
      message("Log2 Normalizing...")
      values <- log2(values + 1)
    }
  }
  
  rownames(values) <- name
  
  return(values)
  
}


#
.getFeatureDF <- function(ArrowFiles = NULL, subGroup = "TileMatrix", threads = getArchRThreads()){
  
  threads <- min(threads, length(ArrowFiles))
  
  .helpFeatureDF <- function(ArrowFile = NULL, subGroup = NULL){
    o <- h5closeAll()
    featureDF <- DataFrame(h5read(ArrowFile, paste0(subGroup,"/Info/FeatureDF")))
    featureDF$seqnames <- Rle(as.character(featureDF$seqnames))
    o <- h5closeAll()
    return(featureDF)
  }
  
  fdf <- .helpFeatureDF(ArrowFiles[1], subGroup = subGroup)
  
  if(length(ArrowFiles) > 1){
    ArrowFiles <- ArrowFiles[-1]
    checkIdentical <- .safelapply(seq_along(ArrowFiles), function(x){
      fdfx <- .helpFeatureDF(ArrowFiles[x], subGroup = subGroup)
      identical(fdfx, fdf)
    }, threads = threads) %>% unlist %>% all
    if(!checkIdentical){
      stop("Error not all FeatureDF for asssay is the same!")
    }
  }
  
  #Re-Order for Split Check!
  newOrder <- split(seq_len(nrow(fdf)), fdf$seqnames) %>% {lapply(seq_along(.), function(x) .[[x]])} %>% Reduce("c", .)
  fdf[newOrder,]
  
}


#
.safelapply <- function(..., threads = 1, preschedule = FALSE){
  
  if(tolower(.Platform$OS.type) == "windows"){
    threads <- 1
  }
  
  if(threads > 1){
    
    o <- mclapply(..., mc.cores = threads, mc.preschedule = preschedule)
    
    errorMsg <- list()
    
    for(i in seq_along(o)){ #Make Sure this doesnt explode!
      if(inherits(o[[i]], "try-error")){
        capOut <- utils::capture.output(o[[i]])
        capOut <- capOut[!grepl("attr\\(\\,|try-error", capOut)]
        capOut <- head(capOut, 10)
        capOut <- unlist(lapply(capOut, function(x) substr(x, 1, 250)))
        capOut <- paste0("\t", capOut)
        errorMsg[[length(errorMsg) + 1]] <- paste0(c(paste0("Error Found Iteration ", i, " : "), capOut), "\n")
      }
    }
    
    if(length(errorMsg) != 0){
      
      errorMsg <- unlist(errorMsg)
      errorMsg <- head(errorMsg, 50)
      errorMsg[1] <- paste0("\n", errorMsg[1])
      stop(errorMsg)
      
    }
    
  }else{
    
    o <- lapply(...)
    
  }
  
  o
  
}


#
.getMatFromArrow <- function(
  ArrowFile = NULL, 
  featureDF = NULL, 
  binarize = NULL, 
  cellNames = NULL,
  useMatrix = "TileMatrix", 
  useIndex = FALSE,
  threads = 1
){
  
  if(is.null(featureDF)){
    featureDF <- .getFeatureDF(ArrowFile, useMatrix)
  }
  
  if(any(c("seqnames","idx") %ni% colnames(featureDF))){
    stop("Need to provide featureDF with columns seqnames and idx!")
  }
  
  #Add RowNames for Check at the end
  rownames(featureDF) <- paste0("f", seq_len(nrow(featureDF)))
  
  o <- h5closeAll()
  
  matClass <- h5read(ArrowFile, paste0(useMatrix,"/Info/Class"))
  if(matClass %ni% c("Sparse.Binary.Matrix", "Sparse.Integer.Matrix", "Sparse.Double.Matrix", "Sparse.Assays.Matrix")){
    stop("Arrow Mat is not a valid Sparse Matrix!")
  }
  if(is.null(binarize)){
    if(matClass == "Sparse.Binary.Matrix"){
      binarize <- TRUE
    }else{
      binarize <- FALSE
    }
  }
  if(matClass == "Sparse.Binary.Matrix"){
    if(!binarize){
      stop("Sparse Matrix in Arrow is Binarized! Set binarize = TRUE to use matrix!")
    }
  }
  
  matColNames <- paste0(.sampleName(ArrowFile), "#", h5read(ArrowFile, paste0(useMatrix,"/Info/CellNames")))
  if(!is.null(cellNames)){
    idxCols <- which(matColNames %in% cellNames)
  }else{
    idxCols <- seq_along(matColNames)
  }
  
  seqnames <- unique(featureDF$seqnames)
  
  mat <- .safelapply(seq_along(seqnames), function(x){
    
    seqnamex <- seqnames[x]
    featureDFx <- featureDF[BiocGenerics::which(featureDF$seqnames %bcin% seqnamex),]
    idxRows <- featureDFx$idx
    
    j <- Rle(
      values = h5read(ArrowFile, paste0(useMatrix,"/",seqnamex,"/jValues")), 
      lengths = h5read(ArrowFile, paste0(useMatrix,"/",seqnamex,"/jLengths"))
    )
    
    #Match J
    matchJ <- S4Vectors::match(j, idxCols, nomatch = 0)
    idxJ <- BiocGenerics::which(matchJ > 0)
    if(useIndex){
      i <- h5read(ArrowFile, paste0(useMatrix,"/",seqnamex,"/i"), index = list(idxJ, 1))
    }else{
      i <- h5read(ArrowFile, paste0(useMatrix,"/",seqnamex,"/i"))[idxJ]
    }
    j <- matchJ[idxJ]
    
    #Match I
    matchI <- match(i, idxRows, nomatch = 0)
    idxI <- which(matchI > 0)
    i <- i[idxI]
    j <- j[idxI]
    i <- matchI[idxI]
    
    if(!binarize){
      x <- h5read(ArrowFile, paste0(useMatrix,"/",seqnamex,"/x"))[idxJ][idxI]
    }else{
      x <- rep(1, length(j))
    }
    
    mat <- Matrix::sparseMatrix(
      i=as.vector(i),
      j=j,
      x=x,
      dims = c(length(idxRows), length(idxCols))
    )
    rownames(mat) <- rownames(featureDFx)
    
    rm(matchI, idxI, matchJ, idxJ, featureDFx, idxRows)
    
    return(mat)
    
  }, threads = threads) %>% Reduce("rbind", .)
  
  o <- h5closeAll()
  
  colnames(mat) <- matColNames[idxCols]
  
  #Double Check Order!
  mat <- mat[rownames(featureDF), , drop = FALSE]
  rownames(mat) <- NULL
  
  if(!is.null(cellNames)){
    mat <- mat[,cellNames,drop=FALSE]
  }
  
  return(mat)
  
}


#
.sampleName <- function(ArrowFile = NULL){
  o <- h5closeAll()
  sampleName <- h5read(ArrowFile, paste0("Metadata/Sample"))
  o <- h5closeAll()
  return(sampleName)
}


####################################
# Log Tools
####################################

#' Set ArchR Logging
#' 
#' This function will set ArchR logging
#'
#' @param useLogs A boolean describing whether to use logging with ArchR.
#' @export
addArchRLogging <- function(useLogs = TRUE){
  .validInput(input = useLogs, name = "useLogs", valid = "boolean")
  message("Setting ArchRLogging = ", useLogs)
  options(ArchR.logging = useLogs)
  return(invisible(0))
}

#' Get ArchR Logging
#' 
#' This function will get ArchR logging
#'
#' @export
getArchRLogging <- function(){
  ArchRLogging <- options()[["ArchR.logging"]]
  if(!is.logical(ArchRLogging)){
    options(ArchR.logging = TRUE)
    return(TRUE)
  }
  ArchRLogging
}

#' Set ArchR Debugging
#' 
#' This function will set ArchR Debugging which will save an RDS if an error is encountered.
#'
#' @param debug A boolean describing whether to use logging with ArchR.
#' @export
addArchRDebugging <- function(debug = FALSE){
  .validInput(input = debug, name = "debug", valid = "boolean")
  message("Setting ArchRDebugging = ", debug)
  options(ArchR.logging = debug)
  return(invisible(0))
}

#' Get ArchR Debugging
#' 
#' This function will get ArchR Debugging which will save an RDS if an error is encountered.
#'
#' @export
getArchRDebugging <- function(){
  ArchRDebugging <- options()[["ArchR.debugging"]]
  if(!is.logical(ArchRDebugging)){
    options(ArchR.debugging = FALSE)
    return(FALSE)
  }
  ArchRDebugging
}

#' Create a Log File for ArchR
#' 
#' This function will create a log file for ArchR functions. If ArchRLogging is not TRUE
#' this function will return NULL.
#'
#' @param name A character string to add a more descriptive name in log file.
#' @param logDir The path to a directory where log files should be written.
#' @export
createLogFile <- function(
  name = NULL,
  logDir = "ArchRLogs",
  useLogs = getArchRLogging()
){
  
  .validInput(input = name, name = "name", valid = "character")
  .validInput(input = logDir, name = "logDir", valid = "character")
  .validInput(input = useLogs, name = "useLogs", valid = "boolean")
  
  if(!useLogs){
    return(NULL)
  }
  dir.create(logDir, showWarnings = FALSE)
  if(is.null(name)){
    logFile <- .tempfile(pattern = "ArchR", fileext = ".log", tmpdir = logDir)
  }else{
    logFile <- .tempfile(pattern = paste0("ArchR-", name), fileext = ".log", tmpdir = logDir)
  }
  logFile
}

.messageDiffTime <- function(...){ #Deprecated
  .logDiffTime(...)
}

.logDiffTime <- function(
  main = "",
  t1 = NULL,
  verbose = TRUE,
  addHeader = FALSE,
  t2 = Sys.time(),
  units = "mins",
  header = "###########",
  tail = "elapsed.",
  precision = 3,
  logFile = NULL,
  useLogs = getArchRLogging()
){
  
  if(verbose){
    
    timeStamp <- tryCatch({
      dt <- abs(round(difftime(t2, t1, units = units),precision))
      if(addHeader){
        msg <- sprintf("%s\n%s : %s, %s %s %s\n%s", header, Sys.time(), main, dt, units, tail, header)
      }else{
        msg <- sprintf("%s : %s, %s %s %s", Sys.time(), main, dt, units, tail)
      }
      message(msg)
    }, error = function(x){
      message("Time Error : ", x)
    })
    
  }
  
  if(!useLogs){
    return(invisible(0))
  }
  
  if(!is.null(logFile)){
    if(file.exists(logFile)){
      logStamp <- tryCatch({
        dt <- abs(round(difftime(t2, t1, units = units),precision))
        if(addHeader){
          msg <- sprintf("%s\n%s : %s, %s %s %s\n%s", header, Sys.time(), main, dt, units, tail, header)
        }else{
          msg <- sprintf("%s : %s, %s %s %s", Sys.time(), main, dt, units, tail)
        }
        cat(paste0(msg,"\n"), file = logFile, append = TRUE)
      }, error = function(x){
        0
      })
    }
  }
  
  return(invisible(0))
  
}

.startLogging <- function(
  logFile = NULL, 
  useLogs = getArchRLogging()
){
  
  if(!useLogs){
    return(invisible(0))
  }
  
  if(is.null(logFile)){
    return(invisible(0))
  }
  
  if(file.exists(logFile)){
    return(invisible(0))
  }
  
  .getRam <- function(OS = .Platform$OS.type){
    if(grepl("linux", OS, ignore.case = TRUE)){
      ram <- paste0("Linux : ", as.numeric(system("awk '/MemTotal/ {print $2}' /proc/meminfo", intern = TRUE)))
    }else if(grepl("unix", OS, ignore.case = TRUE)){
      ram <- system("/usr/sbin/system_profiler SPHardwareDataType", intern = TRUE)
      ram <- paste0("MAC : ", gsub("Memory:","",gsub(" ","", grep("Memory", ram, value = TRUE))))
    }else{
      ram <- NA
    }
  }
  
  message("ArchR logging to : ", logFile, 
          "\nIf there is an issue, please report to github with logFile!")
  
  #Begin With
  cat(.ArchRLogo(ascii = "Package", messageLogo = FALSE), file = logFile, append = FALSE) 
  cat("\nLogging With ArchR!\n\n", file = logFile, append = TRUE) 
  cat(paste0("Start Time : ",Sys.time(),"\n\n"), file = logFile, append = TRUE)
  
  #ArchR Info
  cat("------- ArchR Info\n\n", file = logFile, append = TRUE)
  cat(paste0("ArchRThreads = ", getArchRThreads()), file = logFile, append = TRUE)
  tryCatch({
    if(!is.null(getArchRGenome())){
      cat(paste0("\nArchRGenome = ", getArchRGenome()), file = logFile, append = TRUE)
    }
  }, error = function(x){
  })
  cat("\n\n", file = logFile, append = TRUE)
  
  #Add Info
  cat("------- System Info\n\n", file = logFile, append = TRUE)
  cat(paste0("Computer OS = ", .Platform$OS.type), file = logFile, append = TRUE)
  tryCatch({
    cat(paste0("\nTotal Cores = ", detectCores()), file = logFile, append = TRUE)
  }, error = function(x){
  })
  # tryCatch({
  #     cat(paste0("\nTotal RAM = ", .getRam()), file = logFile, append = TRUE)
  # }, error = function(x){
  # })
  cat("\n\n", file = logFile, append = TRUE)
  
  #Session Info
  cat("------- Session Info\n\n", file = logFile, append = TRUE)
  utils::capture.output(sessionInfo(), file = logFile, append = TRUE)
  cat("\n\n------- Log Info\n\n", file = logFile, append = TRUE)
  
  return(invisible(0))
  
}

.logMessage <- function(
  ..., 
  logFile = NULL,
  verbose = FALSE,   
  useLogs = getArchRLogging()
){
  
  msg <- utils::capture.output(message(...), type = "message")
  msg <- paste0(msg, collapse = "\n")
  
  if(is.null(msg)){
    stop("Message must be provided when logging!")
  }
  
  if(verbose){
    message(sprintf("%s : %s", Sys.time(), msg))
  }
  
  if(!useLogs){
    return(invisible(0))
  }
  
  if(is.null(logFile)){
    return(invisible(0))
  }
  
  cat(sprintf("\n%s : %s\n", Sys.time(), msg), file = logFile, append = TRUE)
  
  return(invisible(0))
  
}

.logHeader <- function(
  name = NULL, 
  logFile = NULL,   
  useLogs = getArchRLogging()
){
  
  if(!useLogs){
    return(invisible(0))
  }
  
  if(is.null(logFile)){
    return(invisible(0))
  }
  
  if(is.null(name)){
    stop("Name must be provided when logging!")
  }
  
  header <- "###########"
  cat(sprintf("\n%s\n%s : %s\n%s\n\n", header, Sys.time(), name, header), file = logFile, append = TRUE)
  
  return(invisible(0))
}

.logStop <- function(
  ..., 
  logFile = NULL,
  useLogs = getArchRLogging()
){
  
  msg <- utils::capture.output(message(...), type = "message")
  msg <- paste0(msg, collapse = "\n")
  
  if(is.null(msg)){
    stop("Message must be provided when logging!")
  }
  
  if(useLogs){
    if(!is.null(logFile)){
      cat(sprintf("\n%s : %s\n", Sys.time(), msg), file = logFile, append = TRUE)
    }
  }
  
  stop(sprintf("%s\n", msg), call. = FALSE)
  
  return(invisible(0))
  
}

.logError <- function(
  e = NULL,
  fn = NULL,
  info = NULL, 
  errorList = NULL,
  logFile = NULL,   
  useLogs = getArchRLogging(),
  throwError = TRUE,
  debug = getArchRDebugging()
){
  
  header <- "************************************************************"
  
  if(is.null(logFile)){
    useLogs <- FALSE
  }
  
  if(useLogs){
    #To Log File
    cat(sprintf("\n%s\n%s : ERROR Found in %s for %s \nLogFile = %s\n\n", header, Sys.time(), fn, info, logFile), file = logFile, append = TRUE)
    
    utils::capture.output(print(e), file = logFile, append = TRUE)
    
    if(!is.null(errorList)){
      tryCatch({
        #.safeSaveRDS(errorList, "Save-Error.rds")
        .logThis(errorList, name = "errorList", logFile)
      }, error = function(e){
        cat("Error recording errorList", file = logFile, append = TRUE)
      })
    }
    
    cat(sprintf("\n%s\n\n", header), file = logFile, append = TRUE)
  }
  
  #To Console
  cat(sprintf("\n%s\n%s : ERROR Found in %s for %s \nLogFile = %s\n\n", header, Sys.time(), fn, info, logFile))
  
  if(debug){
    if(!is.null(errorList)){
      debugFile <- paste0(gsub("\\.log", "", logFile), "-debug.rds")
      cat(sprintf("\n%s : ArchRDebugging is set to TRUE, DebugFile = %s\n", Sys.time(), debugFile))
      .safeSaveRDS(errorList, debugFile)
    }
  }
  
  print(e)
  
  cat(sprintf("\n%s\n\n", header))
  
  if(throwError) stop("Exiting See Error Above")
  
  return(invisible(0))
  
}


.logThis <- function(
  x = NULL, 
  name = NULL, 
  logFile = NULL, 
  useLogs = getArchRLogging()
){
  
  if(!useLogs){
    return(invisible(0))
  }
  
  if(is.null(logFile)){
    return(invisible(0))
  }
  
  if(!file.exists(logFile)){
    stop("logFile does not exist! Something may have deleted this file! Exiting...")
  }
  if(is.null(name)){
    stop("Name must be provided when logging!")
  }
  cat(paste0("\n", Sys.time(), " : ", name, ", Class = ", class(x), "\n"), file = logFile, append = TRUE)
  
  if(missing(x)){
    cat("Data is Missing\n\n", file = logFile, append = TRUE)
    return(invisible(0))
  }
  
  if(is.matrix(x)){
    
    px <- x[head(seq_len(nrow(x)), 5), head(seq_len(ncol(x)), 5), drop = FALSE]
    suppressMessages(utils::capture.output(print(px), file = logFile, append = TRUE))
    cat("\n", file = logFile, append = TRUE)
    cat(paste0(name, ": nRows = ", nrow(x), ", nCols = ", ncol(x), "\n"), file = logFile, append = TRUE)
    
  }else if(is.data.frame(x)){
    
    cat(paste0(name, ": nRows = ", nrow(x), ", nCols = ", ncol(x), "\n"), file = logFile, append = TRUE)
    suppressMessages(utils::capture.output(print(head(x)), file = logFile, append = TRUE))
    cat("\n", file = logFile, append = TRUE)
    
  }else if(is(x, "dgCMatrix") | is(x, "dgeMatrix")){
    
    cat(paste0(name, ": nRows = ", nrow(x), ", nCols = ", ncol(x), "\n"), file = logFile, append = TRUE)
    cat(paste0(name, ": NonZeroEntries = ", length(x@x), ", EntryRange = [ ", paste0(range(x@x), collapse=" , "), " ]\n"), file = logFile, append = TRUE)    
    px <- x[head(seq_len(nrow(x)), 5), head(seq_len(ncol(x)), 5), drop = FALSE]
    suppressMessages(utils::capture.output(print(px), file = logFile, append = TRUE))
    cat("\n", file = logFile, append = TRUE)
    
  }else if(is(x, "GRanges")){
    
    suppressMessages(utils::capture.output(print(x), file = logFile, append = TRUE))
    
  }else if(is(x, "SummarizedExperiment")){
    
    suppressMessages(utils::capture.output(print(x), file = logFile, append = TRUE))
    
  }else if(is(x, "DataFrame")){
    
    suppressMessages(utils::capture.output(print(x), file = logFile, append = TRUE))
    
  }else if(is(x, "ArchRProj")){
    
    suppressMessages(utils::capture.output(print(proj), file = logFile, append = TRUE))
    
  }else if(is(x, "SimpleList") | is(x, "list")){
    
    for(i in seq_along(x)){
      
      y <- x[[i]]
      
      if(missing(y)){
        next
      }
      
      if(is.matrix(y)){
        
        cat("\n", file = logFile, append = TRUE)
        cat(paste0(paste0(name,"$", names(x[i])), ": nRows = ", nrow(y), ", nCols = ", ncol(y), "\n"), file = logFile, append = TRUE)        
        px <- y[head(seq_len(nrow(y)), 5), head(seq_len(ncol(y)), 5), drop = FALSE]
        suppressMessages(utils::capture.output(print(px), file = logFile, append = TRUE))
        cat("\n", file = logFile, append = TRUE)
        
      }else if(is.data.frame(y)){
        
        cat("\n", file = logFile, append = TRUE)
        cat(paste0(paste0(name,"$", names(x[i])), ": nRows = ", nrow(y), ", nCols = ", ncol(y), "\n"), file = logFile, append = TRUE)
        suppressMessages(utils::capture.output(print(head(y)), file = logFile, append = TRUE))
        cat("\n", file = logFile, append = TRUE)
        
      }else if(is(y, "dgCMatrix")){
        
        cat("\n", file = logFile, append = TRUE)
        cat(paste0(paste0(name,"$", names(x[i])), ": nRows = ", nrow(y), ", nCols = ", ncol(y), "\n"), file = logFile, append = TRUE)
        cat(paste0(paste0(name,"$", names(x[i])), ": NonZeroEntries = ", length(y@x), ", EntryRange = [ ", paste0(range(y@x), collapse=" , "), " ]\n"), file = logFile, append = TRUE)        
        px <- y[head(seq_len(nrow(y)), 5), head(seq_len(ncol(y)), 5), drop = FALSE]
        suppressMessages(utils::capture.output(print(px), file = logFile, append = TRUE))
        cat("\n", file = logFile, append = TRUE)
        
      }else if(is(y, "SimpleList") | is(y, "list")){
        
        for(j in seq_along(y)){
          
          z <- y[[j]]
          
          if(missing(z)){
            next
          }
          
          if(is.matrix(z)){
            
            cat("\n", file = logFile, append = TRUE)
            cat(paste0(paste0(name,"$", names(y[j])), ": nRows = ", nrow(z), ", nCols = ", ncol(z), "\n"), file = logFile, append = TRUE)            
            px <- z[head(seq_len(nrow(z)), 5), head(seq_len(ncol(z)), 5), drop = FALSE]
            suppressMessages(utils::capture.output(print(px), file = logFile, append = TRUE))
            cat("\n", file = logFile, append = TRUE)
            
          }else if(is.data.frame(z)){
            
            cat("\n", file = logFile, append = TRUE)
            cat(paste0(paste0(name,"$", names(y[j])), ": nRows = ", nrow(z), ", nCols = ", ncol(z), "\n"), file = logFile, append = TRUE)
            suppressMessages(utils::capture.output(print(head(z)), file = logFile, append = TRUE))
            cat("\n", file = logFile, append = TRUE)
            
          }else if(is(z, "dgCMatrix")){
            
            cat("\n", file = logFile, append = TRUE)
            cat(paste0(paste0(name,"$", names(y[j])), ": nRows = ", nrow(z), ", nCols = ", ncol(z), "\n"), file = logFile, append = TRUE)
            cat(paste0(paste0(name,"$", names(y[j])), ": NonZeroEntries = ", length(z@x), ", EntrzRange = [ ", paste0(range(z@x), collapse=" , "), " ]\n"), file = logFile, append = TRUE)            
            px <- z[head(seq_len(nrow(z)), 5), head(seq_len(ncol(z)), 5), drop = FALSE]
            suppressMessages(utils::capture.output(print(px), file = logFile, append = TRUE))
            cat("\n", file = logFile, append = TRUE)
            
          }else if(is(y, "SimpleList") | is(y, "list")){
            
            #Only print 2x nested lists
            
          }else{
            
            tryCatch({
              cat("\n", file = logFile, append = TRUE)
              cat(paste0(paste0(name,"$", names(y[j])), ": length = ", length(z), "\n"), file = logFile, append = TRUE)
              suppressMessages(utils::capture.output(print(head(z)), file = logFile, append = TRUE))
              cat("\n", file = logFile, append = TRUE)
            }, error = function(q){
            })
            
          }
          
        }
        
      }else{
        
        tryCatch({
          cat("\n", file = logFile, append = TRUE)
          cat(paste0(paste0(name,"$", names(x[i])), ": length = ", length(y), "\n"), file = logFile, append = TRUE)
          suppressMessages(utils::capture.output(print(head(y)), file = logFile, append = TRUE))
          cat("\n", file = logFile, append = TRUE)
        }, error = function(q){
        })
        
      }
      
    }
    
  }else{
    
    tryCatch({
      cat("\n", file = logFile, append = TRUE)
      cat(paste0(name, ": length = ", length(x), "\n"), file = logFile, append = TRUE)
      suppressMessages(utils::capture.output(print(head(x)), file = logFile, append = TRUE))
      cat("\n", file = logFile, append = TRUE)
    }, error = function(q){
    })
    
  }
  
  cat("\n", file = logFile, append = TRUE)
  return(invisible(0))
  
}

.endLogging <- function(
  logFile = NULL, 
  useLogs = getArchRLogging()
){
  
  if(!useLogs){
    return(invisible(0))
  }
  
  if(is.null(logFile)){
    return(invisible(0))
  }
  
  rL <- readLines(logFile)
  o <- tryCatch({
    t1 <- gsub("Start Time : ","", grep("Start Time", rL, ignore.case = TRUE, value = TRUE))
    mn <- as.numeric(difftime(Sys.time(), as.POSIXct(t1), units = "mins"))
    hr <- as.numeric(difftime(Sys.time(), as.POSIXct(t1), units = "hours"))
    cat("\n------- Completed\n\n", file = logFile, append = TRUE)
    cat(paste0("End Time : ",Sys.time(),"\n"), file = logFile, append = TRUE)
    cat(paste0("Elapsed Time Minutes = ", mn), file = logFile, append = TRUE)
    cat(paste0("\nElapsed Time Hours = ", hr), file = logFile, append = TRUE)
    cat("\n\n-------\n\n\n\n", file = logFile, append = TRUE)
    message("ArchR logging successful to : ", logFile)
  }, error = function(x){
  })
  
  # tryCatch({
  #   R.utils::gzip(logFile, paste0(logFile, ".gz"))
  #   message("ArchR logging successful to : ", paste0(logFile, ".gz"))
  # }, error = function(x){
  # })
  
  return(invisible(0))
  
}


##########################################################################################
# Validation Methods
##########################################################################################

.validInput <- function(input = NULL, name = NULL, valid = NULL){
  
  valid <- unique(valid)
  
  if(is.character(valid)){
    valid <- tolower(valid)
  }else{
    stop("Validator must be a character!")
  }
  
  if(!is.character(name)){
    stop("name must be a character!")
  }
  
  if("null" %in% tolower(valid)){
    valid <- c("null", valid[which(tolower(valid) != "null")])
  }
  
  av <- FALSE
  
  for(i in seq_along(valid)){
    
    vi <- valid[i]
    
    if(vi == "integer" | vi == "wholenumber"){
      
      if(all(is.numeric(input))){
        #https://stackoverflow.com/questions/3476782/check-if-the-number-is-integer
        cv <- min(abs(c(input%%1, input%%1-1)), na.rm = TRUE) < .Machine$double.eps^0.5
      }else{
        cv <- FALSE
      }
      
    }else if(vi == "null"){
      
      cv <- is.null(input)
      
    }else if(vi == "bool" | vi == "boolean" | vi == "logical"){
      
      cv <- is.logical(input)
      
    }else if(vi == "numeric"){
      
      cv <- is.numeric(input)
      
    }else if(vi == "vector"){
      
      cv <- is.vector(input)
      
    }else if(vi == "matrix"){
      
      cv <- is.matrix(input)
      
    }else if(vi == "sparsematrix"){
      
      cv <- is(input, "dgCMatrix")
      
    }else if(vi == "character"){
      
      cv <- is.character(input)
      
    }else if(vi == "factor"){
      
      cv <- is.factor(input)
      
    }else if(vi == "rlecharacter"){
      
      cv1 <- is(input, "Rle")
      if(cv1){
        cv <- is(input@values, "factor") || is(input@values, "character")
      }else{
        cv <- FALSE
      }
      
    }else if(vi == "palette"){
      
      cv <- all(.isColor(input))
      
    }else if(vi == "timestamp"){
      
      cv <- is(input, "POSIXct")
      
    }else if(vi == "dataframe" | vi == "data.frame" | vi == "df"){
      
      cv1 <- is.data.frame(input)
      cv2 <- is(input, "DataFrame")
      cv <- any(cv1, cv2)
      
    }else if(vi == "fileexists"){
      
      cv <- all(file.exists(input))
      
    }else if(vi == "direxists"){
      
      cv <- all(dir.exists(input))
      
    }else if(vi == "granges" | vi == "gr"){
      
      cv <- is(input, "GRanges")
      
    }else if(vi == "grangeslist" | vi == "grlist"){
      
      cv <- .isGRList(input)
      
    }else if(vi == "list" | vi == "simplelist"){
      
      cv1 <- is.list(input)
      cv2 <- is(input, "SimpleList")
      cv <- any(cv1, cv2)
      
    }else if(vi == "bsgenome"){
      
      cv1 <- is(input, "BSgenome")
      cv2 <- tryCatch({
        library(input)
        eval(parse(text=input))
      }, error = function(e){
        FALSE
      })
      cv <- any(cv1, cv2)
      
    }else if(vi == "se" | vi == "summarizedexperiment"){
      
      cv <- is(input, "SummarizedExperiment")
      
    }else if(vi == "seurat" | vi == "seuratobject"){
      
      cv <- is(input, "Seurat")
      
    }else if(vi == "txdb"){
      
      cv <- is(input, "TxDb")
      
    }else if(vi == "orgdb"){
      
      cv <- is(input, "OrgDb")
      
    }else if(vi == "bsgenome"){
      
      cv <- is(input, "BSgenome")
      
    }else if(vi == "parallelparam"){
      
      cv <- is(input, "BatchtoolsParam")
      
    }else if(vi == "archrproj" | vi == "archrproject"){
      
      cv <- is(input, "ArchRProject")
      ###validObject(input) check this doesnt break anything if we
      ###add it. Useful to make sure all ArrowFiles exist! QQQ
      
    }else{
      
      stop("Validator is not currently supported by ArchR!")
      
    }
    
    if(cv){
      av <- TRUE
      break
    }   
    
  }
  
  if(av){
    
    return(invisible(TRUE))
    
  }else{
    
    stop("Input value for '", name,"' is not a ", paste(valid, collapse="," ), ", (",name," = ",class(input),") please supply valid input!")
    
  }
  
}


#
.tempfile <- function(pattern = "tmp", tmpdir = "tmp", fileext = "", addDOC = TRUE){
  
  dir.create(tmpdir, showWarnings = FALSE)
  
  if(addDOC){
    doc <- paste0("-Date-", Sys.Date(), "_Time-", gsub(":","-", stringr::str_split(Sys.time(), pattern=" ",simplify=TRUE)[1,2]))
  }else{
    doc <- ""
  }
  
  tempfile(pattern = paste0(pattern, "-"), tmpdir = tmpdir, fileext = paste0(doc, fileext))
  
}


#
.requirePackage <- function(x = NULL, load = TRUE, installInfo = NULL, source = NULL){
  if(x %in% rownames(installed.packages())){
    if(load){
      suppressPackageStartupMessages(require(x, character.only = TRUE))
    }else{
      return(0)
    }
  }else{
    if(!is.null(source) & is.null(installInfo)){
      if(tolower(source) == "cran"){
        installInfo <- paste0('install.packages("',x,'")')
      }else if(tolower(source) == "bioc"){
        installInfo <- paste0('BiocManager::install("',x,'")')
      }else{
        stop("Unrecognized package source, available are cran/bioc!")
      }
    }
    if(!is.null(installInfo)){
      stop(paste0("Required package : ", x, " is not installed/found!\n  Package Can Be Installed : ", installInfo))
    }else{
      stop(paste0("Required package : ", x, " is not installed/found!"))
    }
  }
}


#
.ArchRLogo <- function(ascii = "Logo", messageLogo = TRUE){
  Ascii <- list(
    Package = c("
           ___      .______        ______  __    __  .______      
          /   \\\     |   _  \\\      /      ||  |  |  | |   _  \\\     
         /  ^  \\\    |  |_)  |    |  ,----'|  |__|  | |  |_)  |    
        /  /_\\\  \\\   |      /     |  |     |   __   | |      /     
       /  _____  \\\  |  |\\\  \\\\___ |  `----.|  |  |  | |  |\\\  \\\\___.
      /__/     \\__\\ | _| `._____| \\______||__|  |__| | _| `._____|
    "),
    
    #modified from cyu@athena.mit.edu
    Logo = c("
                                                   / |
                                                 /    \\\
            .                                  /      |.
            \\\\\\                              /        |.
              \\\\\\                          /           `|.
                \\\\\\                      /              |.
                  \\\                    /                |\\\
                  \\\\#####\\\           /                  ||
                ==###########>      /                   ||
                 \\\\##==......\\\    /                     ||
            ______ =       =|__ /__                     ||      \\\\\\\
        ,--' ,----`-,__ ___/'  --,-`-===================##========>
       \\\               '        ##_______ _____ ,--,__,=##,__   ///
        ,    __==    ___,-,__,--'#'  ==='      `-'    | ##,-/
        -,____,---'       \\\\####\\\\________________,--\\\\_##,/
           ___      .______        ______  __    __  .______      
          /   \\\     |   _  \\\      /      ||  |  |  | |   _  \\\     
         /  ^  \\\    |  |_)  |    |  ,----'|  |__|  | |  |_)  |    
        /  /_\\\  \\\   |      /     |  |     |   __   | |      /     
       /  _____  \\\  |  |\\\  \\\\___ |  `----.|  |  |  | |  |\\\  \\\\___.
      /__/     \\__\\ | _| `._____| \\______||__|  |__| | _| `._____|
    ")
  )
  
  if(messageLogo){
    message(Ascii[[ascii]])
  }else{
    Ascii[[ascii]]
  }
  
}