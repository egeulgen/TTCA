#' TTCA: Transcript Time Course Analysis
#'
#' Background: The analysis of microarray time series promises a deeper insight into the dynamics of the cellular response following stimulation. A common observation in this type of data is that some genes respond with quick, transient dynamics, while other genes change their expression slowly over time. The existing methods for detecting significant expression dynamics often fail when the expression dynamics show a large heterogeneity. Moreover, these methods often cannot cope with irregular and sparse measurements.
#' Results: The method proposed here is specifically designed for the analysis of perturbation responses. It combines different scores to capture fast and transient dynamics as well as slow expression changes, and performs well in the presence of low replicate numbers and irregular sampling times. The results are given in the form of tables including links to figures showing the expression dynamics of the respective transcript. These allow to quickly recognise the relevance of detection, to identify possible false positives and to discriminate early and late changes in gene expression. An extension of the method allows the analysis of the expression dynamics of functional groups of genes, providing a quick overview of the cellular response. The performance of this package was tested on microarray data derived from lung cancer cells stimulated with epidermal growth factor (EGF).
#' Paper: Albrecht, Marco, et al. (2017)<DOI:10.1186/s12859-016-1440-8>.
#'
#' The package has not be applied to Hi-Seq data yet. The problem is the huge variety in the read counts. An additional transformation of normalized Hi-Seq data might be an option to scale the values between two values like 0 and 1 (Simple idea: Datalog<-log(data, base = max(data))). Not tested. IF you are interested to adjust my package to sequence data, feel free to contact me.
#'
#' @param grp1        Data set with longitudinal sampled data (data.frame)
#' @param grp2        Data set with longitudinal sampled data for comparison (data.frame)
#' @param grp1.time   Time points for data set 1 (vector like: c(0,0,0.5,1,2,4,6,8,12,12)
#' @param grp2.time   Time points for data set 2 (vector like: c(0,0,0.5,3,2,4,6,8,12,12,24)
#' @param lambda      Smoothing parameter in penalty term of quantil regression (default: lambda=0.6 ). Adjust, if fit is too strict or too flexible.
#' @param annot       Annotation for pictures and result (Data.frame with 2 columns with ID and GeneName). (Default: annot=NA)
#' @param pVal        P-value for the local hypothesis test        (default: 0.05).
#' @param codetest    Reduces the data set to 200 features for a quick run of the program. (default: codetest=FALSE)
#' @param PeakMode    Peakmode "norm" uses variance between replicates. If changed to another character value, a normal hypothesis test will be conducted (default: PeakMode="norm")
#' @param timeInt     Defines early, middle and late time period. Defines the middle time period between 4 h and 12 h with timeInt<-c(4,12). (default: timeInt=NULL)
#' @param file        Result folder will be saved at this location   (default: file=getwd() ).
#' @param MaxPics     Limits the number of plots (default: MaxPics=10000)
#' @param Stimulus1   Searches this term together with the gene name in PubMed. Stimulus1="Insulin+like+growth+factor"   ( default: Stimulus1="")
#' @param Stimulus2   Searches this term together with the gene name in PubMed. Stimulus2="epidermal+growth+factor"      ( default: Stimulus2="")
#' @param S           Defines mode. S =="GO" changes programm to gene ontology mode  (default:  S="gene")
#' @param mapGO       Link genes to Gene Ontology terms  (default: mapGO="")
#' @param annotation  Merges the TTCA by rowname with a table of your wish. Example: annotation<-annotation[,c("probeset_id", "gene_name","transkript_id","GO_BP","GO_CC","GO_mf")]  (default: annotation="annotation")
#'
#' @return The R-package delivers a table with different significance values, rankings, p-values. Moreover, it will plot the most important time courses and quality control images.
#'
#' @examples
#' \dontrun{
#'
#' ##########################################
#' #### Gene-ANALYSE
#' ##########################################
#' require(quantreg);require(VennDiagram);require(tcltk2); require(tcltk);
#' require(RISmed);require(Matrix)
#' data(EGF,Control,annot,annotation)
#'
#' S="gene"
#' Control.time <-  c(0,0,0.5,1,4,6,24,24,48,48,48)
#' EGF.time     <-  c(0,0,0.5,0.5,1,2,4,6,8,12,18,24,24,48,48,48)
#' file         =   paste0(getwd(),"/TTCA_Gene")
#' dir.create(file)
#' ######
#' TTCAresult<-TTCA(grp1=EGF, grp1.time=EGF.time, grp2=Control, grp2.time=Control.time,S="gene",
#'                  lambda=0.6, annot=annot, annotation=annotation,pVal=0.05,codetest=FALSE,
#'                  file=file, Stimulus1="epidermal+growth+factor", timeInt=c(4,12), MaxPics =10000)
#'}
#'
#'
#'
#'
#' \dontrun{
#' ##########################################
#' #### GO-ANALYSE
#' ##########################################
#' require(quantreg);require(VennDiagram);require(tcltk2); require(tcltk);
#' require(RISmed);require(Matrix)
#' #source("https://bioconductor.org/biocLite.R")
#' #biocLite("biomaRt")
#' library(biomaRt)
#' data(EGF,Control,annot,annotation)
#'
#' require(biomaRt)
#' ensembl <-  useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")
#' mapGO <- getBM(attributes=c("go_id","name_1006",'affy_hugene_2_0_st_v1'),
#'                filters = 'affy_hugene_2_0_st_v1', values=rownames(annot), mart =ensembl)
#' colnames(mapGO)<-c("go_id","GO_name","probeset_id")
#'
#' S="GO"
#' Control.time <-  c(0,0,0.5,1,4,6,24,24,48,48,48)
#' EGF.time     <-  c(0,0,0.5,0.5,1,2,4,6,8,12,18,24,24,48,48,48)
#' file         =   paste0(getwd(),"/TTCA_GO")
#' dir.create(file)
#'
#' TTCAresult<-TTCA(grp1=EGF, grp1.time=EGF.time, grp2=Control, grp2.time=Control.time,
#'                  S="GO", pVal=0.05,lambda=0.6,codetest=FALSE, file=file,
#'                  Stimulus1="epidermal+growth+factor", timeInt=c(4,12),
#'                  MaxPics=10000, mapGO=mapGO)
#' }
#'
#' @importFrom grDevices dev.off
#' @importFrom grDevices png
#' @importFrom grDevices rgb
#'
#' @importFrom graphics abline
#' @importFrom graphics axis
#' @importFrom graphics hist
#' @importFrom graphics legend
#' @importFrom graphics lines
#' @importFrom graphics par
#' @importFrom graphics plot
#' @importFrom graphics points
#' @importFrom graphics polygon
#' @importFrom graphics title
#'
#' @importFrom stats dcauchy
#' @importFrom stats density
#' @importFrom stats dgamma
#' @importFrom stats dlnorm
#' @importFrom stats dlogis
#' @importFrom stats dnorm
#' @importFrom stats dweibull
#' @importFrom stats median
#' @importFrom stats optimize
#' @importFrom stats pcauchy
#' @importFrom stats pgamma
#' @importFrom stats plnorm
#' @importFrom stats plogis
#' @importFrom stats pnorm
#' @importFrom stats pweibull
#' @importFrom stats quantile
#' @importFrom stats sd
#'
#' @importFrom utils sessionInfo
#' @importFrom utils write.csv
#' @importFrom utils write.table
#'
#'
#' @import MASS
#' @import Matrix
#' @import RISmed
#' @import quantreg
#' @import VennDiagram
#' @import tcltk2
#' @import tcltk
#' @export
TTCA = function (grp1, grp1.time, grp2, grp2.time, lambda = 0.6, annot = NA, 
                 annotation = "annotation", timeInt = NULL, pVal = 0.05, codetest = FALSE, 
                 file = getwd(), MaxPics = 10000, Stimulus1 = "", Stimulus2 = "", 
                 S = "gene", mapGO = "", PeakMode = "norm") 
{
  if (is.data.frame(annot) & !identical(rownames(annot), rownames(grp1)) & 
      !identical(rownames(grp1), rownames(grp2))) {
    stop("Annot must have the same rownames like grp1 and grp2! Not identical rowname vector.")
  }
  if (!is.data.frame(annot) & S != "GO") {
    stop("You forgot the annot file. You can neglect this in GO mode.")
  }
  if (is.null(colnames(mapGO)[1]) & S == "GO") {
    stop("You need a mapGO file. See example in Documentation (?TTCA).")
  }
  print("--------------------------------------------------------------------------------------------------")
  print(paste0("                                                                                                  "))
  print(paste0("                                                                                                  "))
  print(paste0("--                                                                                              --"))
  print(paste0("-----                                                                                        -----"))
  print(paste0("--------                                WELCOME TO TTCA                                   --------"))
  print(paste0("-----                                                                                        -----"))
  print(paste0("--                              transcript time course analysis                                 --"))
  print(paste0("                                                                                                  "))
  print(paste0("                                                                                                  "))
  print(paste0("                                                                                                  "))
  print("--------------------------------------------------------------------------------------------------")
  Funtext <- sample(1:10, 1)
  if (Funtext == 1) {
    print(paste0(format(Sys.time(), "%H:%M:%S"), ": The work begins. Make a coffee or tea break. You deserve it ;D.                           "))
  }
  if (Funtext == 2) {
    print(paste0(format(Sys.time(), "%H:%M:%S"), ": Perfect. TTCA works like crazy.                                                          "))
  }
  if (Funtext == 3) {
    print(paste0(format(Sys.time(), "%H:%M:%S"), ": TTCA begins to enslave your computer. Let us burn the processor.                        "))
  }
  if (Funtext == 4) {
    print(paste0(format(Sys.time(), "%H:%M:%S"), ": TCCA works. Time for day dreaming.                                                      "))
  }
  if (Funtext == 5) {
    print(paste0(format(Sys.time(), "%H:%M:%S"), ": Absolute fantastic. You started TTCA!                                                   "))
  }
  if (Funtext == 6) {
    print(paste0(format(Sys.time(), "%H:%M:%S"), ": TTCA cares now about your data. Do not worry, everything will be fine.                  "))
  }
  if (Funtext == 7) {
    print(paste0(format(Sys.time(), "%H:%M:%S"), ": TTCA does something nice for you. Say something nice to a collegue.                     "))
  }
  if (Funtext == 8) {
    print(paste0(format(Sys.time(), "%H:%M:%S"), ": Aye, aye sir. TTCA has leaved the harbour and comes hopefully back with                 "))
    print(paste0(format(Sys.time(), "%H:%M:%S"), ": magnificient data. Pray to God that no pirates attack us with bugs.                     "))
  }
  if (Funtext == 9) {
    print(paste0(format(Sys.time(), "%H:%M:%S"), ": Aye, aye sir. TTCA has leaved the harbour and comes hopefully back with                 "))
    print(paste0(format(Sys.time(), "%H:%M:%S"), ": magnificient data. Pray to God that no pirates attack us with bugs.                     "))
  }
  if (Funtext == 10) {
    print(paste0(format(Sys.time(), "%H:%M:%S"), ": TTCA has started to wash your data. Hopefully, we will obtain clear data.               "))
  }
  print(paste0("                                                                                                  "))
  Sys.setlocale("LC_TIME", "C")
  Sys.setlocale("LC_ALL", "English")
  Sys.setlocale("LC_TIME", "English")
  print(paste0(format(Sys.time(), "%H:%M:%S"), ": System locale is temporary setted to international standard with Sys.setlocale"))
  print(paste0("          Original: ", strsplit(Sys.setlocale(category = "LC_TIME", 
                                                              locale = ""), ";")[[1]]))
  grp1n <- deparse(substitute(grp1))
  grp2n <- deparse(substitute(grp2))
  startTimeFunction <- Sys.time()
  global <- tkProgressBar(title = paste0("global progress bar ", 
                                         grp1n, " vs. ", grp2n, " (Start time: ", format(startTimeFunction, 
                                                                                         "%H:%M"), ")"), min = 0, max = 100, width = 1000)
  setTkProgressBar(global, 1, label = paste("Creates new folder"))
  folder1 <- paste0(file, "/contrast_", grp1n, "-vs-", grp2n)
  if (exists("folder")) {
    rm(folder)
  }
  for (i in 1:100) {
    folder2 <- paste0(folder1, "_version_", i)
    if (!file.exists(folder2) & !exists("folder")) {
      dir.create(folder2)
      folder <- folder2
    }
  }
  dir.create(paste0(folder, "/Information"))
  dir.create(paste0(folder, "/Quality_Control"))
  dir.create(paste0(folder, "/Significance_And_Pval"))
  dir.create(paste0(folder, "/Saved_Intermediate_Results"))
  dir.create(paste0(folder, "/Significance_And_Pval/PeakScore"))
  dir.create(paste0(folder, "/Venn_and_other_characteristics"))
  dir.create(paste0(folder, "/Resultfigure"))
  rm(folder1, folder2, file)
  Infotext <- list()
  Infotext$header = "Thanks for using method TTCA to extract genes with strongest activity over time."
  Infotext$starttime <- startTimeFunction
  Infotext$grp1.time <- grp1.time
  Infotext$grp2.time <- grp2.time
  Infotext$grp1n <- grp1n
  Infotext$grp2n <- grp2n
  Infotext$lambda <- lambda
  Infotext$pVal <- pVal
  Infotext$codetest = codetest
  Infotext$Stimulus1 = Stimulus1
  Infotext$Stimulus2 = Stimulus2
  Infotext$MaxPics = MaxPics
  Infotext$timeInt = timeInt
  Infotext$annotsize = dim(annot)
  Infotext$annot = head(annot, 10)
  Infotext$grp1size <- dim(grp1)
  Infotext$grp1 <- head(grp1, 10)
  Infotext$grp2size <- dim(grp2)
  Infotext$grp2 <- head(grp2, 10)
  Infotext$sessionInfo <- sessionInfo()
  sink(file = paste0(folder, "/Information/Function_Input.txt"), 
       append = T)
  print(Infotext)
  sink(file = NULL)
  rm(Infotext)
  setTkProgressBar(global, 2, label = paste("New folders have been created. Removes features that are not annotated or do not have positive values."))
  grp1 <- grp1[, order(grp1.time)]
  grp1.time <- grp1.time[order(grp1.time)]
  grp2 <- grp2[, order(grp2.time)]
  grp2.time <- grp2.time[order(grp2.time)]
  if (S == "GO") {
    print(paste0(format(Sys.time(), "%H:%M:%S"), ": GO-Modus.                                           "))
    GO <- ChangeToGO(grp1 = grp1, grp2 = grp2, grp1.time = grp1.time, 
                     mapGO = mapGO)
    save(GO, file = paste0(folder, "/Saved_Intermediate_Results/RawGene_to_GO.RData"))
    grp1 <- GO$grp1
    grp2 <- GO$grp2
    annot = GO$annot
  }
  beginSize = nrow(grp1)
  grp1.grp2 <- merge(grp1, grp2, by = "row.names")
  rownames(grp1.grp2) <- grp1.grp2[, "Row.names"]
  grp1.grp2 <- grp1.grp2[, -1]
  grp1.grp2.time <- c(grp1.time, grp2.time)
  grp <- c(rep(1, length(grp1.time)), rep(2, length(grp2.time)))
  if (S != "GO") {
    rownames(annot) = annot[, "probeset_id"]
    NotAnnot <- rownames((annot[is.na(annot[, "gene_name"]), 
    ]))
    grp1 <- grp1[!(rownames(grp1) %in% NotAnnot), ]
    grp2 <- grp2[!(rownames(grp2) %in% NotAnnot), ]
    grp1.grp2 <- grp1.grp2[!(rownames(grp1.grp2) %in% NotAnnot), 
    ]
    annot <- annot[!(rownames(annot) %in% NotAnnot), ]
    print(paste0(format(Sys.time(), "%H:%M:%S"), ": ", beginSize - 
                   dim(grp1)[1], " not annotated features have been removed (absolute ", 
                 round((1 - dim(grp1)[1]/beginSize) * 100, 1), " % reduction)"))
  }
  info <- nrow(grp1)
  Z1 <- apply(grp1.grp2, 1, max)
  if (min(Z1, na.rm = TRUE) < 0) {
    Z2 <- names(Z1[Z1 < 0])
    grp1 <- grp1[!(rownames(grp1) %in% Z2), ]
    grp2 <- grp2[!(rownames(grp2) %in% Z2), ]
    grp1.grp2 <- grp1.grp2[!(rownames(grp1.grp2) %in% Z2), 
    ]
    print(paste0(format(Sys.time(), "%H:%M:%S"), ": ", info - 
                   dim(grp1)[1], " features have been removed because they have no positive intensity values (absolute ", 
                 round((1 - dim(grp1)[1]/beginSize) * 100, 1), " % reduction)"))
    rm(Z2)
  }
  rm(Z1, info)
  if (codetest == TRUE) {
    grp1.grp2 <- grp1.grp2[sort(sample(1:nrow(grp1.grp2), 
                                       1100)), ]
    MaxPics = 10
    setTkProgressBar(global, 2.5, label = paste("Code testing modus. Works with 1100 features only"))
    print(paste0(format(Sys.time(), "%H:%M:%S"), ": Code testing modus. Works with 1100 features only                                       "))
  }
  setTkProgressBar(global, 3, label = paste("Non-annotated features were removed.Tests for significant distances between the dynamics"))
  Sumreplicates <- sum(duplicated(grp1.time)) + sum(duplicated(grp2.time))
  if (Sumreplicates > 4 & PeakMode == "norm") {
    RESULT <- PeakScore(dat = grp1.grp2, grp = grp, tme = grp1.grp2.time, 
                        folder = folder, global = global)
  }
  else {
    RESULT <- PeakScoreHOH1(dat = grp1.grp2, grp = grp, tme = grp1.grp2.time, 
                            folder = folder, global = global)
  }
  if (nrow(RESULT[RESULT[, "MaxDist"] == 0, ]) > 0) {
    identDynamics <- rownames(RESULT[RESULT[, "MaxDist"] == 
                                       0, ])
    RESULT = RESULT[rownames(RESULT) != identDynamics, ]
    grp1.grp2 = grp1.grp2[rownames(grp1.grp2) != identDynamics, 
    ]
    grp1 = grp1[rownames(grp1) != identDynamics, ]
    grp2 = grp2[rownames(grp2) != identDynamics, ]
  }
  tresP <- DistplotPval(x1 <- RESULT[, "MaxDist", drop = FALSE], 
                        main = "Hypothesistest for peak score", xlab = "Maximum Distance", 
                        pVal = pVal, folder = paste0(folder, "/Significance_And_Pval/"), 
                        ForTest <- "Significant_Peak_Score_alternativ")
  RESULT[, "SignPeakH0"] <- RESULT[, "MaxDist"] > tresP$tres
  RESULT <- merge(RESULT, tresP$Pval, by = "row.names")
  rownames(RESULT) <- RESULT[, 1]
  RESULT <- RESULT[, -1]
  setTkProgressBar(global, 9.25, label = paste("Found significant distances between the dynamics.Program checks whether the dynamics are significant"))
  # save(RESULT, file = paste0(folder, "/Saved_Intermediate_Results/RESULT_afterPeakScore.RData"))
  if (S == "GO") {
    setTkProgressBar(global, 9.25, label = paste("Found significant distances between the dynamics.Program calculates reverse overlap score now"))
    ROS <- revOverlapScore(go = GO, tm1 = grp1.time, tm2 = grp2.time, 
                           folder = folder)
    RESULT <- merge(RESULT, ROS, by = "row.names")
    rownames(RESULT) <- RESULT[, 1]
    RESULT <- RESULT[, -1]
    tresRos <- DistplotPval(x1 <- RESULT[, "RevOverlapScore", 
                                         drop = FALSE], main = "Reverse Overlap Score", xlab = "overlap", 
                            pVal = pVal, folder = paste0(folder, "/Significance_And_Pval/"), 
                            ForTest <- "Reverse_Overlap_Score")
    RESULT[, "SignROS"] <- RESULT[, "RevOverlapScore"] > 
      tresRos$tres
    RESULT <- merge(RESULT, tresRos$Pval, by = "row.names")
    rownames(RESULT) <- RESULT[, 1]
    RESULT <- RESULT[, -1]
    RESULT <- merge(RESULT, GO$groupsize, by = "row.names")
    rownames(RESULT) <- RESULT[, 1]
    RESULT <- RESULT[, -1]
    setTkProgressBar(global, 9.75, label = paste("Calculation of reverse overlap score is done. Program checks whether the dynamics are significant"))
    # save(RESULT, file = paste0(folder, "/Saved_Intermediate_Results/RESULT_after_ROS_and_groupsize.RData"))
  }
  dynScor <- tkProgressBar(title = paste0("progress bar (Start time: ", 
                                          format(Sys.time(), "%H:%M"), ")"), min = 0, max = 100, 
                           width = 1000)
  setTkProgressBar(dynScor, 0, label = paste("Calculating time for dynamic test."))
  timeusage <- proc.time()
  FF1time <- FF.median(dat = grp1.grp2[1:50, ], grp = grp, 
                       tme = grp1.grp2.time)
  tt <- proc.time() - timeusage
  timePerFeature <- sum(as.numeric(tt)[1:2])/50
  timeforSample <- timePerFeature * dim(grp1.grp2)[1]
  setTkProgressBar(global, 10, label = paste0("Program checks whether the dynamics are significant. This step is expected to be completed by ", 
                                              ProgTime(timeforSample), "."))
  setTkProgressBar(dynScor, 1, label = paste0("Program checks whether the dynamics are significant."))
  rm(FF1time, timeusage, timeforSample, timePerFeature, tt)
  Dynamic <- FF.median(dat = grp1.grp2, grp = grp, tme = grp1.grp2.time, 
                       dynScor = dynScor, global = global)
  SaveDynamic <- Dynamic$FF
  # save(SaveDynamic, file = paste0(folder, "/Saved_Intermediate_Results/RESULT-before-DynamicScore_Dynamic.RData"))
  RESULT[, "Dynamic"] <- as.data.frame(Dynamic$FF)
  tresDyn <- DistplotPval(x1 <- RESULT[, "Dynamic", drop = FALSE], 
                          main = "Null hypothesis versus  alternative hypothesis", 
                          xlab = "(h0:residuals)/(h1:residuals)", pVal = pVal, 
                          folder = paste0(folder, "/Significance_And_Pval/"), ForTest <- "Significant_Dynamic_Score")
  RESULT <- merge(RESULT, tresDyn$Pval, by = "row.names")
  rownames(RESULT) <- RESULT[, 1]
  RESULT <- RESULT[, -1]
  RESULT[, "SignDyn"] <- RESULT[, "Dynamic"] > tresDyn$tres
  # save(RESULT, file = paste0(folder, "/Saved_Intermediate_Results/RESULT-with-DynamicScore.RData"))
  close(dynScor)
  rm(tresDyn, dynScor, SaveDynamic)
  setTkProgressBar(global, 40, label = paste("Found significant dynamics. Program checks whether the integrals between dynamics are significant."))
  if (length(timeInt) != 2) {
    T1 <- (max(grp1.grp2.time) - min(grp1.grp2.time))/3
    timeInt <- c(min(grp1.grp2.time) + T1, min(grp1.grp2.time) + 
                   2 * T1)
    print(paste0(format(Sys.time(), "%H:%M:%S"), ": Set early response period before ", 
                 round(timeInt[1], 1), "h and late response period after ", 
                 round(timeInt[2], 1), "h."))
  }
  qr1.qr2 <- cbind(Dynamic$qr1, Dynamic$qr2)
  rownames(qr1.qr2) <- rownames(grp1.grp2)
  colnames(qr1.qr2) <- c(unique(grp1.time), unique(grp2.time))
  grpqr1.qr2 <- c(rep(1, ncol(Dynamic$qr1)), rep(2, ncol(Dynamic$qr2)))
  globMaxDist = max(qr1.qr2) - min(qr1.qr2)
  int <- apply(X = qr1.qr2, MARGIN = 1, FUN = Integral, grp = grpqr1.qr2, 
               grp1.time = unique(grp1.time), grp2.time = unique(grp2.time), 
               timeInt = timeInt, globMaxDist = globMaxDist)
  int = t(int)
  colnames(int) <- c("earlyResponse", "middleResponse", "lateResponse", 
                     "completeResponse")
  for (i in 1:4) {
    tresI <- DistplotPval(x1 <- as.data.frame(int[, i, drop = FALSE]), 
                          xlab = "Integral beetween dynamics", main = paste0("Integral for ", 
                                                                             colnames(int)[i], " significant in dynamic"), 
                          pVal = pVal, folder = paste0(folder, "/Significance_And_Pval/"), 
                          ForTest = colnames(int)[i])
    int <- as.data.frame(int)
    if (i == 1) {
      int[, "signifearlyRes"] = as.numeric(int[, "earlyResponse"]) > 
        tresI$tres
      int <- merge(int, tresI$Pval, by = "row.names")
    }
    if (i == 2) {
      int[, "signifmiddleRes"] = as.numeric(int[, "middleResponse"]) > 
        tresI$tres
      int <- merge(int, tresI$Pval, by = "row.names")
    }
    if (i == 3) {
      int[, "signiflateRes"] = as.numeric(int[, "lateResponse"]) > 
        tresI$tres
      int <- merge(int, tresI$Pval, by = "row.names")
    }
    if (i == 4) {
      int[, "signifResponse"] = as.numeric(int[, "completeResponse"]) > 
        tresI$tres
      int <- merge(int, tresI$Pval, by = "row.names")
    }
    rownames(int) <- int[, 1]
    int <- int[, -1]
  }
  RESULT <- merge(RESULT, int, by = "row.names")
  rownames(RESULT) <- RESULT[, 1]
  RESULT <- RESULT[, -1]
  # save(RESULT, file = paste0(folder, "/Saved_Intermediate_Results/RESULT-after-IntegralScore.RData"))
  if (S == "GO") {
    subRES <- RESULT[, c("SignPeak", "SignDyn", "signifearlyRes", 
                         "signifmiddleRes", "signiflateRes", "SignROS")]
  }
  else {
    subRES <- RESULT[, c("SignPeak", "SignDyn", "signifearlyRes", 
                         "signifmiddleRes", "signiflateRes")]
  }
  sig <- NULL
  for (i in 1:nrow(RESULT)) {
    sig[i] = sum(as.logical(subRES[i, ])) >= 2
  }
  RESULT <- RESULT[sig, ]
  rm(sig, tresI, subRES)
  setTkProgressBar(global, 45, label = paste("Found integrals between dynamics.", 
                                             nrow(RESULT), " features are in at least one criterion significant. Programm is looking for features in Pubmed "))
  # save(RESULT, file = paste0(folder, "/Saved_Intermediate_Results/RESULT-before-PubMedScore.RData"))
  RESULT[, "PubMed"] = 0
  if (suppressWarnings(tryCatch(QueryCount(EUtilsSummary("TTCA")), 
                                error = function(e) {
                                  -1
                                })) == -1) {
    print(paste0(format(Sys.time(), "%H:%M:%S"), ": No access to webpage. Searching Pubmed for publication count was unsuccessful.           "))
  }
  else {
    process <- round(1:10 * (nrow(RESULT))/10)
    if (Stimulus1 != "") {
      if (Stimulus2 == "") {
        Stimulus <- Stimulus1
      }
      else {
        Stimulus <- paste0(Stimulus1, " AND ", Stimulus2)
      }
    }
    print("--------------------------------------------------------------------------------------------------")
    print(" Hint:                                                                                            ")
    print("                                                                                                  ")
    print(" In order not to overload the E-utility servers, NCBI recommends that users limit large jobs to   ")
    print(" either weekends or between 9:00 PM and 5:00 AM Eastern time during weekdays. Failure to comply   ")
    print(" with this policy may result in an IP address being blocked from accessing NCBI.                  ")
    print("                                                                                                  ")
    original <- Sys.timezone()
    Sys.setenv(TZ = "EST")
    InTime <- any(c((as.numeric(format(Sys.time(), "%H")) < 
                       5), (as.numeric(format(Sys.time(), "%H")) > 20), 
                    (weekdays(Sys.time()) %in% c("Saturday", "Sunday"))))
    if (InTime) {
      InTimeText = "(GOOD)                        "
    }
    else {
      InTimeText = "(IP adress could be blocked)  "
    }
    print(paste0(format(Sys.time(), " At the moment it is %A at %I:%M %p Eastern time "), 
                 InTimeText, "             "))
    Sys.setenv(TZ = original)
    print(format(Sys.time(), " and                 %A at %I:%M %p your time                                              "))
    format(Sys.time(), "%I:%M %p")
    print("                                                                         Source: R-package RISmed ")
    print("--------------------------------------------------------------------------------------------------")
    print("                                                                                                  ")
    SlowDownPubmedCall <- seq(1, nrow(RESULT), 3)
    ti <- proc.time()
    for (i in 1:nrow(RESULT)) {
      RESULT[i, "PubMed"] <- suppressWarnings(tryCatch(QueryCount(EUtilsSummary(paste0("(", 
                                                                                       annot[annot[, "probeset_id"] == rownames(RESULT[i, 
                                                                                       ]), "gene_name"], ") AND ", Stimulus))), error = function(e) {
                                                                                         NA
                                                                                       }))
      if (i %in% process) {
        setTkProgressBar(global, (45 + 2 * which(process == 
                                                   i)), label = paste("Found integrals between dynamics.", 
                                                                      nrow(RESULT), " features are in at least one criterion significant. Programm is looking for features in Pubmed (", 
                                                                      which(process == i) * 10, " % ready)"))
      }
      if (i %in% SlowDownPubmedCall) {
        ti2 <- proc.time() - ti
        Sys.sleep(max(0, 1.001 - as.numeric(ti2[3])))
        ti <- proc.time()
      }
    }
    rm(process, ti, ti2, SlowDownPubmedCall, original, InTime, 
       InTimeText)
  }
  RESULT[, c("PubMedscore")] <- logb(RESULT[, c("PubMed")] + 
                                       1, base = (max(c(RESULT[, c("PubMed")], 1), na.rm = TRUE) + 
                                                    1))
  RESULT[is.na(RESULT[, "PubMedscore"]), "PubMedscore"] = 0
  # save(RESULT, file = paste0(folder, "/Saved_Intermediate_Results/RESULT-after-PubMedScore.RData"))
  setTkProgressBar(global, 65, label = paste("Literature search in PubMed is complete. Calculates the Consensus Score."))
  if (S == "GO") {
    RESULTCS <- RESULT[, c("earlyResponse", "middleResponse", 
                           "lateResponse", "MaxDist", "Dynamic", "RevOverlapScore")]
  }
  else {
    if (all(RESULT[, "PubMedscore"] == RESULT[1, "PubMedscore"])) {
      RESULTCS <- RESULT[, c("earlyResponse", "middleResponse", 
                             "lateResponse", "MaxDist", "Dynamic")]
    }
    else {
      RESULTCS <- RESULT[, c("earlyResponse", "middleResponse", 
                             "lateResponse", "MaxDist", "Dynamic", "PubMedscore")]
    }
  }
  RESULTCS <- apply(RESULTCS, 2, function(x) (x - mean(x))/sd(x[x > 
                                                                  0]))
  CS1 <- as.data.frame(rowMeans(RESULTCS[, c("earlyResponse", 
                                             "middleResponse", "lateResponse")]))
  if (S == "GO") {
    CS2 <- rowMeans(cbind(RESULTCS[, c("MaxDist", "Dynamic", 
                                       "RevOverlapScore")], CS1))
  }
  else {
    if (all(RESULT[, "PubMedscore"] == RESULT[1, "PubMedscore"])) {
      CS2 <- rowMeans(cbind(RESULTCS[, c("MaxDist", "Dynamic")], 
                            CS1))
    }
    else {
      CS2 <- rowMeans(cbind(RESULTCS[, c("MaxDist", "Dynamic", 
                                         "PubMedscore")], CS1))
    }
  }
  setTkProgressBar(global, 67.5, label = paste("The consensus score has been calculated. Testing significance and producing pics."))
  png(filename = paste0(folder, "/Venn_and_other_characteristics/Source_for_Consensus_Score.png"), 
      width = 3.25, height = 3.25, units = "in", res = 1200, 
      pointsize = 6)
  suppressWarnings(tryCatch((plot(density(CS2, na.rm = TRUE), 
                                  xlim = c(min(RESULTCS), -min(RESULTCS)), xlab = paste0("standardized distributions with ", 
                                                                                         nrow(RESULTCS), " elements."), col = "black", main = "Source for Consensus Score")), 
                            error = function(e) {
                              NA
                            }))
  suppressWarnings(tryCatch((lines(density(RESULTCS[, c("Dynamic")], 
                                           na.rm = TRUE), col = "blue", lwd = 2)), error = function(e) {
                                             NA
                                           }))
  suppressWarnings(tryCatch((lines(density(RESULTCS[, c("MaxDist")], 
                                           na.rm = TRUE), col = "red", lwd = 2)), error = function(e) {
                                             NA
                                           }))
  suppressWarnings(tryCatch((lines(density(rowMeans(RESULTCS[, 
                                                             c("earlyResponse", "middleResponse", "lateResponse")]), 
                                           na.rm = TRUE), col = "green", lwd = 2)), error = function(e) {
                                             NA
                                           }))
  suppressWarnings(tryCatch((if (S == "GO") {
    lines(density(RESULTCS[, c("RevOverlapScore")], na.rm = TRUE), 
          col = "pink", lwd = 2)
  }
  else {
    lines(density(RESULTCS[, c("PubMedscore")], na.rm = TRUE), 
          col = "pink", lwd = 2)
  }), error = function(e) {
    NA
  }))
  suppressWarnings(tryCatch((lines(density(CS2, na.rm = TRUE), 
                                   col = "black", lwd = 2)), error = function(e) {
                                     NA
                                   }))
  suppressWarnings(tryCatch((if (S == "GO") {
    legend("topright", inset = 0.05, c("Consensus Score", 
                                       "Rev Overlap Score"), fill = c("black", "pink"), 
           horiz = FALSE, bty = "n")
  }
  else {
    legend("topright", inset = 0.05, c("Consensus Score", 
                                       "Pubmed Score"), fill = c("black", "pink"), horiz = FALSE, 
           bty = "n")
  }), error = function(e) {
    NA
  }))
  suppressWarnings(tryCatch((legend("topleft", inset = 0.05, 
                                    c("Distance Score", "Dynamic Score", "Integral Score"), 
                                    fill = c("red", "blue", "green"), horiz = FALSE, bty = "n")), 
                            error = function(e) {
                              NA
                            }))
  dev.off()
  CS3 <- CS2 - min(CS2)
  if (sum(CS3) == 0) {
    CS3 <- CS3 + 0.01
  }
  RESULT[, "ConsensusScore"] <- CS3/max(CS3, na.rm = TRUE)
  # save(RESULT, file = paste0(folder, "/Saved_Intermediate_Results/RESULT-after-ConsensusScore1.RData"))
  tresC <- DistplotPval(x1 <- RESULT[, "ConsensusScore", drop = FALSE], 
                        xlab = "Averarage Score by standarizied Scores", main = paste0("Significance in Consensus Score"), 
                        pVal = pVal, folder = paste0(folder, "/Significance_And_Pval/"), 
                        ForTest = "Significant_Consensus_Score")
  RESULT <- merge(RESULT, tresC$Pval, by = "row.names")
  rownames(RESULT) <- RESULT[, 1]
  RESULT <- RESULT[, -1]
  RESULT[, "signifConsensus"] = as.numeric(RESULT[, "ConsensusScore"]) > 
    tresC$tres
  rm(CS2, CS3, RESULTCS, tresC)
  RESULT <- merge(RESULT, annot, by = "row.names")
  rownames(RESULT) <- RESULT[, 1]
  RESULT <- RESULT[, -1]
  # save(RESULT, file = paste0(folder, "/Saved_Intermediate_Results/RESULT-after-ConsensusScore2.RData"))
  setTkProgressBar(global, 70, label = paste("The Consensus Score is finished. Produce ", 
                                             min(nrow(RESULT), MaxPics), " pics to show the dynamic."))
  prg <- unique(round(seq(min(nrow(RESULT), MaxPics)/20, min(nrow(RESULT), 
                                                             MaxPics), min(nrow(RESULT), MaxPics)/20)))
  RESULT[, "Link"] = NA
  RESULT <- RESULT[rev(order(RESULT[, "ConsensusScore"])), 
  ]
  timeusage <- proc.time()
  if (S == "GO") {
    grp1SDu <- GO$grp1SDu
    grp2SDu <- GO$grp2SDu
    grp1SDd <- GO$grp1SDd
    grp2SDd <- GO$grp2SDd
    GOsz <- GO$groupsize
  }
  for (j in 1:min(nrow(RESULT), MaxPics)) {
    if (S == "GO") {
      GOsd1u = grp1SDu[rownames(grp1SDu) == RESULT[j, "probeset_id"], 
      ]
      GOsd2u = grp2SDu[rownames(grp2SDu) == RESULT[j, "probeset_id"], 
      ]
      GOsd1d = grp1SDd[rownames(grp1SDd) == RESULT[j, "probeset_id"], 
      ]
      GOsd2d = grp2SDd[rownames(grp2SDd) == RESULT[j, "probeset_id"], 
      ]
      GOszj = GOsz[rownames(GOsz) == RESULT[j, "probeset_id"], 
                   1]
    }
    name = RESULT[j, "gene_name"]
    name1 = gsub("[[:punct:]]", " ", name)
    if (nchar(name) > 40) {
      cex1 = 40/nchar(name) * 1.2
    }
    else {
      cex1 = 1.2
    }
    pos = which(rownames(grp1) == rownames(RESULT[j, ]))
    Poi1 = as.numeric(grp1[pos, ])
    Poi2 = as.numeric(grp2[pos, ])
    if (!is.na(sum(Poi1, Poi2))) {
      tmp_df <- data.frame(Poi1 = as.numeric(Poi1),
                           grp1.time = grp1.time)
      p1 <- rqss(Poi1 ~ qss(grp1.time, lambda = 0.6), tau = 0.5, data = tmp_df)
      p1$coef[2:length(p1$coef)] = p1$coef[2:length(p1$coef)] + 
        p1$coef[1]
      p1 <- p1$coef
      
      tmp_df <- data.frame(Poi2 = as.numeric(Poi2),
                           grp2.time = grp2.time)
      p2 <- rqss(Poi2 ~ qss(grp2.time, lambda = 0.6), tau = 0.5, data = tmp_df)
      p2$coef[2:length(p2$coef)] = p2$coef[2:length(p2$coef)] + 
        p2$coef[1]
      p2 <- p2$coef
      gl <- c(Poi1, Poi2, p1, p2)
      if (S == "GO") {
        tmp_df <- data.frame(response = as.numeric(Poi1 + GOsd1u),
                             grp1.time = grp1.time)
        
        p1u <- rqss(response ~ qss(grp1.time,lambda = 0.6), tau = 0.5, data = tmp_df)
        p1u$coef[2:length(p1u$coef)] = p1u$coef[2:length(p1u$coef)] + 
          p1u$coef[1]
        p1u <- p1u$coef
        
        tmp_df <- data.frame(response = as.numeric(Poi1 - GOsd1d),
                             grp1.time = grp1.time)
        p1d <- rqss(response ~ qss(grp1.time, lambda = 0.6), tau = 0.5, data = tmp_df)
        p1d$coef[2:length(p1d$coef)] = p1d$coef[2:length(p1d$coef)] + 
          p1d$coef[1]
        p1d <- p1d$coef
        
        tmp_df <- data.frame(response = as.numeric(Poi2 + GOsd2u),
                             grp2.time = grp2.time)
        p2u <- rqss(response ~ qss(grp2.time, lambda = 0.6), tau = 0.5, data = tmp_df)
        p2u$coef[2:length(p2u$coef)] = p2u$coef[2:length(p2u$coef)] + 
          p2u$coef[1]
        p2u <- p2u$coef
        
        tmp_df <- data.frame(response = as.numeric(Poi2 - GOsd2d),
                             grp2.time = grp2.time)
        p2d <- rqss(response ~ qss(grp2.time, lambda = 0.6), tau = 0.5, data = tmp_df)
        p2d$coef[2:length(p2d$coef)] = p2d$coef[2:length(p2d$coef)] + 
          p2d$coef[1]
        p2d <- p2d$coef
        gl <- c(gl, p1u, p1d, p2u, p2d)
      }
      png(filename = paste0(folder, "/Resultfigure/CS_", 
                            sprintf("%01.3f", RESULT[j, "ConsensusScore"], 
                                    3), "_", gsub(" ", "_", substr(name1, start = 1, 
                                                                   stop = min(40, nchar(name1)))), ".png"), width = 3.25, 
          height = 3.25, units = "in", res = 1200, pointsize = 6)
      gl2 <- max(gl, na.rm = TRUE) - min(gl, na.rm = TRUE)
      if (S == "GO") {
        plot(unique(grp1.time), p1, ylim = c(min(gl) - 
                                               (0.1 * gl2), max(gl, na.rm = TRUE)), type = "l", 
             col = "red", lwd = 2, ylab = paste0("average intensity of ", 
                                                 GOszj, " genes"), xlab = " time [h]", xaxt = "n")
      }
      else {
        plot(unique(grp1.time), p1, ylim = c(min(gl) - 
                                               (0.1 * gl2), max(gl, na.rm = TRUE)), type = "l", 
             col = "red", lwd = 2, ylab = "Intensity", xlab = " time [h]", 
             xaxt = "n")
      }
      points(grp1.time, Poi1, col = "red")
      lines(unique(grp2.time), p2, col = "blue", lwd = 2)
      points(grp2.time, Poi2, col = "blue", pch = 4)
      axis(side = 1, at = round(unique(grp1.time, grp2.time)))
      title(main = list(paste0(name, " (PubMed:", RESULT[j, 
                                                         "PubMed"], ")"), cex = cex1))
      legend("bottom", "groups", c(grp1n, "measurement", 
                                   grp2n, "measurement"), pch = c(15, 1, 15, 4), 
             col = c("red", "red", "blue", "blue"), ncol = 4, 
             bty = "n")
      if (S == "GO") {
        polygon(c(unique(grp1.time), rev(unique(grp1.time))), 
                c(p1u, rev(p1d)), col = rgb(1, 0, 0, 0.1), 
                xpd = FALSE, border = NA)
        polygon(c(unique(grp2.time), rev(unique(grp2.time))), 
                c(p2u, rev(p2d)), col = rgb(0, 0, 1, 0.1), 
                xpd = FALSE, border = NA)
      }
      dev.off()
      rm(cex1)
    }
    if (S == "GO") {
      RESULT[j, "Link"] = paste0("=HYPERLINK(TranslateQuotes", 
                                 paste0("Resultfigure/CS_", sprintf("%01.3f", 
                                                                    RESULT[j, "ConsensusScore"], 3), "_", gsub(" ", 
                                                                                                               "_", substr(name1, start = 1, stop = min(40, 
                                                                                                                                                        nchar(name1)))), ".png"), "TranslateQuotes,TranslateQuotes figure TranslateQuotes)")
    }
    else {
      RESULT[j, "Link"] = paste0("=HYPERLINK(TranslateQuotes", 
                                 paste0("Resultfigure/CS_", sprintf("%01.3f", 
                                                                    RESULT[j, "ConsensusScore"], 3), "_", gsub(" ", 
                                                                                                               "_", substr(name1, start = 1, stop = min(40, 
                                                                                                                                                        nchar(name1)))), ".png"), "TranslateQuotes,TranslateQuotes", 
                                 name1, "TranslateQuotes)")
    }
    if (j %in% prg) {
      idx = 1:length(prg)
      timeUse <- sum(as.numeric(proc.time() - timeusage)[1:2]) * 
        (20 - idx[j == prg])
      setTkProgressBar(global, (70 + idx[j == prg]), label = paste("Produce ", 
                                                                   min(nrow(RESULT), MaxPics), " pics to show the dynamic (", 
                                                                   round(idx[j == prg] * 5), " % done ). This step is expected to be completed by ", 
                                                                   ProgTime(timeUse), "."))
      timeusage <- proc.time()
    }
  }
  if (S == "GO") {
    rm(GOsd1u, GOsd1d, GOsd2u, GOsd2d, Poi1, Poi2, p1u, p1d, 
       p2u, p2d, gl)
  }
  # save(RESULT, file = paste0(folder, "/Saved_Intermediate_Results/RESULT-after-PicsAndLink.RData"))
  selection <- c("signifearlyRes", "signifmiddleRes", "signiflateRes", 
                 "signifResponse")
  category <- c("Early", "Middle", "Late", "Full")
  VENNfunc(result = RESULT, R0 = selection, category = category, 
           SAVEplot = paste0(folder, "/Venn_and_other_characteristics/integralVENN.png"))
  tresCS1 <- DistplotPval(x1 = CS1[, 1, drop = FALSE], xlab = "Combinated Integral Score", 
                          main = paste0("Significance in Integral Score"), pVal = pVal, 
                          folder = paste0(folder, "/Significance_And_Pval/"), ForTest = "Significant_Integral_Score")
  RESULT[, "IntegralScore"] = CS1[, 1] > tresCS1$tres
  selection <- c("SignPeak", "SignDyn", "IntegralScore", "signifConsensus")
  category <- c("Peak", "Dynamic", "Integral", "Consensus")
  VENNfunc(result = RESULT, R0 = selection, category = category, 
           SAVEplot = paste0(folder, "/Venn_and_other_characteristics/ConsensusScore1VENN.png"))
  if (S == "GO") {
    selection <- c("SignPeak", "SignDyn", "IntegralScore", 
                   "SignROS")
    category <- c("Peak", "Dynamic", "Integral", "revOverlap")
  }
  else {
    RESULT[, "sigPubMed"] <- RESULT[, "PubMed"] > 0
    selection <- c("SignPeak", "SignDyn", "IntegralScore", 
                   "sigPubMed")
    category <- c("Peak", "Dynamic", "Integral", "PubMed")
  }
  VENNfunc(result = RESULT, R0 = selection, category = category, 
           SAVEplot = paste0(folder, "/Venn_and_other_characteristics/ConsensusScore2VENN.png"))
  RESULT[, "RankConsensus"] <- rank(-(RESULT[, "ConsensusScore"]), 
                                    ties.method = "min", na.last = TRUE)
  RESULT[, "RankMaxDist"] <- rank(-(RESULT[, "MaxDist"]), ties.method = "min", 
                                  na.last = TRUE)
  RESULT[, "RankDynamic"] <- rank(-(RESULT[, "Dynamic"]), ties.method = "min", 
                                  na.last = TRUE)
  RESULT[, "RankearlyResponse"] <- rank(-(RESULT[, "earlyResponse"]), 
                                        ties.method = "min", na.last = TRUE)
  RESULT[, "RankmiddleResponse"] <- rank(-(RESULT[, "middleResponse"]), 
                                         ties.method = "min", na.last = TRUE)
  RESULT[, "RanklateResponse"] <- rank(-(RESULT[, "lateResponse"]), 
                                       ties.method = "min", na.last = TRUE)
  RESULT[, "RankcompleteResponse"] <- rank(-(RESULT[, "completeResponse"]), 
                                           ties.method = "min", na.last = TRUE)
  if (S == "GO") {
    RESULT[, "RankRevOverlapScore"] <- rank(-(RESULT[, "RevOverlapScore"]), 
                                            ties.method = "min", na.last = TRUE)
    RESULT[, "RankMean"] <- round(rowMeans(cbind(rowMeans(RESULT[, 
                                                                 c("RankearlyResponse", "RankmiddleResponse", "RanklateResponse")]), 
                                                 RESULT[, "RankDynamic"], RESULT[, "RankMaxDist"], 
                                                 RESULT[, "RankRevOverlapScore"])))
    RESULT[, "RankMedian"] <- apply(RESULT[, c("RankearlyResponse", 
                                               "RankmiddleResponse", "RanklateResponse", "RankDynamic", 
                                               "RankMaxDist", "RankRevOverlapScore")], 1, median)
  }
  else {
    RESULT[, "RankMean"] <- round(rowMeans(cbind(rowMeans(RESULT[, 
                                                                 c("RankearlyResponse", "RankmiddleResponse", "RanklateResponse")]), 
                                                 RESULT[, "RankDynamic"], RESULT[, "RankMaxDist"])))
    RESULT[, "RankMedian"] <- apply(RESULT[, c("RankearlyResponse", 
                                               "RankmiddleResponse", "RanklateResponse", "RankDynamic", 
                                               "RankMaxDist")], 1, median)
  }
  # save(RESULT, file = paste0(folder, "/Saved_Intermediate_Results/RESULT-after-Rankbased.RData"))
  header <- paste0("Start time: (", startTimeFunction, ")   End Time: (", 
                   Sys.time(), ")   Systeminfo: (", paste(as.character(unique(paste0(Sys.info()))), 
                                                          collapse = " "), ")")
  if (S == "GO") {
    RESULT_all <- RESULT[, c("gene_name", "ConsensusScore", 
                             "RankMean", "RankMedian", "PubMed", "groupsize", 
                             "Link", "signifConsensus", "SignPeak", "SignPeakH0", 
                             "SignInstability", "SignROS", "SignDyn", "signifearlyRes", 
                             "signifmiddleRes", "signiflateRes", "signifResponse", 
                             "PubMedscore", "MaxDist", "Instability", "RevOverlapScore", 
                             "Dynamic", "earlyResponse", "middleResponse", "lateResponse", 
                             "completeResponse", "Pval_MaxDist", "Pval_RevOverlapScore", 
                             "Pval_Dynamic", "Pval_earlyResponse", "Pval_middleResponse", 
                             "Pval_lateResponse", "Pval_completeResponse", "Pval_ConsensusScore", 
                             "RankConsensus", "RankMaxDist", "RankRevOverlapScore", 
                             "RankDynamic", "RankearlyResponse", "RankmiddleResponse", 
                             "RanklateResponse", "RankcompleteResponse", "probeset_id")]
    colnames(RESULT_all)[1] = "GO_name"
    colnames(RESULT_all)[length(colnames(RESULT_all))] = "GO_id"
  }
  else {
    RESULT_all <- RESULT[, c("gene_name", "ConsensusScore", 
                             "RankMean", "RankMedian", "PubMed", "Link", "signifConsensus", 
                             "SignPeak", "SignPeakH0", "SignInstability", "SignDyn", 
                             "signifearlyRes", "signifmiddleRes", "signiflateRes", 
                             "signifResponse", "PubMedscore", "MaxDist", "Instability", 
                             "Dynamic", "earlyResponse", "middleResponse", "lateResponse", 
                             "completeResponse", "Pval_MaxDist", "Pval_Dynamic", 
                             "Pval_earlyResponse", "Pval_middleResponse", "Pval_lateResponse", 
                             "Pval_completeResponse", "Pval_ConsensusScore", "RankConsensus", 
                             "RankMaxDist", "RankDynamic", "RankearlyResponse", 
                             "RankmiddleResponse", "RanklateResponse", "RankcompleteResponse", 
                             "probeset_id")]
  }
  if (colnames(as.data.frame(annotation))[1] != "annotation") {
    RESULT_all <- merge(RESULT_all, annotation, by = "row.names")
    rownames(RESULT_all) <- RESULT_all[, 1]
    RESULT_all <- RESULT_all[, -1]
  }
  file_name <- paste0(folder, "/RESULT_", grp1n, "-vs-", grp2n, 
                      "_all.tsv")
  my.write(x = RESULT_all, file = file_name, header = header, 
           f = write.table, sep = "\t", row.names = FALSE, quote = F)
  # save(RESULT_all, file = paste0(folder, "/Saved_Intermediate_Results/RESULT_END_all.RData"))
  if (S == "GO") {
    RESULT_rank <- RESULT[, c("gene_name", "RankMean", "RankMedian", 
                              "PubMed", "groupsize", "Link", "MaxDist", "Instability", 
                              "RevOverlapScore", "Dynamic", "earlyResponse", "middleResponse", 
                              "lateResponse", "completeResponse", "SignInstability", 
                              "RankConsensus", "RankMaxDist", "RankRevOverlapScore", 
                              "RankDynamic", "RankearlyResponse", "RankmiddleResponse", 
                              "RanklateResponse", "RankcompleteResponse", "probeset_id")]
    colnames(RESULT_rank)[1] = "GO_name"
    colnames(RESULT_rank)[length(colnames(RESULT_rank))] = "GO_id"
  }
  else {
    RESULT_rank <- RESULT[, c("gene_name", "RankMean", "RankMedian", 
                              "PubMed", "Link", "MaxDist", "Instability", "Dynamic", 
                              "earlyResponse", "middleResponse", "lateResponse", 
                              "completeResponse", "SignInstability", "RankConsensus", 
                              "RankMaxDist", "RankDynamic", "RankearlyResponse", 
                              "RankmiddleResponse", "RanklateResponse", "RankcompleteResponse", 
                              "probeset_id")]
  }
  if (colnames(as.data.frame(annotation))[1] != "annotation") {
    RESULT_rank <- merge(RESULT_rank, annotation, by = "row.names")
    rownames(RESULT_rank) <- RESULT_rank[, 1]
    RESULT_rank <- RESULT_rank[, -1]
  }
  file_name <- paste0(folder, "/RESULT_", grp1n, "-vs-", grp2n, 
                      "_rank.tsv")
  my.write(x = RESULT_rank, file = file_name, header = header, 
           f = write.table, sep = "\t", row.names = FALSE, quote = F)
  # save(RESULT_rank, file = paste0(folder, "/Saved_Intermediate_Results/RESULT_END_rank.RData"))
  if (S == "GO") {
    RESULT_Pval <- RESULT[, c("gene_name", "ConsensusScore", 
                              "Pval_ConsensusScore", "PubMed", "groupsize", "Link", 
                              "MaxDist", "Instability", "RevOverlapScore", "Dynamic", 
                              "earlyResponse", "middleResponse", "lateResponse", 
                              "completeResponse", "SignInstability", "Pval_MaxDist", 
                              "Pval_RevOverlapScore", "Pval_Dynamic", "Pval_earlyResponse", 
                              "Pval_middleResponse", "Pval_lateResponse", "Pval_completeResponse", 
                              "Pval_ConsensusScore", "probeset_id")]
    colnames(RESULT_Pval)[1] = "GO_name"
    colnames(RESULT_Pval)[length(colnames(RESULT_Pval))] = "GO_id"
  }
  else {
    RESULT_Pval <- RESULT[, c("gene_name", "ConsensusScore", 
                              "Pval_ConsensusScore", "PubMed", "Link", "MaxDist", 
                              "Instability", "Dynamic", "earlyResponse", "middleResponse", 
                              "lateResponse", "completeResponse", "SignInstability", 
                              "Pval_MaxDist", "Pval_Dynamic", "Pval_earlyResponse", 
                              "Pval_middleResponse", "Pval_lateResponse", "Pval_completeResponse", 
                              "probeset_id")]
  }
  if (colnames(as.data.frame(annotation))[1] != "annotation") {
    RESULT_Pval <- merge(RESULT_Pval, annotation, by = "row.names")
    rownames(RESULT_Pval) <- RESULT_Pval[, 1]
    RESULT_Pval <- RESULT_Pval[, -1]
  }
  file_name <- paste0(folder, "/RESULT_", grp1n, "-vs-", grp2n, 
                      "_Pval.tsv")
  my.write(x = RESULT_Pval, file = file_name, header = header, 
           f = write.table, sep = "\t", row.names = FALSE, quote = F)
  # save(RESULT_Pval, file = paste0(folder, "/Saved_Intermediate_Results/RESULT_END_Pval.RData"))
  fileConn <- file(paste0(folder, "/If_the_Link_in_table_doesnt_work.txt"))
  writeLines(c("Hi User", " ", "The link does not work automatically. Excel uses symbols for the hyperlink", 
               "that can not be printed directly from R. To make the link available you will", 
               "need to replace certain parts. Mark in Excel the column with the links and", 
               "use then the Excel function [find and replace all].", 
               "-First, replace all [TranslateQuotes] to correct quotes.", 
               "-Second, change the orientation from the slash (/) to a backslash.", 
               "Now, the links should work. Note that the table only together with the folder", 
               "[Resultfigure] can be moved, because otherwise the reference is lost. If some-", 
               "thing is still not work, look at the hyperlink function in Excel and replace", 
               "something if necessary.", " ", "With best regards", 
               " ", "Marco Albrecht ", "(the developer)"), fileConn)
  close(fileConn)
  Infotext <- list()
  Infotext$sessionInfo2 <- sessionInfo()
  Infotext$CalculationTime <- Sys.time() - startTimeFunction
  sink(file = paste0(folder, "/Information/Function_Input.txt"), 
       append = T)
  print(Infotext)
  sink(file = NULL)
  close(global)
  print("                                                                                                  ")
  print("                                                                                                  ")
  print("                  ---  Thanks for analysing with TTCA. Goodbye and take care  ---                 ")
  print("                                                                                                  ")
  print("                                                                                                  ")
  print("                                                                                                  ")
  print("--------------------------------------------------------------------------------------------------")
  print(paste0("Start: ", startTimeFunction, "                                                                        "))
  print(paste0("End:   ", Sys.time(), "                                        Kind Regards: Marco Albrecht    "))
  print("--------------------------------------------------------------------------------------------------")
  return(RESULT_all)
}
#end TTCA
#######################################################################################################################################################################################
#######################################################################################################################################################################################
#######################################################################################################################################################################################
#######################################################################################################################################################################################
#######################################################################################################################################################################################
#######################################################################################################################################################################################
PeakScore<-function(dat,grp, tme,folder,global){
  grp1      = dat[, grp==1]
  grp2      = dat[, grp==2]
  grp1.time = tme[grp==1]
  grp2.time = tme[grp==2]
  ##################################################
  info<-nrow(grp1)
  a1<-floor(nrow(dat)/1000)   # steps
  a2<-ceiling(nrow(dat)/a1)   # more as 1000 per step for dynamic threshold
  A1<-apply(dat, 1,mean)
  A2<-A1[order(A1)]
  tresVec2<-tresVec1<-f11<-f1<-list(NA)
  dvth<-dvSd<-NULL
  inter<-intersect(grp1.time,grp2.time)
  tresForgroupwise=7
  ##################################################
  if(a1>tresForgroupwise){ #only groupwise if database is big enough/ IF CHANGE than change also by mean-SD plot
    for (k in 1:a1){ # for group wise Calculation
      A3<-names(A2[((k-1)*a2+1):(k*a2)])
      grp1a     <-     grp1[ rownames(grp1)  %in%  A3, ]
      grp2a     <-     grp2[ rownames(grp2)  %in%  A3, ]
      L<-rep(NA,max(length(grp1a),length(grp2a),na.rm = TRUE)) #for thresholddetermination
      j1=1
      ##################################################
      grp1e<-grp1a[ ,1:length(inter)]*NA
      grp2e<-grp2a[ ,1:length(inter)]*NA
      colnames(grp1e)<-colnames(grp2e)<-as.character(inter)
      for (i in 1:length(inter)){ # merge columns
        if( sum(grp1.time==inter[i])>1 ){L[j1]<- as.numeric(quantile(apply(grp1a[,grp1.time==inter[i]], 1, sd), probs =0.95,na.rm=TRUE)) # threshold were 95% is below for each multiple replicated timepoint
        j1=j1+1
        grp1e[,i]<- apply(grp1a[,grp1.time==inter[i]], 1, mean) # mean of replicates
        }else{                           grp1e[,i]<- grp1a[,grp1.time==inter[i]]                 # replicate direct
        } # if more than one samples per timepoint are available, than take the mean itensity for this time point.
        if( sum(grp2.time==inter[i])>1 ){L[j1]<- as.numeric(quantile(apply(grp2a[,grp2.time==inter[i]], 1, sd), probs =0.95,na.rm=TRUE)) # threshold were 95% is below for each multiple replicated timepoint
        j1=j1+1
        grp2e[,i]<- apply(grp2a[,grp2.time==inter[i]], 1, mean) # mean of replicates
        }else{                           grp2e[,i]<- grp2a[,grp2.time==inter[i]]                 # replicate direct
        }
      }
      d<-apply(abs(grp1e-grp2e), 1, max)
      ##################################################
      tres<-mean(L[1:(j1-1)])*2 # threshold: mean of all 95til sd ; twofold
      dt<-d           # for groupthreshold
      dt[]=tres       # for groupthreshold
      dvSd=c(dvSd,d)  # standarddeviatien. in each loop add one value to the vector
      dvth=c(dvth,dt) # groupthreshold
      ##################################################
      ## figure
      png(filename = paste0(folder,"/Significance_And_Pval/PeakScore/DistanceDistribution_fromLowerMean",((k-1)*a2+1),"toUpperMean",(k*a2),"_threshold_",round(tres,2),".png"),width = 3.25, height    = 3.25,units     = "in",res       = 1200, pointsize = 6)
      hist(d,breaks=(0:100)/100*max(d,na.rm=TRUE),main="Maximum distance beetween both groups", sub="Everything to the left of the red line was removed (",tres,")")
      abline(v=tres,col="red",lty=1)
      dev.off()
      ## end figure
      ##################################################
      # for threshold mean plot
      tresVec1[[k]]<-mean(A2[((k-1)*a2+1):(k*a2)],na.rm=TRUE)
      tresVec2[[k]]<-tres
      #progress
      setTkProgressBar(global,(3+(k/a1*6)), label=paste("Tests for significant distances between the dynamics (",round(k/a1*100),"% done)"))
      ##################################################
    }
  }else{
    L<-rep(NA,max(length(grp1),length(grp2),na.rm = TRUE)) #for thresholddetermination
    j1=1
    ################################################## k
    grp1e<-grp1[ ,1:length(inter)]*NA
    grp2e<-grp2[ ,1:length(inter)]*NA
    colnames(grp1e)<-colnames(grp2e)<-as.character(inter)
    for (i in 1:length(inter)){ # merge columns
      if( sum(grp1.time==inter[i])>1 ){L[j1]<- as.numeric(quantile(apply(grp1[,grp1.time==inter[i]], 1, sd), probs =0.95,na.rm=TRUE)) # threshold were 95% is below for each multiple replicated timepoint
      j1=j1+1
      grp1e[,i]<- apply(grp1[,grp1.time==inter[i]], 1, mean) # mean of replicates
      }else{                           grp1e[,i]<- grp1[,grp1.time==inter[i]]                 # replicate direct
      } # if more than one samples per timepoint are available, than take the mean itensity for this time point.
      if( sum(grp2.time==inter[i])>1 ){L[j1]<- as.numeric(quantile(apply(grp2[,grp2.time==inter[i]], 1, sd), probs =0.95,na.rm=TRUE)) # threshold were 95% is below for each multiple replicated timepoint
      j1=j1+1
      grp2e[,i]<- apply(grp2[,grp2.time==inter[i]], 1, mean) # mean of replicates
      }else{                           grp2e[,i]<- grp2[,grp2.time==inter[i]]                 # replicate direct
      }
    }
    d<-apply(abs(grp1e-grp2e), 1, max)
    ##################################################
    tres<-mean(L[1:(j1-1)])*2 # threshold: mean of all 95til sd ; twofold
    dt<-d           # for groupthreshold
    dt[]=tres       # for groupthreshold
    dvSd=c(dvSd,d)  # standarddeviation. in each loop add one value to the vector
    dvth=c(dvth,dt) # groupthreshold
    ##################################################
    ## figure
    png(filename = paste0(folder,"/Significance_And_Pval/PeakScore/DistanceDistribution_AllinOne_threshold_",round(tres,2),".png"),width = 3.25, height    = 3.25,units     = "in",res       = 1200, pointsize = 6)
    hist(d,breaks=(0:100)/100*max(d,na.rm=TRUE),main="Maximum distance beetween both groups", sub="Everything to the left of the red line was removed (",tres,")")
    abline(v=tres,col="red",lty=1)
    dev.off()
    ## end figure
  }
  ##################################################
  ## For RESULT
  dvSd <- as.data.frame(dvSd)
  dvth <- as.data.frame(dvth)
  RESULT<-merge(dvth,dvSd ,by="row.names", all = TRUE)
  rownames(RESULT)<-RESULT[,1]
  RESULT<-RESULT[,-1]
  colnames(RESULT)<-c("MaxDistThres","MaxDist")
  RESULT[,"SignPeak"]<-RESULT[,"MaxDist"]> RESULT[,"MaxDistThres"]

  ##################################################
  ## Instability detection:  Median standarddeviation across all replicated timepoints
  Duplic<-unique(grp1.time[duplicated(grp1.time)],grp2.time[duplicated(grp2.time)])
  Ls<-matrix(NA,nrow(dat),length(Duplic))
  for (i in 1:length(Duplic)){
    Ls[,i]<-as.numeric(apply(grp1[  ,grp1.time==Duplic[i]], 1, sd))

  }
  sL<-t(Ls) # remove duplicated columns if any
  Ls<-t(sL[!duplicated(sL),])
  RESULT[,"Instability"]    <-apply(Ls,1,median)
  RESULT[,"SignInstability"]<-RESULT[,"Instability"]> RESULT[,"MaxDistThres"]
  ###############
  ## output preparation
  RESULT[,"MaxDist"] <-RESULT[,"MaxDist"]/(max(dat,na.rm=TRUE)-min(dat,na.rm=TRUE)) #normalize to complete intensity bandbrid of the array
  RESULT[,"Instability"] <-RESULT[,"Instability"]/(max(dat,na.rm=TRUE)-min(dat,na.rm=TRUE)) #normalize to complete intensity bandbrid of the array
  ##################################################
  ## plot
  if(a1>tresForgroupwise){
    tresVec1<-unlist(tresVec1)
    tresVec2<-unlist(tresVec2)
    dynamicMean<-as.numeric(rowMeans(dat))
    #
    png(filename = paste0(folder,"/Venn_and_other_characteristics/a_Mean_Standard_Deviation_Plot.png"),width = 3.25, height    = 3.25,units     = "in",res       = 1200, pointsize = 6)
    plot(tresVec1,tresVec2,type = "p", bty = "n",pch=16,col="red", xlab = "Mean of a dynamic", ylab = "Standard deviation threshold",xlim=c(min(dynamicMean),max(dynamicMean,na.rm = TRUE)),main=paste0("Mean of the dynamics and their related standard derivation"))
    par(new=TRUE)
    hist(dynamicMean, breaks = 100,  prob=TRUE,xlim=c(min(dynamicMean),max(dynamicMean,na.rm = TRUE)), axes = FALSE, xlab = "", ylab = "",main="")
    dev.off()
  }
  ## end Plot
  ##################################################
  ## Output
  RESULT<-RESULT[,c("MaxDist","SignPeak","Instability","SignInstability")]
  save(RESULT, file=paste0(folder,"/Saved_Intermediate_Results/RESULT-with-PeakScore.RData"))
  return(RESULT)
} # end Peakscore
#######################################################################################################################################################################################
PeakScoreHOH1<-function(dat,grp, tme,folder,global){
  grp1      = dat[, grp==1]
  grp2      = dat[, grp==2]
  grp1.time = tme[grp==1]
  grp2.time = tme[grp==2]
  ##################################################
  info<-nrow(grp1)
  inter<-intersect(grp1.time,grp2.time)
  ################################################## k
  grp1e<-grp1[ ,1:length(inter)]*NA
  grp2e<-grp2[ ,1:length(inter)]*NA
  colnames(grp1e)<-colnames(grp2e)<-as.character(inter)

  for(i in 1:length(inter)){ # merge columns
    if( sum(grp1.time==inter[i])>1 ){grp1e[,i]<- apply(grp1[,grp1.time==inter[i]], 1, mean) # mean of replicates
    }else{                           grp1e[,i]<- grp1[,grp1.time==inter[i]]                 # replicate direct
    } # if more than one samples per timepoint are available, than take the mean itensity for this time point.
    if( sum(grp2.time==inter[i])>1 ){grp2e[,i]<- apply(grp2[,grp2.time==inter[i]], 1, mean) # mean of replicates
    }else{                           grp2e[,i]<- grp2[,grp2.time==inter[i]]                 # replicate direct
    }
  }

  PeakHOH1 <- apply(abs(grp1e-grp2e), 1, max)

  RESULT<- as.data.frame(PeakHOH1)
  colnames(RESULT)<-"MaxDist"
  RESULT["Instability"]<-NA
  RESULT["SignInstability"]<-NA
  RESULT["SignPeak"]<-RESULT["SignInstability"] < NA
  RESULT<-RESULT[,c("MaxDist","SignPeak","Instability","SignInstability")]

  save(RESULT, file=paste0(folder,"/Saved_Intermediate_Results/RESULT-with-PeakScore-Without-use-of-replicates.RData"))
  return(RESULT)
} # end PeakscoreHOH1
#######################################################################################################################################################################################
VENNfunc<-function(result,R0,category,SAVEplot){
  R1<-result[,R0]
  R2<-rownames(result)

  B1 <-unique(R2[  R1[,1]])
  E2 <-unique(R2[  R1[,2]])
  G3 <-unique(R2[  R1[,3]])
  M4 <-unique(R2[  R1[,4]])
  str1 <- R0[1]
  str2 <- R0[2]
  str3 <- R0[3]
  str4 <- R0[4]
  png(filename = SAVEplot,width = 3.5,height    = 3.25,units     = "in",res       = 1200, pointsize = 6)

  draw.quad.venn(
    area1=length(B1),
    area2=length(E2),
    area3=length(G3),
    area4=length(M4),
    n12  =length(intersect( B1,E2 )),
    n13  =length(intersect( B1,G3 )),
    n14  =length(intersect( B1,M4 )),
    n23  =length(intersect( E2,G3 )),
    n24  =length(intersect( E2,M4 )),
    n34  =length(intersect( G3,M4 )),
    n123 =length(intersect(intersect( B1,E2 ),G3 )),
    n124 =length(intersect(intersect( B1,E2 ),M4 )),
    n134 =length(intersect(intersect( B1,G3 ),M4 )),
    n234 =length(intersect(intersect( E2,G3 ),M4 )),
    n1234 =length(intersect(intersect(intersect( B1,E2 ),G3 ),M4)),
    category = category,
    fill = c("yellow1", "purple4", "turquoise1", "chartreuse2"),
    label.col = "gray20",
    alpha = 0.5,
    cex = 1,
    col="aliceblue",
    cat.cex = 2,
    cat.col =  c("yellow1", "purple4", "turquoise1", "chartreuse2")
  )
  dev.off()
}
#######################################################################################################################################################################################
Integral<-function(dat="grp1.grp2", grp="grp", timeInt= "timeInt",grp1.time="grp1.time",grp2.time="grp2.time",globMaxDist ){

  # available data for each time course
  dat1 = as.numeric(dat[ grp==1])
  dat2 = as.numeric(dat[ grp==2])
  #Identify time gaps
  time=unique(c(grp1.time,grp2.time, timeInt))
  tim1<-time[!(time %in% grp1.time)]
  tim2<-time[!(time %in% grp2.time)]
  # Gives a table: 1 row: measurements and missing values (NA)
  #                2 row: time points
  RB1<- RB2<-rbind(rep(NA,length(time)),time)
  RB1[1,which(time %in% grp1.time)]<-dat1
  RB2[1,which(time %in% grp2.time)]<-dat2
  # Interpolate missing values of time course 1
  mis1<-which(is.na(RB1[1,])) #where are missing values?
  exis1<-which(!(is.na(RB1[1,]))) #where are existing values?
  if (length(mis1)>0){
    for (i in 1:length(mis1)){
      v1=mis1[i]
      RB1[1,v1]<-RB1[1,max(exis1[exis1<v1])]+((RB1[1,min(exis1[exis1>v1])]-RB1[1,max(exis1[exis1<v1])])/(as.numeric(RB1[2,min(exis1[exis1>v1])])-as.numeric(RB1[2,max(exis1[exis1<v1])])))*as.numeric(RB1[2,v1]-RB1[2,max(exis1[exis1<v1])])
    }}
  # Interpolate missing values of time course 2
  mis2<-which(is.na(RB2[1,])) #where are missing values?
  exis2<-which(!(is.na(RB2[1,]))) #where are existing values?
  if (length(mis2)>0){
    for (i in 1:length(mis2)){
      v2=mis2[i]
      RB2[1,v2]<-RB2[1,max(exis2[exis2<v2])]+((RB2[1,min(exis2[exis2>v2])]-RB2[1,max(exis2[exis2<v2])])/(as.numeric(RB2[2,min(exis2[exis2>v2])])-as.numeric(RB2[2,max(exis2[exis2<v2])])))*as.numeric(RB2[2,v2]-RB2[2,max(exis2[exis2<v2])])
    }}
  # combine both curves in one matrix
  RB<-rbind(RB1[1,],RB2,rep(NA,ncol(RB1)),rep(NA,ncol(RB1)))
  # search curve crossing
  RB[4,]<-c(0, diff(0.5*sign(RB[1,]-RB[2,])))
  sp<-which(RB[4,]!=0)
  if(length(sp)>0){
    # Add columns for the new values at the timepoints for curve crossing
    RB<-cbind(RB,matrix(rep(NA,length(sp)*5),5,length(sp)))
    # calculate Intercept point and add values in new columns
    for (i in 1:length(sp)){
      r=sp[i]
      RB[3,(length(time)+i)]<-RB[3,(r-1)]+((RB[1,(r-1)]-RB[2,(r-1)])*(RB[3,r]-RB[3,(r-1)]))/(RB[2,r]-RB[2,(r-1)]-RB[1,r]+RB[1,(r-1)])
      RB[1,(length(time)+i)]<-RB[1,(r-1)]+((RB[1,r]-RB[1,(r-1)])/(RB[3,r]-RB[3,(r-1)]))*(RB[3,(length(time)+i)]-RB[3,(r-1)])
      RB[2,(length(time)+i)]<-RB[2,(r-1)]+((RB[2,r]-RB[2,(r-1)])/(RB[3,r]-RB[3,(r-1)]))*(RB[3,(length(time)+i)]-RB[3,(r-1)])
    }}
  # Order matrix according time sequence
  RB<-RB[,order(RB[3,])]
  # calculate trapezoids
  RB[5,1]=0
  for (i in 2:ncol(RB)){
    RB[5,i]<-abs((  (RB[1,i]+RB[1,(i-1)])/2  -  (RB[2,i]+RB[2,(i-1)])/2) * as.numeric(RB[3,i]-RB[3,(i-1)]))
  }

  earlyRev<- globMaxDist*as.numeric((RB[3,which(RB[3,]==timeInt[1])]-RB[3,1]))
  middleRev<-globMaxDist*as.numeric((RB[3,which(RB[3,]==timeInt[2])]-RB[3,which(RB[3,]==timeInt[1])]))
  lateRev<-  globMaxDist*as.numeric((RB[3,ncol(RB)]                 -RB[3,which(RB[3,]==timeInt[2])]))
  compRev<-  globMaxDist*as.numeric((RB[3,ncol(RB)]                 -RB[3,1]))

  # period specific integral, normalized by maximum value
  earlyResponse    <-sum(RB[5,c((               1               ):which(RB[3,]==timeInt[1]))])/earlyRev
  middleResponse   <-sum(RB[5,c((which.max(RB[3,]==timeInt[1])+1):which(RB[3,]==timeInt[2]))])/middleRev
  lateResponse     <-sum(RB[5,c((which.max(RB[3,]==timeInt[2])+1): ncol(RB)                )])/lateRev
  completeResponse <-sum(RB[5,                                                              ])/compRev
  sum(earlyResponse,middleResponse,lateResponse)
  completeResponse

  I=c(earlyResponse,middleResponse,lateResponse,completeResponse)
  return(I)
}
#######################################################################################################################################################################################
#upper<-DistplotPval(x1=RESULT[,"MaxDist" ,drop=FALSE],xlab=xlab, pVal=pVal,folder=paste0(folder,"/LocalHypothesis/"))
DistplotPval<-function(x1,xlab="value",main="Null hypothesis versus alternative hypothesis", pVal=pVal,folder,ForTest){
  Res=x1
  name1=colnames(x1)
  name2=paste0("Pval_",name1)
  x1=x1[,1,drop=TRUE] #from data.frame to numeric vector
  x1[is.na(x1)]=0.00001 #ok if one expect only one or two NA
  nbox=150 #number of boxes in destribution
  xfit<-seq(min(x1,na.rm = TRUE),max(x1,na.rm = TRUE),length=2000)
  ## find optimal Method and Parameter for fit
  d1   <-c(seq(0,(max(x1)+(max(x1)-min(x1)))+0.0001,(max(x1)-min(x1))/8 + 1e-5),0.0001,0.001,0.01)# with focus around zero and small trend in positiv direction
  d1   <-d1[order(d1)]
  dmax <-max(abs(d1))
  d3<-d2   <-1:length(d1)*NA
  for (i in 1:length(d1)){
    A<-optimalMethodeDistr( delta=d1[i],x1=x1,dmax=dmax)
    d2[i]<-A$opt
    d3[i]<-A$method
    rm(A)
  }
  bestMethod<-d3[which(d2==max(d2,na.rm = TRUE))[1]]
  delta<-d1[which(d2==max(d2,na.rm = TRUE))[1]]
  ## transformation
  x1 = x1 + delta
  xfit<-seq(min(x1,na.rm = TRUE),max(x1,na.rm = TRUE),length=2000)
  #
  parG<-suppressWarnings(fitdistr(x1,bestMethod))
  if(bestMethod=="cauchy")   {yfit<-dcauchy(xfit,location=as.numeric(parG[]$estimate[1]),scale   = as.numeric(parG[]$estimate[2]))
  Res[,name2]= apply(X=Res[,name1,drop=FALSE],MARGIN=1, function(x) pcauchy((x+delta),location=as.numeric(parG[]$estimate[1]),scale   = as.numeric(parG[]$estimate[2]),lower.tail = FALSE)) #Pval=P[X > x]
  fminforPval<-function(x, parG,pVal){
    a=abs(pVal-pcauchy(x,location=as.numeric(parG[]$estimate[1]),scale   = as.numeric(parG[]$estimate[2]),lower.tail = FALSE))
  }
  upper= optimize(f = fminforPval , lower = min(x1),upper = max(x1), maximum = FALSE,tol = .0000001,parG=parG,pVal=pVal   )$minimum
  }
  if(bestMethod=="gamma")    {yfit<-dgamma(xfit,shape=as.numeric(parG[]$estimate[1]), rate = as.numeric(parG[]$estimate[2]))
  Res[,name2]= apply(X=Res[,name1,drop=FALSE],MARGIN=1, function(x) pgamma((x+delta),shape=as.numeric(parG[]$estimate[1]), rate = as.numeric(parG[]$estimate[2]),lower.tail = FALSE)) #Pval=P[X > x]
  fminforPval<-function(x, parG,pVal){
    a=abs(pVal-pgamma(x,shape=as.numeric(parG[]$estimate[1]), rate = as.numeric(parG[]$estimate[2]),lower.tail = FALSE))
  }
  upper= optimize(f = fminforPval , lower = min(x1),upper = max(x1), maximum = FALSE,tol = .0000001,parG=parG,pVal=pVal   )$minimum
  }
  if(bestMethod=="lognormal"){yfit<-dlnorm(xfit,meanlog=as.numeric(parG[]$estimate[1]), sdlog = as.numeric(parG[]$estimate[2]))
  Res[,name2]= apply(X=Res[,name1,drop=FALSE],MARGIN=1, function(x) plnorm((x+delta),meanlog=as.numeric(parG[]$estimate[1]), sdlog = as.numeric(parG[]$estimate[2]),lower.tail = FALSE)) #Pval=P[X > x]
  fminforPval<-function(x,parG,pVal){
    a=abs(pVal-plnorm(x,meanlog=as.numeric(parG[]$estimate[1]), sdlog = as.numeric(parG[]$estimate[2]),lower.tail = FALSE))
  }
  upper= optimize(f = fminforPval , lower = min(x1),upper = max(x1), maximum = FALSE,tol = .0000001,  parG=parG,pVal=pVal   )$minimum
  }
  if(bestMethod=="logistic") {yfit<-dlogis(xfit,location=as.numeric(parG[]$estimate[1]),scale   = as.numeric(parG[]$estimate[2]))
  Res[,name2]= apply(X=Res[,name1,drop=FALSE],MARGIN=1, function(x) plogis((x+delta),location=as.numeric(parG[]$estimate[1]),scale   = as.numeric(parG[]$estimate[2]),lower.tail = FALSE)) #Pval=P[X > x]
  fminforPval<-function(x, parG,pVal){
    a=abs(pVal - plogis(x,location=as.numeric(parG[]$estimate[1]),scale   = as.numeric(parG[]$estimate[2]),lower.tail = FALSE))
  }
  upper= optimize(f = fminforPval , lower = min(x1),upper = max(x1), maximum = FALSE,tol = .0000001,parG=parG,pVal=pVal   )$minimum
  }
  if(bestMethod=="normal")   {yfit<-dnorm(xfit,mean=as.numeric(parG[]$estimate[1]), sd = as.numeric(parG[]$estimate[2]))
  Res[,name2]= apply(X=Res[,name1,drop=FALSE],MARGIN=1, function(x) pnorm((x+delta),mean=as.numeric(parG[]$estimate[1]), sd = as.numeric(parG[]$estimate[2]),lower.tail = FALSE)) #Pval=P[X > x]
  fminforPval<-function(x, parG,pVal){
    a=abs(pVal - pnorm(x,mean=as.numeric(parG[]$estimate[1]), sd = as.numeric(parG[]$estimate[2]),lower.tail = FALSE))
  }
  upper= optimize(f = fminforPval , lower = min(x1),upper = max(x1), maximum = FALSE,tol = .0000001, parG=parG,pVal=pVal   )$minimum
  }
  if(bestMethod=="weibull")  {yfit<-dweibull(xfit,shape = as.numeric(parG[]$estimate[1]),scale   = as.numeric(parG[]$estimate[2]))
  Res[,name2]= apply(X=Res[,name1,drop=FALSE],MARGIN=1, function(x) pweibull((x+delta),shape = as.numeric(parG[]$estimate[1]),scale   = as.numeric(parG[]$estimate[2]),lower.tail = FALSE)) #Pval=P[X > x]
  fminforPval<-function(x, parG,pVal){
    a=abs(pVal - pweibull(x,shape = as.numeric(parG[]$estimate[1]),scale   = as.numeric(parG[]$estimate[2]),lower.tail = FALSE))
  }
  upper= optimize(f = fminforPval , lower = min(x1),upper = max(x1), maximum = FALSE,tol = .0000001,  parG=parG,pVal=pVal   )$minimum
  }

  yfit   =yfit*diff(hist(x1,breaks=nbox,plot = FALSE)$mids[1:2])*length(x1)
  pthresu=yfit[which(xfit>=upper)[1]]
  # end calculate y-values. push x-Achses back
  x1 = x1  - delta
  xfit=xfit- delta
  upper=upper-delta

  # Pic
  png(filename = paste0(folder,ForTest,"_PVal-threshold_",round(upper,3),".png"),width = 3.25,
      height    = 3.25,units     = "in",res       = 1200, pointsize = 6)
  hist(x1,breaks=nbox, main=main,xlab=xlab)#,xaxt="n")
  lines(xfit, yfit, col="blue", lwd=2)
  polygon(c(upper,upper+max(xfit)/150,upper+max(xfit)/150,upper),c(pthresu*1.5,pthresu*1.2,0,0),col="red4",border="red4")
  legend("topright", inset=.05,c(paste0("Fit with ",bestMethod," distribution"), paste0("Significant dynamics with p-value: ",round(pVal,3)),paste0("threshold: ",round(upper,3))), fill=c("blue","lightcoral","red4"), horiz=FALSE,bty = "n")
  polygon(c(xfit[xfit>upper],rev(xfit[xfit>upper])),c(yfit[xfit>upper],yfit[xfit>upper]*0),col=rgb(1, 0, 0,0.5),xpd=FALSE)#
  dev.off()
  #End Pic
  A=list()
  A$tres=upper
  A$Pval=Res[,name2,drop=FALSE]
  return(A)
}
#######################################################################################################################################################################################
optimalMethodeDistr<-function(delta, x1,dmax ){
  x1=  x1+delta
  xfit<-seq(min(x1,na.rm = TRUE),max(x1,na.rm = TRUE),length=2000)
  # search destribution with bestk fit
  cmet<-c( "cauchy", "gamma", "lognormal", "logistic","normal","weibull")
  loglik<-NULL
  for (i in 1:length(cmet)){
    loglik[i]<-suppressWarnings(tryCatch(fitdistr(x1,cmet[i])$loglik, error=function(e)  {return(-Inf)}))  #test godness of fit, if somthing produce error, return -inf
  }
  A=list()
  A$opt<-max(loglik,na.rm = TRUE)-(abs(delta)/dmax)*0.01*max(loglik,na.rm = TRUE)# with penaltyterm to have minimum delta if the result is equal
  A$method<-cmet[which(loglik==max(loglik,na.rm = TRUE))]
  return(A)
}

#######################################################################################################################################################################################
my.write <- function(x, file, header, f = write.csv, ...){
  # create and open the file connection
  datafile <- file(file, open = 'wt')
  # close on exit
  on.exit(close(datafile))
  # if a header is defined, write it to the file
  if(!missing(header)) writeLines(header,con=datafile)
  # write the file using the defined function and required addition arguments
  f(x, datafile,...)
}
#######################################################################################################################################################################################
#FF.median(dat=grp1.grp2, grp=grp, tme=grp1.grp2.time,dynScor=dynScor,global=global)

FF.median = function(dat, grp, tme,dynScor=0,global=0){
  dat1 = dat[, grp==1]
  dat2 = dat[, grp==2]
  tme1 = tme[grp==1]
  tme2 = tme[grp==2]

  weight<-weightvec(grp1.time=tme1,grp2.time=tme2)
  prg<-unique(round((seq(nrow(dat)/32,nrow(dat),  nrow(dat)/32))))
  timeusage<-proc.time()
  rss0 = sapply(1:nrow(dat) , medianr, dat=dat , xx=tme  ,w=weight$grp1and2 , k=2,prg=prg,dynScor=dynScor)
  timeUse<-sum(as.numeric(proc.time()-timeusage)[1:2])*2
  if (as.character(global )[1]!="0"){setTkProgressBar(global, 20, label=paste0( "Program checks whether the dynamics are significant. This step is expected to be completed by ",ProgTime(timeUse),"."))}
  if (as.character(dynScor)[1]!="0"){setTkProgressBar(dynScor,34, label=paste0( "Program checks whether the dynamics are significant. (",34," % done)"))}
  weight<-weightvec(grp1.time=tme1,grp2.time=tme2)
  prg<-unique(round((seq(nrow(dat)/32,nrow(dat),  nrow(dat)/32))))
  timeusage<-proc.time()
  # fill up a list in list here. first list indicate gene. second list indicate error value a and fitted values res.
  Output1 = lapply(1:nrow(dat1), medianrS, dat=dat1, xx=tme1 ,w=weight$grp1 ,k=34,prg=prg,dynScor=dynScor)
  timeUse<-sum(as.numeric(proc.time()-timeusage)[1:2])
  if (as.character(global )[1]!="0"){setTkProgressBar(global, 30, label=paste0( "Program checks whether the dynamics are significant. This step is expected to be completed by ",ProgTime(timeUse),"."))}
  if (as.character(dynScor)[1]!="0"){setTkProgressBar(dynScor,66, label=paste0( "Program checks whether the dynamics are significant. (",66," % done)"))}
  Output2 = lapply(1:nrow(dat2), medianrS, dat=dat2, xx=tme2 ,w=weight$grp2 ,k=66,prg=prg,dynScor=dynScor)
  if (as.character(global )[1]!="0"){setTkProgressBar(global, 39, label=paste0( "Program checks whether the dynamics are significant. This step will be finished soon."))}
  if (as.character(dynScor)[1]!="0"){setTkProgressBar(dynScor,99, label=paste0( "Program checks whether the dynamics are significant. (",99," % done)"))}
  # unlist list in list.
  RolfA1<-1:length(Output1)*NA
  RolfRes1<- matrix( NA, length(Output1),length(Output1[[1]]$res))
  for (i in 1:length(Output1) ){
    RolfA1[i]<-Output1[[i]]$a
    RolfRes1[i,]<-Output1[[i]]$res
  }
  RolfA2<-1:length(Output2)*NA
  RolfRes2<- matrix( NA, length(Output2),length(Output2[[1]]$res))
  for (i in 1:length(Output2) ){
    RolfA2[i]<-Output2[[i]]$a
    RolfRes2[i,]<-Output2[[i]]$res
  }

  rss1  = colSums(rbind(RolfA1, RolfA2))
  # Make output list
  Dynamic=list()
  Dynamic$FF <- rss0/rss1
  Dynamic$qr1<-RolfRes1
  Dynamic$qr2<-RolfRes2
  return(Dynamic)

}
#######################################################################################################################################################################################
ProgTime<-function(timeforSample)
{ tt1<-as.numeric(unlist(strsplit( format(Sys.time(), "%H:%M:%S"),":")[1]))
  Sincemidnight<-tt1[1]*3600+tt1[2]*60+tt1[3]
  prognostTime<-timeforSample+Sincemidnight
  min<-h<-d<-0
  d<-floor(prognostTime/(24*3600))
  prognostTime<-prognostTime-24*3600*d
  h<-floor(prognostTime/(3600))
  prognostTime<-prognostTime-3600*h
  min<-floor(prognostTime/(60))
  if (h<12){ uhrzeit<-paste0(sprintf("%02.0f", h),":",sprintf("%02.0f", min)," AM")} else { uhrzeit<-paste0(sprintf("%02.0f", h-12),":",sprintf("%02.0f", min)," PM")}
  if (d==0){ day<-"today" }
  if (d==1){ day<-"tomorrow" }
  if (d>1){ day<- paste0("in ",d," days") }
  time<-paste0( uhrzeit," ",day )
  return(time)
}
#######################################################################################################################################################################################
# grp1.time <- c(0,0,0.5,0.5,1,2,3,3,3,4,5,6,8,12,18,24,24,24,48,48,48,50,50,50)
# grp2.time <- c(0,0,0,0.5,0.5,2,4,6,8,12,24,24,24,48,48,48,48,48)
# weight<-weightvec(grp1.time=grp1.time,grp2.time=grp2.time)

weightvec<-function(grp1.time,grp2.time){

  d<-c(grp1.time,grp2.time)
  e<-rep(1, length(d))
  a<-unique(c(grp1.time,grp2.time))
  b<-unique(grp1.time)
  f<-unique(grp2.time)
  # Compare only common time points
  e[!(d %in% intersect(b,f))]<-0
  e1<-e[1:length(grp1.time)]
  e2<-e[(length(grp1.time)+1):length(e)]
  e3<-e

  g1<-as.data.frame(table(grp1.time))
  for (i in 1:length(grp1.time)){
    e1[i]<-e1[i]*(1/g1[g1$grp1.time==grp1.time[i],2] )} # each common time-point weight have the same value 1
  e1<-e1*(length(grp2.time)/length(grp1.time))  # models with more datapoints have more residuels, but perhaps a better fit.
  # all residuels weighted to residuel per datapoint
  # otherwise a model with less datapoints have a benefit, because there more flexible in fitting
  #sum(e1/(length(grp1.time)))
  #sum(e1*(length(grp2.time)/length(grp1.time)))
  #sum(e1*(length(grp2.time)/(length(grp1.time)+length(grp2.time))))

  g2<-as.data.frame(table(grp2.time))
  for (i in 1:length(grp2.time)){
    e2[i]<-e2[i]*(1/g2[g2$grp2.time==grp2.time[i],2] )}
  e2<-e2*(length(grp1.time)/length(grp2.time))
  #sum(e2/(length(grp2.time)))
  #sum(e2*(length(grp1.time)/length(grp2.time)))
  #sum(e2*(length(grp1.time)/(length(grp1.time)+length(grp2.time))))

  sum(c(e1,e2))
  g3<-as.data.frame(table(d))
  for (i in 1:length(d)){
    e3[i]<-e3[i]*(1/g3[g3$d==d[i],2] )}
  e3<-e3*sum(c(e1,e2))/sum(e3)                 # balance the local null hypotheses out

  #sum(e3/length(e3))

  A<-list()
  A$grp1<-e1
  A$grp2<-e2
  A$grp1and2<-e3

  return(A)
}
#######################################################################################################################################################################################
# res2 = sapply(1:nrow(dat2), medianr, dat=dat2, xx=tme2 ,w=weight$grp2 ,k=66,prg=prg,dynScor=dynScor)
medianr = function(i, dat, xx ,w,k,prg,dynScor)
{
  tmp_df <- data.frame(response = as.numeric(dat[i, ]),
                       xx = xx)
  res <- rqss(response ~ qss(xx, lambda = 0.6), tau = 0.5, data = tmp_df) #knots/estimated Intensity for each time point
  res$coef[2:length(res$coef)]=res$coef[2:length(res$coef)]+  res$coef[1]
  res<-res$coef
  z1<-rle(xx[order(xx)])$lengths   # replikates per timepoint
  z2<-cumsum(c(1,z1)) # cumative sum
  H<-matrix(NA,1,sum(z1))
  # Intensity for each timepoint multiplicated with amount of replicates:
  for (j in 1:length(z1)){
    H[z2[j]:(z2[j]+z1[j]-1)]=rep(res[j],z1[j])
  }
  l<-abs(H-dat[i,]) # residual= estimated Intensity minus rawdata
  a<-sum(l*w)   # each value have a particular weight
  # progressbar
  if ((i %in% prg) & (as.character(dynScor)[1]!="0") ){
    idx=1:length(prg)
    setTkProgressBar(dynScor, (idx[i==prg]+k), label=paste0( "Program checks whether the dynamics are significant. (",idx[i==prg]+k," % done)"))
  }
  return(a)
}
######################################################################
# like medianr but with fitted values as additional feedback
medianrS = function(i, dat, xx ,w,k,prg,dynScor)
{
  tmp_df <- data.frame(response = as.numeric(dat[i, ]),
                       xx = xx)
  res <- rqss(response ~ qss(xx, lambda = 0.6), tau = 0.5, data = tmp_df) #knots/estimated Intensity for each time point
  res$coef[2:length(res$coef)]=res$coef[2:length(res$coef)]+  res$coef[1]
  res<-res$coef
  z1<-rle(xx[order(xx)])$lengths   # replikates per timepoint
  z2<-cumsum(c(1,z1)) # cumative sum
  H<-matrix(NA,1,sum(z1)) #prepare an empty vector
  # Intensity for each timepoint multiplicated with amount of replicates:
  for (j in 1:length(z1)){
    H[z2[j]:(z2[j]+z1[j]-1)]=rep(res[j],z1[j]) # New fitted value times number of replikates to obtain a vector with the same size as original measurement vector is.
  }
  l<-abs(H-dat[i,]) # residual= estimated Intensity minus rawdata
  a<-sum(l*w)   # each value have a particular weight
  # progressbar
  if ((i %in% prg) & (as.character(dynScor)[1]!="0") ){
    idx=1:length(prg)
    setTkProgressBar(dynScor, (idx[i==prg]+k), label=paste0( "Program checks whether the dynamics are significant. (",idx[i==prg]+k," % done)"))
  }
  Output=list()
  Output$a<-a
  Output$res<-res
  return(Output)
}
#################################################################################################################################
ChangeToGO<-function(grp1=grp1,grp2=grp2,grp1.time=grp1.time, mapGO){

  #
  dat=merge(grp1,grp2,by="row.names")
  rownames(dat) <- dat[,"Row.names"]
  dat           <- dat[,-1]
  grp           <- c(rep(1, ncol(grp1)),rep(2, ncol(grp2)))# data groups


  A3<-unique(mapGO[,"go_id"])
  B<-matrix(NA,length(A3),ncol(dat))
  colnames(B)=colnames(dat)
  rownames(B)=A3
  Sd<-Su<-B #S for standarddeviation, B for mean Itensity
  groupsize<-B[ ,1,drop=FALSE]
  colnames(groupsize)<-"groupsize"
  prggo<-unique(round((seq(nrow(B)/80,nrow(B),  nrow(B)/80)))) #for progress
  idxgo=1:length(prggo)
  nstart=which(duplicated(grp1.time)==TRUE)[1] #number of replicates at starttime (0h)
  for (i in 1:nrow(B)){
    A4 = mapGO[mapGO[,"go_id"]==rownames(B)[i],"probeset_id"] #probeset Ids for a given GO group
    A5 = dat[rownames(dat) %in% A4 ,] # expression data for each probeset_ID for a given group
    A6 = t(apply(A5,1,function(x) x- mean(x[1:nstart])))# zeropoint is mean centered for all genes in this group
    if(nrow(A6)>2){
      B[i,]=apply(A6,2,mean)
      A7<-A6*NA
      groupsize[i,1]<-nrow(A6)
      for (k in 1:nrow(A6)){
        A7[k,] = A6[k,] > B[i,]
      }
      for (j in 1:ncol(A6)){
        su      <- A6[A7[ ,j]==1,j] # over the mean
        Su[i,j] <- mean(su-B[i,j])
        sd      <- A6[A7[ ,j]==0,j] # under the mean
        Sd[i,j] <- -mean(sd-B[i,j])
      }}}
  B<-B[!is.na(B[,1]),]
  Su<-Su[!is.na(Su[,1]),]
  Sd<-Sd[!is.na(Sd[,1]),]
  foundall<-intersect(intersect(rownames(B),rownames(Su)),rownames(Sd)) # get sure that all have the same Go groups
  B<-B[rownames(B) %in% foundall,]
  Su<-Su[rownames(Su) %in% foundall,]
  Sd<-Sd[rownames(Sd) %in% foundall,]
  print(paste0(format(Sys.time(), "%H:%M:%S"),": ",nrow(B)," Gene Ontology groups can be analyzed                                               "))
  ####
  A3<-mapGO[,c("go_id","GO_name")]
  colnames(A3)=c("probeset_id","gene_name") #only for the programm
  A4<-A3[A3[,"probeset_id"] %in% rownames(B),]
  A4<-A4[!duplicated(A4[,"probeset_id"]),] #ID in row have to be unique
  rownames(A4)<-A4[,"probeset_id"]
  ###
  GO<-list()
  GO$grp1<-B[,grp==1]
  GO$grp2<-B[,grp==2]
  GO$annot<-A4
  GO$groupsize<-groupsize
  GO$grp1SDu<-Su[,grp==1]
  GO$grp2SDu<-Su[,grp==2]
  GO$grp1SDd<-Sd[,grp==1]
  GO$grp2SDd<-Sd[,grp==2]
  return(GO)
}
#################################################################################################################################
revOverlapScore<-function(go=GO,tm1=grp1.time,tm2=grp2.time,folder=folder){
  GO<-NULL
  grp1.time<-NULL
  grp2.time<-NULL
  tmi<-intersect(tm1,tm2)
  grp1<-go$grp1
  grp2<-go$grp2
  grp1SDu <-go$grp1SDu
  grp2SDu <-go$grp2SDu
  grp1SDd <-go$grp1SDd
  grp2SDd <-go$grp2SDd
  GOsz    <-go$groupsize
  ROS<-grp1[,1,drop=FALSE]*NA
  colnames(ROS)<-"RevOverlapScore"
  for(n in 1:nrow(ROS)){
    ovl=1:length(tmi)*NA
    for(m in 1:length(tmi)){
      m1<-mean(   grp1[n, tm1==tmi[m]])
      m2<-mean(   grp2[n, tm2==tmi[m]])
      sdu1<-mean(grp1SDu[n, tm1==tmi[m]])
      sdu2<-mean(grp2SDu[n, tm2==tmi[m]])
      sdd1<-mean(grp1SDd[n, tm1==tmi[m]])
      sdd2<-mean(grp2SDd[n, tm2==tmi[m]])
      if(all(m1==m2,sdu1==sdu2,sdd1==sdd2)){ ovl[m]=NA }else{ #ignore identical columns. Example 0 control
        ovl[m]= max((m1-sdd1)-(m2+sdu2),(m2-sdd2)-(m1+sdu1))/(0.5*(sdu1+sdu2+sdd1+sdd2))
      }}
    ROS[n,]<-max(ovl, na.rm = TRUE)
  }
  # for Pic
  png(filename = paste0(folder,"/Venn_and_other_characteristics/ReverseOverlapscore_MintoMax.png"),width = 5.25,height    = 6.25,units     = "in",res       = 1200, pointsize = 6)
  par(mfrow=c(2,1))
  M<-which(ROS==max(ROS, na.rm = TRUE))[1]
  Mr<-c(grp1[M,],grp2[M,],grp1[M,]+grp1SDu[M,],grp2[M,]+grp2SDu[M,],rev(grp1[M,]-grp1SDd[M,]),grp2[M,]-grp2SDd[M,])
  plot(tm1,grp1[M,], type="l" ,col="red",ylab="Intensity change",xlab="time",main=paste0("maximum reverse overlap score: ",round(max(ROS),2)," (",-round(max(ROS)*100,1),"% minimum overlap)"), ylim=range(Mr))
  lines(tm2,grp2[M,],col="blue")
  polygon(c((tm1),rev((tm1))),c(grp1[M,]+grp1SDu[M,],rev(grp1[M,]-grp1SDd[M,])),col=rgb(1, 0, 0,0.1),xpd=FALSE,border = NA)
  polygon(c((tm2),rev((tm2))),c(grp2[M,]+grp2SDu[M,],rev(grp2[M,]-grp2SDd[M,])),col=rgb(0, 0, 1,0.1),xpd=FALSE,border = NA)
  M<-which(ROS==min(ROS, na.rm = TRUE))[1]
  Mr<-c(grp1[M,],grp2[M,],grp1[M,]+grp1SDu[M,],grp2[M,]+grp2SDu[M,],rev(grp1[M,]-grp1SDd[M,]),grp2[M,]-grp2SDd[M,])
  plot(tm1,grp1[M,], type="l" ,col="red",ylab="Intensity change",xlab="time",main=paste0("Minimum reverse overlap score: ",round(min(ROS),2)," (",round(-min(ROS)*100,1),"% minimum overlap)"),ylim=range(Mr))
  lines(tm2,grp2[M,],col="blue")
  polygon(c((tm1),rev((tm1))),c(grp1[M,]+grp1SDu[M,],rev(grp1[M,]-grp1SDd[M,])),col=rgb(1, 0, 0,0.1),xpd=FALSE,border = NA)
  polygon(c((tm2),rev((tm2))),c(grp2[M,]+grp2SDu[M,],rev(grp2[M,]-grp2SDd[M,])),col=rgb(0, 0, 1,0.1),xpd=FALSE,border = NA)
  dev.off()
  # end Pic
  return(ROS)
} #end function
