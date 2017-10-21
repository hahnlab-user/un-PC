#' Unbundled PCA (un-PC) uses PCA results together with geographic sampling information to infer
#' past patterns of migration on the landscape
#' @param inputToProcess is either a single file or a directory containing multiple PCA result files to process.
#' @param outputPrefix is the name prefix used to generate the output .r file from data import, the .rds file from
#' pairwiseCalc, and the .pdf files from the unPC plotting (depending on which of these tasks are activated in the
#' function flags)
#' @param runDataImport boolean TRUE/FALSE. TRUE Runs the data import module; FALSE does not.
#' @param runPairwiseCalc boolean TRUE/FALSE. TRUE Runs the pairwise calculation module (which must follow the import module); FALSE does not.
#' @param geogrCoords specifies the path to the file containing the geographic coordinates for each individual
#' represented in the input file(s).
#' If the input PCA results were calculated with smartPCA following the msLandscape
#' pipeline, the geogrCoords are the automatically generated file of focal population locations and are then used
#' to calculate the pairwise unPC values between populations.
#' If the input PCA results are not from msLandscape, then the geogrCoords need to be the coordinates for each individual (i.e.
#' the number of rows in geogrCoords must be the same as the number of rows in the input PCA data). These geogrCoords
#' are then used to generate population labels for each individual and these reduced population-level coordinates
#' are then used to calculate the pairwise unPC values between populations. Currently, this functionality exists only when visualizing
#' results from a single file, not all files in a directory.
#' @param roundEarth boolean TRUE/FALSE. TRUE uses the Haversine formula to calculate the distance between pairs of populations on a globe;
#' FALSE uses Cartesian distance on a plane instead.
#' @param firstPC the number of the first principal component to use in calculating the unPC values; this is used as a column index
#' into the inputToProcess file(s).
#' @param secondPC the number of the second principal component to use in calculating the unPC values; this is used as a column index
#' into the inputToProcess file(s).
#' @param runPlotting boolean TRUE/FALSE. TRUE Runs the unPC plotting module (which must follow the pairwise calculation module); FALSE does not.
#' @param geogrCoordsForPlotting specifies the path to the file containing the geographic coordinates for each individual
#' represented in the input file(s). These geographic coordinates are used in generating the unPC plot ONLY (not calculating the unPC
#' scores). If not specified, the geogrCoords are used for both unPC calculation and plotting.
#' @param plotWithMap boolean TRUE/FALSE. TRUE includes the portion of the world map specified by the geographic coordinates used for plotting;
#' FALSE does not include any map
#' @param applyManualColor boolean TRUE/FALSE (this may be removed later). TRUE is the manual coloring to use for the manuscript; FALSE is dynamic
#' coloring based on the range of unPC values for the given dataset.
#' @param colorBrewerPalette is a string specifying a color palette in the RColorBrewer package. This is only used
#' if applyManualColor is FALSE, otherwise it is ignored.
#' @param ellipseWidth if specified is used as the width of the plotted ellipses instead of the default value. This can help to fix cases where
#' ellipses are wider than the length of the ellipse (especially a problem for small ellipses),
#' which causes them to appear like they are connecting populations other than the populutations
#' they are actually connecting. Setting this takes some trial and error based on the geographic range of the data being processed. Try 1 or 0.5
#' as a start.
#' @param populationPointNormalization float. This is a normalization (division) factor for the size of the points representing the number of individuals
#' sampled from each population. By default the normalization is 2, which halves the size of all points (all sizes will be divided by 2).
#' @param runAggregated boolean TRUE/FALSE. If inputToProcess is a directory containing multiple files, this controls whether the output
#' from all the files is averaged and then plotting (TRUE), or whether the output from each file is plotted individually (FALSE).
#' @param savePlotsToPdf boolean TRUE/FALSE. Whether to automatically save the plots to pdf (the default; TRUE), or to display the plots on the screen (FALSE)
#' @return None
#' @name unPC
#' @export


# CURRENTLY REQUIRES DEFINED FILE EXTENSION (.EVEC) FOR SEARCHING IN DIRECTORIES. MAY NEED TO RE-EVALUATE. 031717

# If input is a single file, currently the PCA results need 1 more row than the sample coords (because 1st row is skipped as a header for results straight from smartPCA).
# It also requires a leading column of sample identifiers (again for results from smartPCA)
# Can get around this by putting adding a dummy header as the first row in the PCA results file and a column of dummy sample names as the first column. May re-evaluate this import method. 042017.

# Need to implement a switch for mapping US states (requires the 'state' database instead of the 'world' database that is otherwise used - 'world' does not have state boundaries)? 042017

# This is un-PC standardized for msLandscape.

# Improved in 2 major ways -
# Get to cleanly select the geographic coords to use in un-PC calculation (required), and can optionally select different geographic coords to use
# with plotting (otherwise the coords used in the un-PC calculation are used for plotting as well)

# Can toggle whether to display a map over the un-PC visualizations (the portion of the world map is chosen based on the long/lat range of values in the
# geographic coords used for plotting). Can also toggle whether geographic dists used in un-PC calculation are Cartesian or on the globe using the Haversine formula.
# 031417


# abind makes it nice to bind matrices to a current 3-D array
# library(abind)
# library(ggplot2)
# library(plyr)
# library(shape)


# Need to clean up the dim==3 etc flags in pairwise calc.
# Prob if want to run PCA results from single file, use the coords for pairwise calc as the IMPORT coords, then
# use the popn aggregated coords for the pairwise calc coords, and can choose a different coord set for plotting.
# If want to use different...I don't know, need to think about more.

unPC <- function(inputToProcess, outputPrefix = "unPC_visualization", runDataImport = TRUE,
                 runPairwiseCalc = TRUE, geogrCoords, roundEarth = FALSE, firstPC = 1, secondPC = 2,
                 runPlotting = TRUE, geogrCoordsForPlotting = NULL, plotWithMap = FALSE, applyManualColor = FALSE,
                 colorBrewerPalette = NULL, ellipseWidth = NULL, populationPointNormalization = 2,
                 runAggregated = TRUE, savePlotsToPdf = TRUE){

start.time <- Sys.time()
#print(class(populationPointNormalization))
# ****************
# Setting up the modularity (ability to run any/all of the 3 contained modules independently)
# This is set up to be modular with running any of the three modules in it (modules 2 and 3 will only
# run if the whole pipeline has been run first, and will currently error out if not)

# These are boolean 0/1 flags for whether to run the modules or not.
#runDataImport <- 1
#runPairwiseCalc <- 1
#runMapping <- 1

# Boolean for whether to run the pairwiseCalc and runMapping to aggregate the results of each PCA run saved in the
# current folder, or whether to preserve the results for each run separately to generate a 3-D output matrix after
# pairwiseCalc (each layer is the pairwise dist calc for one PCA run), and create separate plots for each layer.
#runAggregated <- 1

#THIS NEEDS TO BE MADE REQUIRED FOR THE UN-PC CALCULATION; OPTIONALLY CAN USE A DIFF GEOGR COORD SET FOR THE PLOTTING.

# Boolean for whether to read the geographic coordinates for the different populations from an external file (1) - specified below, or
# to calculate them for the fullHexGrid completely within this script (0). Used to specify DCT sampling locations
#useGeogrCoordsFromFile <- 1

# Boolean for whether to use the Haversine formula to calculate the geographic distance between populations (1) or Cartesian distance (0)
#roundEarth <- 0

# Can change the PCAs to use in the pairwise comparisons here (PCs 1-10) - for the unPC pairwise diffs calc.
#firstPC <- 1
#secondPC <- 2

# Boolean for whether to use the geographic coords from file (above) to only calculate the unPC distance and use the regular fullHexGrid
# for the plotting, or whether to use the geographic coords for the plotting as well. This is only used if the useGeogrCoordsFromFile flag
# above is set to 1 for true.

#useGeogrCoordsFromFileForPlotting <- 1

# Should a map be plotted with the un-PC results
#plotWithMap <- 0

# To make the modularity work, need to specify the file names to use here:
# For import - need to specify the working dir. Then the .evec input files are found automatically
#               need to specify name of .r output file (will become input to pairwise distance calc)

# For pairwise dist calc - need to specify the name of the .r input file (which is the output from the evec import)
#                           need to specify the name of the .rds output file (which is the input to the mapper)

# For the mapper - need to specify the name of the .rds input file (which is the output from the pairwise dist calc)
#                   need to specify the name of the automatically saved plot files

# Path to the folder that contains the .evec files to process (and also can be used to modify where the various output
# files - described below) are saved.
#filePathToUse <- "~/Box Sync/ghostLands/ghostLands_PCAResultsForUnPC/ghostLands_fullHexGridWithBufferPopns_10IndPerPop_migrBarrierSimToEEMS/"
#filePathToUse <- "~/Box Sync/ghostLands/ghostLands_PCAResultsForUnPC/ghostLandsFullHexGrid_parallelSimAndAnalyze_091216_constMigr_m3/allPCAOutputConvertedForUnPC/"
#filePathToUse <- "~/Box Sync/ghostLands/ghostLands_PCAResultsForUnPC/ghostLands_fullHexGridWithBufferPopns_samplingLikeEEMSFig2CBottom_migrBarrierSimToEEMS_msRecode_030617/"


# This is fairly hacky seeming to tell whether the input to process is a directory of files to process or a single file

# Used to get the re to use for searching against the file names (will need to update for the msLandscape output format;
# don't worry too much about details; handled with more detail in the data import module below.)

reToUse <- glob2rx("*.evec")
#reToUse <- glob2rx("*.txt")

# For POPRES, and maybe all smartPCA output
# reToUse <- glob2rx("*.eigs")

fileListOutput <- list.files(path = inputToProcess, pattern = reToUse)

# The input is a dir with files to process
if(length(fileListOutput > 0)){
  #print(paste0("Length is:", length(fileListOutput)))
  setwd(inputToProcess)

  inputIsADir <- 1

} else{

  # Used below in file import separate from the list of file import if the input was a directory.
  inputIsADir <- 0

  # This would be a single file and needs to be handled here.
  #print(inputToProcess)

  # Parse the path from the input file name, and use that as the output path.
  # Returns a list
  initialSplit <- strsplit(inputToProcess, "/")

  # This removes the file name from the list of path elements and returns a char vector
  secondSplit <- initialSplit[[1]][1:lengths(initialSplit) - 1]

  # Construct the path - the collapse is essential
  pathToUse <- paste(secondSplit, collapse = "/")

  # Generate the file name of the population aggregated coordinates that will be created and saved
  # during the file import process. If the user has specified running all 3 modules of unPC at the same
  # time (and has input a single file, meaning that the input popn coords are individual-based instead of
  # popn-based like for dirs from the msLandscape pipeline), then after the individual coords are aggregated
  # into popn-level coords, these popn-aggregated coords are automatically
  # (and silently) used for pairwise calc. If the user only specifies running the pairwise calc module (with or without
  # graphing), then they have to enter this popn-aggregated coords to use, other

  # Split on '.' using escape characters
  splitPopnCoordsName <- strsplit(geogrCoords, "\\.")

  popnCoordsFileStem <- splitPopnCoordsName[[1]][1]

  popnAggregatedCoordsFileName <- paste0(popnCoordsFileStem, "_popnAggregated.txt")

  setwd(pathToUse)

}


# Make the various output file names. This also allows seamless toggling of the different
# modules (assuming the other modules have been run at least once with the currently specified outputPrefix.)
rOutputFileName <- paste0(outputPrefix, "_evecArray.r")

rdsRelPathPlusFileName <- paste0(outputPrefix, "_pairwiseDistCalc_unPC.rds")

plotTitleRoot <- paste0(outputPrefix, "_unPC_outputPlot_")


# If this is a single PCA results file and the data import has been set to run, then
# use the user-entered (individual-level) geogr coords for the data import, but pass the
# popn-level geogr coords (that are calculated in the data import module below, but already
# specified by a var defined above) to the pairwise
# unPC calc instead (with a message)
if(inputIsADir == 0 & isTRUE(runDataImport)){
  geogrCoordsForImport <- geogrCoords
  geogrCoords <- popnAggregatedCoordsFileName
  print("Running data import using the specified individual-level geographic coordinates file,")
  print("and running the pairwise difference calculation using the same coordinates after they have been aggregated by population.")
}

# setwd(filePathToUse)

# For saving the output evec and eval merged arrays (after the data import and before the diffs calc)
#rOutputFileName="./ghostLandsFullHexGrid_10IndPerPop_migrBarrierSimToEEMS_121416_evecsForUnPC_usingDynamicFlipping.r"
#rOutputFileName="./ghostLandsFullHexGrid_10IndPerPop_m3_121416_evecsForUnPC_usingDynamicFlipping.r"
# rOutputFileName="./ghostLands_fullHexGridWithBufferPopns_samplingLikeEEMSFig2CBottom_mmigrBarrierSimToEEMS_msRecode_030617_v2_evecsForUnPC_usingDynamicFlipping.r"


# This is the file name for saving the pairwise differences matrix (after the diffs calc and before the graphing).
# This needs to be the relative path + the filename
#rdsRelPathPlusFileName <- "./10IndPerPop_migrBarrierSimToEEMS_121516_DCTCorrected_popnSamplingTableForPlot.rds"
#rdsRelPathPlusFileName <- "./10IndPerPop_m3EvenSampling_121516_DCTCorrected_popnSamplingTableForPlot.rds"
# rdsRelPathPlusFileName <- "./ghostLands_fullHexGridWithBufferPopns_samplingLikeEEMSFig2CBottom_mmigrBarrierSimToEEMS_msRecode_030617_v2_DCTCorrected_popnSamplingTableForPlot_Aggregated.rds"

# This is a file of the geographic coordinates to use for each population (instead of the default calculated fullHexGrid popn points that are generated within this script)
# The file should have the x coords in col 1, y coords in col 2, and the population number (left to right, top to bottom) in col3 (with headers), and space delimited
# NOTE - CURRENTLY THIS FILE NEEDS TO HAVE HEADER 'x' for Longitude, and 'y' for Latitude and 'popNum_fromTopLeft' for the population numbers
# because later indexing depends on those names for the imported data.frame columns. This is also used almost directly (only x and y are interchanged) for plotting the location
# of the populations in the ggplot output.
# This makes the code more readable, but less flexible in terms of the input file formatting.
# NOTE - THE APPEARANCE OF THE GRAPHS FOR DIFFERENT COORDINATES IS CRITICALLY DETERMINED BY THE WIDTH OF THE ELLIPSES USED (THE 'ry' IN THE ELLIPSE MAKE CALL). IF THIS IS TOO
# LARGE THE GRAPHS CAN OBFUSCATE PATTERNS AND LOOK BIZZARE WITH THEIR GEOMETRY. GOOD RULE OF THUMB SO FAR - WIDTH OF OVALS SHOULD BE SET TO BE ~1/40 OF THE RANGE OF THE GEOGRAPHIC COORDS
# USED. This is now changed automatically based on the coordinate range of the data 10/12/16.

# For DCT Use below. (for now, need to manually edit DCT coords side by side the focal population coords that are successfully screened using my Python screen seqs script.)
# geogrCoordsFileToUseFor_unPCCalc <- "~/Box Sync/ghostLands/ghostLands_PCAResultsForUnPC/ghostLands_fullHexGridWithBufferPopns_samplingLikeEEMSFig2CBottom_migrBarrierSimToEEMS_msRecode_030617/ghostLands_EEMSFig2CBottom_DCTCoords_for_unPC.txt"

# geogrCoordsFileToUseFor_plotting <- "~/Box Sync/ghostLands/ghostLands_PCAResultsForUnPC/ghostLands_fullHexGridWithBufferPopns_samplingLikeEEMSFig2CBottom_migrBarrierSimToEEMS_msRecode_030617/ghostLands_Blank_10by10__haploidSamples_ms_nsam_1000_screened_1_time_popnCoordsFile_EEMSFig2CBot.txt"

# For PCA-based correction (using mean locations for the 100 iterations of m3) with the dynamic flipping sorting of PCA results to ensure
# they are all in proper register (ie popn 1 in upper left corner) 121416
#geogrCoordsFileToUse <- "~/Box Sync/ghostLands/ghostLands_PCAResultsForUnPC/ghostLands_meanPC1PC2Coords_from_m3_100iters_forDCTLike_unPCLocationCorrection_revised_121416.txt"

# This is the rel path plus the STEM ONLY of the output file name for the maps. This is also used to make the plot title, so it
# needs to be descriptive. Because
# both maps (most diff on top, and least diff on top) are generated together, the mapping
# script appends a '_mostDiffOnTop.pdf' or '_leastDiffOnTop.pdf' tag to the end of the file stem
# before saving the pdf.
#pdfMapRelPathPlusOutputFileName <- "./ghostLands_fullHexGridWithBufferPopns_10IndPerPop_migrBarrierSimToEEMS_121616_DCTCorr_with3StDevWideColorBins_"
#pdfMapRelPathPlusOutputFileName <- "./ghostLands_fullHexGridWithBufferPopns_10IndPerPop_m3_121616_DCTCorr_with3StDevWideColorBins_"
# pdfMapRelPathPlusOutputFileName <- "./ghostLands_fullHexGridWithBufferPopns_samplingLikeEEMSFig2CBottom_mmigrBarrierSimToEEMS_msRecode_030617_v2_DCTCorrected_with3StDevWideColorBins_Aggregated_DCTpointsPlotted"


# Finished setting up modularity; starting modules.
#*****************************

# --------------------
# Data importing

if(isTRUE(runDataImport)){

    # This script imports the results of multiple PCA runs (multiple .evec files; all the evec
    # files in the current directory) and parses them into a 3-D array (with the third dim = the
    # iteration number of the ms simulation, and constisting of the the .evec results matrix from
    # that iteration) to save as an .rds file which can then be used for easily plotting either all the results
    # or a specified subset.

    # This also parses the eigenvalues written in the first line of each .evec file, and
    # concatenates them into a matrix (each row is a different simulation iteration number)

    # This has to exclude any characters (ie the iteration numbers in the evec holder) when
    # making the combined array in order for it to be a numeric matrix and not a character matrix
    # (which ggplot can't handle with the plotting). This now (3/22/16) creates output that the plotter script
    # can use for plotting.

    # Merging these results with those from other independent runs (and runs of this script):
    # use:
    # require("abind")
    # holder <- abind(result1,result2,along=3)
    # iterNameHolder <- ""
    # for(iteration in seq(1,length(GaussianRF_EnvDiffResultsTabulate[,1,1]),1)){
    # iterNameHolder <- c(iterNameHolder,paste("Iter.#",iteration,sep=""))
    #}
    # dimnames(GaussianRF_EnvDiffResultsTabulate)[[1]] <- iterNameHolder[-1]
    # dimnames(GaussianRF_EnvDiffResultsTabulate)[[2]] <- c("Avg.Env.Value_Species1","Avg.Env.Value_Species2","Diff_in_Avg_Env_Value(Sp1-Sp2)")
    # dimnames(GaussianRF_EnvDiffResultsTabulate)[[3]] - names of different input files here.

    # 3/24/16 This is working correctly for the small lattice ghostLand tests for parameter tuning

    # This is a bit ugly. On the import, need to only import .evec files (not others in the dir), and need to import them in numeric
    # order based on their iteration # (so the layers of the 3-D array will be made in numeric order, too). This is currently implemented with 2 loops - 1 to
    # find just the .evec files, and parse their iteration numbers, and 1 to use the iteration numbers to import the data into
    # the 3-D array in numeric order of their iterations.

    # If the specified input was a single file (instead of a dir of files, and therefore is presumably a file of PCA results from outside
    # the msLandscape simulation toolbox). Import it almost identical to PCA results from msLandscape toolbox, except we cannot count on
    # and population identifying information in these PCA results like we can for the msLandscape simulation toolbox. Therefore we use the
    # population coordinates file (which has to list the coords for each individual, NOT population in this case (ie num rows == num indivs)
    # to generate population labels based on unique lat/lon combinations and then assign each individual one of the population labels).
    # The coords are then condensed to the population-level for use in the pairwise calculation and plotting (if another set of coords is not
    # specified).
    if(inputIsADir == 0){

      rawEvecTable <- read.table(inputToProcess, sep="", skip=1)
      # Drop the first column (the sample identifiers)
      PCAEvec <- rawEvecTable[,-1]

      # Read the top line of the PCA results (the eigenvalues) separately
      rawEvalVec <- read.table(inputToProcess,sep="\t",nrows=1)

      # These need to be lat lon pairs for each individual
      indivCoords <- read.table(geogrCoordsForImport, sep = "", header = FALSE)

      # Confirmed trigger 031617.
      if (nrow(indivCoords) != nrow(PCAEvec)){
        stop("ERROR. The number of rows (i.e. individuals) in the specified coordinate file does not match the number of rows in the PCA results. Exiting.")
      }

      names(indivCoords) <- c("Lat","Lon")
      # The unique populations are the unique combinations of lat and lon coordinates among the individuals.
      uniquePopnCoords <- unique(indivCoords[c("Lat","Lon")])

      # Sort the uniquePopnCoords descending for lat and ascending for lon (i.e. listing coords from top left
      # to bottom right, matching the ordering of the popn coords files from the msLandscape toolbox)
      sortedUniquePopnCoords <- uniquePopnCoords[order(-uniquePopnCoords[,1], uniquePopnCoords[,2]),]

      # Generate a file of population level aggregated coordinates (these are in the same order as the population numbers are assigned
      # e.g. the coordinates on the first line are the coordinates for the popn assigned to #1)
      write.table(x = sortedUniquePopnCoords, file = popnAggregatedCoordsFileName, sep = "\t",row.names = FALSE, col.names = FALSE)

      # Now we need to find which entry (row) of the sortedUniquePopnCoords corresponds to each individuals coordinates -
      # This row number matching in the sortedUniquePopnCoords is how the population numbers are assigned to each individual.
      # This matching is based on unique lat/lon combinations.
      popnNumbersForIndividuals <- match(paste(indivCoords[,1],indivCoords[,2]), paste(sortedUniquePopnCoords[,1], sortedUniquePopnCoords[,2]))

      # Now append this vector of populaton numbers to column 11 of the PCAEvec. NOTE - this popn vect needs
      # to be in column 11 for correct parsing for the pairwise calculation below. This still allows columns for the
      # top 10 eigenvectors to be included (since the first col in the output is dropped above.). This will replace
      # any existing column 11.

      # Pad the PCAEvec with columns of zero if it has less than 11 columns. This assumes that the PCs chosen to use are not outside of the PCs
      # that exist (as non-zero pads) in the input data (this is currently not a tested condition.)
      if (ncol(PCAEvec) < 11){
        for( colNum in seq(ncol(PCAEvec) + 1, 11, 1)){
          columnToAdd <- vector(mode = "numeric", length = nrow(PCAEvec))
          PCAEvec <- cbind(PCAEvec, columnToAdd)
        }
      }

      # Already checked in the 'stop' above to make sure this is safe.
      PCAEvec[,11] <- popnNumbersForIndividuals

      # Processing done. Now get them ready for saving below (using common var names).
      allEvecArray <- PCAEvec[,1:11]
      allEvalMatrix <- rawEvalVec
      #print(allEvecArray)


    }

    if(inputIsADir == 1){
      # generate a list of all the files in the folder
      fileList <- list.files(path=getwd(), full.names=F, recursive=FALSE)

      # Loop through the files, and process the ones that are evec files (end with "tabDelimForRPlotting-evec.txt" or any other
      # string ending defined below).

      # This is the original file suffix to search
      #fileSuffixToSearch <- "tabDelimForRPlotting-evec.txt"

      # This is more general
      fileSuffixToSearch <- ".evec"

      lengthFileSuffixToSearch <- nchar(fileSuffixToSearch)

      # This will be a named list acting like a dict in Python to store the parsed iter number (used below
      # for adding files in the correct, sorted order to the output holder) and the file name. This makes
      # the script more flexible to differences in file names.
      screenedFileList <- list()

      # Still need to vectorize with lapply. 8/15/16
      for(file in fileList){

        print(paste("The file name is:", file, sep = ""))

        # If the name of the current file ends with the correct file suffix for an evec file, then add its values
        # to the 3-D array of results.

        if( substr(file,(nchar(file) - lengthFileSuffixToSearch + 1),nchar(file)) == fileSuffixToSearch){
          #print("evec File")

          # This is the original parsing.
          #splitFile1 <- strsplit(file, "-")

          # This is more general
          splitFile1 <- strsplit(file, "_")

          # Indexing into the split file name string (split with both '-' and '_' to get the iteration number). This
          # makes it flexible to get the iter num regardless of how many digits it has.

          # This is the original parsing
          #iterNum <- strsplit(splitFile1[[1]][2], "-")[[1]][1]

          # This is more general - splits on the file extension (and assumes the iteration
          # number in the file name follows format XXX_XXXX_##.XXX where ## is the iteration number
          # and occurs immediately before the file extension (regardless of how many "_" characters
          # occur in the file name). "." has to be escaped.
          splitFile2 <- strsplit(splitFile1[[1]][length(splitFile1[[1]])], "\\.")
          iterNum <- splitFile2[[1]][1]

          # This evaluates the value of iterNum and puts it as a name in the list (this named list is functioning like a Python dict)
          screenedFileList[[iterNum]] <- file

        } else{
          print(paste0("WARNING - Skipping file: ", file, " because it does not appear to be an evec file (based on its file name suffix + extension)."))
        }
      }

      iterNumHolder <- names(screenedFileList)

      # Sort the iteration numbers of the .evec files in the folder. For this to sort in full numeric order, need to convert
      # to numeric for the sort, then convert back to character for use in the second for loop below.
      sortedIterNumHolder <- as.character(sort(as.numeric(iterNumHolder)))

      # Ideally I would do all this avoiding for loops using like lapply, but a for loop makes
      # it so much simpler to put the output into the holder array in numerical order of the ms iterations
      # (as opposed to the order the files are encountered by R in the directory, which is not numeric)

      # Still need to vectorize (prob with apply). 8/15/16
      for(num in sortedIterNumHolder){
        #print(num)
        # NEED TO HARD-CODE THIS PASTE COMMAND TO MATCH THE FILENAMES IN THE DIR AND TO INSERT
        # THE NUMBERS IN THE CORRECT SPOTS.

        # This is for the previous file names (before 7/16) with 2 "_Iter"s
        #currFileName <- paste("ghostLands_msSimConvertedForEigenstratPCA_Iter-",num,"_Iter_",num,"-tabDelimForRPlotting-evec.txt",sep="")

        # Another (fun) variation after re-running existing ms sim files through the revised ms output coder 022617
        #currFileName <- paste("ghostLands_migrBarr_evenSmpl_forPCA_Iter_Iter-",num,"-tabDelimForRPlotting-evec.txt",sep="")

        # This is for the high throughput simulated files (starting 7/16) that have 1 "_Iter"
        #currFileName <- paste("ghostLands_msSimConvertedForEigenstratPCA_Iter-",num,"-tabDelimForRPlotting-evec.txt",sep="")

        currFileName <- screenedFileList[[num]]

        #print(currFileName)
        # These are the eigenvectors only (skipping the first line of eigenvalues)
        rawEvecTable <- read.table(currFileName, sep="", skip=1)
        #print(rawEvecTable[1,])

        # Drop the first column (the sample identifiers)
        PCAEvec=as.matrix(rawEvecTable[,-1])

        # Reflect the PC1, and PC2 values as necessary so the points for population 1 are in the upper left
        # corner of the plot (ie negative x (PC1) values and positive y (PC2) values). 121416

        # Sometimes the first row has NA entries, so need to find the first row in both PCA 1 and PCA 2 (they really should
        # be the same, but this does not suppose that), and use that to test the values for whether transposing is necessary.
        firstNonNARowPC1 <- which(!is.na(PCAEvec[,1]))[1]
        firstNonNARowPC2 <- which(!is.na(PCAEvec[,2]))[1]

        rowNumForReflectionTest <- ifelse(firstNonNARowPC1 >= firstNonNARowPC2, firstNonNARowPC1, firstNonNARowPC2)

        if(PCAEvec[rowNumForReflectionTest,1] < 0 & PCAEvec[rowNumForReflectionTest,2] > 0){
          processedPCAEvec <- PCAEvec
        }

        if(PCAEvec[rowNumForReflectionTest,1] < 0 & PCAEvec[rowNumForReflectionTest,2] < 0){
          processedPCAEvec <- cbind(PCAEvec[,1], -1*PCAEvec[,2], PCAEvec[,3:11])
          #print("Testing.")
        }

        if(PCAEvec[rowNumForReflectionTest,1] > 0 & PCAEvec[rowNumForReflectionTest,2] > 0){

          processedPCAEvec <- cbind(-1*PCAEvec[,1], PCAEvec[,2:11])
        }

        if(PCAEvec[rowNumForReflectionTest,1] > 0 & PCAEvec[rowNumForReflectionTest,2] < 0){
          processedPCAEvec <- cbind(-1*PCAEvec[,1], -1*(PCAEvec[,2]), PCAEvec[,3:11])
          #print("Test.")
        }

        # Create the array
        if(num == 1){
          #print(paste(num,nrow(processedPCAEvec)))
          allEvecArray <- processedPCAEvec
        } else{
          #print(paste(num,nrow(processedPCAEvec)))
          allEvecArray <- abind::abind(allEvecArray, processedPCAEvec, along=3)
          print("appended.")
        }

        # Read in the eigenvals. from the first row of the evec file (these are the
        # variances explained by each eigenvect.) The comment.char = "" lets it ignore
        # The leading '#' output by smartPCA.
        rawEvalVec <- read.table(currFileName,sep="",nrows=1, comment.char = "")
        PCAEval <- as.matrix(rawEvalVec[-1])

        if(num == 1){
          allEvalMatrix <- PCAEval
        } else{

          allEvalMatrix <- rbind(allEvalMatrix,PCAEval)
        }

      }

      # Add the row, column, and 3rd dimensional names for the 3-D array of evec values.
      dimnames(allEvecArray)[[1]] <- paste("Individual#", seq(1,nrow(allEvecArray[,,1])))
      dimnames(allEvecArray)[[2]] <- c(paste("Evec#", seq(1,10)),"Population#")
      print(dimnames(allEvecArray)[[3]])
      dimnames(allEvecArray)[[3]] <- paste("ms Sim. Iteration#", sortedIterNumHolder)
    }

    save(allEvecArray,allEvalMatrix, file = rOutputFileName)

}
# End Data read, parse, and save
#--------------------------------------

#--------------------------------------
# Start the pairwise geographic and PCA distance calculation module

# This is adapted from my Novembre et al 2008 PCA dist pairwise calculator to work with the
# ghostLands fullHexGrid (10 x 10 hexagonal grid of sampled populations each separated by unsampled ghost
# populations in 6 directions and the whole grid bordered by 2 fully ghost 'buffered' rings of ghost populations.)

# This calculates the average pairwise differences in both PCA dist and 'geographic' dist for each sampled population
# in the grid, which will be used to generate the map of genetic discontinuities using the ratio of PCA dist to geographic
# dist (just like for my Novembre et al 2008 PCA genetic diff mapper)

# This calculates the correct 'geographic' coordinates of the tiles of the simulation based on a 10 x 10 grid of focal populations
# (ie the fullHexGrid). These coordinates are defined in the same way as they are defined in my ggplot plotter for this
# grid (which is how I will be mapping the PCAs onto the grid, so I can use this output directly for the mapping)

# All for loops replaced with apply statments

# This seems to be working right with the runAggregated flag choices 7/22/16

if(isTRUE(runPairwiseCalc)){

    # This is the PCA data (from Eigenstrat run as part of my simulation/analysis pipeline on Karst). This is the
    # parsed PCA results data after running my script: ghostLands_EigenstratPCAEvecResultsParse.R (which also prepares
    # the PCA results for direct plotting using the separate script). This loads variables 'allEvalMatrix' and 'allEvecArray'.
    # The 'allEvalMatrix' is actually a row vector of the eigenvectors for the PCA.
    # The 'allEvecArray' is either a matrix (if 1 run of ms simulation used), or a 3-D array (if >1 run of the ms simulation used)
    load(rOutputFileName)

    # Plotting the raw PCA1 PCA2 values for all iterations with very low alpha value
    #plot(allEvecArray[,1],allEvecArray[,2], xlab = "PC1", ylab = "PC2", main = "Raw PCA results from constant migration scenario (all 100 iterations)", col = rgb(0,0,0,0.1), cex.lab = 1.2, cex.main = 1.2)

    #Import the geographic coordinates to use for calculating the unPC scores. This should have no header and be in 2 cols (whitespace delim) with the first
    # col lat, and the 2nd col lon.
    popnAggregateGeogrCoords <- read.table(geogrCoords, sep = "", header = FALSE)

    names(popnAggregateGeogrCoords) <- c("y","x")

    # Plotting the transformed coordinates
    #plot(popnAggregateGeogrCoords$x, popnAggregateGeogrCoords$y, xlab = "Expected PC1", ylab = "Expected PC2", main = "PCA with constant migr and even sampling expected locations for PC1 PC2 points", col = rgb(0,0,0), cex.lab = 1.2, cex.main = 1.2)


    # This is the function that calculates the pairwise comparison values for each of the
    # comparisons. It is called using apply() below on the pairwiseCompareHolder

    # 7/18/16 Indexing directly from popAveregedPCAEvecs like it currently does is MUCH faster
    # than creating a copy of the current 2 rows being considered, and then indexing in those.

    # This is the function used to tabulate the different calculations needed for each pairwise comparison (this is called by apply()). This is the same
    # function is used to tabulate the pairwise diffs regardless of whether the tabulations should be aggregated across different PCA results (specified by the runAggregated flag)
    # at the start of the script. It does this by being called in different contexts (different apply calls) depending on whether they should be aggregated or not.
    tabulateMetrics <- function(array, popAggregatedData){
        #print(paste(array[1], array[2]))
        popNum1 <- array[1]
        popNum2 <- array[2]

        #print(paste(popNum1,popNum2))
        # These lat/lons are from the external file with geographic coordinates to use for the unPC calculations. Get the lat/lon for the
        # two populations currently being compared.
        lat1 <- popnAggregateGeogrCoords$y[popNum1]
        lat2 <- popnAggregateGeogrCoords$y[popNum2]

        lon1 <- popnAggregateGeogrCoords$x[popNum1]
        lon2 <- popnAggregateGeogrCoords$x[popNum2]

        # Calculates the geodesic distance between two points specified by radian latitude/longitude using the
        # Haversine formula (hf)
        haversineDist <- function(lat1, lat2, lon1, lon2) {
            # Earth radium (km)
            R <- 6371
            radLon1 <- lon1*(pi/180)
            radLon2 <- lon2*(pi/180)
            radLat1 <- lat1*(pi/180)
            radLat2 <- lat2*(pi/180)
            #print('**********')
            #print(paste(lon1,lon2,lat1,lat2))
            longDiff <- (radLon2 - radLon1)
            latDiff <- (radLat2 - radLat1)
            #print(paste(delta.long,delta.lat))
            a <- sin(latDiff/2)^2 + (cos(radLat1) * cos(radLat2) * sin(longDiff/2)^2)
            #print(a)
            c <- 2 * asin(min(1,sqrt(a)))
            #print(c)
            d = R * c
            #print(d)
            return(d) # Distance in km
        }


        # Calculate the geographic distance between the current pair of points either on a globe or a plane.
        if(isTRUE(roundEarth)){
            geogrDist <- haversineDist(lat1,lat2,lon1,lon2)
        } else{
            geogrDist <- sqrt((lat1 - lat2)^2 + (lon1 - lon2)^2)
        }


        #first and second PCs (usually 1 and 2) for the first population in the current comparison.
        firstPCA1 <- popAggregatedData[popNum1,firstPC]
        firstPCA2 <- popAggregatedData[popNum1,secondPC]

        #first and second PCs (usually 1 and 2) for the second population in the current comparison.
        secondPCA1 <- popAggregatedData[popNum2,firstPC]
        secondPCA2 <- popAggregatedData[popNum2,secondPC]

        PCADist <- sqrt((firstPCA1 - secondPCA1)^2 + (firstPCA2 - secondPCA2)^2)


        # This ratio of the PCA dist to the geographic dist for the 2 populations currently being compared is their un-PC score
        unPCValue <- PCADist / geogrDist

        # Build a vector to hold the relevant information (especially the unPCValue) for the pair of populations currently being compared.
        pairwiseDiffHolder <- cbind(popNum1,popNum2,lon1,lat1,lon2,lat2,geogrDist,firstPCA1,firstPCA2,secondPCA1,secondPCA2,PCADist,unPCValue)

        # Return the pairwiseDiffHolder upon each function call; apply automatically concatenates the tempHolder vectors it gets returned to it.
        return(pairwiseDiffHolder)


    }



    #---------------------
    # This is for aggregating all results together (by taking avg) before calculating the pairwise diffs

    if(isTRUE(runAggregated)){

        # For testing
        #allEvecArray <- allEvecArray[1:10,]

        # For 3-D input arrays (ie from multiple ms runs), could turn them into a single very long matrix here, then do the
        # aggregate step on that to effectively average all the runs.

        # * Try to change to dplyr later. 7/18/16
        #require("plyr")

        # This converts the x by y by z 3-D array from multiple ms runs into an x*z by y matrix for use in all subsequent steps
        # The .id = NULL is necessary to keep plyr from adding a column keeping track of which z dimension matrix each row
        # in the re-shaped matrix came from. This otherwise prevents it from being able to convert into a numerical matrix.
        # This is otherwise working correctly.
        if(length(dim(allEvecArray)) == 3){

            # Make a vector of counts of each population sampled in the first matrix of the 3-D 'stack'.
            # This will be used for plotting population circles with sizes proportional to the number of individuals sampled in each population.

            # Although all sampling should be the same for each matrix in the 3-D array,
            # for very uneven sampling with strong migration patterns, smartPCA may drop PCA scores for all individuals from some of the sampled populations.
            # This looks through the summary tables (vectors) of counts of individuals in each population from each iteration run and catalogs the length (ie number of entries in each).
            # The assumption is that
            # the longest of the summary tables represents the true population sampling with no missing populations. This appears to work fairly well for
            # large numbers of iterations.  NOTE : If not all populations are represented even from the longest table though, this will give an error when plotting.

            calcTableLength <- function(inputMatrix){
                tableLength <- length(table(inputMatrix[,11]))
                return(tableLength)
            }

            tableLengths <- apply(X = allEvecArray, MARGIN = 3, calcTableLength)

            representativeTableIndex <- which(tableLengths == max(tableLengths))[1]

            samplingTableForPlotting <- table(allEvecArray[,11,representativeTableIndex])
            allEvecArray <- as.matrix(plyr::adply(allEvecArray,3,.fun = NULL, .id = NULL))

        } else{
            samplingTableForPlotting <- table(allEvecArray[,11])
        }

        # Average each column of the allEvecArray by population (ie clump it on a population level, not the individual level
        # that it initially has.) popAveragedPCAEvecs is a dataframe with 12 columns (aggregate adds a groups column as column 1)
        popAveragedPCAEvecs <- aggregate(allEvecArray[,1:11], by = list(allEvecArray[,11]),FUN=mean, na.rm = TRUE)

        # The aggregate call above adds a 'groups' column before any of the data. In that case, subset the aggregated Evecs
        # back to the original number of cols (remove the groups column)
        if (ncol(popAveragedPCAEvecs) == 12){
            popAveragedPCAEvecs <- popAveragedPCAEvecs[,2:12]
        }

        # The number of populations is the number of rows in the popAveragedPCAEvecs
        numPopns <- nrow(popAveragedPCAEvecs)


        # Revised code with no for loops (1 apply() used with custom function instead)

        # Start by creating a 2-column matrix of the different pairwise comparisons that need to be made (popn1 in col1, popn2 in col2).
        # NOTE - this is made using the actual population numbers used as entries in the column, NOT the length of the column, and therefore
        # is robust to missing population numbers

        rawCompareHolder <- cbind(rep(popAveragedPCAEvecs[,11],each=numPopns),rep(popAveragedPCAEvecs[,11] ,numPopns))

        # Now keep only the unique pairwise comparisons (non-self comparisons too, so like lower triangular part of the
        # comparison matrix). These are the comparisons where the entry in col2 is > the entry in col1.
        pairwiseCompareHolder <- rawCompareHolder[rawCompareHolder[,1] < rawCompareHolder[,2],]



        # Apply the tabulateMetrics function over margin 1 (rows) of the pairwiseCompareHolder.
        # This is taking the place of a for loop; because it is working row-wise, it passes both
        # values for the current row (ie the numbers of the two populations that are currently
        # being compared) to the tabulateMetrics function, which is then able to unpack the values
        # using indexing to then use the values in turn as row indexes to define the two rows to
        # use for tabulating the pairwise comparison.

        # Running the tabulations using apply() instead of even a single for loop over the pairwiseCompareHolder
        # instead (not nested for loops), runs about 3x faster (~6secs with apply vs. ~ 17 secs with the for loop)


        # apply automatically concatenates (cbinds, not rbinds) the individual vector of values it has returned to it by the tabulateMetrics function each time it is called into a matrix
        # of output. Need to pass the popAveragedPCAEvecs as additional input because the function is defined above
        # (and before popAveragedPCAEvecs exists.)
        rawOutput <- apply(pairwiseCompareHolder, 1, tabulateMetrics, popAggregatedData = popAveragedPCAEvecs)


        pairwiseDiffsOutput <- as.data.frame(t(rawOutput))

        colnames(pairwiseDiffsOutput) <- c("population_1", "population_2", "x_1", "y_1", "x_2", "y_2",
                                            "geogrDist", "PC1_1", "PC2_1", "PC1_2", "PC2_2", "PCDist",
                                            "ratioPCToGeogrDist")

        saveRDS(pairwiseDiffsOutput,rdsRelPathPlusFileName)
        # Strip the ".rds" extension off the rds saving file name to build the name for saving the sampling vector for plotting
        rdsFileStem <- substr(x = rdsRelPathPlusFileName, start = 1, stop = nchar(rdsRelPathPlusFileName) - 4)
        saveRDS(samplingTableForPlotting, paste(rdsFileStem, "_popnSamplingTableForPlot.rds", sep = ""))
    }
    #----------------------------------------------------

    #----------------------------------------------------
    # This is for calculating the pairwise diffs for each set of output results separately (not aggregated). Each set of results
    # is represented by a 'layer' in the 3rd dimension of the allEvecArray

    if(!isTRUE(runAggregated)){
        print("Running without aggregation.")
        #numPopns <- nrow(allEvecArray[,,1])

        aggregateIterationCount <- 0


        #For testing ONLY -
        # subset the allEvecArray to it has 10 entries in the 3rd dimension instead of 100 (for speed of debugging)
        # allEvecArray <- allEvecArray[,,1:3]

        # This tests for the presence of rows with NAs in them across all iterations and returns the count of NA entries for each entry (to be turned into a boolean flag and in very special cases where all iterations contain some NA rows, the counts
        # will choose the best iteration to try for making the samplingTable. For each iteration, some rows containing NAs are fine as long as there are other rows from the same population
        # that do not have any NAs and therefore are possible to aggregate across. But if all inidividuals from the same population only have NA entries, then the aggregation will fail for that
        # population and will return a shorter output than expected. Use this simple testMatrix function to determine the number of populations to expect from the population aggregation
        # based on one of the iterations where there are no NA rows and therfore the max number of represented populations.
        testMatrixForNA <- function(inputMatrix) {
            sumNAEntries <- sum(is.na(inputMatrix[,11]))
            return(sumNAEntries)
        }

        #print(length(dim(allEvecArray)))


        # If the sum of the hasNA '1' flags is less than the length, then at least 1 iteration has all of its individuals without NAs, and can use that to find
        # the number of represented populations

        if(length(dim(allEvecArray)) == 3){

          numIterations <- dim(allEvecArray)[3]

          sumNAPerIterationHolder <- apply(allEvecArray, 3, testMatrixForNA)
          hasNAHolder <- ifelse(sumNAPerIterationHolder > 0, 1, 0)

          if(sum(hasNAHolder) < length(hasNAHolder)){
              #print("In sum.")
              noNARepresentativeIter <- which(hasNAHolder == 0)[1]
              numRepresentedPopns <- nrow(aggregate(allEvecArray[,1:11, noNARepresentativeIter], by = list(allEvecArray[,11, noNARepresentativeIter]),FUN=mean, na.rm = TRUE))

              # This is a vector of counts of each population sampled.
              # This will be used for plotting population circles with sizes proportional to the number of individuals sampled in each population.
              samplingTableForPlotting <- table(allEvecArray[,11,noNARepresentativeIter])

          } else{
              print("WARNING - ALL OF THE ITERATIONS HAVE AT LEAST 1 INDIVIDUAL (ROW) IN THEIR RESULTS WITH ALL NA ENTRIES.")
              print("This will fail if 1 of the populations in 1 of the iterations has all NA values for every individual.")
              print("If this runs, it means that the NA entries only affect isolated individuals from each population in each iteration.")

              iterationWithMinNumNAs <- which(sumNAPerIterationHolder == min(sumNAPerIterationHolder))[1]

              numRepresentedPopns <- nrow(aggregate(allEvecArray[,1:11, iterationWithMinNumNAs], by = list(allEvecArray[,11, iterationWithMinNumNAs]),FUN=mean, na.rm = TRUE))

              # This is a vector of counts of each population sampled.
              # This tries to make the sampling table with the best chance of success by using the iteration number that contains the fewest number of NA
              # row entries; this still may fail despite this (see print message above).
              samplingTableForPlotting <- table(allEvecArray[,11,iterationWithMinNumNAs])

          }


        } else if (length(dim(allEvecArray)) == 2){

          numIterations <- 1

          numRepresentedPopns <- nrow(aggregate(allEvecArray[,1:11], by = list(allEvecArray[,11]),FUN=mean, na.rm = TRUE))

          # This is a vector of counts of each population sampled.
          samplingTableForPlotting <- table(allEvecArray[,11])

        }



        keptMatrixList = list()

        # This aggregation step is purposefully a for loop and not an apply call (there is an apply call to tabulateMetrics for each iteration of the for loop).
        # This allows the aggregation to skip aggregating the PCA output file for any iteration where the PCA results for all inidividuals of a given population have entries of 'NA',
        # which otherwise causes an error, and allows easier counting of the progress of processing the different iterations.

        for(iteration in seq(1,numIterations,1)){

            aggregateIterationCount <- aggregateIterationCount + 1
            print(paste0("Processing individual: ", aggregateIterationCount))

            if(length(dim(allEvecArray)) == 3){
              inputArrayToAggregate <- allEvecArray[,,iteration]
            } else if (length(dim(allEvecArray)) == 2){
              inputArrayToAggregate <- allEvecArray
            }

            # Average each column of the allEvecArray by population (ie clump it on a population level, not the individual level
            # that it initially has.) popAveragedPCAEvecs is a dataframe with 12 columns (aggregate adds a groups column as column 1)
            popAveragedPCAEvecs <- aggregate(inputArrayToAggregate[,1:11], by = list(inputArrayToAggregate[,11]),FUN=mean, na.rm = TRUE)

            # The aggregate call above adds a 'groups' column before any of the data. In that case, subset the aggregated Evecs
            # back to the original number of cols (remove the groups column)
            if (ncol(popAveragedPCAEvecs) == 12){
                popAveragedPCAEvecs <- popAveragedPCAEvecs[,2:12]
            }

            #print(popAveragedPCAEvecs)

            # The number of populations is the number of rows in the popAveragedPCAEvecs
            numPopns <- nrow(popAveragedPCAEvecs)

            # Check if the number of populations in the current iteration matches the expected number (from an iteration with no populations that are only represented by 'NA' entries for all individuals)
            # Process it if true. Give a warning message about the iterations that do not meet this criteria and are therefore not processed further (they would cause an error for the whole batch of files otherwise)
            if(numPopns == numRepresentedPopns){

                # Now generate the indices of the pairwise comparisons using the numPopns

                # Start by creating a 2-column matrix of the different pairwise comparisons that need to be made (popn1 in col1, popn2 in col2).
                # NOTE - this is made using the actual population numbers used as entries in the column, NOT the length of the column, and therefore
                # is robust to missing population numbers

                rawCompareHolder <- cbind(rep(popAveragedPCAEvecs[,11],each=numPopns),rep(popAveragedPCAEvecs[,11] ,numPopns))

                # Now keep only the unique pairwise comparisons (non-self comparisons too, so like lower triangular part of the
                # comparison matrix). These are the comparisons where the entry in col2 is > the entry in col1.
                pairwiseCompareHolder <- rawCompareHolder[rawCompareHolder[,1] < rawCompareHolder[,2],]

                # apply automatically concatenates (cbinds, not rbinds) the individual vector of values it has returned to it by the tabulateMetrics function each time it is called into a matrix
                # of output. Need to pass the popAveragedPCAEvecs to it as an additional parameter. This is looping through the pairwise comparisons and tabulating them to return
                # to the calling apply() call.
                rawOutput <- apply(pairwiseCompareHolder, 1, tabulateMetrics, popAggregatedData = popAveragedPCAEvecs)

                # Add the output matrix to the list.
                keptMatrixList[[iteration]] <- t(rawOutput)

            } else{
                print(paste0("WARNING: Iteration #", iteration, " has been dropped from the pairwise differences calculation (and graphing if selected) because it contains one or more populations with 'NA' entries for every individual in the population."))

            }

        }


        # Else the matrix will have 13 columns instead; This works fine if there is only 1 iteration to process.
        pairwiseDiffs3DArrayOutput <- abind::abind(keptMatrixList, along=3)

        dimnames(pairwiseDiffs3DArrayOutput) <- list(as.numeric(seq(1,nrow(pairwiseDiffs3DArrayOutput[,,1]),1)),
                                                        c("population_1", "population_2", "x_1", "y_1", "x_2", "y_2",
                                                          "geogrDist", "PC1_1", "PC2_1", "PC1_2", "PC2_2", "PCDist",
                                                          "ratioPCToGeogrDist"),
                                                         paste("resultsFrom_PCAIteration#",seq(length(pairwiseDiffs3DArrayOutput[1,1,])),sep=""))



        saveRDS(pairwiseDiffs3DArrayOutput, rdsRelPathPlusFileName)
        # Strip the ".rds" extension off the rds saving file name to build the name for saving the sampling vector for plotting
        rdsFileStem <- substr(x = rdsRelPathPlusFileName, start = 1, stop = nchar(rdsRelPathPlusFileName) - 4)
        saveRDS(samplingTableForPlotting, paste(rdsFileStem, "_popnSamplingTableForPlot.rds", sep = ""))

        #print(pairwiseDiffs3DArrayOutput)
        #print(class(pairwiseDiffs3DArrayOutput))
        #print(dim(pairwiseDiffs3DArrayOutput))

    }

}
# End of the pairwise distance calculation module
#-----------------------------------------

#-----------------------------------------
# Start of the graphing module

if(isTRUE(runPlotting)){

    #library(ggplot2)
    # the shape library is for the ellipse generation
    #library(shape)

    # This takes matrix output from the ghostLands_fullHexGrid_PCA_Geographic_DiffsCalc.R script.
    # The matrix is the pairwise geographic and PCA dists between each 'tile' on the fullHexGrid.
    # The coordinates in the matrix are set up for direct plotting here. This layers the connecting ellipses on
    # the plot and colors them by (PCA dist/ geogr dist) values, like in my plotter for the Novembre et al 2008 data.

    # Use the variable specifying the path to save the .rds file of pairwise distances, specified in the previous module,
    # to load the .rds file for plotting here

    pairwiseTileDiffHolder <- readRDS(rdsRelPathPlusFileName)

    # Load in the sampling vector of number of individuals sampled from each population, which will be used to make the size of the plotted point for
    # each population proportional to the number of individuals sampled.
    rdsFileStem <- substr(x = rdsRelPathPlusFileName, start = 1, stop = nchar(rdsRelPathPlusFileName) - 4)
    samplingTableForPlotting_df <- as.data.frame(readRDS(paste(rdsFileStem, "_popnSamplingTableForPlot.rds", sep = "")))

    # If the input data is in a 3-D array (signalling that it is the result of calculating the pairwise dists separately for each PCA results iteration),
    # then set runAggregated to 0 for false, regardless of whether the user set it to 1 (give a warning if this is the case),
    # because the aggregation has to happen in the calculation module,
    # and cannot happen here in the plotter.
    if(length(dimnames(pairwiseTileDiffHolder)) == 3 & isTRUE(runAggregated)){

        print("WARNING - runAggregated is set to TRUE for the plotting, but the results were")
        print("not aggregated across all simulation iterations during the pairwise calculation step for this dataset which is")
        print("required for aggregated plotting. Plotting the non-aggregated data instead (i.e. one set of plots per iteration.")
        runAggregated = 0

    }


    # Making the name for the plot title using the root of the filename (Iter# will be appended if runAggregated == 0)
    # plotTitleRoot <- strsplit(pdfMapRelPathPlusOutputFileName, "./")[[1]][2]

    #plotTitle <- "ghostLands_6216_fullHexGridWithBufferPopns_verticalRiverConnecter_sampledLikeEEMSFig2CBotRotated90DegCW_PC1_PC2"

    #-----------------------------
    generateUnPCPlot <- function(currentIterationNumber, pairwiseTileHolder){

        # If this is part of a 3-D array (ie multiple PCA outputs from independent simulations), then
        # extract the current layer of the array (passed by the lapply call to this function) and convert
        # from matrix to a data frame (needs to be data frame for correct graphing.)
        if(length(dimnames(pairwiseTileHolder)) == 3){

            pairwiseTileDiffs <- as.data.frame(pairwiseTileHolder[,,currentIterationNumber])

        } else{

            pairwiseTileDiffs <- pairwiseTileHolder
        }

        # If there is an entry in the file to use specifically for the geographic coordinates for plotting, then do that; otherwise use the
        # geographic coordinates from the file used to calculate the pairwise diffs.
        if (!is.null(geogrCoordsForPlotting)){
          # Use the population numbers for each comparison (popn 1 and popn2 separately) as indices into the
          # vectors of x and y coords for each population from the external file (max popn# === the number of coords lines in the external file)
          # The coords in the external file have to be (y,x) as in (long, lat)

          coordsForPlotting <- read.table(geogrCoordsForPlotting, sep = "", header = FALSE)

          yCoordPopn1 <- coordsForPlotting[,1][pairwiseTileDiffs$population_1]
          xCoordPopn1 <- coordsForPlotting[,2][pairwiseTileDiffs$population_1]
          yCoordPopn2 <- coordsForPlotting[,1][pairwiseTileDiffs$population_2]
          xCoordPopn2 <- coordsForPlotting[,2][pairwiseTileDiffs$population_2]

          pairwiseCoordsForPlotting <- data.frame("y_1" = yCoordPopn1, "x_1" = xCoordPopn1, "y_2" = yCoordPopn2, "x_2" = xCoordPopn2)

          # Need to assert that the number of rows in coords ForPlotting matches the number of rows in the samplingTableForPlotting (which was made previously
          # in the pairwise calc step)
          if(nrow(coordsForPlotting) != nrow(samplingTableForPlotting_df)){
            stop(paste0("ERROR. The number of populations (i.e. rows) in the specified coordinates file for plotting does not match the number of populations tabulated in the pairwise calculation step, which was: ", nrow(samplingTableForPlotting_df),". Exiting."))
          }

          populationCoordHolder <- data.frame("Lat" = coordsForPlotting[,1], "Lon" = coordsForPlotting[,2],"numIndivsSampled" = samplingTableForPlotting_df$Freq)

          # And generate a very basic boolean color vector to make the colors for the unsampled and the sampled populations different.
          populationCoordHolder$popColoring <- ifelse(populationCoordHolder$numIndivsSampled == 0,0,1)


        } else{

          #Import the geographic coordinates that were also used for calculating the unPC scores. This should have no header and be in 2 cols (whitespace delim) with the first
          # col lat, and the 2nd col lon.
          coordsForPlotting <- read.table(geogrCoords, sep = "", header = FALSE)

          populationCoordHolder <- data.frame("Lat" = coordsForPlotting[,1], "Lon" = coordsForPlotting[,2],"numIndivsSampled" = samplingTableForPlotting_df$Freq)

          # And generate a very basic boolean color vector to make the colors for the unsampled and the sampled populations different.
          populationCoordHolder$popColoring <- ifelse(populationCoordHolder$numIndivsSampled == 0,0,1)

          # Because using the same coords for plotting as were used for the pairwise calculations, can use the coordinates directly from the pairwise calculation output.
          pairwiseCoordsForPlotting <- data.frame("y_1" = pairwiseTileDiffs$y_1, "x_1" = pairwiseTileDiffs$x_1, "y_2" = pairwiseTileDiffs$y_2, "x_2" = pairwiseTileDiffs$x_2)

        }


        #numTileRows <- max(pairwiseTileDiffs$y_1)
        # Because of the staggering hex tiles, there are twice as many indices as actual tiles in any given column
        #numTileCols <- max(pairwiseTileDiffs$x_1)/2

        minPCADist <- min(pairwiseTileDiffs$PCDist)
        maxPCADist <- max(pairwiseTileDiffs$PCDist)

        PCAGeogrRatioDist <- pairwiseTileDiffs$ratioPCToGeogrDist

        # The metric to use to calculate the bins for the different colors. The actual
        # metric to use to determine the correct color for each pairwise comparison is
        # defined in the for loop below. If want to change the coloring metric (ie. PCA dist only
        # or PCADist:GeogrDist, then need to change it both here and in the colorDist variable in the
        # for loop)
        minColorDist <- min(PCAGeogrRatioDist)
        maxColorDist <- max(PCAGeogrRatioDist)
        #
        #spanColorDist <- maxPCADist - minPCADist
        spanColorDist <- maxColorDist - minColorDist

        # Choose the number of distinct color classes to use for plotting the
        # ellipses
        #       numColorGroups <- 30

        # The PCA distance for each of the equally spaced color groups in the
        # range from min to max PCA dist.
        #spanColorGroups <- spanColorDist / numColorGroups


        if(isTRUE(applyManualColor)){
          # -----------------------

          # This is for making any pairwise connection with an unPC value within x stdevs of the mean unPC value
          # for DCT corrected values from the m3 scenario white
          # and then defining 5 x stdev bins on either side of that distribution for coloring progressively
          # darker to red and green. For the DCT used as of 121516 (Made with DCT of double the actual num of popn rows)
          # a width of 3 stdevs seems good.
          # These are the mean and stdev unPC values from the m3 data analyzed using the DCT correction.

          # These are mean and stdev values to use for coloring BEFORE ms file re-coding (before ~022517)
          #stDevUnPCVal <- 0.003830205
          #meanUnPCVal <- 0.04473611

          # These are mean and stdev values to use for coloring AFTER ms file re-coding 022717
          meanUnPCVal <- 0.04437096
          stDevUnPCVal <- 0.0049017


          # The number of color groups on either side of the middle (white) group
          numColorGroups <- 5

          spanColorGroups <- stDevUnPCVal

          # it's pretty neat that this sort of dynamic addition/subtraction of an arbitrary number of values in a sequence is possible. This puts together the vector from highest to lowest
          # value, but findIntervals needs it increasing, so reverse it below
          # Can change the width of the bins (in units of standard deviations of the unPC value distribution from the DCT correction
          # using the m3 migration scenario) below. This assembles the intervals on both sides of the mean, so half the desired interval
          # width should be added/subtracted to the mean value below, then the full width added/subtracted to all subsequent bins.
          maxBreaksVect <- c(meanUnPCVal + (1.5 * stDevUnPCVal), meanUnPCVal + (1.5 * stDevUnPCVal) + (seq(1,numColorGroups,1) * (3*stDevUnPCVal)))
          minBreaksVect <- c(meanUnPCVal - (1.5 * stDevUnPCVal), meanUnPCVal - (1.5 * stDevUnPCVal) - (seq(1,numColorGroups,1) * (3*stDevUnPCVal)))
          colorBreaksVect <- c(rev(minBreaksVect), maxBreaksVect)

          colorBrewerColors <- c("#8E0152","#C51B7D","#DE77AE","#F1B6DA","#FDE0EF","#F7F7F7","#E6F5D0","#B8E186","#7FBC41","#4D9221","#276419")

          # Add orange for the highest interval - Not necessary for DCT correction, only empirical m3 corrections.
          #colorBrewerColors <- c("#FFA500", colorBrewerColors)

          # Need to reverse them to put the red as most differentiated
          colorBrewerColors <- rev(colorBrewerColors)

          # ----------------------------------------
        } else{

          if(!is.null(colorBrewerPalette)){

            # Check that RColorBrewer is installed. Verified this triggers 031717
            if (!requireNamespace("RColorBrewer", quietly = TRUE)){
              stop("The package 'RColorBrewer' is required to change the color palettes used for plotting. Please install it and run unPC again.", call. = FALSE)
            }

            availPalettes <- RColorBrewer::brewer.pal.info

            if(colorBrewerPalette %in% rownames(availPalettes)){
              paletteToUse <- colorBrewerPalette
            } else{
              # Tested. fires correctly 031717
              print(paste0("Did not recognize the specified color palette: ", colorBrewerPalette, " in the RColorBrewer package."))
              print("Using palette 'PiYG' by default.")
              paletteToUse <- "PiYG"
            }

          } else{
            # Default color palette to use
            paletteToUse <- "PiYG"
          }

          # Use the brewer.pal.info to get the number of colors avail for each color scheme, and then use the max number.
          availPalettes <- RColorBrewer::brewer.pal.info
          numColors <- availPalettes[which(rownames(availPalettes) == paletteToUse), 1]

          minUnPCDist <- min(pairwiseTileDiffs$ratioPCToGeogrDist)
          maxUnPCDist <- max(pairwiseTileDiffs$ratioPCToGeogrDist)
          meanUnPCDist <- mean(pairwiseTileDiffs$ratioPCToGeogrDist)

          maxToMeanUnPCDist <- maxUnPCDist - meanUnPCDist

          meanToMinUnPCDist <- meanUnPCDist - minUnPCDist

          # Set up the number of bins to be the same as the maximum number colors for the specified palette,
          # with the mean value centered in the bin with middle color (generally close to white)
          # This breaks the side of the unPC value distribution that has the greatest range from the mean into #colors/2 bins.
          # Note - for an odd number of colors, this will make the bin with the mean value the middle bin, but for an even number of
          # bins, the bin containing the mean will be offset by 1 bin, making it a color other than white.
          # The same bin edges are enforced on the other side of the distribution, but this means that if the distribution
          # is highly skewed, the colors on only one side of the distribution get into the darker shades, and the colors on the opposite end (other plot) will
          # be lighter shades.

          # If the distr. is right skewed, make the colorBreaksVect by subtracting from the max value, else make it by adding to the min value.
          # NOTE - the colorBreaksVect can be neg if the 'if' condition is met.
          if(maxToMeanUnPCDist > meanToMinUnPCDist){
            binWidth <- maxToMeanUnPCDist / (numColors / 2)
            colorBreaksVect <- c(maxUnPCDist, maxUnPCDist - seq(1,numColors,1) * binWidth)
            # The colorBreaksVect needs to be non-decreasing in the findInterval call below, so fix that here.
            colorBreaksVect <- rev(colorBreaksVect)

          } else{
            binWidth <- meanToMinUnPCDist / (numColors / 2)
            colorBreaksVect <- c(minUnPCDist, minUnPCDist + seq(1,numColors,1) * binWidth)
          }

          colorBrewerColors <- RColorBrewer::brewer.pal(n = numColors, name = paletteToUse)

          # flipping the colors makes them consistent with the manually defined colors above.
          colorBrewerColors <- rev(colorBrewerColors)

        }

        # This is the dr parameter in the ellipse call to determine how many points
        # to use to define the border of the ellipse (this is in units of radians per point)
        ellipseDrawResolution <- 0.25

        spanCoordSize_x <- max(pairwiseCoordsForPlotting$x_1) - min(pairwiseCoordsForPlotting$x_1)
        spanCoordSize_y <- max(pairwiseCoordsForPlotting$y_1) - min(pairwiseCoordsForPlotting$y_1)

        # Set the largest span (in the x or y direction) as the one to use.
        spanCoordSize <- ifelse(spanCoordSize_x > spanCoordSize_y, spanCoordSize_x, spanCoordSize_y)

        # This is the width of the ellipses to use based on the span of the coordinates.
        # THIS TESTING IS NOT WORKING WELL RIGHT NOW - DOES NOT PICK OUT BAD ENTRIES. ENTRIES LIKE 2 OR '2' ARE FINE. 031717
        if(!is.null(ellipseWidth)){
          #print("In not null.")
          # First need to make sure the entered ellipseWidth can be coerced to numeric. Using this trycatch
          # give more friendly and more descriptive warning/error messages than default R
          # The testResult will be NA if not able to coerce to numeric, in which case set it manually to a
          # default value.
          testResult <- tryCatch({
            as.numeric(ellipseWidth)
          }, warning = function(w){
            print("WARNING: Not able to convert the entered ellipse width into a number. Setting to default value instead.")
            return(w)
          }, error = function(e){
            print("ERROR: Could not parse the entered ellipse width. Setting to default value instead.")
            return(e)
          })

          if (!is.na(testResult)) {
            calcWidth <- as.numeric(ellipseWidth)
          } else{
            # if fails the numeric test above, then set to default value
            calcWidth <- spanCoordSize / 40
          }

        } else{
          # If the ellipse width was not specified, set it to the default value
          calcWidth <- spanCoordSize / 40
        }

        # Get the number of ellipses that will be drawn (1 per row in the pairwiseTileDiffs)
        numEllipses <- length(pairwiseCoordsForPlotting[,1])

        # Use the drawing resolution to calculate the number of points that will be
        # generated for each ellipse
        numPointsPerEllipse <- ceiling((2*pi)/ellipseDrawResolution) + 1

        # pre-allocate a matrix to hold all of the ellipse coordinates
        ellipseHolder <- matrix(nrow=numPointsPerEllipse*numEllipses,ncol=2)

        # pre-allocate a vector to hold the associated group number for each ellipse (important for
        # plotting them)
        ellipseGroupHolder <- vector(mode="numeric", length=numEllipses)

        # Pre-allocate vector for the color groups of the different ellipses
        colorGroupHolder <- vector(mode="numeric", length=numEllipses)

        lat1 <- as.numeric(pairwiseCoordsForPlotting$y_1)
        lon1 <- as.numeric(pairwiseCoordsForPlotting$x_1)
        lat2 <- as.numeric(pairwiseCoordsForPlotting$y_2)
        lon2 <- as.numeric(pairwiseCoordsForPlotting$x_2)

        latDiff <- lat2 - lat1
        lonDiff <- lon2 - lon1

        # For testing only; use the geogrDist below for the unPC oval construction
        lat_lonDist <- sqrt((latDiff)^2 + (lonDiff)^2)

        # Because these are staggered hexagonal tiles, the rows(y) are good to directly calculate the distance
        # between tiles, but the columns (x) are not, because adjacent tiles in the same row have indices that
        # change in steps of 2 (to accomodate the indices of the staggered tiles that link between them on
        # the row above and the row below them), so need to divide this column number difference by 2 to get an
        # accurate distance (See calc module above for more details)
        #geogrDist <- pairwiseTileDiffs$geogrDist


        ellipseAngle <- atan(latDiff/lonDiff)*(180/pi)

        minLat <- pmin(lat1, lat2)
        minLon <- pmin(lon1, lon2)

        ellipseCenterLat <- minLat + (0.5*abs(latDiff))
        ellipseCenterLon <- minLon + (0.5*abs(lonDiff))

        ellipseDrawResolutionVect <- rep(ellipseDrawResolution,length(pairwiseTileDiffs[,1]))

        ellipseValueMatrix <- cbind(lat_lonDist,ellipseAngle,ellipseCenterLon,ellipseCenterLat)


        # Define the function for calculating each ellipse (this will be called by the apply() call below in order
        # to calculate as many ellipses as needed). The resolution for each ellipse is fixed and passed in to this function
        # as a second (constant) variable in the apply() call
        ellipseCalc <- function(input,resolution, ellipseWidth){

            geogrDist <- input[1]
            ellipseAngle <- input[2]
            ellipseDrawResolution <- resolution
            ellipseCenterLon <- input[3]
            ellipseCenterLat <- input[4]


            #print(geogrDist)
            #print(ellipseAngle)
            #print(ellipseCenterLon)
            #print(ellipseCenterLat)

            calcEllipse <- shape::getellipse(rx = geogrDist/2, ry = ellipseWidth, angle = ellipseAngle, dr = ellipseDrawResolution, mid = c(ellipseCenterLon,ellipseCenterLat))
            #print(calcEllipse)
            return(calcEllipse)
        }

        # For testing only
        # ellipseTester <- ellipseValueMatrix[1:10,]

        # Use apply here instead of a for loop. This goes through the input data from the matrix 1 row at a time, calling ellipseCalc function and also
        # passing the fixed parameter value for the drawing resolution (in the var resolution). The apply() call passes the row of values to the called function
        # which unpacks the values and returns the calculated ellipse. The returned ellipse coordinates are originally a 2 col matrix
        # (X in col1 and Y in col2), but when returned to apply, then get converted from a matrix to a vector (with X as the first half of the rows, and Y as the
        # other half.) Apply then also cbinds the output from each ellipseCalc call to form a matrix of results (with ncol = number of calcEllipse calls, and therefore the number of ellipses)
        generatedEllipses <- apply(ellipseValueMatrix,1,ellipseCalc,resolution=ellipseDrawResolution, ellipseWidth = calcWidth)

        # This starts conversion by first splitting the generatedEllipses output matrix (rowwise) in half; the first half (the first 50% of rows on all cols) contains the x coords, the second half
        # contains the y coords
        generatedEllipsesXCoordsMat <- generatedEllipses[1:(nrow(generatedEllipses)/2),]
        generatedEllipsesYCoordsMat <- generatedEllipses[((nrow(generatedEllipses)/2) + 1) : nrow(generatedEllipses),]

        # Need to get numEllipse/numPointsPerEllipse information in order to create the ellipseGroup vector that identifies all the entries that correspond to each ellipse.
        # Calc the number of ellipses created (which is the number of cols in generatedEllipses)
        numEllipses <- ncol(generatedEllipses)

        # The number of points defining each ellipse is the number of rows (entries) in the XCoordsMat (or the YCoordsMat)
        numPointsPerEllipse <- nrow(generatedEllipsesXCoordsMat)

        # This converts the matrices of X and Y values to vectors (column-wise transformation indices going from the matrix to the vector)
        generatedEllipsesXCoordsVect <- as.vector(generatedEllipsesXCoordsMat)
        generatedEllipsesYCoordsVect <- as.vector(generatedEllipsesYCoordsMat)

        # This makes a grouping vect for the ellipses with length = numEllipses * numPointsPerEllipse
        generatedEllipsesGroupVect <- rep(seq(1,numEllipses,1), each = numPointsPerEllipse)

        # This is a great function. findIntervals essentially does a histogram bin of the PCA to Geographic distance ratio data given the binning sequence, then forces the top and bottom points to be in the
        # first and last bins. It returns a vect of the bin # (starting at 1) for the data value for each point. This gives a vect the length of the number of ellipses.
        generatedEllipsesColorIntervalVect <- findInterval(PCAGeogrRatioDist, colorBreaksVect, all.inside = TRUE)

        # Make a color group vector for the ellipses based on their value of PC:Geographic distance. Like the grouping vect above, this has length = numEllipses * numPointsPerEllipse
        generatedEllipsesColorGroupVect <- rep(generatedEllipsesColorIntervalVect, each = numPointsPerEllipse)

        # Create a color vector using the grouping vector as an index
        generatedEllipsesColorsVect <- colorBrewerColors[generatedEllipsesColorGroupVect]

        #Merge the ellipse data (x,y coords.) and group values into a single matrix
        calcEllipseDf <- data.frame("Longitude" = generatedEllipsesXCoordsVect, "Latitude" = generatedEllipsesYCoordsVect,
                                    "Group" = generatedEllipsesGroupVect, "ColorGroup" = generatedEllipsesColorGroupVect, "Color" = generatedEllipsesColorsVect, stringsAsFactors = FALSE)

        # This sorts the ellipses by color group (starting with 1, the color group
        # representing the smallest PCA distances.

        # Change the way the sorting works here to change the ordering of the ellipses
        # in the plot - use: 'order(ColorGroup)' to sort ellipses in ascending order of their
        # color group so
        # that the ellipses with the largest PCADist/GeogrDist values are plotted last (ie on top
        # of all others). Use: 'order(-ColorGroup)' to sort ellipses in descending of their
        # color group so that ellipses with the smallest PCADist/GeogrDist values are plotted last
        # (ie on top of the others)

        # Sort so the least diff are at the top
        sortedCalcEllipseDf_leastDiffOnTop <- calcEllipseDf[with(calcEllipseDf, order(-ColorGroup)), ]

        # Sort so the most diff are at the top
        sortedCalcEllipseDf_mostDiffOnTop <- calcEllipseDf[with(calcEllipseDf, order(ColorGroup)), ]

        # ggplot plots the ellipses in the order of the grouping variable, with '1' being
        # plotted first, and so on. Since the original
        # ellipse grouping was agnostic to the color groups, ellipses are plotted
        # without regard to their color group. Here set up another grouping vector
        # that codes the ellipses correctly after sorting (ellipses with the smallest
        # PCA distances are coded with the lower group numbers and will therefore
        # be plotted under the ellipses representing higher PCA distances.)
        sortedGroupHolder <- rep(seq(1,numEllipses,1), each = numPointsPerEllipse)

        sortedCalcEllipseDf_mostDiffOnTop$GroupsAfterColorSorting <- sortedGroupHolder
        sortedCalcEllipseDf_leastDiffOnTop$GroupsAfterColorSorting <- sortedGroupHolder


        # All of this plots in 1 call to ggplot - the ellipses first, then the country outlines overlaid on them.
        # The ellipses plot automatically (each one separately) based on the groups from the ellipseGroupHolder, and the colors in the colorGroupHolder

        # This creates the output visualizations for both the most and least differentiated ellipses stacked on top.
        create_unPCMap <- function(unPCDataFrameForPlot, populationCoordHolder, plotTitle, outputFileName){

          if(isTRUE(plotWithMap)) {

            mappingData <- ggplot2::map_data("world")

            # For plotting with US states
            #mappingData <- ggplot2::map_data("state")

            xMin <- min(unPCDataFrameForPlot$Lon)
            xMax <- max(unPCDataFrameForPlot$Lon)

            yMin <- min(unPCDataFrameForPlot$Lat)
            yMax <- max(unPCDataFrameForPlot$Lat)

            xlims <- c(xMin - 2, xMax + 2)
            ylims <- c(yMin - 2, yMax + 2)

            # Check that PBSmapping (which is a suggested package and only used for map clipping) is installed. Verified this triggers 031717
            if (!requireNamespace("PBSmapping", quietly = TRUE)){
              stop("The package 'PBSmapping' is required for creating the maps. Please install it and run unPC again.", call. = FALSE)
            }

            # Need to change the names to match those expected by PBSmapping::clipPolys
            names(mappingData) <- c("X", "Y", "PID", "POS", "region", "subregion")

            clippedMap <- PBSmapping::clipPolys(mappingData, xlim = xlims, ylim = ylims)

            print(class(clippedMap))

            # To get the point scaling right, need to put the size outside the aes; Now using ln(num Indivs sampled per population) for the point scaling, and it is working well.
            unPCMap <- ggplot2::ggplot(data = unPCDataFrameForPlot, ggplot2::aes(x=Longitude, y=Latitude)) +
              ggplot2::geom_polygon(group=unPCDataFrameForPlot$GroupsAfterColorSorting, fill=unPCDataFrameForPlot$Color) +
              ggplot2::geom_polygon(data = clippedMap, ggplot2::aes(x = X, y = Y, group = PID), fill = NA, linetype = 1, size = 0.2, color = "black") +
              ggplot2::geom_point(data=populationCoordHolder, ggplot2::aes(x=Lon, y=Lat, color = factor(popColoring)), size = populationCoordHolder$numIndivsSampled / populationPointNormalization) +
              ggplot2::scale_color_manual(values = c("black", "gray")) +
              ggplot2::theme_bw() + ggplot2::theme(legend.position="none", axis.text.x = ggplot2::element_blank(), axis.title.x = ggplot2::element_blank(),
                                axis.text.y = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank(),
                                panel.grid.major = ggplot2::element_blank(),
                                axis.ticks.x = ggplot2::element_blank(),
                                axis.ticks.y = ggplot2::element_blank()) +
              ggplot2::coord_map(xlim = xlims, ylim = ylims)

            if(isTRUE(savePlotsToPdf)){
              pdf(outputFileName, width = 8, height = 8)
              print(unPCMap)
              dev.off()
            } else{
              print(unPCMap)
            }

          } else {

            xMin <- min(unPCDataFrameForPlot$Lon)
            xMax <- max(unPCDataFrameForPlot$Lon)

            yMin <- min(unPCDataFrameForPlot$Lat)
            yMax <- max(unPCDataFrameForPlot$Lat)

            unPCMap <- ggplot2::ggplot(data = unPCDataFrameForPlot, ggplot2::aes(x=Longitude, y=Latitude)) +
              ggplot2::geom_polygon(group=unPCDataFrameForPlot$GroupsAfterColorSorting, fill=unPCDataFrameForPlot$Color) +
              ggplot2::geom_point(data=populationCoordHolder, ggplot2::aes(x=Lon, y=Lat, color = factor(popColoring)), size = populationCoordHolder$numIndivsSampled / populationPointNormalization) +
              ggplot2::scale_color_manual(values = c("black", "gray")) +
              ggplot2::theme(legend.position="none", panel.grid = ggplot2::element_blank(), panel.background = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank(),
                                 axis.text.x = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank(), axis.title.x = ggplot2::element_blank(),
                                 axis.ticks.y = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank()) +
              ggplot2::coord_cartesian(xlim = c(xMin, xMax), ylim = c(yMin, yMax))

            #print(populationCoordHolder$numIndivsSampled / populationPointNormalization)

            if(isTRUE(savePlotsToPdf)){
              pdf(outputFileName, width = 8, height = 8)
              print(unPCMap)
              dev.off()
            } else{
              print(unPCMap)
            }

          }

        }


        # Call both plotters (most and least differentiated ellipses stacked on top) sequentially here using the
        # function above.

        # If the currentIterationNumber is set to -1, this means the function was called using population-level mean aggregated data for the plotting, so name these specially.
        # The the number is anything else, this is a plot of the results from a specific PCA analysis of a simulation iteration, so add the name of the iteration to the graph file
        # and title.
        if(currentIterationNumber == -1){

          plotTitle <- paste("Population-level mean aggregated results_mostDiffOnTop",sep="")
          pdfFileName <- paste(plotTitleRoot,"_popnLevelMeanAggregated_mostDiffOnTop.pdf",sep="")

          create_unPCMap(sortedCalcEllipseDf_mostDiffOnTop, populationCoordHolder, plotTitle, pdfFileName)

          plotTitle <- paste("Population-level mean aggregated results_leastDiffOnTop",sep="")
          pdfFileName <- paste(plotTitleRoot,"_popnLevelMeanAggregated_leastDiffOnTop.pdf",sep="")

          create_unPCMap(sortedCalcEllipseDf_leastDiffOnTop, populationCoordHolder, plotTitle, pdfFileName)

        } else{

          plotTitle <- paste("Iteration_#", as.character(currentIterationNumber), "_mostDiffOnTop.pdf", sep="")
          pdfFileName <- paste(plotTitleRoot,"_Iteration_#", as.character(currentIterationNumber), "_mostDiffOnTop.pdf", sep = "")

          create_unPCMap(sortedCalcEllipseDf_mostDiffOnTop, populationCoordHolder, plotTitle, pdfFileName)

          plotTitle <- paste("Iteration_#", as.character(currentIterationNumber), "_leastDiffOnTop.pdf", sep="")
          pdfFileName <- paste(plotTitleRoot,"_Iteration_#", as.character(currentIterationNumber), "_leastDiffOnTop.pdf", sep = "")

          create_unPCMap(sortedCalcEllipseDf_leastDiffOnTop, populationCoordHolder, plotTitle, pdfFileName)
        }


    }


    # If the input data is a matrix, then call generateUnPCPlot directly. This passes -1 as the iteration number to trigger different file name
    # and plot title construction for this aggregated graph
    if(length(dimnames(pairwiseTileDiffHolder)) == 2) {
        generateUnPCPlot(-1, pairwiseTileDiffHolder)
    }

    # If the input data is a 3-D array, then need to call generateUnPCPlot on each
    # layer of it (representing the output from PCA of a single simulation iteration) in turn
    # to generate the plots for each. Using lapply like a for loop to be able to also pass in the iteration (layer) number
    # of the current layer being graphed. This is important for making the file name and the plot title.
    # This passes the index value (made with the seq() call) to generateUnPCPlot by default, and also passes another
    # argument, which is the current layer of the 3-D array to graph, which is converted to a data.frame so it
    # can be plotted correctly.
    # NOTE - If iterations were dropped from the pairwise diffs calculation in the above module due to missing data, they are dropped here too, but the iteration numbers that
    # were kept (non-consecutive) will not match the iteration numbers printed in the plots (consecutive).
    if(length(dimnames(pairwiseTileDiffHolder)) == 3) {
        lapply(seq(1:length(pairwiseTileDiffHolder[1,1,])), generateUnPCPlot, pairwiseTileHolder = pairwiseTileDiffHolder)

    }
}

# End pairwise distance mapper
#-------------------------------------------



end.time <- Sys.time()
time.taken <- difftime(end.time, start.time, units = "secs")
print(paste("Took: ",time.taken, "seconds"))

}

