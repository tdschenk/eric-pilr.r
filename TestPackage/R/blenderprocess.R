## Blender data import and processing script
# Eric Ruud, MEI Research Ltd. eruud@meinergy.com

## Load libraries
library(jsonlite)

## Define Functions
derivative <- function(x, data_interval, derivative_window = 8,
                       average_points = 1) {
  ## Calculate the number of data points needed to cover the desired
  ## interval time
  seconds_in_minute <- 60
  resolution <- data_interval / seconds_in_minute
  data_points <- derivative_window / resolution

  ## Create vector to sum the future points
  f <- rep(0, data_points + 1)
  ## what if average_points is large?
  f[1:average_points] <- 1
  f[(length(f) - average_points + 1) : length(f)] <- -1

  ## Filter the selected Header row using the filter vector
  ## determined above

  dVector <- stats::filter(x, f) / derivative_window
  dVector[is.na(dVector)] <- 0
  dVector
}

## Define variables
# Filetype JSON
filetype_json <- fromJSON("C:\\Users\\Eric\\Dropbox\\kfactor\\blender-data-file-handler.json")

# Get the path to the blender file
# This single select code may be required for non-windows platforms
# blender_file <- file.choose(new = TRUE)

# select multiple files (may be windows only)
# blender_files <- choose.files("*.csv")
# blender_file <- blender_files[1]
blender_file <- "C:\\Users\\Eric\\Google Drive\\Eric_PILR\\blender\\Data\\UNC_TELEDYNE_200LPM-1-Infusion-160823.csv"
blender_files <- c(blender_file)

# Create a filename for export
export_file <- paste(strsplit(basename(blender_file),".csv"),"_processed.tsv",collapse = "", sep = "")

# Get the column names to the blender file
colnames <- read.csv(blender_file, check.names=FALSE, nrows = 1, skip = 1, header = FALSE, stringsAsFactors = FALSE)

## Build the colnames vector to import only columns we want with readcsv
# We'll loop over each element in the vector of column names to see if we want it
# If we don't, we'll change it to NULL so that we don't import it
importclass <- vector(, length(colnames))
importname <- vector(, length(colnames))
for (i in 1:length(colnames))
{
  for (colname in filetype_json$parser$data_mappings$from_name)
  {
    # Check the JSON for the column names
    if (colname == colnames[i])
    {
      importclass[i] <- NA
      importname[i] <- colnames[i]
      break
    }
    else
    {
      importclass[i] <- "NULL"
      importname[i] <- "NULL"
    }
  }
}

# import valid columns
import <- read.csv(blender_file, check.names=FALSE, skip = 2, stringsAsFactors = FALSE, col.names = importname, colClasses = importclass)

# if multiple files are selected, append to first
# assume same format for now
if (length(blender_files)>1)
{
  for (i in 2:length(blender_files))
  {
    import2 <- read.csv(blender_files[i], check.names=FALSE, skip = 3, stringsAsFactors = FALSE, colClasses = importclass)
    import[1:nrow(import2)+nrow(import),] <- import2
  }
}


# Change the date/time to POSIX format
# Assume format and time column name (CalRQ files)
import$Time <- as.POSIXct(import$Time, format = "%m/%d/%Y %H:%M:%S")

# Calc a data interval for the derivative
data_interval <- median(diff(as.numeric(import$Time)))

# Remove data with BIOSFlow as a NaN since we can't use it
import <- import[complete.cases(import$BIOSFlow),]

# We'll try to figure out where the data intervals are now
# Take the derivative, assume values above the average of the
# derivative mean the value has switched.
biosderiv <- abs(derivative(import$BIOSFlow, data_interval, derivative_window = 2))

counter <- 0
breaks <- 0
for (i in 1:length(biosderiv))
{
  # Check if deriv is above mean
  if (biosderiv[i] > mean(biosderiv))
  {
    # if it is, add to the counter vector
    counter <- counter + 1
  }
  else
  {
    # if it's not above the average, check to see if it was previously, for at least 2 points.
    if (counter > 1)
    {
      # if so, note this breakpoint and zero the counter vector to find the next break
      breaks <- append(breaks,i-2)
      counter <- 0
    }
  }
}

# In case we have a negative value in breaks
breaks <- breaks[breaks>=0]

# Loop over the data with the breaks we found.
# Initalize variables
processed = as.data.frame(matrix(ncol=16, nrow=length(breaks) - 1))
names(processed) = c("BIOSF1","BIOSF2","BIOSF3","BIOSF4", "MFC1", "MFC2", "MFC3", "MFC4", "MFC1Error","MFC2Error","MFC3Error","MFC4Error","Index1","Index2","Index3","Index4")

exportdata = as.data.frame(matrix(ncol=6))
names(exportdata) = c("BIOSF","MFC1", "MFC2", "MFC3", "MFC4", "Index")

for (i in 1:(length(breaks) - 1))
{
  # One subset of data
  # 4 MFCs may be installed. Find out which was used for this set
  for (j in 1:4)
  {
    # loop over all 4 mfcs, if they != null
    if (!(is.null(eval(parse(text = paste("import$MFCFlow_",j,sep = ""))))))
    {
      # Vector exists
      # Check to see if there was flow for this break
      # We will see if the values in this section are less than 2% of full scale

      limit <- max(eval(parse(text = paste("import$MFCFlow_",j,sep = ""))))*.02
      bioslimit <- max(import$BIOSFlow)*.02

      section <- eval(parse(text = paste("import$MFCFlow_",j,"[",breaks[i],":",breaks[i+1],"]",sep = "")))
      biossection <- eval(parse(text = paste("import$BIOSFlow[",breaks[i],":",breaks[i+1],"]",sep = "")))

      # Check if this break is above that limit for this MFC
      # We don't want to analyze the MFC if it's not used for this part
      # Ignore if mode is under 2% of fullscale
      if (median(section) > limit)
      {
        # this MFC was used for this section
        # (not used)
        # MFCName <- paste("import$MFCFlow_",j,sep = "",collapse = "")

        # Filter to values within 2% of fullscale of median
        upperlimit <- median(section) + limit
        lowerlimit <- upperlimit - 2*limit
        biosupperlimit <- median(biossection) + bioslimit
        bioslowerlimit <- biosupperlimit - 2*bioslimit

        selection <- (section < upperlimit & section > lowerlimit & biossection < biosupperlimit & biossection > bioslowerlimit)

        # attempt to detect biosmodifier
        # depending on the instrument biosf might be 1000x mfc read value
        # compare error and choose what has less error
        if (abs(mean((section-biossection)/biossection)) < abs(mean((section*1000-biossection)/biossection)))
        {
          biosmodifier <- 1
        } else {
          biosmodifier <- 1000
        }

        # correct bios read and assign selected vars
        biossection <- biossection[selection]/biosmodifier
        section <- section[selection]

        ## calculate some values
        # mean of MFC data
        processed[i,4+j] <- mean(section)

        # error vs. BIOS
        processed[i,8+j] <- abs(mean((section-biossection)/biossection))

        # mean of BIOSF data that passed our filter
        processed[i,j] <- mean(biossection)

        offset <- breaks[i]
        # Build out an array of the data used in calculations
        if (length(section)>0)
        {
        exportdata[1:length(section)+offset,1] <- biossection
        exportdata[1:length(section)+offset,j+1] <- section
        exportdata[1:length(section)+offset,6] <- which(selection == TRUE) + offset
        }
      }
    }
  }
}

# Build a CSV for export
# Check to see which MFCs have data then export them

# save avg temp and pressure data once
tempavg <- mean(import$BIOSTemp, na.rm = TRUE)
pressavg <- mean(import$BIOSBP, na.rm = TRUE)
cat(c("Cal Temp (C)\tCal Pressure (mmHg)\tCal Date\tFilename\n",tempavg,"\t",pressavg,"\t",strftime(import$Time[1],"%D"),"\t",basename(blender_file),"\n"), file=export_file, append=FALSE)

# in case other files were used
if (length(blender_files>1))
{
  for (i in 2:length(blender_files))
  cat(c("\t\t\t",basename(blender_files[i]),"\n"), file=export_file, append=TRUE)
} else
{
  cat(c("\n"), file=export_file, append=TRUE)
}

# Write our data to a CSV/TSV file
for (i in 1:4)
{
  # Is there at least one element that is not NA?
  navalues <- !(is.na(processed[4+i]))
  # If so, process this MFC
  if (any(navalues))
  {
    # find which values are not NA
    bios <- subset(processed[i],navalues)
    mfcflow <- subset(processed[4+i],navalues)
    error <- subset(processed[8+i],navalues)

    # set rownames to null so ordering references row nums
    rownames(bios) <- NULL
    rownames(mfcflow) <- NULL
    rownames(error) <- NULL

    # sort based on MFCFlow
    ordering <- unlist(sort(mfcflow[,1], index.return = TRUE)[2])
    bios <- bios[ordering,1]
    mfcflow <- mfcflow[ordering,1]
    error <- error[ordering,1]

    # build a character vector for CalRQ
    # 6 decmal places, format is MFC1,BIOS1,MFC2,BIOS2, etc
    for (j in 1:length(mfcflow))
    {
      if (j==1)
      {
        calstring <- c(sprintf(fmt = "%.6f",mfcflow[j]),sprintf(fmt = "%.6f",bios[j]))
      } else
      {
        calstring <- append(calstring,c(sprintf(fmt = "%.6f",mfcflow[j]),sprintf(fmt = "%.6f",bios[j])))
      }

    }

    # Write out the averages we just calc'd with a header
    write.table(list(bios,mfcflow,error),file = export_file, row.names = FALSE, append = TRUE, sep = "\t", col.names = c("BIOSF",paste("MFCFlow",i,sep=""),"Error"))

    # Write out avg error and the cal string we calc'd
    cat(c("\tAvg Error (%):\t",mean(error),"\tCal String:\t",paste(calstring,collapse=","),"\t\n\n"), file=export_file, append=TRUE)

    ## For CalRQ calibration, we use a 10-point curve
    # Try to make a 10-point curve if we have more than 10 points
    if (length(mfcflow) > 10)
    {
      targets <- seq(from=min(mfcflow), to=max(mfcflow), length.out = 10)

      # initalize vectors
      calindex <- vector(,10)
      for (k in 1:length(targets))
      {
        calindex[k] <- which.min(abs(mfcflow-targets[k]))
      }
      # Store the data
      calbios <- bios[calindex]
      calmfc <- mfcflow[calindex]
      calerror <- error[calindex]

      # Note that this is a subset of data
      cat(c("10-Point Calibration Subset\n"), file=export_file, append=TRUE)

      # Write out the averages we just calc'd with a header
      write.table(list(calbios,calmfc,calerror),file = export_file, row.names = FALSE, append = TRUE, sep = "\t", col.names = c("BIOSF",paste("MFCFlow",i,sep=""),"Error"))

      # build a character vector for CalRQ
      # 6 decmal places, format is MFC1,BIOS1,MFC2,BIOS2, etc
      for (j in 1:10)
      {
        if (j==1)
        {
          calstring <- c(sprintf(fmt = "%.6f",calmfc[j]),sprintf(fmt = "%.6f",calbios[j]))
        } else
        {
          calstring <- append(calstring,c(sprintf(fmt = "%.6f",calmfc[j]),sprintf(fmt = "%.6f",calbios[j])))
        }

      }

      # Write out avg error and the cal string we calc'd
      cat(c("\tAvg Error (%):\t",mean(calerror),"\tCal String:\t",paste(calstring,collapse=","),"\t\n\n"), file=export_file, append=TRUE)
    }

  }
}

# Append raw data to TSV
cat(c("Data Used\n"), file=export_file, append=TRUE)
write.table(exportdata[!is.na(exportdata$BIOSF),],file = export_file, row.names = FALSE, append = TRUE, sep = "\t", col.names = c("BIOSF","MFC1", "MFC2", "MFC3", "MFC4", "Index"))
