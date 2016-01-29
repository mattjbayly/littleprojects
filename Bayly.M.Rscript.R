#---------------------------------------------------------
# Set key user parameters for analysis 
#---------------------------------------------------------
	# Use existing dataset from hard drive or download from Water Survey Canada website?
	Set_Source <- "Hard drive" # or "Re-download"
	StationID <- "08HB006"
	StartDate <- "1992-01-01"
	EndDate <- "2011-12-31"
	Time_Unit <- "days"
	directory <- "C:/Users/DW/Desktop/Matts_misc/Spring 2016/4.Submitted Applications/EcoFish_data_analyst/Evaluation project"
	filename <- "Daily__Dec-06-2013_10_08_10PM__ddf.csv"
	myPARAM <- 1 # Type of data: 1=discharge, 2=water level
	
#---------------------------------------------------------
# Load in data and check structure 
#---------------------------------------------------------
# Set file directory:
	setwd(directory); 
	#dir() # view files in directory 
# Read in csv file:
	dat <- read.csv(filename)
# Check data structure:
	#head(dat, 3) # check data structure
	#tail(dat)
	# Make corrections: exclude "disclaimer" and file notes (from 'id' column)
		dat <- dat[which(dat$ID==StationID),]
		dat$ID <- droplevels(dat$ID)
		summary(dat, 3) # check data structure: str(dat)
		# - looks good
# Reclassify date & PARAM column  		
	dat$Date <- as.Date(dat$Date) # Reformat to date
	dat$PARAM <- as.numeric(dat$PARAM) # Reformat to date
	range(dat$Date)
	
#---------------------------------------------------------
# Subset data to specified date range:
	# After and equal to start date and before or equal to end date 
		dat2 = subset(dat, Date >= as.Date(StartDate) & Date <= as.Date(EndDate) & PARAM == myPARAM)

#---------------------------------------------------------
# Find and characterize missing data: 
#---------------------------------------------------------
	DescribeMissing <- function(mydat, MissingRows, MissingDates, printMissAll){
		if (MissingRows == TRUE) {
			# Find any difference time units in series & count them 
				expected_length <- seq(as.Date(StartDate), as.Date(EndDate), Time_Unit)
				print(paste0(length(expected_length) - length(mydat[,c("Date")]),
					" ", Time_Unit, " row # difference from expected sequence in time series"))
			# Which time units (if any) are different from expected time sequence?	
				these_match <- match(as.character(expected_length), as.character(mydat[,c("Date")]))
				if (printMissAll == TRUE) {
					print(paste0("These missing row are: ", mydat[,c("Date")][-c(these_match)]))
			} }
		if (MissingDates == TRUE) {
			# Find missing date-data values in sequence: 
				mymissing <- mydat[!complete.cases(mydat[,c("Value")]),]
				print(paste0(length(mymissing$Date), " days have missing NA values"))
				result <- rle(diff(julian(mymissing$Date)))
				print(paste0("Warning there were up to ", max(result$lengths),
					" sequential dates with missing values!"))
				print(rev(sort(result$lengths)))
					
				if (printMissAll == TRUE) {
					print(mymissing)
				}
			}
		}

# Run Missing data report function on dataset:		
	DescribeMissing(mydat=dat2, MissingRows=TRUE, MissingDates=TRUE, printMissAll=FALSE)
	
#---------------------------------------------------------------------
# Interpolate over missing data using linear interpolation and splines: 
#---------------------------------------------------------------------
	# Simple linear interpolation over missing values
		dat2$LinInt <- approx(dat2$Date, dat2$Value, 
			n=length(seq(as.Date(StartDate), as.Date(EndDate), Time_Unit)))[[2]]		
	# Interpolation with spline, fmm
		dat2$SplineInt <- spline(dat2$Date, dat2$Value, 
			n=length(seq(as.Date(StartDate), as.Date(EndDate), Time_Unit)))[[2]]		
	# Interpolation with spline, PERIODIC
		dat2$SplineIntPeri <- spline(dat2$Date, dat2$Value, 
			n=length(seq(as.Date(StartDate), as.Date(EndDate), Time_Unit)), method="periodic")[[2]]
	
#---------------------------------------------------------------------
# Generate Summary Table: 
#---------------------------------------------------------------------
	# Make re-usable function to generate customized discharge summary statistics:
DischargeStatistics <- function(mydat, DischareColumn){	
	
	# Discharge Summary Statistics:
		df <- data.frame(Statistic=c("Mean Annual Discharge (MAD)", 
					"Median", "Min daily", "Max daily", "10th percentile", "20th percentile"), Discharge=NA)
		df$Discharge[1] <- summary(mydat[,c(DischareColumn)])["Mean"]	
		df$Discharge[2] <- summary(mydat[,c(DischareColumn)])["Median"]
		df$Discharge[3] <- summary(mydat[,c(DischareColumn)])["Min."]
		df$Discharge[4] <- summary(mydat[,c(DischareColumn)])["Max."]
		df$Discharge[5] <- quantile(mydat[,c(DischareColumn)], c(.10), na.rm=TRUE) 
		df$Discharge[6] <- quantile(mydat[,c(DischareColumn)], c(.20), na.rm=TRUE) 
		df$Discharge <- round(df$Discharge, 1)
	# Monthly average discharge and percent MAD:
		monthName <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
		df2 <- data.frame(Month=monthName, AverageDischarge=NA, PercentMAD=NA)
			mean_monthly <- tapply(mydat[,c(DischareColumn)],format(mydat[,c("Date")],"%m"),mean, na.rm=TRUE)
			monthly_percent <- as.numeric(mean_monthly)/summary(mydat[,c(DischareColumn)])["Mean"]
		df2$AverageDischarge <- round(mean_monthly, 1)
		df2$PercentMAD <- round((monthly_percent*100), 0)

	# Merge output to list object
		FinalOutput <- list("Total" = df, "Monthly" = df2)
		return(FinalOutput)
}

	# Basic summary table from linear interpolation:
	SummaryTable <- DischargeStatistics(mydat=dat2, DischareColumn="LinInt")

#---------------------------------------------------------------------
# Test sensitivity of interpolation method over missing data
	# e.g. linear vs. spline_fmm vs. period spline
	# Missing Values filled with spline interpolation vs. linear interpolation  
		mapply("/",DischargeStatistics(mydat=dat2, DischareColumn="SplineInt"),DischargeStatistics(mydat=dat2, DischareColumn="LinInt"),SIMPLIFY = FALSE)
	# Missing Values filled with period spline interpolation vs. linear interpolation  
		mapply("/",DischargeStatistics(mydat=dat2, DischareColumn="SplineIntPeri"),DischargeStatistics(mydat=dat2, DischareColumn="LinInt"),SIMPLIFY = FALSE)
	# Missing Values left as NA vs. linear interpolation  
		mapply("/",DischargeStatistics(mydat=dat2, DischareColumn="Value"),DischargeStatistics(mydat=dat2, DischareColumn="LinInt"),SIMPLIFY = FALSE)
	
#---------------------------------------------------------------------
# Create 366X20 matrix of days of year & 20 year period: 
#---------------------------------------------------------------------
		DatSeq <- seq(as.Date(StartDate), as.Date(EndDate), "years")
		YearSeq <- as.numeric(format(DatSeq,"%Y"))
		dat2$YearSeq <- as.numeric(format(dat2$Date,"%Y"))
	# Create empty matrix to fill:
		bigmatrix <- matrix(NA, nrow = 366, ncol = length(YearSeq))
	# fill in yearly values: 
		for(i in 1:length(YearSeq)){
			# loop through years individually
			YearS <- dat2[which(dat2$YearSeq == YearSeq[i]), ] # single year data
			yearlength <- dim(YearS)[1] # adjust fill for number of days in year
			bigmatrix[1:yearlength,i] <- YearS$LinInt # fill matrix column for year
		}
		colnames(bigmatrix) <- YearSeq # better labels
		rownames(bigmatrix)<- 1:366
		head(bigmatrix, 2); dim(bigmatrix)
	
#---------------------------------------------------------------------
# Create a report-quality plot that displays the annual trends in discharge: 
#---------------------------------------------------------------------
		# function to add transparency to plot
			makeTransparent = function(..., alpha=0.5) {
					  if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")
						  alpha = floor(255*alpha)  
						  newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
						  .makeTransparent = function(col, alpha) {
							rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
					  }
					  newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)
					  return(newColor)
					}
	
	# Build figure that shows all years as semi-transparent grey lines first 
		plot(bigmatrix[,1], type="l", ylab=expression("Discharge ("~m^{3} ~ "/ Second)"), xlab="",
				ylim=c(0, max(dat2$LinInt) + 2), xaxt="n", lwd=2, col=makeTransparent("lightgrey", alpha=0.5))

	# Add separate line-year to plot
		for(i in 2:length(YearSeq)){ # first year has already been added to the base plot.
				# loop through years individually
				lines(bigmatrix[,i], lwd=2, col=makeTransparent("lightgrey", alpha=0.5))
			}
	# Add line for median, max & min 
		lines(apply(bigmatrix, 1, max, na.rm=TRUE), lwd=2, col=makeTransparent("red", alpha=0.5))
		lines(apply(bigmatrix, 1, min, na.rm=TRUE), lwd=2, col=makeTransparent("red", alpha=0.5))
		lines(apply(bigmatrix, 1, mean, na.rm=TRUE), lwd=2, col=makeTransparent("blue", alpha=0.5))
		lines(apply(bigmatrix, 1, median, na.rm=TRUE), lwd=2, col=makeTransparent("green", alpha=0.5))

	# Axis label unique for each month: 
		# the first of each month in the 366 day sequence
			my_at <- which(format(strptime(1:366, format="%j"), format="%d")=="01")
			my_labels <- format(strptime(1:366, format="%j"), format="%d-%b")[my_at]
			# add axis to plot
			axis(1, at = my_at, labels=my_labels)	
	
#---------------------------------------------------------------------
# Save/write final objects for report: 
#---------------------------------------------------------------------	
	# Summary Statistics Table:
		write.csv(SummaryTable, row.names=FALSE,
			file=paste0("SummaryTable_", StationID, ".", min(YearSeq), "-", max(YearSeq),".csv"))
	
	#---------------------------------------------------------------------
	# Save Big matrix 366*20:
		write.csv(bigmatrix, row.names=FALSE,
			file=paste0("Daily_", StationID, ".", min(YearSeq), "-", max(YearSeq),".csv"))
	
	#---------------------------------------------------------------------
	# Save final plot:
		# EcoFish Template settings
			pdf(width=8.5-(1.22+0.90), height=5,
				file=paste0("Annual_", StationID, ".", min(YearSeq), "-", max(YearSeq),".pdf"), family="Times")
			
		par(mar=c(5,4.2,2.1,2.1))
			Log_bigmatrix <- log(bigmatrix)
			Log_bigmatrix <- Log_bigmatrix[1:365,]

			# Log scale plot
			myseq <- seq(0, 300, by=25); Logmyseq <- log(myseq); Logmyseq[1] <- 0
			plot(Log_bigmatrix[,1], type="n", ylab=expression("Discharge ("~m^{3} ~ "/ Second)"), xlab="",
				ylim=c(min(Log_bigmatrix, na.rm=TRUE), max(Log_bigmatrix, na.rm=TRUE)), xaxt="n", yaxt="n", lwd=2, col=makeTransparent("lightgrey", alpha=0.5))
			# add truncated y-axis ticks 
				for(i in 1:length(Logmyseq)){
					abline(h=Logmyseq[i], col=makeTransparent("darkgrey", alpha=0.7))
				}
			# Add line for median, max & min 
				lines(apply(Log_bigmatrix, 1, max, na.rm=TRUE), lwd=1.2, col=makeTransparent("green", alpha=0.8))
				lines(apply(Log_bigmatrix, 1, min, na.rm=TRUE), lwd=1.2, col=makeTransparent("blue", alpha=0.8))
				#lines(apply(Log_bigmatrix, 1, mean, na.rm=TRUE), lwd=2, col=makeTransparent("blue", alpha=0.5))
				lines(apply(Log_bigmatrix, 1, median, na.rm=TRUE), lwd=1.2, col=makeTransparent("red", alpha=0.8))
			# add axis to plot
				axis(1, at = my_at, labels=my_labels, cex.axis=0.8, las=2)
				axis(2, at=Logmyseq, labels=as.character(myseq), cex.axis=0.8)

			#legend
				# setup for no margins on the legend
				legend('top','groups',c("Max","Min","Median"), lwd=c(1.2, 1.2, 1.2),
					col=c(makeTransparent("green", alpha=0.5),
					makeTransparent("blue", alpha=0.5),makeTransparent("red", alpha=0.5)),
					ncol=3,bty ="n")
		
			# close file connection for plot
			dev.off()
			
#---------------------------------------------------------------------
# END