### Simple peak detection
### Benjamin Bell 
### https://github.com/benbell95/peak-detection
### https://www.benjaminbell.co.uk
### Work in progress - use at own risk
### Noisy spectra will affect results of peak detection

### Function to remove peaks within close proximity (min distance between peaks)
# Peaks will be grouped that fall within distance of each other
# All but the highest peak in each group will be removed.
# peaks = peak height table (As produced by spd_peaks function)
# d = minimum distance between peaks
spd_dist <- function(peaks, d=8, suppress=FALSE) {
	# Check d is greater than 0
	if(d<1) stop("Distance value should be greater than 0")
	# Check object class to ensure correct input data
	if(class(peaks)[1]=="list") stop("\n(!) Use lapply() to use function with a list.\nFor example: lapply(peaks, function(x) spd_dist(x, d=8)")
	if(length(class(peaks)) < 2 || class(peaks)[2]!="spd_peaks") stop("\n(!) Peak/trough data required before running spd_dist() function. \nPlease run spd_peaks() or spd_troughs() function first.")
	### Get heights
	h <- peaks$height
	# Which peaks are within distance
	d_d <- ifelse(diff(peaks[,1]) < d, TRUE, FALSE)
	# Get TRUE peaks // +1 to correct index
	d_p <- which(d_d==TRUE) + 1
	# Check whether any are within distance, and only run code if there are
	if(length(d_p)>0) {
		# Split these peaks into groups
		d_gr <- split(d_p, cumsum(c(1, diff(d_p) != 1)))
		# Add the preceeding index value to each group
		d1 <- unlist(lapply(d_gr, "[", 1)) -1
		d_gr <- mapply(c, d1, d_gr, SIMPLIFY=FALSE)
		### Get the biggest peak heights from each group
		mxp_in <- lapply(d_gr, function(x) which.max(h[x]))
		mxp <- mapply("[", d_gr, mxp_in, SIMPLIFY=FALSE)
		### Combine peaks (biggest peaks < d, and peaks > d)
		d_op <- (1:length(peaks[,1]))[-as.numeric(unlist(d_gr))]
		p <- c(d_op, as.numeric(mxp)) |> sort()
		# Subset original peak data
		peaks <- peaks[p,]
		# Message how many peaks were removed
		r <- length(h) - length(peaks[,2])
		mppm <- ifelse(r==1, "peak", "peaks")
		if(suppress==FALSE) {
			message(paste0(r, mppm, " removed using the specified distance (d=", d, ").")) }
		return(peaks)
	} else {
		if(suppress==FALSE) {
		message(paste0("0 peaks removed using the specified distance (d=", d, ")."))
		}
		return(peaks)
	}
}


### Input data for spd_dist comes from the spd_peaks function ###
### Function can be called from spd_peaks directly, so no need to run seperately, but there is the option to do so ###


### Simple peak detection
# Requires R version 4.1 or higher (due to use of pipe)
# x should be numeric/matrix/data.frame of x/y values
# n = the number of neighbours to check either side of peak -
# mp = (optional) minimum peak height (e.g. peaks above the mean, or sd etc.)
# d = (optional) minimum distance between peak heights 
# suppress = logical. Change to true to suppress messages.
#test <- 2
#ifelse(test == 1, "peak", "peaks")
# type shouldn't really be changed - used internally

spd_peaks <- function(x, n=3, mp, d, type="peak", suppress=FALSE) {  
	# Check for type variable (from trough function)
	#if(missing(type)==TRUE) {
	#	type <- "peak"
	#}
    # Create indexes
    l <- length(x)
    x1 <- 1:l
    x0 <- x1 - 1
    x2 <- x1 + 1
    # Remove first and last values
    x1 <- x1[-c(1, l)]
    x0 <- x0[-c(1, l)]
    x2 <- x2[-c(1, l)]
    # Initial peak finder - find all peaks    
    p_all_test <- ifelse(x[x1] > x[x0] & x[x1] > x[x2], TRUE, FALSE)
    # Which peak (+1 to correct index)
    p_all <- which(p_all_test==TRUE) + 1
	### Check if peaks above min peak height
    if(!missing(mp)) {
        p_all <- p_all[which(as.numeric(x[p_all]) >= mp)]
		if(suppress==FALSE) {
			mpr <- length(which(as.numeric(x[p_all]) >= mp))
			mppm <- ifelse(mpr==1, "peak", "peaks")
			message(paste0(mpr, mppm, " removed which were below minimum peak height (mp=", mp, ")."))
		}
    }
    ### Check neighbours to confirm peak
    if(n > 1) {
        # Check index - if peaks fall within start-n and n-end, remove them to avoid errors
        # Consequence that peaks occuring close to the start and end of the dataset will not be caught with high n values
        p_n <- as.numeric(p_all)
        c_n <- ifelse(p_n <= n | p_n >= l-n, TRUE, FALSE)
        # Remove
        p_n <- p_n[which(c_n==FALSE)]
        # Data store
        pn_d <- matrix(nrow=length(p_n), ncol=n)
        # Tests
        j <- 0:(n-1)
        for(i in 1:n) {
                pn_d[,i] <- ifelse(x[p_n-j[i]] > x[p_n-i] & x[p_n+j[i]] > x[p_n+i], TRUE, FALSE)
        }
        # Check which peaks passed all tests
		pn_all <- apply(pn_d, 1, table) |> lapply("[", "TRUE") |> as.matrix() |> (\(x) p_n[which(x[,1]==n)])()
        # Get original index + peak height
		peaks <- data.frame(index=pn_all, peak=names(x)[pn_all], height=as.numeric(x[pn_all]), type=type)
    } else {
        # Get original index + peak height
        peaks <- data.frame(index=as.numeric(p_all), peak=names(p_all), height=as.numeric(x[p_all]), type=type)
    }
	class(peaks) <- c("data.frame", "spd_peaks")
    # Return peak data
	# Check minimum distance between peaks
	if(!missing(d)&&d>0) {
		peaks <- spd_dist(peaks=peaks, d=d, suppress=suppress)
	} else {
    	return(peaks)
	}
}


### Function to calculate slope angle of the peaks
# x should be the original data (used in the spd_peaks function)
# peaks = the output from the spd_peaks function
spd_slope <- function(x, peaks, n=3, abs=TRUE, r=2) {
	# Slope Angle in degrees
	# m = (y2 – y1)/(x2 – x1) 
	x1 <- peaks$index
	x0 <- peaks$index - n
	x2 <- peaks$index + n
	y1 <- x[peaks$index]
	y0 <- x[peaks$index - n]
	y2 <- x[peaks$index + n]
	# Do this for left and right slopes
	m0 <- atan((y0 - y1)/(x0 - x1)) * 180 / pi
	m2 <- atan((y2 - y1)/(x2 - x1)) * 180 / pi
	# Convert both angles to positives?
	if(abs==TRUE) {
		m0 <- abs(m0)
		m2 <- abs(m2)
	}
	# Mean slope angle
	mm <- rowMeans(cbind(m0, m2))
	# Percent difference between the angles
	# Greater difference means more irregular peak shape
	md <- abs(m0 - m2) / mm * 100
	# Combine
	df <- data.frame(peaks, angle_left=m0, angle_right=m2, angle_mean=round(mm, r), diff_percent=round(md, r))
	rownames(df) <- NULL
	return(df)

}

### Trough detection
# Simply inverts the data and then runs it through the peak detection function
# Set to rename output data from "peaks" to "troughs" - but this can be disabled using rename=FALSE

spd_troughs <- function(x, n=3, mp, d, rename=TRUE) {
	# Invert data
	x1 <- x * -1
	# Run through peak detection
	t <- spd_peaks(x=x1, n=n, mp=mp, d=d, type="trough")
	# Get actual trough heights from original data
	t$height <- x[t$index]
	if(rename==TRUE) {
		colnames(t)[2] <- "trough"
	}
	return(t)

}


### Plot peak height labels (onto a spectra that has already been plotted)
# Can also plot trough height labels
# x = peak height table generated by spd_peak function
# at = peak position, either "peak" (occurs at the peak height) or "top" (occurs at the top of the plot). For troughs, use either "trough", or "bottom"
# xat = whether to use index or peak height for x position
# stagger = whether to stagger the peak labels position (on the y axis)
# offset = offsets the peak label on the y axis (converts to percent of plot coordinates)
# labels = can supply your own to overwrite defaults

spd_plot_peaks <- function(x, at="peak", xat="index", offset=5, stagger=TRUE, labels, col="red", arrows=TRUE, type=1, ...) {
	### Labels
	if(!missing(labels)) {
		labels <- labels
	} else {
		labels <- x[,2]
	}
	# Get coordinates from plot
	co <- par("usr")
	y1 <- diff(co[3:4]) / 100
	# x axis
	if(xat=="index") {
		xp <- x[,1]
	} else {
		xp <- x[,2]
	}
	### Guess x axis	
	#xp <- x[,1]
	#if(abs(co[1]) > 5) {
	#	xp <- as.numeric(x[,2])
	#}
	### Label position 
	if(at=="peak") {
		lab_at <- x$height + (offset * y1)
		s <- 6
		ay <- y1
	} else if(at=="top") {
		lab_at <- rep(y1 + par("usr")[4], length(x$height))
		s <- -2
		ay <- y1
	} else if(at=="trough") {
		lab_at <- x$height + (offset * - y1)
		s <- -6
		ay <- y1 * -1
	} else if(at=="bottom") {
		lab_at <- rep(y1 + par("usr")[3], length(x$height))
		s <- -2
		ay <- y1 * -1
	}
	if(stagger==TRUE) {
		suppressWarnings(lab_at <- lab_at + c(0, (y1 * s)))
	}
	### Arrows
	if(arrows==TRUE) {
		arrows(x0=xp, y0=x$height, x1=xp, y1=lab_at - ay, code=type, length=0.15, col=col, lwd=2, xpd=NA)
	}
	### Add labels to plot
	text(x=xp, y=lab_at, labels=labels, col=col, xpd=NA, ...)
}