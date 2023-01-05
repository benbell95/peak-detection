### Peak area calculation functions
### Benjamin Bell 
### https://github.com/benbell95/peak-detection
### https://www.benjaminbell.co.uk
### Work in progress - use at own risk
### Noisy spectra will affect results of peak detection

###### Requires peak_detection.r functions to work ######

### Shoelace function
# This is for calculation the area
spd_shoelace <- function(x, y) {
    # Check lengths of x and y match
    if(length(x) != length(y)) stop("Lengths of x and y must match")
    # Check for NA's and return NA if true
    if(any(is.na(x))|any(is.na(y))==TRUE) {
        a <- NA
        return(a)
        break
    }
    # Get length
    n <- length(x)
    # Create empty vector for results
    r <- vector("numeric", n)
    # Create secondary index
    j <- c(2:n, 1)
    # Calculate area
    for(i in 1:n) {
        r[i] <- x[i] * y[j[i]] - x[j[i]] * y[i]
    }
    # Calculate sum (as absolute value), then halve
    a <- abs(sum(r)) * 0.5
    return(a)
}

### Refine function for peak area detection
# This uses linear regression to refine the boundaries of the peak
# x = pxy from spd_area function
# pxy is a list of peak area x and y values
spd_area_refine <- function(x, r) {
    # Get first and last values
    pxy_h <- lapply(x, "[", 1,)
    pxy_t <- lapply(x, \(x) x[lengths(x[1]),])
    # Combine values and perform regression
    pxy_lm <- Map(function(x, y) rbind(x,y), x=pxy_h, y=pxy_t) |> lapply(\(x) if(class(x)=="data.frame") lm(y ~ x_in, data=x))
    # Predict new values from lm
    nd <- lapply(x, "[", 1)
    pxy_p <- Map(function(x, y) if(class(x)=="lm") predict(x, newdata=y), x=pxy_lm, y=nd)
    # Check for negative y values from real to predicted.
    # Must round down the y values, otherwise get errors - 1st and last value should both equal 0
    # Get y values
    pxy_y <- lapply(x, "[", 3)
     # Calculate differences
    pxy_d <- Map(function(x, y) if(class(x)=="data.frame") round(x, r) - round(y, r), x=pxy_y, y=pxy_p)
    # Get length of data
    pxy_dl <- lapply(pxy_d, \(x) if(class(x)=="data.frame") length(x[,1]))
    # Split data in two
    pxy_d2 <- lapply(pxy_d, function(x) if(class(x)=="data.frame") split(x, cut(1:nrow(x), 2, labels=FALSE)))
    # Get length of split data
    pxy_d2l <- lapply(pxy_d2, \(x) lapply(x, \(x) if(class(x)=="data.frame") length(x[,1])))
    # ...and find negative values
    pxy_n <- lapply(pxy_d2, \(x) lapply(x, \(x) if(class(x)=="data.frame") which(x[,1] < 0)))
    # Add the length of split data 1 to split data 2
    for(i in seq_along(pxy_n)) {
        if(length(pxy_n[[i]]$'2')==0)
            next
        # Otherwise add
        pxy_n[[i]]$'2' <- pxy_n[[i]]$'2' + pxy_d2l[[i]]$'1'
    }
    # Remove values 
    for(i in seq_along(x)) {
        if(class(x[[i]])!="data.frame")
            next
        if(length(pxy_n[[i]]$'2')>0) {
            # 1st pass - remove end data points first
            x[[i]] <- x[[i]][-c((min(pxy_n[[i]]$'2')):pxy_dl[[i]]),]
        }
        if(length(pxy_n[[i]]$'1')>0) {
            # 2nd pass - remove starting data points
            x[[i]] <- x[[i]][-c(1:((max(pxy_n[[i]]$'1') -1))),]
        }
    }
    return(x)
}

### Function to detect and calculate peak area
# Requires shoelace function
# Requires functions from peak detection script
# Get peak area for "single" peak - chosen by user
# Uses the peak detection function
# Uses the shoelace function to calculate area
# Input data = matrix of spectra (e.g. individual grains in a sample)
# Use lapply to run it on a list of matrices
# Peak should be the peak that you want the area for
# Can also include adjacent peaks e.g. "1512|1516|1520"
# Must use grep syntax for multiple peaks
### n, mp, d - from spd_peaks function, use where ...
# area - calulate on "in" (index) or "wn" (wavenumber)
# Refine - if TRUE, perform linear regression to refine peak boundaries
# r - round peak heights for linear regression (refine=TRUE)

spd_area <- function(x, peak, n=3, area="in", refine=TRUE, r=8, ...) {
    # Get wavenumber resolution
    wnr <- as.numeric(row.names(x)[1]) - as.numeric(row.names(x)[2]) 
    # Get peaks first - user defined peak criteria
    p <- apply(x, 2, \(x) spd_peaks(x, n=n, suppress=TRUE, ...))
    # Get troughs - get all possible troughs!
    t <- apply(x, 2, \(x) spd_troughs(x, n=1, rename=FALSE))
    # Combine peak and trough data, then reorder and reset rownames
    pt <- mapply(rbind, p, t, SIMPLIFY=FALSE) |> lapply(\(x) "[" (x[order(x$index),])) |> lapply(\(x) {row.names(x) <- NULL;x})
    # Get the row index for the peak of interest 
    pt_p <- lapply(pt, \(x) "["(x, grep(peak, x$peak),)) |> lapply(\(x) "["(x, grep("peak", x$type),)) |> lapply(\(x) as.numeric(row.names(x)))
    # Get adjacent troughs
    l <- length(pt)
    # In a loop
    pt_pp <- vector("list", l)
    for(i in 1:l) {
       pt_pp[[i]] <- pt[[i]][c(pt_p[[i]]-1, pt_p[[i]]+1),] 
    }
    # Get index values
    pt_ppi <- lapply(pt_pp, \(x) "[" (x$index))
    # Get peak values
    pt_ppp <- lapply(pt_pp, \(x) "[" (as.numeric(x$peak)))
    # For peak area, get x and y values
    # In a loop
    pxy <- vector("list", l)
    for(i in 1:l) {
        # Skip when no data
        if(length(pt_ppi[[i]])==0)
            next
        x_in <- pt_ppi[[i]][1]:pt_ppi[[i]][2]
        x_p <- seq(pt_ppp[[i]][1], pt_ppp[[i]][2], -wnr)
        pxy[[i]] <- data.frame(x_in=x_in, x_p=x_p, y=x[,i][x_in])
    }
    # Refine peak boundaries using regression
    if(refine==TRUE) {
        pxy <- spd_area_refine(x=pxy, r=r)
    }
    # Get area
    ax <- ifelse(area=="in", 1, 2)
    area <- lapply(pxy, \(x) if(length(x)>0) spd_shoelace(x=x[,ax], y=x[,3]))
    # Return data in a list
    data <- list(peak_xy=pxy, area=area, peaks=p, troughs=t)
    return(data)
}

#### Manual peak area calculation (calculate area between two wavenumbers)
# Bands should be two peak numbers
# Peaks should be in same order as the data, e.g. ascending or decending
spd_area_bands <- function(x, bands, area="in", refine=TRUE, r=8) {
    # Check length of bands (must be 2)
    if(length(bands)!=2) stop("Must supply two wavenumbers")
    # Get wavenumber resolution
    wnr <- as.numeric(row.names(x)[1]) - as.numeric(row.names(x)[2]) 
    # Get index value from peak numbers (bands)
    in1 <- which(row.names(x)==bands[1])
    in2 <- which(row.names(x)==bands[2])
    # subset data
    x <- x[in1:in2,]
    # For peak area, get x and y values
    l <- ncol(x)
    x_in <- in1:in2
    x_p <- seq(bands[1], bands[2], -wnr)
    # In a loop
    pxy <- vector("list", l)
    for(i in 1:l) {
        pxy[[i]] <- data.frame(x_in=x_in, x_p=x_p, y=x[,i])
    }
    # Refine peak boundaries using regression
    if(refine==TRUE) {
        pxy <- spd_area_refine(x=pxy, r=r)
    }
    # Get area
    ax <- ifelse(area=="in", 1, 2)
    area <- lapply(pxy, \(x) if(length(x)>0) spd_shoelace(x=x[,ax], y=x[,3]))
    # list
    area <- list(peak_xy=pxy, area=area)
    return(area)
}

### Peak area stats
# outlier=TRUE - remove outliers (using IQR)
# m = outlier removal IQR multiplier
spd_area_stats <- function(x, outlier=TRUE, m=1.5) {
    a <- unlist(x)
    # Remove inf values (occur with ratio data)
    a <- a[which(a!=Inf)]
    # Remove outliers
    if(outlier==TRUE) {
        q <- quantile(a)
        iqr <- IQR(a)
        w <- which(a > q[4] + (iqr * m) | a < q[2] - (iqr * m))
        # remove outliers if present
        if(length(w)>0) {
            a <- a[-c(w)]
        }
    }
    # Stats
    me <- mean(a)
    med <- median(a)
    sd <- sd(a)
    se <- sd / sqrt(length(a))
    # Combine
    data <- list(area=a, mean=me, median=med, std_deviation=sd, std_error=se)
    return(data)
}

### Peak area ratio function
# x and y should be the peak area data (a list)
# log - either TRUE or FALSE to log the results
# stats=TRUE - group ratios and provide stats

spd_area_ratio <- function(x, y, log=FALSE, stats=TRUE, outlier=TRUE, m=1.5) {
    if(log==TRUE) {
        rat <- Map(function(x, y) {abs(log(x / y))}, x=x$area, y=y$area)
    } 
    else {
        rat <- Map(function(x, y) {abs(x / y)}, x=x$area, y=y$area)
    }
    if(stats==TRUE) {
        rat <- spd_area_stats(rat, outlier=outlier, m=m)
        names(rat)[1] <- "ratio"
    }
    return(rat)
}