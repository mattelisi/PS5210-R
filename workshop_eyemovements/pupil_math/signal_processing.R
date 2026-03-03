#' Signal processing
#'
#' Gap fill-in (cubic-spline / linear interpolation).
#' Used to interpolate blinks in pupil-size recordings
#' (e.g., Lisi, Bonato & Zorzi, Biological Psychology, 2015).
#'
#' @param x       Numeric vector (signal). NAs mark the gap(s) to fill.
#' @param max     Maximum gap length in milliseconds to interpolate (default 50 ms).
#' @param sp      Sampling period in seconds (default 1/60 s).
#' @param type    Interpolation type: "linear" (default) or "cubic".
#' @return        The signal with eligible gaps filled in (invisibly).
#' @export
fillGap <- function(x, sp = 1/60, max = 50, type = "linear") {
  
  if (!type %in% c("linear", "cubic")) {
    warning("Unknown gap interpolation 'type' parameter!")
    return(invisible(x))
  }
  
  n        <- length(x)
  ID       <- seq_len(n)
  x_interp <- x
  
  # --- identify contiguous runs of NAs via run-length encoding ---------------
  r       <- rle(is.na(x))
  run_end <- cumsum(r$lengths)       # last index of each run
  run_beg <- run_end - r$lengths + 1 # first index of each run
  
  for (k in seq_along(r$lengths)) {
    if (!r$values[k]) next           # not a NA run, skip
    
    beg  <- run_beg[k]
    end  <- run_end[k]
    size <- r$lengths[k]
    gap_ms <- size * (sp * 1000)     # gap duration in milliseconds
    
    # ---- linear interpolation ------------------------------------------------
    if (type == "linear") {
      if (gap_ms > max)  next        # gap too long — leave as NA
      if (beg == 1 || end == n) next # gap at signal boundary — can't interpolate
      
      lvs <- x_interp[beg - 1]      # last valid sample before gap
      fvs <- x_interp[end + 1]      # first valid sample after gap
      
      if (is.na(lvs) || is.na(fvs)) next   # neighbours themselves NA
      
      # denominator = size+1 so interpolation spans [lvs, fvs] exclusive
      m <- (fvs - lvs) / (size + 1)
      x_interp[beg:end] <- lvs + m * seq_len(size)
    }
    
    # ---- cubic-spline interpolation ------------------------------------------
    # For cubic, large gaps are simply left as NA (not excluded via sentinel).
    # The spline is fitted once after the loop using all remaining valid points.
    # We don't need per-gap logic here; gaps > max will remain NA and are
    # naturally excluded from the knot set below.
  }
  
  # ---- cubic: one global spline pass over all still-NA positions -------------
  if (type == "cubic") {
    remaining_na <- which(is.na(x_interp))
    
    if (length(remaining_na) > 0) {
      # Determine which NA runs are within the allowed gap length
      r2       <- rle(is.na(x_interp))   # re-run on (possibly updated) signal
      run_end2 <- cumsum(r2$lengths)
      run_beg2 <- run_end2 - r2$lengths + 1
      fill_idx <- integer(0)
      
      for (k in seq_along(r2$lengths)) {
        if (!r2$values[k]) next
        beg2   <- run_beg2[k]
        end2   <- run_end2[k]
        size2  <- r2$lengths[k]
        gap_ms <- size2 * (sp * 1000)
        if (gap_ms <= max && beg2 > 1 && end2 < n) {
          fill_idx <- c(fill_idx, beg2:end2)
        }
      }
      
      if (length(fill_idx) > 0) {
        ok_idx <- which(!is.na(x))   # use original valid samples as knots
        x_interp[fill_idx] <- spline(ID[ok_idx], x[ok_idx],
                                     xout = fill_idx)$y
      }
    }
  }
  
  invisible(x_interp)
}

#' Miscellaneous helpers
#'
#' Determine if a number is even
#' @param x number to test
#' @return TRUE if x is even, FALSE if odd
#' @export
is.even <- function(x) { x %% 2 == 0 }


#' Signal processing
#'
#' Infinite-impulse-response (IIR) single-pole lowpass filter.
#' @param x  Numeric vector (signal). NAs are passed through unchanged.
#' @param fc Cut-off frequency (Hz).
#' @param sp Sampling period in seconds (default 1/60 s).
#' @return   Filtered signal (same length as x), returned invisibly.
#' @export
filtLP <- function(x, fc, sp = 1/60) {
  RC    <- 1 / (2 * pi * fc)
  alpha <- sp / (sp + RC)       # smoothing factor
  
  y <- x                        # write to separate output, preserve input
  
  for (i in seq_along(x)) {
    if (is.na(x[i])) next       # leave NAs untouched
    
    if (i == 1 || is.na(x[i - 1])) {
      y[i] <- x[i]              # no previous valid sample: pass through
    } else {
      y[i] <- alpha * x[i] + (1 - alpha) * y[i - 1]  # use filtered previous
    }
  }
  
  invisible(y)
}


#' Signal processing
#'
#' Central (symmetric) non-weighted moving-average filter.
#' Each output sample is the mean of the \code{2*n+1} samples centred on it.
#' The first and last \code{n} samples are left unfiltered (copied as-is).
#'
#' @param x Numeric vector (signal).
#' @param n Half-window size (must be a positive integer; total window = 2n+1).
#' @return  Filtered signal (same length as x), returned invisibly.
#' @export
filtCMA <- function(x, n) {
  if (is.even(n)) {
    warning("Filter half-window size 'n' should be odd so the full window (2n+1) is symmetric!")
    return(invisible(NULL))
  }
  
  len <- length(x)
  if (is.null(x) || len == 0) return(invisible(NULL))
  
  y <- x                          # pre-fill with original (handles edges)
  
  # only compute where the full symmetric window fits inside the signal
  for (i in seq_len(len)) {
    if (is.na(x[i])) next
    if (i > n && i <= len - n) {
      y[i] <- mean(x[(i - n):(i + n)], na.rm = TRUE)
    }
  }
  
  invisible(y)
}


#' Signal processing
#'
#' Central (symmetric) moving-median filter.
#' Each output sample is the median of the \code{2*n+1} samples centred on it.
#' The first and last \code{n} samples are left unfiltered (copied as-is).
#'
#' @param x Numeric vector (signal).
#' @param n Half-window size (must be a positive integer; total window = 2n+1).
#' @return  Filtered signal (same length as x), returned invisibly.
#' @export
filtMM <- function(x, n) {
  if (is.even(n)) {
    warning("Filter half-window size 'n' should be odd so the full window (2n+1) is symmetric!")
    return(invisible(NULL))
  }
  
  len <- length(x)
  if (is.null(x) || len == 0) return(invisible(NULL))
  
  y <- x                          # pre-fill with original (handles edges)
  
  # only compute where the full symmetric window fits inside the signal
  for (i in seq_len(len)) {
    if (is.na(x[i])) next
    if (i > n && i <= len - n) {
      y[i] <- median(x[(i - n):(i + n)], na.rm = TRUE)
    }
  }
  
  invisible(y)
}

