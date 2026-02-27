va2pix <- function(va, scr) {
  # Convert visual angle (deg) to pixels.
  scr$subDist * tan(va * pi / 180) / (scr$width / (10 * scr$xres))
}

movmean_partial <- function(x, k) {
  n <- length(x)
  left <- floor((k - 1) / 2)
  right <- k - 1 - left
  out <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    lo <- max(1L, i - left)
    hi <- min(n, i + right)
    w <- x[lo:hi]
    out[i] <- if (all(is.na(w))) NA_real_ else mean(w, na.rm = TRUE)
  }
  out
}

movmean2d <- function(x, k) {
  cbind(movmean_partial(x[, 1], k), movmean_partial(x[, 2], k))
}

vecvel <- function(xx, sampling, type = 2L) {
  if (!is.matrix(xx)) {
    xx <- as.matrix(xx)
  }
  n <- nrow(xx)
  v <- matrix(0, n, 2)
  if (type == 1L) {
    if (n >= 3L) {
      v[2:(n - 1), ] <- (sampling / 2) * (xx[3:n, ] - xx[1:(n - 2), ])
    }
  } else if (type == 2L) {
    if (n >= 5L) {
      v[3:(n - 2), ] <- (sampling / 6) * (
        xx[5:n, ] + xx[4:(n - 1), ] - xx[2:(n - 3), ] - xx[1:(n - 4), ]
      )
    }
    if (n >= 3L) {
      v[2, ] <- (sampling / 2) * (xx[3, ] - xx[1, ])
      v[n - 1, ] <- (sampling / 2) * (xx[n, ] - xx[n - 2, ])
    }
  } else {
    stop("type must be 1 or 2")
  }
  v
}

microsacc_merge <- function(x, vel, vfac = 5, mindur = 8, merge_interval = 10) {
  vx <- vel[, 1]
  vy <- vel[, 2]
  
  msdx <- sqrt(stats::median(vx^2, na.rm = TRUE) - stats::median(vx, na.rm = TRUE)^2)
  msdy <- sqrt(stats::median(vy^2, na.rm = TRUE) - stats::median(vy, na.rm = TRUE)^2)
  
  if (!is.finite(msdx) || msdx < .Machine$double.xmin) {
    msdx <- sqrt(mean(vx^2, na.rm = TRUE) - mean(vx, na.rm = TRUE)^2)
  }
  if (!is.finite(msdy) || msdy < .Machine$double.xmin) {
    msdy <- sqrt(mean(vy^2, na.rm = TRUE) - mean(vy, na.rm = TRUE)^2)
  }
  
  radiusx <- vfac * msdx
  radiusy <- vfac * msdy
  radius <- c(radiusx, radiusy)
  
  if (!is.finite(radiusx) || !is.finite(radiusy) || radiusx == 0 || radiusy == 0) {
    return(list(msac = matrix(numeric(0), ncol = 7), radius = radius))
  }
  
  test <- (vx / radiusx)^2 + (vy / radiusy)^2
  indx <- which(test > 1)
  
  if (length(indx) == 0) {
    return(list(msac = matrix(numeric(0), ncol = 7), radius = radius))
  }
  
  starts <- c(indx[1], indx[which(diff(indx) > 1) + 1])
  ends <- c(indx[which(diff(indx) > 1)], indx[length(indx)])
  durs <- ends - starts + 1
  keep <- durs >= mindur
  
  sac <- cbind(starts[keep], ends[keep])
  if (nrow(sac) == 0) {
    return(list(msac = matrix(numeric(0), ncol = 7), radius = radius))
  }
  
  merged <- sac[1, , drop = FALSE]
  if (nrow(sac) > 1) {
    for (i in 2:nrow(sac)) {
      last <- nrow(merged)
      if ((sac[i, 1] - merged[last, 2]) <= merge_interval) {
        merged[last, 2] <- sac[i, 2]
      } else {
        merged <- rbind(merged, sac[i, , drop = FALSE])
      }
    }
  }
  
  signed_range <- function(v) {
    idx <- which(!is.na(v))
    if (length(idx) == 0) return(NA_real_)
    vv <- v[idx]
    i_min <- idx[which.min(vv)[1]]
    i_max <- idx[which.max(vv)[1]]
    sign(i_max - i_min) * (max(vv) - min(vv))
  }
  
  msac <- cbind(merged, matrix(NA_real_, nrow(merged), 5))
  for (s in seq_len(nrow(merged))) {
    a <- merged[s, 1]
    b <- merged[s, 2]
    vm <- sqrt(vel[a:b, 1]^2 + vel[a:b, 2]^2)
    msac[s, 3] <- if (all(is.na(vm))) NA_real_ else max(vm, na.rm = TRUE)
    
    dx <- x[b, 1] - x[a, 1]
    dy <- x[b, 2] - x[a, 2]
    msac[s, 4] <- dx
    msac[s, 5] <- dy
    
    msac[s, 6] <- signed_range(x[a:b, 1])
    msac[s, 7] <- signed_range(x[a:b, 2])
  }
  
  list(msac = msac, radius = radius)
}

saccpar <- function(sac) {
  if (is.null(sac) || nrow(sac) == 0) {
    return(matrix(numeric(0), ncol = 8))
  }
  a <- sac[, 1]
  b <- sac[, 2]
  d <- b - a + 1
  vpeak <- sac[, 3]
  dist <- sqrt(sac[, 4]^2 + sac[, 5]^2)
  angd <- atan2(sac[, 5], sac[, 4])
  ampl <- sqrt(sac[, 6]^2 + sac[, 7]^2)
  anga <- atan2(sac[, 7], sac[, 6])
  cbind(a, b, d, vpeak, dist, angd, ampl, anga)
}

closest_index <- function(t, time_vec) {
  m <- match(t, time_vec)
  if (!is.na(m)) return(m)
  which.min(abs(time_vec - t))
}
