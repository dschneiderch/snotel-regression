#function to calculate wateryear from date column
water_yr <- function(dates, start_month=10) {
     # Convert dates into POSIXlt
     dates.posix = as.POSIXlt(dates)
     # Year offset
     offset = ifelse(dates.posix$mon >= start_month - 1, 1, 0)
     # Water year
     adj.year = dates.posix$year + 1900 + offset
     # Return the water year
     adj.year
}
