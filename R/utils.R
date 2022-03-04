library(lubridate)

#' @title Return some initial settings
get_settings <- function(h = 4, upswing_thresh = 0.5,
                         min_thresh = 20,
                         start_date = ymd("2020-05-01"),
                         end_date = ymd("2022-03-01")){
    return(list(
        h = h,
        upswing_thresh = upswing_thresh,
        min_thresh = min_thresh,
        start_date = start_date,
        end_date = end_date
    ))

}
