---
title: "Evaluating the ensemble point forecast"
author: "Quang Nguyen"
date: "Last compiled on 2022-05-26"
output: 
    #github_document: default
    html_document:
        keep_md: yes
        toc: true
        toc_float: true
        code_folding: hide
---



# Data pre-loading and processing


```r
library(here)
library(covidcast)
library(epiprocess)
library(zoltr)
library(tidyverse)
library(ggsci)
library(tsibble)
library(covidHubUtils)
library(lubridate)
library(rlang)
library(patchwork)
library(pROC)
here::i_am(path = "notebooks/ensemble-evaluation.Rmd")
source(here("R", "utils.R"))
# (settings <- get_settings(start_date = "2020-06-01", end_date = "2022-03-01"))
theme_set(theme_bw())
```

First let's set some parameters


```r
d_range <- seq(ymd("2021-01-01"), ymd("2022-01-01"), by = 7)
curr_date <- "2022-01-01"
h <- 4
wk_ahead <- 1
inc_case_targets <- paste(1:h, "wk ahead inc case")
surge_thresh <- 0.5
min_inc <- 20
```

# Loading forecasts and underlying data

Using `covidHubUtils` and `zoltar` we load underlying data as well as forecasts


```r
pred_case <- load_forecasts(
    models = "COVIDhub-ensemble", 
    dates = d_range, 
    date_window_size = 6,
    locations = "Massachusetts", 
    types = "point", 
    targets = inc_case_targets, 
    source = "zoltar", 
    verbose = FALSE, 
    as_of = curr_date, 
    hub = c("US")
)
```

```
## get_token(): POST: https://zoltardata.com/api-token-auth/
```

```
## get_resource(): GET: https://zoltardata.com/api/projects/
```

```
## get_resource(): GET: https://zoltardata.com/api/project/44/models/
```

```
## get_resource(): GET: https://zoltardata.com/api/project/44/timezeros/
```

```r
truth_data <- load_truth(
    truth_source = "JHU", 
    target_variable = "inc case", 
    locations = "Massachusetts"
)

true_range <- pred_case %>% pull(target_end_date) %>% 
    unique() %>% 
    lubridate::as_date()

truth_epidf <- truth_data %>% 
    select(-c(model, location, target_variable, location_name, 
              abbreviation, full_location_name)) %>% 
    dplyr::rename("time_value" = "target_end_date") %>%
    filter(time_value %in% true_range) %>%
    as_epi_df(geo_type = "state")
```

# Surge classification using relative change growth rate formulation

We define surge for a given date using relative change growth rate formulation times the bandwidth using the implementation from `epiprocess`: 

$$\frac{1}{h} * \left(\frac{\bar{B}}{\bar{A}} - 1\right) = \frac{1}{h} * \left(\frac{\bar{B} - \bar{A}}{\bar{A}}\right) = \\ \frac{1}{h} * \left(\frac{(h)^{-1}\left(\sum_{t = T+1}^{T + h} Y_t - \sum_{t = T+1-h}^{T} Y_t\right)}{(h)^{-1}\sum_{t = T+1-h}^{T} Y_t}\right) = \frac{1}{h} R^{h}_{T + h}$$

A surge is defined for time-point $T$ as the difference in cumulative incident cases between the periods of $T+1$ and $T+h$ and $T$ and $T-h$. As such, an $h$-week ahead forecaster is a nowcaster of whether or not we're currently in a surge.


```r
truth_epidf <- truth_epidf %>% 
    mutate(gr = growth_rate(y = value, method = "rel_change", h = h) * h) %>%
    mutate(surge = case_when(
        gr >= surge_thresh & value >= min_inc ~ TRUE,
        TRUE ~ FALSE
    ))

ggplot(truth_epidf, aes(x = time_value, y = value)) +
    geom_point(aes(col = surge), size = 2.5) + geom_line(alpha = 0.5) + 
    labs(x = "Weekly data", y = "Incident cases by week")
```

![](ensemble-evaluation_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

# Nowcasting surges    

Here we use the following procedures for a focal timepoint $T$ and bandwidth $h$ (for example, 4)  

1. We take a time period from $T+1- h$ to $T+h$. For example at 2021-01-23, we'd be taking the period from 2021-01-02 to 2021-02-20.     
2. We take the truth period to be from $T+1-h$ to $T$ (inclusive) and the forecasting period to be from $T+1$ to $T+h$. The truth period would have real underlying incident cases while the forecasting period has forecast incident cases at times 1-$h$ weeks ahead. Due to the forecasting date being on Monday instead of exactly one week before the proposed target date, we take forecast incident values from the forecast date closest to the time period defined at $T+1$. For example, at 2021-01-23, we would take forecast values for 2021-01-30 onwards from a forecast date of 2021-01-25.  
3. We compute the growth rate at time $T$ using these two periods as per the formula above   
4. We then classify periods as surges using the definition and thresholds defined above.  


```r
# this function combines real case counts from time points t-h to t and forecasted 
# case counts from t+1 to t+h. Growth rate at time t is then estimated using the relative change 
# method
mismatch_slide <- function(slide_df, pred, h){
    query_dates <- slide_df$time_value
    req_len <- h * 2
    # if not enough weeks for prediction
    if (length(query_dates) != req_len){
        return(NA_real_)
    } else {
        t_date <- query_dates[1:(req_len - h)]
        ref <- tail(t_date, n = 1)
        p_date <- query_dates[(req_len - h + 1):req_len]
        t_data <- slide_df %>% filter(time_value %in% t_date) %>% 
            mutate(type = "true")
        # due to weird issues, the ref date + 1 week should be the 1 week ahead forecast
        f_date <- pred %>% filter(target_end_date == head(p_date, n = 1) & horizon == 1) %>% 
            pull(forecast_date) %>% unique()
        p_data <- pred %>% 
            filter(forecast_date == f_date) %>%
            filter(target_end_date %in% p_date) %>% 
            select(geo_value, target_end_date, value, population, geo_type) %>%
            dplyr::rename("time_value" = target_end_date) %>%
            mutate(type = "pred")
        combine <- bind_rows(t_data, p_data)
        gr <- combine %>% mutate(gr_pred = growth_rate(y = value, h = h, method = "rel_change") * h) %>% 
            filter(time_value == ref) %>% pull(gr_pred)
        return(gr)
    }
}
combined_df <- truth_epidf %>% epi_slide(~mismatch_slide(slide_df = .x, pred = pred_case, h = h), 
                                          n = 2 * 7 * h, align = "center", new_col_name = "pred_gr")
combined_df <- combined_df %>% mutate(surge_pred = case_when(
    pred_gr >= surge_thresh & value >= min_inc ~ TRUE,
    is.na(pred_gr) ~ NA,
    TRUE ~ FALSE
))
head(combined_df, n = 20)
```

```
## # A tibble: 20 × 9
##    geo_value time_value value population geo_type      gr surge pred_gr
##    <chr>     <date>     <dbl>      <dbl> <chr>      <dbl> <lgl>   <dbl>
##  1 ma        2021-01-02 34579    6892503 state    -0.0132 FALSE NA     
##  2 ma        2021-01-09 42954    6892503 state    -0.356  FALSE NA     
##  3 ma        2021-01-16 38591    6892503 state    -0.477  FALSE NA     
##  4 ma        2021-01-23 30367    6892503 state    -0.511  FALSE -0.246 
##  5 ma        2021-01-30 25267    6892503 state    -0.575  FALSE -0.415 
##  6 ma        2021-02-06 19467    6892503 state    -0.568  FALSE -0.470 
##  7 ma        2021-02-13 15068    6892503 state    -0.504  FALSE -0.571 
##  8 ma        2021-02-20 11869    6892503 state    -0.371  FALSE -0.562 
##  9 ma        2021-02-27 11916    6892503 state    -0.179  FALSE -0.405 
## 10 ma        2021-03-06 10251    6892503 state     0.0855 FALSE -0.329 
## 11 ma        2021-03-13 10685    6892503 state     0.274  FALSE -0.148 
## 12 ma        2021-03-20 12241    6892503 state     0.295  FALSE  0.108 
## 13 ma        2021-03-27 14699    6892503 state     0.128  FALSE  0.279 
## 14 ma        2021-04-03 15676    6892503 state    -0.117  FALSE  0.254 
## 15 ma        2021-04-10 14346    6892503 state    -0.316  FALSE -0.0160
## 16 ma        2021-04-17 13672    6892503 state    -0.482  FALSE -0.132 
## 17 ma        2021-04-24 10320    6892503 state    -0.570  FALSE -0.463 
## 18 ma        2021-05-01  8709    6892503 state    -0.649  FALSE -0.547 
## 19 ma        2021-05-08  6261    6892503 state    -0.707  FALSE -0.600 
## 20 ma        2021-05-15  4950    6892503 state    -0.757  FALSE -0.619 
## # … with 1 more variable: surge_pred <lgl>
```

Plotting growth rate when estimated 


```r
gr1 <- ggplot(combined_df, aes(x = time_value, y = gr)) + 
    geom_line(aes(col = "Est. using real data")) +
    geom_line(aes(x = time_value, y = pred_gr, 
                  col = "Est. using real data and forecast")) + 
    labs(col ="Growth Rate", x = "Week", y = "Growth Rate (adjusted for bandwidth)") +
    geom_hline(yintercept = 0.5, col = "red")

combined_df %>% filter(month(time_value) %in% c(6,7,8))
```

```
## An `epi_df` object, with metadata:
## * geo_type  = state
## * time_type = day
## * as_of     = 2022-05-26 13:44:45
## 
## # A tibble: 13 × 9
##    geo_value time_value value population geo_type     gr surge pred_gr
##  * <chr>     <date>     <dbl>      <dbl> <chr>     <dbl> <lgl>   <dbl>
##  1 ma        2021-06-05  1157    6892503 state    -0.807 FALSE  -0.779
##  2 ma        2021-06-12   896    6892503 state    -0.729 FALSE  -0.736
##  3 ma        2021-06-19   518    6892503 state    -0.360 FALSE  -0.730
##  4 ma        2021-06-26   384    6892503 state     0.842 TRUE   -0.695
##  5 ma        2021-07-03   400    6892503 state     3.39  TRUE   -0.598
##  6 ma        2021-07-10   692    6892503 state     6.81  TRUE   -0.284
##  7 ma        2021-07-17  1442    6892503 state     6.61  TRUE    0.684
##  8 ma        2021-07-24  2908    6892503 state     4.25  TRUE    1.24 
##  9 ma        2021-07-31  4600    6892503 state     2.53  TRUE    1.19 
## 10 ma        2021-08-07  6615    6892503 state     1.52  TRUE    0.888
## 11 ma        2021-08-14  8091    6892503 state     0.915 TRUE    0.476
## 12 ma        2021-08-21  9249    6892503 state     0.654 TRUE    0.304
## 13 ma        2021-08-28 10060    6892503 state     0.442 FALSE   0.181
## # … with 1 more variable: surge_pred <lgl>
```


Plotting incident cases and surge classification


```r
p1 <- ggplot(combined_df, aes(x = time_value, y = value)) +
    geom_point(aes(col = surge)) + 
    geom_line(alpha = 0.5) + 
    labs(x = "Week", col = "Surge classification from real data", y = "Incident cases")
p2 <- ggplot(combined_df, aes(x = time_value, y = value)) + 
    geom_point(aes(col = surge_pred)) + 
    geom_line(alpha = 0.5) + 
    labs(x = "Week", col = "Nowcasting surges", y = "Incident cases")

p1 /p2
```

![](ensemble-evaluation_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

# Misclassification rate    


```r
combined_df %>% filter(!is.na(surge_pred)) %>% 
    summarise(misclass_ensemble = mean(surge != surge_pred), 
              misclass_never = mean(surge),
              sens = sum(surge & surge_pred)/sum(surge), 
              spec = sum(!surge & !surge_pred)/sum(!surge))
```

```
## # A tibble: 1 × 4
##   misclass_ensemble misclass_never  sens  spec
##               <dbl>          <dbl> <dbl> <dbl>
## 1             0.174          0.326 0.467     1
```

```r
pROC::roc(surge ~ pred_gr, data = combined_df, subset = !is.na(surge_pred))
```

```
## Setting levels: control = FALSE, case = TRUE
```

```
## Setting direction: controls < cases
```

```
## 
## Call:
## roc.formula(formula = surge ~ pred_gr, data = combined_df, subset = !is.na(surge_pred))
## 
## Data: pred_gr in 31 controls (surge FALSE) < 15 cases (surge TRUE).
## Area under the curve: 0.8129
```

