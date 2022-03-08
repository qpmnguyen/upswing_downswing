---
title: "Latest update for upswings"
author: "Quang Nguyen"
date: "2022-03-08"
output:
    html_document:
        keep_md: yes
        toc: true
        toc_float: true
        code_folding: hide
---



# Setting up  

```r
library(here)
library(covidcast)
library(epiprocess)
library(tidyverse)
library(ggsci)
library(tsibble)
library(patchwork)
here::i_am(path = "notebooks/01-latest-upswings.Rmd")
source(here("R", "utils.R"))
(settings <- get_settings(start_date = "2021-01-01", end_date = "2022-03-01"))
```

```
## $h
## [1] 4
## 
## $upswing_thresh
## [1] 0.5
## 
## $min_thresh
## [1] 20
## 
## $start_date
## [1] "2021-01-01"
## 
## $end_date
## [1] "2022-03-01"
```
# Load data   

For cumulative data, using confirmed cumulative numbers and aggregate by epiweek. For incidence data, using the 7-day average numbers.  

```r
print("Cumulative data")
```

```
## [1] "Cumulative data"
```

```r
(df_cum <- covidcast_signal(data_source = "jhu-csse", 
                          signal = "confirmed_cumulative_num", 
                          start_day = settings$start_date, 
                          end_day = settings$end_date, 
                          geo_type = "state",
                          geo_values = c("ma"),
                          as_of = Sys.Date()) %>% 
        as_epi_df(geo_type = "state", time_type = "day", as_of = max(.$issue)) %>%
        select(geo_value, time_value, cases = value) %>% as_tsibble() %>% 
        index_by(epiweek = ~yearweek(., week_start = 7)) %>% 
        group_by(geo_value) %>% 
        summarize(cases = max(cases)) %>% ungroup() %>% 
        as_tibble() %>% dplyr::rename(time_value = epiweek) %>% 
        as_epi_df(geo_type = "state", time_type = "week"))
```

```
## Fetched day 2021-01-01 to 2022-03-01: num_entries = 425
```

```
## An `epi_df` object, with metadata:
## * geo_type  = state
## * time_type = week
## * as_of     = 2022-03-08 13:02:54
## 
## # A tibble: 62 × 3
##    geo_value time_value  cases
##  * <chr>         <week>  <dbl>
##  1 ma          2020 W53 384181
##  2 ma          2021 W01 427135
##  3 ma          2021 W02 465726
##  4 ma          2021 W03 496093
##  5 ma          2021 W04 521360
##  6 ma          2021 W05 540827
##  7 ma          2021 W06 555895
##  8 ma          2021 W07 567764
##  9 ma          2021 W08 579680
## 10 ma          2021 W09 589931
## # … with 52 more rows
```

```r
print("Incidence data")
```

```
## [1] "Incidence data"
```

```r
(df_inc <- covidcast_signal(data_source = "jhu-csse", 
                          signal = "confirmed_7dav_incidence_num", 
                          start_day = settings$start_date, 
                          end_day = settings$end_date, 
                          geo_type = "state",
                          geo_values = c("ma"),
                          as_of = Sys.Date()) %>% 
        as_epi_df(geo_type = "state", time_type = "day", as_of = max(.$issue)) %>% 
        select(geo_value, time_value, cases = value) %>% as_tsibble() %>%
        index_by(epiweek = ~yearweek(., week_start = 7)) %>% 
        group_by(geo_value) %>% 
        summarise(cases = sum(cases)) %>% ungroup() %>% 
        as_tibble() %>% dplyr::rename(time_value = epiweek) %>% 
        as_epi_df(geo_type = "state", time_type = "week")
        
)
```

```
## Fetched day 2021-01-01 to 2022-03-01: num_entries = 425
```

```
## An `epi_df` object, with metadata:
## * geo_type  = state
## * time_type = week
## * as_of     = 2022-03-08 13:02:57
## 
## # A tibble: 62 × 3
##    geo_value time_value  cases
##  * <chr>         <week>  <dbl>
##  1 ma          2020 W53  9690.
##  2 ma          2021 W01 38105.
##  3 ma          2021 W02 42989.
##  4 ma          2021 W03 33430.
##  5 ma          2021 W04 28184.
##  6 ma          2021 W05 21432.
##  7 ma          2021 W06 17714.
##  8 ma          2021 W07 12999.
##  9 ma          2021 W08 11711.
## 10 ma          2021 W09 10940.
## # … with 52 more rows
```

# Compute upswings on real data  

The `rel_change` method in the `growth_rate` function from `epiprocess` defines relative change at focal time $T$ with bandwith $h$ as: 

$$\frac{1}{h} * \left(\frac{\bar{B}}{\bar{A}} - 1\right) = \frac{1}{h} * \left(\frac{\bar{B} - \bar{A}}{\bar{A}}\right) = \\ \frac{1}{h} * \left(\frac{(h/2)^{-1}\left(\sum_{t = T+1}^{T + h/2} Y_t - \sum_{t = T+1-h/2}^{T} Y_t\right)}{(h/2)^{-1}\sum_{t = T+1-h/2}^{T} Y_t}\right) = \frac{1}{h} R^{h/2}_{T + h/2}$$

where $R_{T+h/2}^{h/2}$ is the **the actual $h/2$-epiweek-incidence relative change** as defined in notebook 7 on computing upswings on incidence data. Using the `rel_change` option in the `growth_rate` function from `epiprocess` with bandwidth being `2*h`.    


```r
surge_inc <- df_inc %>% mutate(growth_raw = growth_rate(time_value, cases, method = "rel_change", 
                                                        h = (settings$h * 2)), 
                  growth_adj = growth_raw * settings$h * 2) %>%
    mutate(surge_raw = case_when(
        growth_raw >= settings$upswing_thresh & cases >= settings$min_thresh * settings$h ~ TRUE, 
        growth_raw < settings$upswing_thresh & cases >= settings$min_thresh * settings$h ~ FALSE,
        cases < settings$min_thresh * settings$h ~ NA
    )) %>% 
    mutate(surge_adj = case_when(
        growth_adj >= settings$upswing_thresh & cases >= settings$min_thresh * settings$h ~ TRUE, 
        growth_adj < settings$upswing_thresh & cases >= settings$min_thresh * settings$h ~ FALSE,
        cases < settings$min_thresh * settings$h ~ NA
    ))


p1 <- ggplot(surge_inc, aes(x = time_value, y = cases)) + geom_point(aes(col = surge_raw)) + theme_bw() + 
    labs(x = "Epiweek", y = "Incident cases", col = "Surge", title = "Using growth_rate function")
p2 <- ggplot(surge_inc, aes(x = time_value, y = cases)) + geom_point(aes(col = surge_adj)) + theme_bw() + 
    labs(x = "Epiweek", y = "Incident cases", col = "Surge", title = "growth_rate multiplied by bandwith")
```


```r
surge_cum <- df_cum %>% mutate(prev_cumulative = lag(cases, order_by = time_value, n = settings$h), 
    h_ew_inc = cases - prev_cumulative, 
    prev_h_ew_inc = lag(h_ew_inc, order_by = time_value, n = settings$h * 7), 
    rel_change = (h_ew_inc - prev_h_ew_inc) / prev_h_ew_inc, 
    surge = case_when(
        rel_change >= settings$upswing_thresh & h_ew_inc >= settings$min_thresh * settings$h ~ TRUE,
        rel_change < settings$upswing_thresh & h_ew_inc >= settings$min_thresh * settings$h ~ FALSE,
        h_ew_inc < settings$min_thresh * settings$h ~ NA),
    prev_surge = lag(surge, order_by = time_value, n = settings$h * 7), 
    cumulative = cases,
    inc = cases - lag(cases, order_by = time_value, n = 7),
    target_end_date = time_value)

p3 <- surge_cum %>% 
    ggplot(aes(x = time_value, y = inc, col = surge)) + geom_point() + 
    theme_bw() + labs(x = "Epiweek", y = "Incidence cases", title = "Surge calculated using cumulative data") 
    
p1 / p2 / p3    
```

```
## Warning: Removed 7 rows containing missing values (geom_point).
```

![](01-latest-upswings_files/figure-html/surge_cum-1.png)<!-- -->

# Basic evaluation for the ensemble forecaster  


```r
library(covidHubUtils)
dates <- seq(ymd(settings$start_date), 
             ymd("2022-02-01"), by = 7)
forecast_case <- load_forecasts(
    models = "COVIDhub-ensemble", 
    dates = "2022-03-01", 
    locations = "Massachusetts", 
    types = "point", 
    targets = "1 wk ahead inc case"
)
```




