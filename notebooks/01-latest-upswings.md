---
title: "Latest update for upswings"
author: "Quang Nguyen"
date: "2022-03-03"
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
here::i_am(path = "notebooks/01-latest-upswings.Rmd")
source(here("R", "utils.R"))
(settings <- get_settings(start_date = "2021-01-01", end_date = "2021-12-31"))
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
## [1] "2021-12-31"
```
# Load data 

```r
cases <- covidcast_signal(data_source = "jhu-csse", 
                          signal = "confirmed_incidence_num", 
                          start_day = settings$start_date, 
                          end_day = settings$end_date, 
                          geo_type = "state",
                          geo_values = c("ma"),
                          as_of = Sys.Date())
```

```
## Fetched day 2021-01-01 to 2021-12-31: num_entries = 365
```

```r
(df <- as_epi_df(cases, geo_type = "state", time_type = "day", as_of = max(cases$issue)) %>% 
        select(geo_value, time_value, cases = value) %>% as_tsibble() %>%
        index_by(epiweek = ~yearweek(., week_start = 7)) %>% 
        group_by(geo_value) %>% 
        summarise(cases = sum(cases)) %>% ungroup() %>% 
        as_tibble() %>% rename(time_value = epiweek) %>% 
        as_epi_df(geo_type = "state", time_type = "day", as_of = max(cases$issue))
        
)
```

```
## An `epi_df` object, with metadata:
## * geo_type  = state
## * time_type = day
## * as_of     = 2022-01-01
## 
## # A tibble: 53 x 3
##    geo_value time_value cases
##  * <chr>         <week> <dbl>
##  1 ma          2020 W53  9003
##  2 ma          2021 W01 42954
##  3 ma          2021 W02 38591
##  4 ma          2021 W03 30367
##  5 ma          2021 W04 25267
##  6 ma          2021 W05 19467
##  7 ma          2021 W06 15068
##  8 ma          2021 W07 11869
##  9 ma          2021 W08 11916
## 10 ma          2021 W09 10251
## # ... with 43 more rows
```

# Compute upswings on real data  

The `rel_change` method in the `growth_rate` function from `epiprocess` defines relative change at focal time $T$ with bandwith $h$ as: 

$$\frac{1}{h} * \left(\frac{\bar{B}}{\bar{A}} - 1\right) = \frac{1}{h} * \left(\frac{\bar{B} - \bar{A}}{\bar{A}}\right) = 
\frac{1}{h} * \left(\frac{(h/2)^{-1}\left(\sum_{t = T+1}^{T + h/2} Y_t - \sum_{t = T+1-h/2}^{T} Y_t\right)}{(h/2)^{-1}\sum_{t = T+1-h/2}^{T} Y_t}\right) = \frac{1}{h} R^{h/2}_{T + h/2}$$

where $R_{T+h/2}^{h/2}$ is the **the actual $h/2$-epiweek-incidence relative change** as defined in notebook 7 on computing upswings on incidence data. 

Using the `rel_change` option in the `growth_rate` function from `epiprocess` with `h` being the set threshold   


```r
df %>% mutate(c_gr_raw = growth_rate(time_value, cases, method = "rel_change", h = settings$h * 2), 
              c_gr = c_gr_raw * settings$h * 2)


ggplot(df, aes(x = time_value, y = cases)) + 
           geom_point(aes(col = as.factor(upswing))) + theme_bw() +
    scale_color_aaas() + geom_line() +
    labs(x = "Epiweek", y = "Incident Cases (by Week)", col = "Upswing", title = "Year: 2021 in MA")
```




