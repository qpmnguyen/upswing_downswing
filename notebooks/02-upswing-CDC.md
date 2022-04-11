---
title: "Comparing surge definition to CDC"
author: "Quang Nguyen"
date: "Last compiled on 2022-04-11"
output: 
    html_document:
        keep_md: yes
        toc: true
        toc_float: true
        code_folding: hide
    github_document: default
---



# Data pre-loading and processing  


```r
library(here)
library(ggforce)
library(covidcast)
library(epiprocess)
library(tidyverse)
library(ggsci)
library(tsibble)
library(patchwork)
here::i_am(path = "notebooks/02-upswing-CDC.Rmd")
source(here("R", "utils.R"))
# (settings <- get_settings(start_date = "2020-06-01", end_date = "2022-03-01"))
```

Loading CDC Rt estimates  


```r
cdc_rt <- read_csv(file = here("data", "cdc_rt.csv")) %>% select(-1)
```

```
## New names:
## * `` -> ...1
```

```
## Rows: 48106 Columns: 15
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: ","
## chr   (3): state, type, version
## dbl  (10): ...1, median, mean, sd, lower_20, upper_20, lower_50, upper_50, l...
## lgl   (1): strat
## date  (1): date
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

```r
cdc_gr <- read_csv(file = here("data", "cdc_gr.csv")) %>% select(-1)
```

```
## New names:
## * `` -> ...1
## Rows: 42700 Columns: 15── Column specification ────────────────────────────────────────────────────────
## Delimiter: ","
## chr   (3): state, type, version
## dbl  (10): ...1, median, mean, sd, lower_20, upper_20, lower_50, upper_50, l...
## lgl   (1): strat
## date  (1): date
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

```r
start_date <- as.character(min(cdc_rt$date[cdc_rt$version == "v254"]))
end_date <- as.character(max(cdc_rt$date[cdc_rt$version == "v254"]))
h <- 4
surge_thresh <- 0.5
min_inc <- 20


head(cdc_rt)
```

```
## # A tibble: 6 × 14
##   state  date       type  median  mean    sd lower_20 upper_20 lower_50 upper_50
##   <chr>  <date>     <chr>  <dbl> <dbl> <dbl>    <dbl>    <dbl>    <dbl>    <dbl>
## 1 Alaba… 2020-08-10 esti…   0.96  0.96  0.02     0.96     0.97     0.95     0.97
## 2 Alaba… 2020-08-11 esti…   0.96  0.96  0.02     0.96     0.97     0.95     0.97
## 3 Alaba… 2020-08-12 esti…   0.96  0.96  0.02     0.96     0.97     0.95     0.97
## 4 Alaba… 2020-08-13 esti…   0.96  0.96  0.02     0.96     0.97     0.95     0.97
## 5 Alaba… 2020-08-14 esti…   0.96  0.96  0.02     0.96     0.97     0.95     0.97
## 6 Alaba… 2020-08-15 esti…   0.96  0.97  0.02     0.96     0.97     0.95     0.98
## # … with 4 more variables: lower_90 <dbl>, upper_90 <dbl>, version <chr>,
## #   strat <lgl>
```

```r
head(cdc_gr)
```

```
## # A tibble: 6 × 14
##   state date       type  median   mean    sd lower_20 upper_20 lower_50 upper_50
##   <chr> <date>     <chr>  <dbl>  <dbl> <dbl>    <dbl>    <dbl>    <dbl>    <dbl>
## 1 Alab… 2020-08-10 esti… -0.011 -0.012 0.006   -0.012   -0.01    -0.015   -0.008
## 2 Alab… 2020-08-11 esti… -0.011 -0.011 0.006   -0.011   -0.009   -0.014   -0.008
## 3 Alab… 2020-08-12 esti… -0.011 -0.011 0.006   -0.011   -0.009   -0.013   -0.007
## 4 Alab… 2020-08-13 esti… -0.011 -0.01  0.006   -0.012   -0.009   -0.013   -0.007
## 5 Alab… 2020-08-14 esti… -0.01  -0.01  0.006   -0.011   -0.009   -0.014   -0.008
## 6 Alab… 2020-08-15 esti… -0.01  -0.01  0.006   -0.011   -0.009   -0.013   -0.007
## # … with 4 more variables: lower_90 <dbl>, upper_90 <dbl>, version <chr>,
## #   strat <lgl>
```


```r
print(paste("Computing surge data from", start_date, "to", end_date))
```

```
## [1] "Computing surge data from 2021-03-05 to 2021-06-14"
```

```r
# load data from covid_cast signal
df <- covidcast_signal(data_source = "jhu-csse", 
                          signal = "confirmed_incidence_num", 
                          start_day = start_date, 
                          end_day = end_date, 
                          geo_type = "state",
                          as_of = Sys.Date()) 
```

```
## Fetched day 2021-03-05 to 2021-05-08: num_entries = 3640
```

```
## Fetched day 2021-05-09 to 2021-06-14: num_entries = 2072
```

```r
# convert to epi-df 
# after temporal operations in tsibble have to back-convert to epi_df
df <- df %>% 
        as_epi_df(geo_type = "state", time_type = "day", as_of = max(.$issue)) %>% 
        select(geo_value, time_value, cases = value) %>% as_tsibble() %>%
        index_by(epiweek = ~yearweek(., week_start = 7)) %>% 
        group_by(geo_value) %>% 
        summarise(cases = sum(cases)) %>% ungroup() %>% 
        as_tibble() %>% dplyr::rename(time_value = epiweek) %>% 
        as_epi_df(geo_type = "state", time_type = "week") %>%
        mutate(geo_value = covidcast::abbr_to_name(geo_value, ignore.case = TRUE))

head(df)
```

```
## # A tibble: 6 × 3
##   geo_value time_value cases
##   <chr>         <week> <dbl>
## 1 Alaska      2021 W09   149
## 2 Alaska      2021 W10   928
## 3 Alaska      2021 W11   897
## 4 Alaska      2021 W12  1097
## 5 Alaska      2021 W13  1234
## 6 Alaska      2021 W14  1274
```

# Surge/Downswing computation  

Compute upswing using growth_rate from epiprocess and use approximate labels. For surge classification using growth rate, we use a the weekly percent increase threshold in incidence to be 50% or more assuming that the incident case per week is at least 20. For surge classification using Rt, classified a week as increasing if the Rt had a 95% probability of being greater than or less than 1 (weekly Rt estimates were obtained by averaging daily estimates across the entire week).  


```r
df <- df %>% group_by(geo_value) %>% 
    mutate(gr = growth_rate(y = cases, method = "rel_change", h = h) * h) %>% 
    mutate(surge = case_when(
        gr >= surge_thresh & cases >= min_inc ~ "increasing",
        gr <= -surge_thresh & cases >= min_inc ~ "decreasing", 
        TRUE ~ "unclassified"
    ))

rt <- cdc_rt %>% filter(version == "v254") %>% 
    mutate(time_value = yearweek(date, week_start = 7)) %>% 
    rename("geo_value" = "state") %>% 
    group_by(geo_value, time_value) %>% 
    summarise(rt = mean(mean), rt_upper = mean(upper_90), rt_lower = mean(lower_90)) %>% 
    mutate(surge = case_when(
        rt_lower > 1 ~ "increasing", 
        rt_upper < 1 ~ "decreasing",
        TRUE ~ "unclassified"
    ))
```

```
## `summarise()` has grouped output by 'geo_value'. You can override using the
## `.groups` argument.
```

```r
head(df)
```

```
## # A tibble: 6 × 5
##   geo_value time_value cases     gr surge       
##   <chr>         <week> <dbl>  <dbl> <chr>       
## 1 Alaska      2021 W09   149  9.56  increasing  
## 2 Alaska      2021 W10   928  1.45  increasing  
## 3 Alaska      2021 W11   897  0.910 increasing  
## 4 Alaska      2021 W12  1097  0.548 increasing  
## 5 Alaska      2021 W13  1234  0.105 unclassified
## 6 Alaska      2021 W14  1274 -0.101 unclassified
```

```r
head(rt)
```

```
## # A tibble: 6 × 6
## # Groups:   geo_value [1]
##   geo_value time_value    rt rt_upper rt_lower surge       
##   <chr>         <week> <dbl>    <dbl>    <dbl> <chr>       
## 1 Alabama     2021 W09 0.892    1.03     0.747 unclassified
## 2 Alabama     2021 W10 0.881    0.993    0.755 decreasing  
## 3 Alabama     2021 W11 0.879    0.975    0.765 decreasing  
## 4 Alabama     2021 W12 0.911    0.985    0.819 decreasing  
## 5 Alabama     2021 W13 0.967    1.05     0.896 unclassified
## 6 Alabama     2021 W14 1.01     1.12     0.945 unclassified
```

```r
joint_df <- inner_join(df, rt, by=c("geo_value", "time_value")) %>% 
    pivot_longer(c(surge.x, surge.y), names_to = "source", values_to = "surge") %>%
    mutate(source = if_else(source == "surge.x", "gr", "cdc"))

cols <- c("unclassified" = "#003f5c", "increasing" = "#bc5090", "decreasing" = "#ffa600")
```


# Per state plots {.tabset}


```r
for (i in joint_df %>% pull(geo_value) %>% unique()){
    cat('##', i, '\n')
    p1 <- ggplot(joint_df %>% filter(geo_value == i) %>% filter(source == "gr"), 
                 aes(x = time_value, y = cases)) + 
        geom_line(alpha = 0.5) + 
        geom_point(aes(col = surge)) + 
        geom_rect(aes(ymin = 0, ymax = max(cases) + 10, xmin = time_value - 0.5, 
                      xmax = time_value + 0.5, fill = surge), alpha = 0.3) + 
        scale_color_manual(values = cols) + 
        scale_fill_manual(values = cols) + 
        theme_bw() +
        labs(y = "Cases", x = "Epiweek (Year-Week)", col = "Phase classification", fill = "Phase classification", 
             title = str_wrap("Phase classification via growth rate thresholded at 0.5", width = 30))
    p2 <- ggplot(joint_df %>% filter(geo_value == i) %>% filter(source == "cdc"), 
                 aes(x = time_value, y = rt)) +
        geom_point() + 
        geom_line(alpha = 0.5) + 
        geom_point(aes(col = surge)) + 
        geom_ribbon(aes(ymin = rt_lower, ymax = rt_upper), alpha = 0.4) +
        geom_rect(aes(ymin = 0.5, ymax = max(rt) + 0.3, xmin = time_value - 0.5, 
                      xmax = time_value + 0.5, fill = surge), alpha = 0.3) + 
        scale_color_manual(values = cols) + 
        scale_fill_manual(values = cols) + 
        theme_bw() +
        labs(y = "Rt", x = "Epiweek (Year-Week)", col = "Phase classification", fill = "Phase classification", 
             title = str_wrap("Phase classification via 95% CI estimates of Rt", width = 30))

    
    x <- joint_df$gr[joint_df$geo_value == i]
    y <- joint_df$rt[joint_df$geo_value == i]
    
    y <- y[-which(is.na(x))]
    x <- na.omit(x)
    
    p3 <- ggplot(joint_df %>% filter(geo_value == i), 
                 aes(x = rank(gr), y = rank(rt))) + 
        geom_point() + 
        labs(x = "rank(growth_rate)", y = "rank(Rt)", 
             title = paste("Spearman correlation:", round(cor(x,y, method = "spearman"), 3))) +
        theme_bw()
    
    layout <- "
    AAA#
    AAA#
    AAAC
    BBBC
    BBB#
    BBB#
    "
        
    print(p1 + p2 + p3 + plot_layout(guides = "collect", design = layout))
    cat("\n \n")
}
```

## Alaska 
![](02-upswing-CDC_files/figure-html/unnamed-chunk-4-1.png)<!-- -->
 
## Alabama 
![](02-upswing-CDC_files/figure-html/unnamed-chunk-4-2.png)<!-- -->
 
## Arkansas 
![](02-upswing-CDC_files/figure-html/unnamed-chunk-4-3.png)<!-- -->
 
## Arizona 
![](02-upswing-CDC_files/figure-html/unnamed-chunk-4-4.png)<!-- -->
 
## California 
![](02-upswing-CDC_files/figure-html/unnamed-chunk-4-5.png)<!-- -->
 
## Colorado 
![](02-upswing-CDC_files/figure-html/unnamed-chunk-4-6.png)<!-- -->
 
## Connecticut 
![](02-upswing-CDC_files/figure-html/unnamed-chunk-4-7.png)<!-- -->
 
## District of Columbia 
![](02-upswing-CDC_files/figure-html/unnamed-chunk-4-8.png)<!-- -->
 
## Delaware 
![](02-upswing-CDC_files/figure-html/unnamed-chunk-4-9.png)<!-- -->
 
## Florida 
![](02-upswing-CDC_files/figure-html/unnamed-chunk-4-10.png)<!-- -->
 
## Georgia 
![](02-upswing-CDC_files/figure-html/unnamed-chunk-4-11.png)<!-- -->
 
## Hawaii 
![](02-upswing-CDC_files/figure-html/unnamed-chunk-4-12.png)<!-- -->
 
## Iowa 
![](02-upswing-CDC_files/figure-html/unnamed-chunk-4-13.png)<!-- -->
 
## Idaho 
![](02-upswing-CDC_files/figure-html/unnamed-chunk-4-14.png)<!-- -->
 
## Illinois 
![](02-upswing-CDC_files/figure-html/unnamed-chunk-4-15.png)<!-- -->
 
## Indiana 
![](02-upswing-CDC_files/figure-html/unnamed-chunk-4-16.png)<!-- -->
 
## Kansas 
![](02-upswing-CDC_files/figure-html/unnamed-chunk-4-17.png)<!-- -->
 
## Kentucky 
![](02-upswing-CDC_files/figure-html/unnamed-chunk-4-18.png)<!-- -->
 
## Louisiana 
![](02-upswing-CDC_files/figure-html/unnamed-chunk-4-19.png)<!-- -->
 
## Massachusetts 
![](02-upswing-CDC_files/figure-html/unnamed-chunk-4-20.png)<!-- -->
 
## Maryland 
![](02-upswing-CDC_files/figure-html/unnamed-chunk-4-21.png)<!-- -->
 
## Maine 
![](02-upswing-CDC_files/figure-html/unnamed-chunk-4-22.png)<!-- -->
 
## Michigan 
![](02-upswing-CDC_files/figure-html/unnamed-chunk-4-23.png)<!-- -->
 
## Minnesota 
![](02-upswing-CDC_files/figure-html/unnamed-chunk-4-24.png)<!-- -->
 
## Missouri 
![](02-upswing-CDC_files/figure-html/unnamed-chunk-4-25.png)<!-- -->
 
## Mississippi 
![](02-upswing-CDC_files/figure-html/unnamed-chunk-4-26.png)<!-- -->
 
## Montana 
![](02-upswing-CDC_files/figure-html/unnamed-chunk-4-27.png)<!-- -->
 
## North Carolina 
![](02-upswing-CDC_files/figure-html/unnamed-chunk-4-28.png)<!-- -->
 
## North Dakota 
![](02-upswing-CDC_files/figure-html/unnamed-chunk-4-29.png)<!-- -->
 
## Nebraska 
![](02-upswing-CDC_files/figure-html/unnamed-chunk-4-30.png)<!-- -->
 
## New Hampshire 
![](02-upswing-CDC_files/figure-html/unnamed-chunk-4-31.png)<!-- -->
 
## New Jersey 
![](02-upswing-CDC_files/figure-html/unnamed-chunk-4-32.png)<!-- -->
 
## New Mexico 
![](02-upswing-CDC_files/figure-html/unnamed-chunk-4-33.png)<!-- -->
 
## Nevada 
![](02-upswing-CDC_files/figure-html/unnamed-chunk-4-34.png)<!-- -->
 
## New York 
![](02-upswing-CDC_files/figure-html/unnamed-chunk-4-35.png)<!-- -->
 
## Ohio 
![](02-upswing-CDC_files/figure-html/unnamed-chunk-4-36.png)<!-- -->
 
## Oklahoma 
![](02-upswing-CDC_files/figure-html/unnamed-chunk-4-37.png)<!-- -->
 
## Oregon 
![](02-upswing-CDC_files/figure-html/unnamed-chunk-4-38.png)<!-- -->
 
## Pennsylvania 
![](02-upswing-CDC_files/figure-html/unnamed-chunk-4-39.png)<!-- -->
 
## Rhode Island 
![](02-upswing-CDC_files/figure-html/unnamed-chunk-4-40.png)<!-- -->
 
## South Carolina 
![](02-upswing-CDC_files/figure-html/unnamed-chunk-4-41.png)<!-- -->
 
## South Dakota 
![](02-upswing-CDC_files/figure-html/unnamed-chunk-4-42.png)<!-- -->
 
## Tennessee 
![](02-upswing-CDC_files/figure-html/unnamed-chunk-4-43.png)<!-- -->
 
## Texas 
![](02-upswing-CDC_files/figure-html/unnamed-chunk-4-44.png)<!-- -->
 
## Utah 
![](02-upswing-CDC_files/figure-html/unnamed-chunk-4-45.png)<!-- -->
 
## Virginia 
![](02-upswing-CDC_files/figure-html/unnamed-chunk-4-46.png)<!-- -->
 
## Vermont 
![](02-upswing-CDC_files/figure-html/unnamed-chunk-4-47.png)<!-- -->
 
## Washington 
![](02-upswing-CDC_files/figure-html/unnamed-chunk-4-48.png)<!-- -->
 
## Wisconsin 
![](02-upswing-CDC_files/figure-html/unnamed-chunk-4-49.png)<!-- -->
 
## West Virginia 
![](02-upswing-CDC_files/figure-html/unnamed-chunk-4-50.png)<!-- -->
 
## Wyoming 
![](02-upswing-CDC_files/figure-html/unnamed-chunk-4-51.png)<!-- -->
 




