---
title: "Comparing surge definition to CDC"
author: "Quang Nguyen"
date: "Last compiled on 2022-06-28"
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
library(tidyverse)
library(ggsci)
library(tsibble)
library(patchwork)
here::i_am(path = "notebooks/upswing-CDC.Rmd")
source(here("R", "utils.R"))
# (settings <- get_settings(start_date = "2020-06-01", end_date = "2022-03-01"))
```

Loading CDC Rt estimates


```r
cdc_rt <- read_csv(file = here("data", "cdc_rt.csv")) %>% select(-1)
```

```
## New names:
## Rows: 48106 Columns: 15
## ── Column specification
## ────────────────────────────────────────────────────────────────────────────────────────────────────────── Delimiter: "," chr
## (3): state, type, version dbl (10): ...1, median, mean, sd, lower_20, upper_20, lower_50, upper_50, lower_90, upper_90 lgl (1):
## strat date (1): date
## ℹ Use `spec()` to retrieve the full column specification for this data. ℹ Specify the column types or set `show_col_types =
## FALSE` to quiet this message.
## • `` -> `...1`
```

```r
cdc_gr <- read_csv(file = here("data", "cdc_gr.csv")) %>% select(-1)
```

```
## New names:
## Rows: 42700 Columns: 15
## ── Column specification
## ────────────────────────────────────────────────────────────────────────────────────────────────────────── Delimiter: "," chr
## (3): state, type, version dbl (10): ...1, median, mean, sd, lower_20, upper_20, lower_50, upper_50, lower_90, upper_90 lgl (1):
## strat date (1): date
## ℹ Use `spec()` to retrieve the full column specification for this data. ℹ Specify the column types or set `show_col_types =
## FALSE` to quiet this message.
## • `` -> `...1`
```

```r
start_date <- as.character(min(cdc_rt$date))
end_date <- as.character(max(cdc_rt$date))
h <- 4
surge_thresh <- 0.5
min_inc <- 20


head(cdc_rt)
```

```
## # A tibble: 6 × 14
##   state   date       type     median  mean    sd lower_20 upper_20 lower_50 upper_50 lower_90 upper_90 version strat
##   <chr>   <date>     <chr>     <dbl> <dbl> <dbl>    <dbl>    <dbl>    <dbl>    <dbl>    <dbl>    <dbl> <chr>   <lgl>
## 1 Alabama 2020-08-10 estimate   0.96  0.96  0.02     0.96     0.97     0.95     0.97     0.93     0.99 v001    NA   
## 2 Alabama 2020-08-11 estimate   0.96  0.96  0.02     0.96     0.97     0.95     0.97     0.93     0.99 v001    NA   
## 3 Alabama 2020-08-12 estimate   0.96  0.96  0.02     0.96     0.97     0.95     0.97     0.92     0.99 v001    NA   
## 4 Alabama 2020-08-13 estimate   0.96  0.96  0.02     0.96     0.97     0.95     0.97     0.93     1    v001    NA   
## 5 Alabama 2020-08-14 estimate   0.96  0.96  0.02     0.96     0.97     0.95     0.97     0.93     1    v001    NA   
## 6 Alabama 2020-08-15 estimate   0.96  0.97  0.02     0.96     0.97     0.95     0.98     0.93     1    v001    NA
```

```r
head(cdc_gr)
```

```
## # A tibble: 6 × 14
##   state   date       type     median   mean    sd lower_20 upper_20 lower_50 upper_50 lower_90 upper_90 version strat
##   <chr>   <date>     <chr>     <dbl>  <dbl> <dbl>    <dbl>    <dbl>    <dbl>    <dbl>    <dbl>    <dbl> <chr>   <lgl>
## 1 Alabama 2020-08-10 estimate -0.011 -0.012 0.006   -0.012   -0.01    -0.015   -0.008   -0.021   -0.004 v001    NA   
## 2 Alabama 2020-08-11 estimate -0.011 -0.011 0.006   -0.011   -0.009   -0.014   -0.008   -0.02    -0.003 v001    NA   
## 3 Alabama 2020-08-12 estimate -0.011 -0.011 0.006   -0.011   -0.009   -0.013   -0.007   -0.02    -0.002 v001    NA   
## 4 Alabama 2020-08-13 estimate -0.011 -0.01  0.006   -0.012   -0.009   -0.013   -0.007   -0.02    -0.001 v001    NA   
## 5 Alabama 2020-08-14 estimate -0.01  -0.01  0.006   -0.011   -0.009   -0.014   -0.008   -0.021   -0.002 v001    NA   
## 6 Alabama 2020-08-15 estimate -0.01  -0.01  0.006   -0.011   -0.009   -0.013   -0.007   -0.02    -0.001 v001    NA
```

Obtain case-count data from JHU using start and end dates as defined by the CDC file which is from 2020-08-10 and 2021-06-14.


```r
print(paste("Computing surge data from", start_date, "to", end_date))
```

```
## [1] "Computing surge data from 2020-08-10 to 2021-06-14"
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
## Fetched day 2020-08-10 to 2020-10-13: num_entries = 3640
```

```
## Fetched day 2020-10-14 to 2020-12-17: num_entries = 3640
```

```
## Fetched day 2020-12-18 to 2021-02-20: num_entries = 3640
```

```
## Fetched day 2021-02-21 to 2021-04-26: num_entries = 3640
```

```
## Fetched day 2021-04-27 to 2021-06-14: num_entries = 2744
```

Here, we convert `covidcast_signal` data into `epi_df` and `tsibble`. We have to back-convert to `epi_df` to work with `growth_rate` function.


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
## 1 Alaska      2020 W33   480
## 2 Alaska      2020 W34   547
## 3 Alaska      2020 W35   533
## 4 Alaska      2020 W36   496
## 5 Alaska      2020 W37   546
## 6 Alaska      2020 W38   552
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

rt <- cdc_rt %>% 
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
## `summarise()` has grouped output by 'geo_value'. You can override using the `.groups` argument.
```

```r
head(df)
```

```
## # A tibble: 6 × 5
##   geo_value time_value cases     gr surge       
##   <chr>         <week> <dbl>  <dbl> <chr>       
## 1 Alaska      2020 W33   480 0.168  unclassified
## 2 Alaska      2020 W34   547 0.0474 unclassified
## 3 Alaska      2020 W35   533 0.120  unclassified
## 4 Alaska      2020 W36   496 0.304  unclassified
## 5 Alaska      2020 W37   546 0.561  increasing  
## 6 Alaska      2020 W38   552 0.938  increasing
```

```r
head(rt)
```

```
## # A tibble: 6 × 6
## # Groups:   geo_value [1]
##   geo_value time_value    rt rt_upper rt_lower surge       
##   <chr>         <week> <dbl>    <dbl>    <dbl> <chr>       
## 1 Alabama     2020 W33 0.962    0.995    0.928 decreasing  
## 2 Alabama     2020 W34 0.976    1.01     0.933 unclassified
## 3 Alabama     2020 W35 0.97     1.02     0.919 unclassified
## 4 Alabama     2020 W36 0.962    1.01     0.916 unclassified
## 5 Alabama     2020 W37 0.976    1.03     0.902 unclassified
## 6 Alabama     2020 W38 0.991    1.05     0.928 unclassified
```

```r
joint_df <- inner_join(df, rt, by=c("geo_value", "time_value")) %>% 
    pivot_longer(c(surge.x, surge.y), names_to = "source", values_to = "surge") %>%
    mutate(source = if_else(source == "surge.x", "gr", "cdc"))

cols <- c("unclassified" = "#003f5c", "increasing" = "#bc5090", "decreasing" = "#ffa600")
```

# Per state plots of increasing/decreasing/unclassified {.tabset}

These plots represent per-week classification based on the schema above. As such, some weeks are considered increasing/decreasing while others are "unclassified"


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
![](upswing-CDC_files/figure-html/unnamed-chunk-5-1.png)<!-- -->
 
## Alabama 
![](upswing-CDC_files/figure-html/unnamed-chunk-5-2.png)<!-- -->
 
## Arkansas 
![](upswing-CDC_files/figure-html/unnamed-chunk-5-3.png)<!-- -->
 
## Arizona 
![](upswing-CDC_files/figure-html/unnamed-chunk-5-4.png)<!-- -->
 
## California 
![](upswing-CDC_files/figure-html/unnamed-chunk-5-5.png)<!-- -->
 
## Colorado 
![](upswing-CDC_files/figure-html/unnamed-chunk-5-6.png)<!-- -->
 
## Connecticut 
![](upswing-CDC_files/figure-html/unnamed-chunk-5-7.png)<!-- -->
 
## District of Columbia 
![](upswing-CDC_files/figure-html/unnamed-chunk-5-8.png)<!-- -->
 
## Delaware 
![](upswing-CDC_files/figure-html/unnamed-chunk-5-9.png)<!-- -->
 
## Florida 
![](upswing-CDC_files/figure-html/unnamed-chunk-5-10.png)<!-- -->
 
## Georgia 
![](upswing-CDC_files/figure-html/unnamed-chunk-5-11.png)<!-- -->
 
## Hawaii 
![](upswing-CDC_files/figure-html/unnamed-chunk-5-12.png)<!-- -->
 
## Iowa 
![](upswing-CDC_files/figure-html/unnamed-chunk-5-13.png)<!-- -->
 
## Idaho 
![](upswing-CDC_files/figure-html/unnamed-chunk-5-14.png)<!-- -->
 
## Illinois 
![](upswing-CDC_files/figure-html/unnamed-chunk-5-15.png)<!-- -->
 
## Indiana 
![](upswing-CDC_files/figure-html/unnamed-chunk-5-16.png)<!-- -->
 
## Kansas 
![](upswing-CDC_files/figure-html/unnamed-chunk-5-17.png)<!-- -->
 
## Kentucky 
![](upswing-CDC_files/figure-html/unnamed-chunk-5-18.png)<!-- -->
 
## Louisiana 
![](upswing-CDC_files/figure-html/unnamed-chunk-5-19.png)<!-- -->
 
## Massachusetts 
![](upswing-CDC_files/figure-html/unnamed-chunk-5-20.png)<!-- -->
 
## Maryland 
![](upswing-CDC_files/figure-html/unnamed-chunk-5-21.png)<!-- -->
 
## Maine 
![](upswing-CDC_files/figure-html/unnamed-chunk-5-22.png)<!-- -->
 
## Michigan 
![](upswing-CDC_files/figure-html/unnamed-chunk-5-23.png)<!-- -->
 
## Minnesota 
![](upswing-CDC_files/figure-html/unnamed-chunk-5-24.png)<!-- -->
 
## Missouri 
![](upswing-CDC_files/figure-html/unnamed-chunk-5-25.png)<!-- -->
 
## Mississippi 
![](upswing-CDC_files/figure-html/unnamed-chunk-5-26.png)<!-- -->
 
## Montana 
![](upswing-CDC_files/figure-html/unnamed-chunk-5-27.png)<!-- -->
 
## North Carolina 
![](upswing-CDC_files/figure-html/unnamed-chunk-5-28.png)<!-- -->
 
## North Dakota 
![](upswing-CDC_files/figure-html/unnamed-chunk-5-29.png)<!-- -->
 
## Nebraska 
![](upswing-CDC_files/figure-html/unnamed-chunk-5-30.png)<!-- -->
 
## New Hampshire 
![](upswing-CDC_files/figure-html/unnamed-chunk-5-31.png)<!-- -->
 
## New Jersey 
![](upswing-CDC_files/figure-html/unnamed-chunk-5-32.png)<!-- -->
 
## New Mexico 
![](upswing-CDC_files/figure-html/unnamed-chunk-5-33.png)<!-- -->
 
## Nevada 
![](upswing-CDC_files/figure-html/unnamed-chunk-5-34.png)<!-- -->
 
## New York 
![](upswing-CDC_files/figure-html/unnamed-chunk-5-35.png)<!-- -->
 
## Ohio 
![](upswing-CDC_files/figure-html/unnamed-chunk-5-36.png)<!-- -->
 
## Oklahoma 
![](upswing-CDC_files/figure-html/unnamed-chunk-5-37.png)<!-- -->
 
## Oregon 
![](upswing-CDC_files/figure-html/unnamed-chunk-5-38.png)<!-- -->
 
## Pennsylvania 
![](upswing-CDC_files/figure-html/unnamed-chunk-5-39.png)<!-- -->
 
## Rhode Island 
![](upswing-CDC_files/figure-html/unnamed-chunk-5-40.png)<!-- -->
 
## South Carolina 
![](upswing-CDC_files/figure-html/unnamed-chunk-5-41.png)<!-- -->
 
## South Dakota 
![](upswing-CDC_files/figure-html/unnamed-chunk-5-42.png)<!-- -->
 
## Tennessee 
![](upswing-CDC_files/figure-html/unnamed-chunk-5-43.png)<!-- -->
 
## Texas 
![](upswing-CDC_files/figure-html/unnamed-chunk-5-44.png)<!-- -->
 
## Utah 
![](upswing-CDC_files/figure-html/unnamed-chunk-5-45.png)<!-- -->
 
## Virginia 
![](upswing-CDC_files/figure-html/unnamed-chunk-5-46.png)<!-- -->
 
## Vermont 
![](upswing-CDC_files/figure-html/unnamed-chunk-5-47.png)<!-- -->
 
## Washington 
![](upswing-CDC_files/figure-html/unnamed-chunk-5-48.png)<!-- -->
 
## Wisconsin 
![](upswing-CDC_files/figure-html/unnamed-chunk-5-49.png)<!-- -->
 
## West Virginia 
![](upswing-CDC_files/figure-html/unnamed-chunk-5-50.png)<!-- -->
 
## Wyoming 
![](upswing-CDC_files/figure-html/unnamed-chunk-5-51.png)<!-- -->
 

# Per-state plots with full phase classification 

These plots represent per-state classifications with include full phase classifications. Specifically:

- For weeks between two increasing phases, we classified them as increasing and for weeks between two decreasing periods, we classified them as decreasing.  
- Weeks between increasing and decreasing phases were classified as peaks, whereas nadirs were defined as periods between decreasing and increasing phases.  
- If peaks or nadirs only contained an interval of 1 week, we included the week before and after in the peak/nadir classification.
- Periods at the beginning or the end of the analysis period were classified as a continuation of whichever phase preceded or followed them.    



```r
# phase classification was done by counting all "unclassified" as missing 
# we then use count_gaps in order to define all complete phases and then assign labels to them  
# and then loop through the main data frame and assign labels based on gaps (if surge is unclassified)
# for some reason using map functions doesn't work when trying to filter time_value >= .from and <= .to

get_lab <- function(df, geo, time){
    return(df %>% filter(geo_value == geo) %>% 
               filter(time_value == time) %>%
               pull(surge))
}

phase_clasification <- function(src){
    data <- joint_df %>% filter(source == src)
    gaps <- data %>% filter(surge != "unclassified") %>% 
        as_tsibble(key = geo_value) %>% tsibble::count_gaps(.full = TRUE) 
    
    
    gaps <- gaps %>% mutate(label = pmap_chr(gaps, function(geo_value, .from, .to, ...){
        before <- get_lab(data, geo_value, .from - 1)
        after <- get_lab(data, geo_value, .to + 1)
        if (length(before) == 0){
            val <- after
        } else if (length(after) == 0){
            val <- before
        } else {
            val <- case_when(
                before == "increasing" & after == "decreasing" ~ "peaks",
                before == "decreasing" & after == "increasing" ~ "nadirs",
                before == "increasing" & after == "increasing" ~ "increasing",
                before == "decreasing" & after == "decreasing" ~ "decreasing",
                TRUE ~ NA_character_
            )
        }
        return(val)
    }))
    
    phases <- vector(length = nrow(data))
    for (i in seq_along(phases)){
        if (data[i,]$surge != "unclassified"){
            lab <- data[i,]$surge
        } else {
             t <- data[i,]$time_value
             g <- data[i,]$geo_value
             lab <- gaps %>% filter(geo_value == g) %>% 
                 filter(.from <= t) %>% 
                 filter(.to >= t) %>% pull(label) %>% as.vector()
             if (length(lab) == 0){
                 lab <- NA_character_
             }
        }
        phases[i] <- lab
    }
    data$phases <- phases
    return(data)
}
```

Applying phase classification to CDC Rt and case count data sets  


```r
gr_phase <- phase_clasification("gr")
```

```
## Using `time_value` as index variable.
```

```r
rt_phase <- phase_clasification("cdc")
```

```
## Using `time_value` as index variable.
```

```r
cols <- c("peaks" = "#003f5c", "increasing" = "#bc5090", "decreasing" = "#ffa600", "nadirs" = "#58508d")
```

## Per-state plots {.tabset}


```r
for (i in joint_df %>% pull(geo_value) %>% unique()){
    cat('###', i, '\n')
    p1 <- ggplot(gr_phase %>% filter(geo_value == i), 
                 aes(x = time_value, y = cases)) + 
        geom_line(alpha = 0.5) + 
        geom_point(aes(col = phases)) + 
        geom_rect(aes(ymin = 0, ymax = max(cases) + 10, xmin = time_value - 0.5, 
                      xmax = time_value + 0.5, fill = phases), alpha = 0.3) + 
        scale_color_manual(values = cols) + 
        scale_fill_manual(values = cols) + 
        theme_bw() +
        labs(y = "Cases", x = "Epiweek (Year-Week)", col = "Phase classification", fill = "Phase classification", 
             title = str_wrap("Phase classification via growth rate thresholded at 0.5", width = 30))
    p2 <- ggplot(rt_phase %>% filter(geo_value == i), 
                 aes(x = time_value, y = rt)) +
        geom_point() + 
        geom_line(alpha = 0.5) + 
        geom_point(aes(col = phases)) + 
        geom_ribbon(aes(ymin = rt_lower, ymax = rt_upper), alpha = 0.4) +
        geom_rect(aes(ymin = 0.5, ymax = max(rt) + 0.3, xmin = time_value - 0.5, 
                      xmax = time_value + 0.5, fill = phases), alpha = 0.3) + 
        scale_color_manual(values = cols) + 
        scale_fill_manual(values = cols) + 
        theme_bw() +
        labs(y = "Rt", x = "Epiweek (Year-Week)", col = "Phase classification", fill = "Phase classification", 
             title = str_wrap("Phase classification via 95% CI estimates of Rt", width = 30))
    print(p1 / p2 + plot_layout(guides = "collect"))
    cat("\n\n")
}
```

### Alaska 
![](upswing-CDC_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

### Alabama 
![](upswing-CDC_files/figure-html/unnamed-chunk-8-2.png)<!-- -->

### Arkansas 
![](upswing-CDC_files/figure-html/unnamed-chunk-8-3.png)<!-- -->

### Arizona 
![](upswing-CDC_files/figure-html/unnamed-chunk-8-4.png)<!-- -->

### California 
![](upswing-CDC_files/figure-html/unnamed-chunk-8-5.png)<!-- -->

### Colorado 
![](upswing-CDC_files/figure-html/unnamed-chunk-8-6.png)<!-- -->

### Connecticut 
![](upswing-CDC_files/figure-html/unnamed-chunk-8-7.png)<!-- -->

### District of Columbia 
![](upswing-CDC_files/figure-html/unnamed-chunk-8-8.png)<!-- -->

### Delaware 
![](upswing-CDC_files/figure-html/unnamed-chunk-8-9.png)<!-- -->

### Florida 
![](upswing-CDC_files/figure-html/unnamed-chunk-8-10.png)<!-- -->

### Georgia 
![](upswing-CDC_files/figure-html/unnamed-chunk-8-11.png)<!-- -->

### Hawaii 
![](upswing-CDC_files/figure-html/unnamed-chunk-8-12.png)<!-- -->

### Iowa 
![](upswing-CDC_files/figure-html/unnamed-chunk-8-13.png)<!-- -->

### Idaho 
![](upswing-CDC_files/figure-html/unnamed-chunk-8-14.png)<!-- -->

### Illinois 
![](upswing-CDC_files/figure-html/unnamed-chunk-8-15.png)<!-- -->

### Indiana 
![](upswing-CDC_files/figure-html/unnamed-chunk-8-16.png)<!-- -->

### Kansas 
![](upswing-CDC_files/figure-html/unnamed-chunk-8-17.png)<!-- -->

### Kentucky 
![](upswing-CDC_files/figure-html/unnamed-chunk-8-18.png)<!-- -->

### Louisiana 
![](upswing-CDC_files/figure-html/unnamed-chunk-8-19.png)<!-- -->

### Massachusetts 
![](upswing-CDC_files/figure-html/unnamed-chunk-8-20.png)<!-- -->

### Maryland 
![](upswing-CDC_files/figure-html/unnamed-chunk-8-21.png)<!-- -->

### Maine 
![](upswing-CDC_files/figure-html/unnamed-chunk-8-22.png)<!-- -->

### Michigan 
![](upswing-CDC_files/figure-html/unnamed-chunk-8-23.png)<!-- -->

### Minnesota 
![](upswing-CDC_files/figure-html/unnamed-chunk-8-24.png)<!-- -->

### Missouri 
![](upswing-CDC_files/figure-html/unnamed-chunk-8-25.png)<!-- -->

### Mississippi 
![](upswing-CDC_files/figure-html/unnamed-chunk-8-26.png)<!-- -->

### Montana 
![](upswing-CDC_files/figure-html/unnamed-chunk-8-27.png)<!-- -->

### North Carolina 
![](upswing-CDC_files/figure-html/unnamed-chunk-8-28.png)<!-- -->

### North Dakota 
![](upswing-CDC_files/figure-html/unnamed-chunk-8-29.png)<!-- -->

### Nebraska 
![](upswing-CDC_files/figure-html/unnamed-chunk-8-30.png)<!-- -->

### New Hampshire 
![](upswing-CDC_files/figure-html/unnamed-chunk-8-31.png)<!-- -->

### New Jersey 
![](upswing-CDC_files/figure-html/unnamed-chunk-8-32.png)<!-- -->

### New Mexico 
![](upswing-CDC_files/figure-html/unnamed-chunk-8-33.png)<!-- -->

### Nevada 
![](upswing-CDC_files/figure-html/unnamed-chunk-8-34.png)<!-- -->

### New York 
![](upswing-CDC_files/figure-html/unnamed-chunk-8-35.png)<!-- -->

### Ohio 
![](upswing-CDC_files/figure-html/unnamed-chunk-8-36.png)<!-- -->

### Oklahoma 
![](upswing-CDC_files/figure-html/unnamed-chunk-8-37.png)<!-- -->

### Oregon 
![](upswing-CDC_files/figure-html/unnamed-chunk-8-38.png)<!-- -->

### Pennsylvania 
![](upswing-CDC_files/figure-html/unnamed-chunk-8-39.png)<!-- -->

### Rhode Island 
![](upswing-CDC_files/figure-html/unnamed-chunk-8-40.png)<!-- -->

### South Carolina 
![](upswing-CDC_files/figure-html/unnamed-chunk-8-41.png)<!-- -->

### South Dakota 
![](upswing-CDC_files/figure-html/unnamed-chunk-8-42.png)<!-- -->

### Tennessee 
![](upswing-CDC_files/figure-html/unnamed-chunk-8-43.png)<!-- -->

### Texas 
![](upswing-CDC_files/figure-html/unnamed-chunk-8-44.png)<!-- -->

### Utah 
![](upswing-CDC_files/figure-html/unnamed-chunk-8-45.png)<!-- -->

### Virginia 
![](upswing-CDC_files/figure-html/unnamed-chunk-8-46.png)<!-- -->

### Vermont 
![](upswing-CDC_files/figure-html/unnamed-chunk-8-47.png)<!-- -->

### Washington 
![](upswing-CDC_files/figure-html/unnamed-chunk-8-48.png)<!-- -->

### Wisconsin 
![](upswing-CDC_files/figure-html/unnamed-chunk-8-49.png)<!-- -->

### West Virginia 
![](upswing-CDC_files/figure-html/unnamed-chunk-8-50.png)<!-- -->

### Wyoming 
![](upswing-CDC_files/figure-html/unnamed-chunk-8-51.png)<!-- -->





