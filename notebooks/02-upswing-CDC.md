Comparing surge definition to CDC
================
Quang Nguyen
2022-03-21

## Data pre-loading and processing

``` r
library(here)
library(covidcast)
library(epiprocess)
library(tidyverse)
library(ggsci)
library(tsibble)
library(patchwork)
here::i_am(path = "notebooks/02-upswing-CDC.Rmd")
source(here("R", "utils.R"))

start_date <- "2020-09-01"
end_date <- "2021-04-02"
h <- 4
surge_thresh <- 0.5
min_inc <- 20
#(settings <- get_settings(start_date = "2020-06-01", end_date = "2022-03-01"))
```

``` r
print(paste("Computing surge data for state of MA from", start_date, "to", end_date))
```

    ## [1] "Computing surge data for state of MA from 2020-09-01 to 2021-04-02"

``` r
df <- covidcast_signal(data_source = "jhu-csse", 
                          signal = "confirmed_incidence_num", 
                          start_day = start_date, 
                          end_day = end_date, 
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
```

    ## Fetched day 2020-09-01 to 2021-04-02: num_entries = 214

## Surge/Downswing computation

Compute upswing using growth_rate from epiprocess and use approximate
labels

``` r
df <- df %>% mutate(gr = growth_rate(y = cases, method = "rel_change", h = h) * h) %>% 
    mutate(base_st = case_when(
        gr >= surge_thresh & cases >= min_inc ~ "increasing",
        gr <= -surge_thresh & cases >= min_inc ~ "decreasing", 
        TRUE ~ "unclassified"
    ))


complete_lab <- df %>% filter(base_st != "unclassified") %>% as_tsibble() %>% tsibble::count_gaps() %>% 
    rowwise() %>%
    mutate(label = case_when(
        df %>% filter(time_value == .from - 1) %>% 
            pull(base_st) == "increasing" & 
        df %>% filter(time_value == .to + 1) %>% 
            pull(base_st) == "decreasing" ~ "peaks",
        
        df %>% filter(time_value == .from - 1) %>% 
            pull(base_st) == "decreasing" & 
        df %>% filter(time_value == .to + 1) %>% 
            pull(base_st) == "increasing" ~ "nadirs",
        
        df %>% filter(time_value == .from - 1) %>% 
            pull(base_st) == "decreasing" & 
        df %>% filter(time_value == .to + 1) %>% 
            pull(base_st) == "decreasing" ~ "decreasing",
        
                df %>% filter(time_value == .from - 1) %>% 
            pull(base_st) == "increasing" & 
        df %>% filter(time_value == .to + 1) %>% 
            pull(base_st) == "increasing" ~ "increasing"))


df <- df %>% mutate(phase = map2_chr(time_value, base_st, ~{
    if (.y != "unclassified"){
        return(.y)
    } else {
        lab <- complete_lab %>% filter(.x >= .from & .x <= .to) %>% 
            pull(label)
        if(length(lab) == 0){
            return(NA_character_)
        } else {
            return(lab)
        }
    }
}))
```

Plotting only increasing/decreasing/unclassified based on point
classification

``` r
ggplot(df, aes(x = time_value, y = cases)) + 
    geom_point(aes(col = base_st), size = 2) + geom_line() +  
    geom_rect(aes(ymin = 0, ymax = max(cases), xmin = time_value - 0.5, 
                  xmax = time_value + 0.5, fill = base_st), alpha = 0.5) +
    theme_bw()
```

![](02-upswing-CDC_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

Plotting with correct phase assignment based on CDC document

``` r
ggplot(df, aes(x = time_value, y = cases)) + 
    geom_point(aes(col = phase), size = 2) + geom_line() +  
    geom_rect(aes(ymin = 0, ymax = max(cases), xmin = time_value - 0.5, 
                  xmax = time_value + 0.5, fill = phase), alpha = 0.5) +
    theme_bw()
```

![](02-upswing-CDC_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

Reference image from CDC

``` r
knitr::include_graphics(path = here("imgs", "cdc_ref_ma.png"))
```

![](/Users/quangnguyen/projects/upswing_downswing/imgs/cdc_ref_ma.png)<!-- -->
