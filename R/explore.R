combined <- combined %>% ungroup()

test_state <- c("ny", "ca", "ma")
test_model <- models[c(1,3)]
subset <- combined %>% filter(model %in% test_model, geo_value %in% test_state)

test <- subset %>% group_by(model, geo_value) %>% epi_slide(~{
    query_dates <- .x %>% pull(time_value) %>% unique()
    req_len <- h * 2
    if (length(query_dates) != req_len){
        out <- NA_real_
    }
    t_date <- query_dates[1:(req_len - h)]
    ref <- tail(t_date, n = 1)
    p_date <- query_dates[(req_len - h + 1):req_len]
    f_date <- .x %>%
        filter(time_value == head(p_date, n = 1) & horizon == 1) %>%
        pull(forecast_date) %>% unique()

    if (length(f_date) == 0){
        out <- NA_real_
    } else {
        pred <- .x %>% filter(forecast_date == f_date) %>%
            dplyr::pull(pred)

        if (length(pred) < h){
            # this exception is when the model does not forecast $h$ weeks in advance
            # (e.g. h = 4 but horizon only extends to 3)
            out <- NA_real_
        } else {
            # get true values
            obs <- .x %>% filter(time_value %in% t_date) %>%
                select(time_value, obs) %>%
                distinct() %>% pull(obs)

            # put everything in an epi_df for posterity
            new_df <- tibble(
                time_value = c(t_date, p_date),
                value = c(obs, pred),
                geo_value = "placeholder"
            ) %>% as_epi_df()

            out <- new_df %>%
                mutate(gr_pred = growth_rate(y = value, h = h, method = "rel_change") * h) %>%
                filter(time_value == ref) %>% pull(gr_pred)
        }
    }
}, n = 2 * 7 * h, align = "center", new_col_name = "pred_gr")

test2 <- test %>% ungroup()  %>%
    select(-c(forecast_date, horizon, pred)) %>% distinct()

test2 <- test2 %>% mutate(surge_pred = case_when(
    pred_gr >= surge_thresh & obs >= min_inc ~ TRUE,
    is.na(pred_gr) ~ NA,
    TRUE ~ FALSE
)) %>% group_by(geo_value, model) %>%
    epi_slide(~{
       bef <- .x$surge_pred[1]
       focal <- .x$surge_pred[2]
       if (is.na(focal) | is.na(bef)){
           out <- NA
       } else {
           if (bef == FALSE & focal == TRUE){
               out <- TRUE
           } else {
               out <- FALSE
           }
       }
       return(out)
    }, n = 2 * 7 * 1, align = "right", new_col_name = "upswing_pred") %>%
    ungroup()


ggplot(test2 %>% filter(model == test_model[1], geo_value == "ca"),
       aes(x = time_value, y = obs)) + geom_point(aes(col = upswing_pred))
test2


ggplot(gradient_classif %>% filter(model == models[9], geo_value == "ca"),
       aes(x = time_value, y = pred_gr)) + geom_point(aes(col = surge_pred))

models[9]
