### Multinomial projections
library(MASS)
library(tidyverse)
library(nnet)
library(cowplot)
library(effects)
library(rjson)
library(arrow)
theme_set(theme_cowplot())


source('R/mnl_pred_ova2_revised.R')

# Get location data -------------------------------------------------------
if(!file.exists('raw-data/locations.csv')){
  locs <- read_csv('https://raw.githubusercontent.com/cdcepi/FluSight-forecast-hub/refs/heads/main/auxiliary-data/locations.csv') |> 
    select(abbreviation:population)
  write_csv(locs, 'raw-data/locations.csv')
} else{
  locs <- read_csv('raw-data/locations.csv')
}

# Download, read-in, and save variant count file --------------------------
dl_url <- "https://data.nextstrain.org/files/workflows/forecasts-ncov/open/nextstrain_clades/usa.tsv.gz"
tmp_file <- tempfile()
download.file(dl_url, tmp_file)
variantcts <- read_tsv(gzfile(tmp_file))
dir.create('raw-data/historic_counts', showWarnings = FALSE)
write_csv(variantcts, paste0('raw-data/historic_counts/', Sys.Date(), '_variantcts.csv'))
rm(tmp_file)


# variantcts <- read_csv(paste0('raw-data/historic_counts/', Sys.Date(), '_variantcts.csv')) ## Use if running manually

# Clean the variant locations in data -------------------------------------
variantcts |> 
  mutate(location = ifelse(location == 'Washington DC', 'District of Columbia', location)) |> 
  mutate(location = ifelse(location == 'Deleware', 'Delaware', location)) |> 
  mutate(location = ifelse(location == 'Louisana', 'Louisiana', location)) |> 
  mutate(location = ifelse(grepl('Louisiana', location), 'Louisiana', location)) |> 
  mutate(location = ifelse(location == 'Washington DC', 'District of Columbia', location)) |> 
  filter(location %in% locs$location_name) -> variantcts


# Plot the variant counts and proportions ---------------------------------
## Total counts over time
variantcts |> group_by(date) |> 
  summarize(cts = sum(sequences)) |> 
  ggplot(aes(date, cts)) + geom_line()

## Proportion of variants
variantcts |> 
  group_by(location, date) |> 
  mutate(prop = sequences/sum(sequences)) |> 
  filter(location == 'California') |> 
  ggplot(aes(date, prop, fill = clade)) +
  geom_col(width = 1)



# Read in clades for this week --------------------------------------------
todays_date <- Sys.Date()
# todays_date <- ymd('2025-07-16') ## Only use this if you want to test out script not on wednesday

json_location <- paste0("https://raw.githubusercontent.com/reichlab/variant-nowcast-hub/refs/heads/main/auxiliary-data/modeled-clades/",
                        todays_date,
                        ".json")
clades <- fromJSON(paste(readLines(json_location), collapse=""))
fcast_clades <- clades$clades

# fcast_clades_original <- fcast_clades
# fcast_clades <- fcast_clades[2:10]

# Do final cleaning for the analysis --------------------------------------
## Currently removing the most recent 30 days, because it's so noisy and messes with predictions (something to test later though)
## Currently only using data back 180 days (something to test later on)
## Parameters used to control the amount of data used and the dates used for model predictions
fcast_days_from_today <- 10
hindcast_days_from_today <- 31
days_to_look_back <- 120
recent_days_to_drop <- 0

## Supposed to forecast 42 total dates, 31 back and 10 forward
forecast_dates <- seq.Date(ymd(todays_date)-days(hindcast_days_from_today), 
                           ymd(todays_date)+days(fcast_days_from_today), by = 'day')

## Create tibble for controlling data and dates later on
tibble(date = seq.Date(ymd(todays_date)-days(days_to_look_back), 
                       ymd(todays_date)+days(fcast_days_from_today), by = 'day'),
       time = seq_along(date)) |> 
  mutate(forecast_date = date%in% forecast_dates) -> time_control_df

## Create tibble to expand to all dates, clades, and locations
locs |> 
  select(location = location_name) |> 
  expand_grid(tibble(clade = fcast_clades),
              tibble(date = unique(c(variantcts$date, time_control_df$date)))) |> 
  mutate(data_dates = date %in%variantcts$date) -> all_possible_data

variantcts |> 
  mutate(clade = ifelse(clade %in% fcast_clades, clade, 'other')) |>  ## If the clade is in the ones we need to forecast, keep, otherwise change to "Other"
  group_by(location, clade, date) |> 
  summarize(sequences = sum(sequences)) -> full_fcast_df

## This removes all dates for each location where there isn't at least one sequence.
all_possible_data |> 
  filter(data_dates) |> 
  left_join(full_fcast_df, by = c('location', 'clade', 'date')) |> 
  group_by(location, date) |> 
  filter(!all(is.na(sequences))) |> 
  mutate(sequences = ifelse(is.na(sequences), 0, sequences)) |> 
  ungroup() -> full_fcast_df

## This code was used to add in sequences that were missing to help with rare variants, but
## It's too arbitrary of an approach so no longer used
# full_fcast_df |> 
#   group_by(location, clade) |> 
#   summarize(n = sum(sequences)) |> 
#   filter(n ==0) |> 
#   ungroup() -> missing_clade_locations
# 
# for(i in 1:nrow(missing_clade_locations)){
#   temp <- missing_clade_locations |> slice(i)
#   
#   full_fcast_df |> 
#     bind_rows(tibble(location = temp$location,
#                      clade = fcast_clades[which(!(fcast_clades %in% c('recombinant', 'other')))],
#                      date = todays_date - days(1),
#                      data_dates = TRUE,
#                      sequences = 1)) -> full_fcast_df
# }


## Now make the actually tibble for forecasting
full_fcast_df |> 
  inner_join(time_control_df, by = 'date') |> 
  uncount(sequences)  -> mod_fcast_df

mod_fcast_df |> 
  inner_join(mod_fcast_df |> 
               count(location) |> 
               filter(n>=100) |> select(location), by = 'location') -> mod_fcast_df

time_center <- mean(mod_fcast_df$time)
time_scale <- sd(mod_fcast_df$time)
mod_fcast_df <- mod_fcast_df |>
  mutate(time_sc = (time - time_center) / time_scale)


# Fit two-stage model and estimate uncertainty ----------------------------
## Stage A: national clade trend over time
## Stage B: regional clade offsets (no regional time slope)

mod_fcast_df <- mod_fcast_df |>
  group_by(location) |>
  mutate(location_total_sequences = n()) |>
  ungroup() |>
  filter(location_total_sequences >= 5) |>
  select(-location_total_sequences)

present_clades <- fcast_clades[fcast_clades %in% unique(mod_fcast_df$clade)]
if (length(present_clades) < 3) {
  stop("Need at least 3 clades with observed counts to fit the two-stage multinomial model.")
}

mod_fcast_df <- mod_fcast_df |>
  filter(clade %in% present_clades) |>
  mutate(clade = factor(clade, levels = present_clades))

stage_a_df <- mod_fcast_df |>
  count(clade, time_sc, name = "weight")

decay_value <- 0.05
# For stronger regularization, try: decay_value <- 0.1

stage_a_mod <- multinom(clade ~ time_sc,
                        data = stage_a_df,
                        weights = weight,
                        maxit = 300,
                        Hess = TRUE,
                        decay = decay_value,
                        trace = FALSE)

simulate_multinom_eta <- function(multinom_model, newdata, nsim = 100) {
  design_mat <- model.matrix(delete.response(terms(multinom_model)), newdata)
  coef_mat <- coef(multinom_model)
  if (is.vector(coef_mat)) {
    coef_mat <- matrix(coef_mat, nrow = 1)
  }
  coef_mu <- as.vector(t(coef_mat))
  coef_cov <- MASS::ginv(multinom_model$Hessian)
  coef_draws <- MASS::mvrnorm(nsim, coef_mu, coef_cov)
  if (is.vector(coef_draws)) {
    coef_draws <- matrix(coef_draws, nrow = 1)
  }
  
  classes <- multinom_model$lev
  nclasses <- length(classes)
  eta <- array(0, dim = c(nrow(newdata), nclasses, nsim))
  
  for (sim in seq_len(nsim)) {
    beta_sim <- matrix(coef_draws[sim, seq_along(coef_mu)],
                       nrow = nclasses - 1,
                       byrow = TRUE)
    eta[, 2:nclasses, sim] <- design_mat %*% t(beta_sim)
  }
  
  list(eta = eta, classes = classes)
}

predict_stage_a_probs <- function(stage_a_model, pred_df, nsim = 100) {
  stage_a <- simulate_multinom_eta(stage_a_model, pred_df, nsim = nsim)
  nobs <- nrow(pred_df)
  nclasses <- length(stage_a$classes)
  probs <- array(NA_real_, dim = c(nsim, nclasses, nobs))
  
  for (sim in seq_len(nsim)) {
    logits <- stage_a$eta[, , sim]
    logits <- logits - apply(logits, 1, max)
    exp_logits <- exp(logits)
    sim_probs <- exp_logits / rowSums(exp_logits)
    probs[sim, , ] <- t(sim_probs)
  }
  
  list(P = probs, classes = stage_a$classes)
}

assign_region <- function(abbrev) {
  case_when(
    abbrev %in% c("CT", "ME", "MA", "NH", "RI", "VT", "NJ", "NY", "PA") ~ "Northeast",
    abbrev %in% c("IL", "IN", "MI", "OH", "WI", "IA", "KS", "MN", "MO", "NE", "ND", "SD") ~ "Midwest",
    abbrev %in% c("DE", "DC", "FL", "GA", "MD", "NC", "SC", "VA", "WV", "AL", "KY", "MS", "TN", "AR", "LA", "OK", "TX") ~ "South",
    abbrev %in% c("AZ", "CO", "ID", "MT", "NV", "NM", "UT", "WY", "AK", "CA", "HI", "OR", "WA") ~ "West",
    TRUE ~ "Other"
  )
}

location_region_map <- locs |>
  filter(location_name != "US") |>
  select(location = location_name, abbreviation) |>
  mutate(region = assign_region(abbreviation)) |>
  select(location, region)

training_counts <- full_fcast_df |>
  filter(clade %in% present_clades, location != "US") |>
  inner_join(time_control_df |> select(date), by = "date") |>
  inner_join(location_region_map, by = "location")

nat_props <- training_counts |>
  group_by(clade) |>
  summarize(n = sum(sequences), .groups = "drop") |>
  mutate(nat_prop = n / sum(n)) |>
  select(clade, nat_prop)

tau_region <- 200
tau_state <- 50
state_pool_k <- 200

region_props <- training_counts |>
  group_by(region, clade) |>
  summarize(n = sum(sequences), .groups = "drop") |>
  tidyr::complete(region, clade = present_clades, fill = list(n = 0)) |>
  group_by(region) |>
  mutate(region_total = sum(n)) |>
  ungroup() |>
  left_join(nat_props, by = "clade") |>
  mutate(region_prop = (n + tau_region * nat_prop) / (region_total + tau_region),
         region_offset = region_prop / nat_prop)

state_props <- training_counts |>
  group_by(location, region, clade) |>
  summarize(n = sum(sequences), .groups = "drop") |>
  tidyr::complete(location, clade = present_clades, fill = list(n = 0)) |>
  left_join(location_region_map, by = "location") |>
  select(location, region = region.y, clade, n) |>
  group_by(location, region) |>
  mutate(state_total = sum(n)) |>
  ungroup() |>
  left_join(region_props |> select(region, clade, region_prop), by = c("region", "clade")) |>
  mutate(state_prop = (n + tau_state * region_prop) / (state_total + tau_state),
         state_offset = state_prop / region_prop) |>
  select(location, clade, state_total, state_offset)

offset_table <- location_region_map |>
  tidyr::crossing(clade = present_clades) |>
  left_join(region_props |> select(region, clade, region_offset), by = c("region", "clade")) |>
  left_join(state_props, by = c("location", "clade")) |>
  mutate(state_total = replace_na(state_total, 0),
         state_offset = replace_na(state_offset, 1),
         region_offset = replace_na(region_offset, 1),
         state_weight = state_total / (state_total + state_pool_k),
         total_offset = region_offset * (state_offset ^ state_weight))

location_offset_mat <- offset_table |>
  select(location, clade, total_offset) |>
  tidyr::pivot_wider(names_from = clade, values_from = total_offset) |>
  arrange(location)

location_offset_matrix <- as.matrix(location_offset_mat |> select(all_of(present_clades)))
rownames(location_offset_matrix) <- location_offset_mat$location

predict_hierarchical_probs <- function(stage_a_model, pred_df, offset_matrix, nsim = 100) {
  stage_a <- predict_stage_a_probs(stage_a_model, pred_df, nsim = nsim)
  classes <- stage_a$classes
  probs <- stage_a$P
  
  for (i in seq_len(nrow(pred_df))) {
    loc_i <- pred_df$location[i]
    if (loc_i %in% rownames(offset_matrix)) {
      offsets <- offset_matrix[loc_i, classes]
    } else {
      offsets <- rep(1, length(classes))
      names(offsets) <- classes
    }
    for (sim in seq_len(nsim)) {
      p <- probs[sim, , i] * offsets
      probs[sim, , i] <- p / sum(p)
    }
  }
  
  list(P = probs, classes = classes)
}

pred_samp_df <- all_possible_data |>
  inner_join(time_control_df, by = "date") |>
  filter(location != "US") |>
  distinct(location, date, time) |>
  inner_join(mod_fcast_df |> distinct(location), by = "location") |>
  arrange(location, time) |>
  mutate(time_sc = (time - time_center) / time_scale)

two_stage_preds <- predict_hierarchical_probs(stage_a_mod, pred_samp_df, location_offset_matrix, nsim = 100)

trajectory_preds <- apply(two_stage_preds$P, 3, identity, simplify = FALSE) |>
  map(as_tibble) |>
  bind_rows(.id = "obs_id") |>
  mutate(obs_id = as.integer(obs_id)) |>
  left_join(pred_samp_df |>
              mutate(obs_id = row_number()) |>
              select(obs_id, location, date, time),
            by = "obs_id") |>
  group_by(obs_id) |>
  mutate(samp = row_number()) |>
  ungroup()

prob_cols <- setdiff(colnames(trajectory_preds), c("obs_id", "location", "date", "time", "samp"))
if (length(prob_cols) != length(two_stage_preds$classes)) {
  stop(paste0("Prediction column mismatch. Found ", length(prob_cols),
              " probability columns but ", length(two_stage_preds$classes), " clades."))
}
colnames(trajectory_preds)[match(prob_cols, colnames(trajectory_preds))] <- two_stage_preds$classes

prediction_trajectories <- trajectory_preds |>
  select(location, date, time, samp, all_of(two_stage_preds$classes)) |>
  gather(clade, prevalence, -location, -date, -time, -samp)

prediction_plotting <- prediction_trajectories |>
  group_by(location, date, time, clade) |>
  summarize(mean = mean(prevalence),
            lower = quantile(prevalence, probs = 0.025),
            upper = quantile(prevalence, probs = 0.975),
            .groups = "drop")

# Plot national stage-A trajectories --------------------------------------
national_pred_df <- time_control_df |>
  distinct(date, time) |>
  arrange(time) |>
  mutate(time_sc = (time - time_center) / time_scale)

national_stage_a_preds <- predict_stage_a_probs(stage_a_mod, national_pred_df, nsim = 200)

national_plot_df <- apply(national_stage_a_preds$P, 3, identity, simplify = FALSE) |>
  map(as_tibble) |>
  bind_rows(.id = "obs_id") |>
  mutate(obs_id = as.integer(obs_id)) |>
  left_join(national_pred_df |>
              mutate(obs_id = row_number()) |>
              select(obs_id, date, time),
            by = "obs_id")

national_prob_cols <- setdiff(colnames(national_plot_df), c("obs_id", "date", "time"))
colnames(national_plot_df)[match(national_prob_cols, colnames(national_plot_df))] <- national_stage_a_preds$classes

national_plot_df <- national_plot_df |>
  gather(clade, prevalence, -obs_id, -date, -time) |>
  group_by(date, time, clade) |>
  summarize(mean = mean(prevalence),
            lower = quantile(prevalence, probs = 0.025),
            upper = quantile(prevalence, probs = 0.975),
            .groups = "drop")

national_obs_df <- full_fcast_df |>
  filter(clade %in% present_clades) |>
  group_by(date, clade) |>
  summarize(sequences = sum(sequences), .groups = "drop") |>
  group_by(date) |>
  mutate(prob = sequences / sum(sequences)) |>
  ungroup() |>
  inner_join(national_pred_df |> select(date), by = "date")

ggplot(national_plot_df, aes(date, mean, color = clade, fill = clade)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .2, color = NA) +
  geom_line() +
  geom_line(data = national_obs_df, aes(y = prob, color = clade), linetype = 2) +
  geom_point(data = national_obs_df, aes(y = prob, color = clade), size = 1) +
  facet_wrap(~clade) +
  labs(title = "National Stage-A Clade Trajectories",
       x = NULL,
       y = "Prevalence",
       color = NULL,
       fill = NULL) +
  theme(legend.position = 'none') -> national_stage_a_plot

dir.create('figures/rt-forecasts', showWarnings = FALSE)
ggsave(filename = paste0('figures/rt-forecasts/', todays_date, '-national-stage-a-trajectories.pdf'),
       plot = national_stage_a_plot,
       width = 12,
       height = 9,
       bg = 'white')




# Plot trajectories for each state ----------------------------------------
get_gam_fit <- function(clade_data){
  ## Function that produces a nice smoother line through the data
  ## The line is useful to compare with forecast projections to make sure nothing is really off
  require(gam)
  if(nrow(clade_data)<4){
    tibble(
      date = clade_data$date,
      time=clade_data$time,
      pred = clade_data$prob
    )
  } else{
    mod <- gam(prob ~ s(time, 3), data = clade_data, weights = clade_data$sequences+1)  
    tibble(
      date = clade_data$date,
      time = clade_data$time,
      pred = mod$fitted.values
    )
  }
}


dir.create('figures/rt-forecasts', showWarnings = FALSE)
pdf(paste0('figures/rt-forecasts/', todays_date, 'trajectory-forecasts.pdf'), width = 12, height = 9, bg='white')
for(loc in unique(prediction_plotting$location)) {
  print(loc)
  
  loc_data <- full_fcast_df |> 
    filter(location == loc) |> 
    group_by(location, date) |> 
    mutate(prob = sequences/sum(sequences)) |> 
    ungroup() |> 
    filter(date > ymd(todays_date)-days(360)) |> 
    mutate(time = as.numeric(difftime(date, min(date), units = 'days')))
  
  loc_fit_probs <- prediction_trajectories |> 
    filter(location == loc) |> 
    mutate(time = as.numeric(time))
  
  if(nrow(loc_data) == 0){
    loc_fit_probs |> 
      ggplot(aes(date, prevalence, color = clade, group = samp), alpha = .2) +
      geom_line() +
      facet_wrap(~clade) +
      labs(title = loc, x = NULL, y = "Prevalence", color = NULL, fill = NULL) +
      scale_x_date(limits = c((ymd(todays_date)-days(360)), ymd(todays_date))) +
      theme(legend.position = 'none') -> loc_plot
  } else{
    loc_data |> 
      nest(data = c('date', 'time', 'prob', 'sequences')) |> 
      mutate(preds = map(.x = data, .f = get_gam_fit)) |> 
      select(-data) |> 
      unnest(preds) |> 
      ggplot(aes(date, pred, color = clade)) +
      geom_point(data = loc_data, aes(x = date, y =prob, color = clade)) +
      geom_line(lty = 2) +
      facet_wrap(~clade) +
      labs(title = loc, x = NULL, y = "Prevalence", color = NULL, fill = NULL) +
      theme(legend.position = 'none') +
      geom_line(data = loc_fit_probs, aes(date, prevalence, group = samp), alpha = .2) -> loc_plot
  }
    
  print(loc_plot)
}
dev.off()



# Produce the files for submission ----------------------------------------

prediction_trajectories |>
  ungroup() |> 
  mutate(time = as.numeric(time)) |> 
  left_join(time_control_df |> select(time, forecast_date), by = 'time') |> 
  left_join(locs |> select(location_name, abbreviation), by = c('location' = 'location_name')) |> 
  filter(forecast_date) |> 
  arrange(time) |> 
  mutate(nowcast_date = todays_date,
         target_date = date,
         location = abbreviation,
         output_type = 'sample',
         output_type_id = paste0(location, samp),
         value = prevalence) |> 
  select(nowcast_date,
         target_date,
         location,
         clade, 
         output_type,
         output_type_id,
         value) -> predictions_for_submission

## Ensure every sample trajectory has all required clades; fill missing with zero prevalence
required_submission_index <- predictions_for_submission |>
  distinct(nowcast_date, target_date, location, output_type, output_type_id)

predictions_for_submission |>
  anti_join(tibble(clade = fcast_clades), by = "clade") |>
  distinct(clade) -> unexpected_clades

if (nrow(unexpected_clades) > 0) {
  warning(paste0("Dropping unexpected clades from submission: ",
                 paste(unexpected_clades$clade, collapse = ", ")))
}

predictions_for_submission |>
  filter(clade %in% fcast_clades) -> predictions_for_submission

required_submission_index |>
  tidyr::crossing(tibble(clade = fcast_clades)) |>
  anti_join(predictions_for_submission,
            by = c("nowcast_date", "target_date", "location",
                   "output_type", "output_type_id", "clade")) |>
  nrow() -> missing_clade_rows_before_fill

if (missing_clade_rows_before_fill > 0) {
  message(paste0("Adding ", missing_clade_rows_before_fill, " missing clade rows with value = 0."))
}

required_submission_index |>
  tidyr::crossing(tibble(clade = fcast_clades)) |>
  left_join(predictions_for_submission,
            by = c("nowcast_date", "target_date", "location",
                   "output_type", "output_type_id", "clade")) |>
  mutate(value = replace_na(value, 0)) |>
  arrange(target_date, location, output_type_id, clade) -> predictions_for_submission

missing_clade_rows_after_fill <- required_submission_index |>
  tidyr::crossing(tibble(clade = fcast_clades)) |>
  anti_join(predictions_for_submission,
            by = c("nowcast_date", "target_date", "location",
                   "output_type", "output_type_id", "clade")) |>
  nrow()

if (missing_clade_rows_after_fill > 0) {
  warning(paste0("Submission still missing ", missing_clade_rows_after_fill, " clade rows after completion."))
}
           
## There should be no rows that show up for this, because that would mean
## trajectories sum about 100% for prevalence
predictions_for_submission |> 
  group_by(location, target_date, output_type_id) |> 
  summarize(prev = sum(value)) |> 
  filter(abs(prev-1)>.001) 


dir.create('processed-data/rt-forecasts', showWarnings = FALSE)
predictions_for_submission |> 
  write_parquet(sink = paste0('processed-data/rt-forecasts/', todays_date, '-UGA-multicast.parquet'))


# Check the file ----------------------------------------------------------
## Move the file to the proper directory and then use hubvalidations to check it
library(hubValidations)


setwd('../variant-nowcast-hub/')
validate_submission(hub_path = '.',
                    file_path = paste0('UGA-multicast/', todays_date, '-UGA-multicast.parquet')) -> sub_validation

# Want all \green checkmarks
sub_validation

## Want to make sure there are no missing required values
sub_validation$req_vals$missing

setwd('../multicast/')

# Trash code not used -----------------------------------------------------

# summary(multi_mod)

# multi_mod <- multinom(clade ~ time + location, 
#                       data=mod_fcast_df, 
#                       maxit = 200,
#                       Hess=T) 
# MASS::ginv(multi_mod$Hessian) -> inv_hess
# 
# multi_mod$coefnames
# 
# coef_samps <- MASS::mvrnorm(2, c(coef(multi_mod)), inv_hess)
# coef_samps
# 
# model.matrix(multi_mod) -> x
# plogis((coef_samps[1,1:7] + coef_samps[1,8:14]*200)/
#          sum(coef_samps[1,1:7] + coef_samps[1,8:14]*200)) 
# 
# t(coef_samps[1,seq(1,14,by=7)]) %*% t(x) -> mat
# 
# mod_fcast_df |> 
#   mutate(var1_pred_prob = exp(c(mat))) |> 
#   count(location, time)
#   

## No longer needed because of mnlpred function
# ## Calculates the lo and hi of 99% confidence interval for the mean prevalence at forecast dates for plotting...
# fit_probs <- Effect(c('time', 'location_num'), se = list(compute = T, level = 0.99, type = 'pointwise'),
#                         multi_mod,
#                         xlevels = list(time = time_control_df$time[time_control_df$forecast_date],
#                                        location_num = unique(mod_fcast_df$location_num)))

# # Format model outputs ----------------------------------------------------
# fit_probs %>% 
#   as_tibble() %>% 
#   select(time, location_num, contains('prob'), -contains('se')) %>% 
#   gather(metric, prob, -time, -location_num) |> 
#   mutate(metric = str_replace(metric, pattern = 'L.prob', replacement = 'lo')) %>% 
#   mutate(metric = str_replace(metric, pattern = 'U.prob', replacement = 'hi')) %>% 
#   separate(metric, into = c('metric', 'clade'), sep = '\\.', extra = 'merge') %>% 
#   mutate(clade = str_remove(clade, pattern = 'X')) |> 
#   spread(metric, prob) %>% 
#   mutate(date = min(mod_fcast_df$date)+days(time)) |> 
#   left_join(dummy_loc_conversion, by = 'location_num') -> fit_probs_clean


# 
# 
# pdf("figures/state-forecasts-summary.pdf", width = 12, height = 9)
# for(loc in unique(prediction_plotting$location)) {
# 
#   loc_data <- full_fcast_df |>
#     filter(location == loc) |>
#     group_by(location, date) |>
#     mutate(prob = sequences/sum(sequences)) |>
#     ungroup() |>
#     filter(date > ymd(todays_date)-days(360)) |>
#     mutate(time = as.numeric(difftime(date, min(date), units = 'days')))
# 
#   loc_fit_probs <- prediction_plotting |>
#     filter(location == loc) |>
#     left_join(time_control_df, by = 'time')
# 
#   print(loc)
#   loc_data |>
#     nest(data = c('date', 'time', 'prob', 'sequences')) |>
#     mutate(preds = map(.x = data, .f = get_gam_fit)) |>
#     select(-data) |>
#     unnest(preds) |>
#     ggplot(aes(date, pred, color = clade)) +
#     geom_ribbon(data = loc_fit_probs, aes(x = date, ymin = lower, ymax = upper, fill = clade), alpha = .2, color = NA, inherit.aes=F) +
#     geom_point(data = loc_data, aes(x = date, y =prob, color = clade)) +
#     geom_line(lty = 2) +
#     facet_wrap(~clade) +
#     labs(title = loc, x = NULL, y = "Prevalence", color = NULL, fill = NULL) +
#     geom_line(data = loc_fit_probs, aes(date, mean)) -> loc_plot
# 
#   print(loc_plot)
# }
# dev.off()
