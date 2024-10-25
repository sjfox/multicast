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
# todays_date <- ymd('2024-10-23') ## Only use this if you want to test out script
todays_date <- Sys.Date()


json_location <- paste0("https://raw.githubusercontent.com/reichlab/variant-nowcast-hub/refs/heads/main/auxiliary-data/modeled-clades/",
                        todays_date,
                        ".json")
clades <- fromJSON(paste(readLines(json_location), collapse=""))
fcast_clades <- clades$clades

# Do final cleaning for the analysis --------------------------------------
## Currently removing the most recent 30 days, because it's so noisy and messes with predictions (something to test later though)
## Currently only using data back 180 days (something to test later on)
## Parameters used to control the amount of data used and the dates used for model predictions
fcast_days_from_today <- 10
hindcast_days_from_today <- 31
days_to_look_back <- 180
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



## Now make the actually tibble for forecasting
full_fcast_df |> 
  inner_join(time_control_df, by = 'date') |> 
  uncount(sequences)  -> mod_fcast_df



# Fit the model and estimate uncertainty ----------------------------------
library(fastDummies)

## This turns the categorical variable into dummy columns with zeroes and ones to indicate which 
## location the data are a part of. multinom() does this automatically, but
## it doesn't work with the mnlpred function I'm using for simulating
## once we do it manually, it works with mnlpred... A bit annoying but not too bad
dummy_cols(mod_fcast_df |> 
             mutate(location = str_replace_all(location, fixed(" "), "")), select_columns = 'location') |> 
  select(-location, ) |> 
  select(clade, time, starts_with('location_')) -> mod_fcast_df

multi_mod <- multinom(clade ~ ., 
                          data=mod_fcast_df, 
                          maxit = 200,
                          Hess=T) 


# Get simulated trajectories and ribbons for plotting ------------------------------------

get_multinom_loc_preds <- function(location_name, 
                                   df, 
                                   multinomial_model,
                                   ntraj = 100){
  # browser()
  ## Function that gets trajectories and summary predictions for a given location
  
  print(paste0("Getting simulated trajectories for: ", location_name))
  
  mnl_pred_ova2_revised(multinomial_model, 
                        data = df,
                        x = 'time', 
                        xvals = min(df$time):max(df$time), ## Can change to whatever "time" values you want to explore trajectories
                        # xvals = seq(100,400, by = 20), 
                        nsim = ntraj) -> preds
  
  ## Munges the prediction df output into proper format for using
  apply(preds$P, 3, identity, simplify=FALSE) |> 
    map(.f = as_tibble) |> 
    bind_rows(.id = 'time') |> 
    group_by(time) |> 
    mutate(samp = seq_along(time)) -> trajectory_preds
  
  colnames(trajectory_preds) <- c('time', sort(unique(df$clade)), 'samp')
  trajectory_preds |> 
    gather(clade, prevalence, -time, -samp) |> 
    mutate(location = location_name) -> trajectory_preds
  
  ## If you want to plot a trajectory for specific clade/location to zoom in
  # trajectory_preds |> 
  #   filter(clade == '24F') |> 
  #   ggplot(aes(as.numeric(time), prevalence, group = samp)) + geom_line()
  
  
  preds$plotdata |> 
    as_tibble() |> 
    mutate(location = location_name) -> summary_preds
  
  list(trajectory_preds, summary_preds)
}

## Creates the "predictor" variable database that is used with mnlpred function to
## Get the prevalence estimates for each location at each timestep
dummy_cols(all_possible_data |> 
             inner_join(time_control_df, by = 'date') |> 
             filter(location !='US') |> 
             mutate(location = str_replace_all(location, fixed(" "), "")), 
           select_columns = 'location') -> pred_samp_df


pred_samp_df |> 
  nest(data = -location) |> 
  mutate(preds = map2(location, data, 
                      get_multinom_loc_preds,
                      multinomial_model = multi_mod)) -> all_predictions
  

all_predictions$preds |> 
  map(.f = function(x) x[[1]]) |> 
  bind_rows() |> 
  left_join(locs |> 
              select(location_name) |> 
              mutate(location = str_replace_all(location_name, fixed(" "), ""))) |> 
  mutate(location = location_name) |> 
  select(-location_name) -> prediction_trajectories

all_predictions$preds |> 
  map(.f = function(x) x[[2]]) |> 
  bind_rows() |> 
  left_join(locs |> 
              select(location_name) |> 
              mutate(location = str_replace_all(location_name, fixed(" "), ""))) |> 
  mutate(location = location_name) |> 
  select(-location_name) -> prediction_plotting




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
  
  loc_data <- full_fcast_df |> 
    filter(location == loc) |> 
    group_by(location, date) |> 
    mutate(prob = sequences/sum(sequences)) |> 
    ungroup() |> 
    filter(date > ymd(todays_date)-days(360)) |> 
    mutate(time = as.numeric(difftime(date, min(date), units = 'days')))
  
  loc_fit_probs <- prediction_trajectories |> 
    filter(location == loc) |> 
    mutate(time = as.numeric(time)) |> 
    left_join(time_control_df, by = 'time')
  
  print(loc)
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
  
  print(loc_plot)
}
dev.off()



# Produce the files for submission ----------------------------------------

prediction_trajectories |>
  ungroup() |> 
  mutate(time = as.numeric(time)) |> 
  left_join(time_control_df, by = 'time') |> 
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
## Something is broken for some reason - need to investigate it later
library(hubValidations)
validate_submission(hub_path = '../variant-nowcast-hub/',
                    file_path = 'UGA-multicast/2024-10-23-UGA-multicast.parquet') -> sub_validation

# Want all \green checkmarks
sub_validation

## Want to make sure there are no missing required values
sub_validation$req_vals$missing

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