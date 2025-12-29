library(mgcv)
library(tidyverse)
library(lubridate)
library(ggplot2)
library(dplyr)
library(data.table)
library(httr2)
library(jsonlite)
setwd("[insert file path to working directory]")
datadir <- "./data/"
outdir <- "./output/"
indir <- "./input/"

# Read in observed deaths data
dat_monthly <- read.csv(paste0(datadir,"monthly data.csv"), sep=',', header=TRUE)%>%
  mutate(date_yr_mnth = as.Date(date_yr_mnth, format="%m/%d/%Y"))

# Read in estimates from the CDC
cdc_excess <- read.csv(paste0(datadir,"CDC Excess Deaths/Excess_Deaths_Associated_with_COVID-19.csv"), sep=',', header=TRUE) %>%
  dplyr::filter((State == "United States") & (Type == "Predicted (weighted)") & (Outcome == "All causes"))%>%
  dplyr::mutate(year = as.numeric(Year)) %>%
  dplyr::select(-c(Year))

# Read in estimates from Karlinsky & Kobak 
wmd_excess <- read.csv(paste0(indir,"WMD/excess-mortality-timeseries.csv"), sep=',', header=TRUE) %>%
  filter(country_name == "United States")

# ------------------------ CDC ------------------------
cdc_yearly <- cdc_excess %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(expected_deaths = sum(Average.Expected.Count,na.rm=T),
                   excess = sum(Excess.Estimate,na.rm=T)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(year>=2020)

write.csv(cdc_yearly,file=paste0(outdir,"cdc_yearly.csv"), row.names = FALSE)

# ------------------------ World Mortality Dataset ------------------------
wmd_yearly <- wmd_excess %>%
  group_by(year) %>%
  summarise(excess = sum(excess.deaths)) %>%
  ungroup()

write.csv(wmd_yearly,file=paste0(outdir,"wmd_yearly.csv"), row.names = FALSE)

# ------------------------ Msemburi et. al  ------------------------ 
# Link to developer's repository: https://github.com/WHOexcessc19/Codebase
# The below code is adapted from spline_model_monthly.R and Final_Sampling.R, for the US case study.

dat_mo_f1 <- dat_monthly %>% rename(observed = outcome)
dat_mo_f1$iso3 <- "USA"
observed2020 <- dat_mo_f1 %>% 
  rename(date = date_yr_mnth) %>%
  filter(date >= "2020-01-01")

## Clean data for use in gam function
exp_obs_bymonthyear <- dat_mo_f1 %>% 
  select(iso3, year, month, observed) %>%
  filter(year < 2020) %>% 
  mutate(country_num = as.numeric(as.factor(iso3)),
         observed = floor(observed)) %>% 
  arrange(iso3, year, month)

## Create data frame to store expected monthly mortality in 2020-2021
# gamma_E and gamma_delta will contain, for each country time period, the 
# parameters of the gamma distribution for the expected mortality
acm_predictions <- data.frame(iso3 = rep(countries, each = 60),
                              year = rep(c(rep(2020, 12), 
                                           rep(2021, 12),
                                           rep(2022, 12),
                                           rep(2023, 12),
                                           rep(2024, 12)), 
                                         times = length(countries)),
                              month = rep(c(1:12, 1:12, 1:12, 1:12, 1:12), 
                                          times = length(countries)),
                              expected_acm = NA,
                              expected_acm_se = NA,
                              expected_log_acm = NA,
                              expected_log_acm_se = NA,
                              gamma_E = NA,
                              gamma_delta = NA,
                              gamma_E_nb = NA,
                              gamma_delta_nb = NA) %>%
  mutate(country_num = as.numeric(as.factor(iso3)))

exp_obs_bymonthyear$expected_acm <- NA
exp_obs_bymonthyear$expected_acm_se <- NA
exp_obs_bymonthyear$expected_log_acm <- NA
exp_obs_bymonthyear$expected_log_acm_se <- NA
exp_obs_bymonthyear$gamma_E <- NA
exp_obs_bymonthyear$gamma_delta <- NA

num_samples <- 10000

## Predict expected monthly mortality in 2020-2021
set.seed(123)
for(i in 1:max(exp_obs_bymonthyear$country_num)){
  whichs <- which(exp_obs_bymonthyear$country_num == i)
  temp <- exp_obs_bymonthyear[whichs, ]
  
  # Fit gam
  # If there are less than 3 years of historical data, use a linear annual 
  # trend, rather than a spline trend
  if(length(unique(temp$year)) < 3){
    print(temp$iso3[1])
    annual_model <- gam(observed ~ year +
                          s(month, bs = "cc", k = length(unique(temp$month))),
                        data = temp, family = nb(theta = NULL, link = "log"))
  } else{
    annual_model <- gam(observed ~ s(year, k = length(unique(temp$year))) +
                          s(month, bs = "cc", k = length(unique(temp$month))),
                        data = temp, family = nb(theta = NULL, link = "log"))
  }
  overd <- exp(annual_model$family$getTheta())
  
  # Get predictions 
  pred <- predict(annual_model,
                  se.fit = TRUE,
                  type = "response",
                  newdata = data.frame(year = c(rep(2020, 12),
                                                rep(2021, 12),
                                                rep(2022, 12),
                                                rep(2023, 12),
                                                rep(2024, 12)),
                                       month = c(1:12, 1:12, 1:12, 1:12, 1:12)))
  whichs_pred <- which(acm_predictions$country_num == i)
  acm_predictions[whichs_pred, "expected_acm"] <- pred$fit
  acm_predictions[whichs_pred, "expected_acm_se"] <- pred$se.fit
  pred_log <- predict(annual_model,
                      se.fit = TRUE,
                      newdata = data.frame(year = c(rep(2020, 12),
                                                    rep(2021, 12),
                                                    rep(2022, 12),
                                                    rep(2023, 12),
                                                    rep(2024, 12)),
                                           month = c(1:12, 1:12, 1:12, 1:12, 1:12)))
  acm_predictions[whichs_pred, "expected_log_acm"] <- pred_log$fit
  acm_predictions[whichs_pred, "expected_log_acm_se"] <- pred_log$se.fit
  
  # Get gamma parameters
  gamma_E <- rep(0, 60)
  gamma_delta <- rep(0, 60)
  gamma_E_nb <- rep(0, 60)
  gamma_delta_nb <- rep(0, 60)
  for(j in 1:60){
    samples <- exp(rnorm(num_samples, mean = pred_log$fit[j], 
                         sd = pred_log$se.fit[j]))
    
    gamma_E[j] <- mean(samples)
    gamma_delta[j] <- ((gamma_E[j]) ^ 2) / var(samples)
    
    samples_nb <- rnbinom(num_samples, size = overd, mu = samples)
    
    gamma_E_nb[j] <- mean(samples_nb)
    gamma_delta_nb[j] <- ((gamma_E_nb[j]) ^ 2) / var(samples_nb)
  }
  acm_predictions[whichs_pred, "gamma_E"] <- gamma_E
  acm_predictions[whichs_pred, "gamma_delta"] <- gamma_delta
  acm_predictions[whichs_pred, "gamma_E_nb"] <- gamma_E_nb
  acm_predictions[whichs_pred, "gamma_delta_nb"] <- gamma_delta_nb
  
  # Get fitted values for historical time periods
  pred_hist <- predict(annual_model, se.fit = TRUE, type = "response")
  exp_obs_bymonthyear[whichs, "expected_acm"] <- pred_hist$fit
  exp_obs_bymonthyear[whichs, "expected_acm_se"] <- pred_hist$se.fit
  pred_log_hist <- predict(annual_model, se.fit = TRUE)
  exp_obs_bymonthyear[whichs, "expected_log_acm"] <- pred_log_hist$fit
  exp_obs_bymonthyear[whichs, "expected_log_acm_se"] <- pred_log_hist$se.fit
  
  num_hist <- length(pred_log_hist$fit)
  gamma_E_hist <- rep(0, num_hist)
  gamma_delta_hist <- rep(0, num_hist)
  for(j in 1:num_hist){
    samples <- exp(rnorm(num_samples, mean = pred_log_hist$fit[j], 
                         sd = pred_log_hist$se.fit[j]))
    
    gamma_E_hist[j] <- mean(samples)
    gamma_delta_hist[j] <- ((gamma_E_hist[j]) ^ 2) / var(samples)
    
  }
  exp_obs_bymonthyear[whichs, "gamma_E"] <- gamma_E_hist
  exp_obs_bymonthyear[whichs, "gamma_delta"] <- gamma_delta_hist
}

## Obtain excess estimate
# Number of samples of estimates
num_samps = 1000
# How many country time points
df.inla <- acm_predictions %>%
  left_join(observed2020[,1:4], by=c("year", "month"))%>%
  mutate(WHO_region = "AMRO",
         Country = iso3,
         months = month)
country.time.n = nrow(df.inla)

## Function for getting expecteds, excess deaths, adapted from the code by WHO
estimates.sampling <- function(country.time.n,num_samps){
  
  # Create dfs of Country, WHO Region, and months to merge in estimates with later
  excess.df.expec <- df.inla %>%
    dplyr::select(Country,WHO_region,year,months)
  excess.df.acm <- df.inla %>%
    dplyr::select(Country,WHO_region,year,months)
  excess.df.excess <- df.inla %>%
    dplyr::select(Country,WHO_region,year,months)
  
  # Loop over samples to replace estimates from subnational and mixed models, and fix benchmarking
  for(i in 1:num_samps){
    
    print(paste("Starting Iteration: ",i))
    
    # Create iteration i covariate model estimates df
    set.seed(42 + i)
    excess.df.i <- data.frame(
      expected.delta=df.inla$gamma_delta,
      expected.E=df.inla$gamma_E,Country=df.inla$Country,
      year=df.inla$year,
      months=df.inla$months,
      observed.ind=ifelse(!is.na(df.inla$observed),TRUE,FALSE),
      observed=df.inla$observed,
      expected.sim=rgamma(country.time.n,shape=df.inla$gamma_delta,rate=df.inla$gamma_delta/df.inla$gamma_E))
    
    # Iteration i vectors for expecteds, acm, and excess
    # Estimated acm are replaced by observed when available
    expec.i <- excess.df.i$expected.sim
    acm.i <- df.inla$observed
    excess.i <- acm.i- excess.df.i$expected.sim 
    
    # Add Sample i excess for each country, time point to df
    excess.df.expec <- cbind(excess.df.expec,as.data.frame(expec.i))
    excess.df.acm <- cbind(excess.df.acm,as.data.frame(unlist(acm.i)))
    excess.df.excess <- cbind(excess.df.excess,as.data.frame(unlist(excess.i)))
  }
  
  final.df <- list(excess.df.acm, excess.df.expec, excess.df.excess)
  return(final.df)
  
}

# Setup and run the country estimates function
country.df <- estimates.sampling(country.time.n,
                                 num_samps)

# Format estimates
expected = country.df[[2]]
acm = country.df[[1]]
excess = country.df[[3]]
colnames(expected) <- c("Country", "WHO_region", "year", "months", 
                        paste0("expec", 1:num_samps))
colnames(acm) <- c("Country", "WHO_region", "year", "months", 
                   paste0("acm", 1:num_samps))
colnames(excess) <- c("Country", "WHO_region", "year", "months", 
                      paste0("excess", 1:num_samps))

excess_summarized = cbind(excess[,1:4],excess=rowMeans(excess[,5:(num_samps+4)]),
                          excess_lower95=apply(excess[,5:(num_samps+4)],1,quantile,0.025),
                          excess_upper95=apply(excess[,5:(num_samps+4)],1,quantile,0.975))
expected_summarized = cbind(expected[,1:4],expected=rowMeans(expected[,5:(num_samps+4)]),
                            expected_lower95=apply(expected[,5:(num_samps+4)],1,quantile,0.025),
                            expected_upper95=apply(expected[,5:(num_samps+4)],1,quantile,0.975))

who_monthly <- merge(expected_summarized,excess_summarized,by=c("year","months", "Country", "WHO_region")) %>%
  select(-c(Country,WHO_region)) %>%
  rename(expected_deaths = expected,
         month = months) %>%
  mutate(month_year=format(as.Date(paste0(year,"-",month,"-","01")), "%Y-%m"))%>%
  relocate(month_year,.before=year)%>%
  arrange(year,month)

# Sum across months and then take quantiles to obtain yearly estimates
expected_long <- expected %>%
  pivot_longer(
    cols = starts_with("expec"),
    names_to = "iteration", 
    values_to = "expected"
  ) %>%
  mutate(iteration = as.numeric(gsub(".*?([0-9]+).*", "\\1", iteration)))

excess_long <- excess %>%
  pivot_longer(
    cols = starts_with("excess"),
    names_to = "iteration",
    values_to = "excess"
  ) %>%
  mutate(iteration = as.numeric(gsub(".*?([0-9]+).*", "\\1", iteration)))

yrly_iter <- merge(expected_long,excess_long,by=c("year","months","iteration", "Country", "WHO_region")) %>%
  group_by(iteration, year)%>%
  summarise(yearly_expected = sum(expected),
            yearly_excess = sum(excess)) %>%
  ungroup()

who_yearly <- yrly_iter %>%
  group_by(year)%>%
  summarise(expected_deaths = mean(yearly_expected),
            expected_lower95 = quantile(yearly_expected,0.025),
            expected_upper95 = quantile(yearly_expected,0.975),
            excess = mean(yearly_excess),
            excess_lower95 = quantile(yearly_excess,0.025),
            excess_upper95 = quantile(yearly_excess,0.975)) %>%
  ungroup()

write.csv(who_yearly,file=paste0(outdir,"who_yearly.csv"), row.names = FALSE)
write.csv(who_monthly,file=paste0(outdir,"who_monthly.csv"), row.names = FALSE)

# ------------------------ The Economist  ------------------------ 
# Link to developer's repository: https://github.com/TheEconomist/covid-19-excess-deaths-tracker
# The below code is adapted from excess_deaths_script.R, for the US case study.

## Function for estimating excess deaths as written by The Economist
get_excess_deaths <- function(df, expected_deaths_model, frequency="weekly", calculate=TRUE, train_model=TRUE){
  # Define formulas and count number of regions
  weekly_formula <- as.formula(total_deaths_per_day ~ year + week)
  weekly_regional_formula <- as.formula(total_deaths_per_day ~ year + week + region +
                                          region:year + region:week)
  monthly_formula <- as.formula(total_deaths_per_day ~ year + month)
  monthly_regional_formula <- as.formula(total_deaths_per_day ~ year + month + region +
                                           region:year + region:month)
  quarterly_formula <- as.formula(total_deaths_per_day ~ year + quarter)
  quarterly_regional_formula <- as.formula(total_deaths_per_day ~ year + quarter + region +
                                             region:year + region:quarter)
  df_regions <- length(unique(df$region))
  
  # Convert weeks and months into fixed effects
  if(frequency == "weekly") {
    df <- df %>% mutate(week = as.character(week))
  } else if (frequency == "monthly") {
    df <- df %>% mutate(month = as.character(month))
  } else if (frequency == "quarterly") {
    df <- df %>% mutate(quarter = as.character(quarter))
  }
  
  # Identify the correct formula for the dataframe
  if(frequency == "weekly" & df_regions == 1) {
    expected_deaths_formula <- weekly_formula
  } else if (frequency == "weekly" & df_regions > 1) {
    expected_deaths_formula <- weekly_regional_formula
  } else if (frequency == "monthly" & df_regions == 1) {
    expected_deaths_formula <- monthly_formula
  } else if (frequency == "monthly" & df_regions > 1)  {
    expected_deaths_formula <- monthly_regional_formula
  } else if (frequency == "quarterly" & df_regions == 1) {
    expected_deaths_formula <- quarterly_formula
  } else if (frequency == "quarterly" & df_regions > 1)  {
    expected_deaths_formula <- quarterly_regional_formula
  }
  
  # Calculate expected deaths
  if(calculate == FALSE) {
    
    # Use pre-existing official model
    expected_deaths <- df %>% filter(year >= 2020)
    
  } else if(train_model == FALSE) {
    
    # Use previously trained Economist model
    expected_deaths <- df %>% filter(year >= 2020) %>%
      mutate(expected_deaths = predict(expected_deaths_model,.) * days)
    
  } else if(frequency == "weekly") {
    
    # Create dataframe of week 53, using week 52 and 53 in previous years
    week_53_df <- df %>%
      filter(week %in% c("52","53")) %>% mutate(week = "53", week_53 = 1)
    
    # Train an Economist weekly model
    train_df <- df %>% 
      filter(week != "53") %>%
      bind_rows(week_53_df) %>%
      filter(end_date <= as.Date("2019-12-30")) %>%
      mutate(total_deaths_per_day = total_deaths / days)
    expected_deaths_model <- lm(expected_deaths_formula,train_df)
    expected_deaths <- df %>% filter(year >= 2020) %>%
      mutate(expected_deaths = predict(expected_deaths_model,newdata=.) * days)
    
  } else if(frequency %in% c("monthly","quarterly")) {
    
    # Train an Economist monthly or quarterly model
    train_df <- df %>% 
      filter(end_date <= as.Date("2019-12-30")) %>%
      mutate(total_deaths_per_day = total_deaths / days)
    expected_deaths_model <- lm(expected_deaths_formula,train_df)
    expected_deaths <- df %>% filter(year >= 2020) %>%
      mutate(expected_deaths = predict(expected_deaths_model,newdata=.,se.fit=TRUE)$fit * days,
             # note added below line in order to obtain SE
             expected_deaths_se = predict(expected_deaths_model,newdata=.,se.fit=TRUE)$se.fit * days)
  }
  
  # Set expected deaths to be non-negative (on implementation, this had no impact, but makes code more robust to the off chance that some country had extremely strong time effects and extremely strong downward trend in deaths over time):
  expected_deaths$expected_deaths[expected_deaths$expected_deaths < 0] <- 0
  
  # Calculate excess deaths
  excess_deaths <- expected_deaths %>%
    mutate(excess_deaths = total_deaths - expected_deaths,
           non_covid_deaths = total_deaths - covid_deaths,
           region_code = as.character(region_code)) %>%
    mutate(covid_deaths_per_100k = covid_deaths / population * 100000,
           excess_deaths_per_100k = excess_deaths / population * 100000,
           excess_deaths_pct_change = 
             ((expected_deaths + excess_deaths) / expected_deaths) - 1)
  
  # Calculate weekly rates for monthly and quarterly data
  if(frequency %in% c("monthly","quarterly")) {
    
    excess_deaths <- excess_deaths %>%
      mutate(total_deaths_per_7_days = total_deaths / days * 7,
             covid_deaths_per_7_days = covid_deaths / days * 7,
             expected_deaths_per_7_days = expected_deaths / days * 7,
             excess_deaths_per_7_days = excess_deaths / days * 7,
             non_covid_deaths_per_7_days = non_covid_deaths / days * 7,
             covid_deaths_per_100k_per_7_days = covid_deaths_per_100k / days * 7,
             excess_deaths_per_100k_per_7_days = excess_deaths_per_100k / days * 7)
    
  }
  
  list(expected_deaths_model,excess_deaths)
}

## Further process data
dat_mo_f2 <- dat_monthly %>%
  mutate(country = "United States",
         region = "United States",
         region_code = 1,
         expected_deaths = "TBC") %>%
  rename(total_deaths = outcome) %>%
  mutate(start_date = c(seq(make_date(2013, 12, 30), make_date(2024, 11, 30), by = "month")),
         end_date = c(seq(make_date(2014, 01, 30), make_date(2024, 12, 30), by = "month")),
         days = lead(as.numeric(interval(ymd(start_date),ymd(end_date)) %/% days(1))),
         days = ifelse(is.na(days),31,days))

res_mo <- get_excess_deaths(dat_mo_f2,
                            expected_deaths_model = NULL,
                            frequency = "monthly")

results_monthly <- res_mo[[2]]
results_monthly$excess <- results_monthly$total_deaths-results_monthly$expected_deaths
economist_monthly <- results_monthly %>%
  mutate(expected_lower95 = expected_deaths - (1.96*expected_deaths_se),
         expected_upper95 = expected_deaths + (1.96*expected_deaths_se),
         excess_lower95 = excess - (1.96*expected_deaths_se),
         excess_upper95 = excess + (1.96*expected_deaths_se))%>%
  rename(observed = total_deaths) %>%
  select(year,month,observed,expected_deaths, expected_deaths_se, expected_lower95, expected_upper95, excess, excess_lower95, excess_upper95)

## Obtain yearly counts and UIs
dat_mo_f2$total_deaths_per_day <- dat_mo_f2$total_deaths / dat_mo_f2$days
dat_mo_f2 <- dat_mo_f2 %>%
  mutate(month = as.character(month))

expected_mod <- res_mo[[1]] 

economist_year_excess <- function(year, alpha = 0.95) {
  
  data_y <- dat_mo_f2 %>% dplyr::filter(year == !!year)
  if (nrow(data_y) == 0L) stop("No monthly rows found for year = ", year)
  
  # deaths per day were used in the regression model
  days_y <- as.numeric(data_y$days)
  
  # predicted per-day deaths + full covariance
  Terms  <- delete.response(terms(expected_mod))
  X_new  <- model.matrix(Terms, data_y)
  
  beta     <- coef(expected_mod)
  Sigma_b  <- vcov(expected_mod)
  
  # predicted per-day means
  fit_day <- as.numeric(X_new %*% beta)
  
  # covariance matrix of daily predictions (M × M)
  Sigma_day <- X_new %*% Sigma_b %*% t(X_new)
  
  # convert per-day -> per-month totals
  D <- diag(days_y)
  
  fit_month   <- as.numeric(D %*% fit_day)
  Sigma_month <- D %*% Sigma_day %*% D
  
  # aggregate months -> yearly expected deaths
  w <- rep(1, length(days_y))
  
  year_mean <- sum(fit_month)
  year_var  <- as.numeric(t(w) %*% Sigma_month %*% w)
  year_se   <- sqrt(year_var)
  
  # CI for expected baseline deaths
  z <- qnorm((1 + alpha) / 2)
  exp_lower <- year_mean - z * year_se
  exp_upper <- year_mean + z * year_se
  
  obs_y <- dat_mo_f2%>%
    dplyr::filter(year == !!year) %>%
    dplyr::summarise(obs = sum(total_deaths), .groups = "drop") %>%
    dplyr::pull(obs)
  
  excess_point <- obs_y - year_mean
  excess_lower <- obs_y - exp_upper   
  excess_upper <- obs_y - exp_lower
  
  tibble::tibble(
    year = year,
    observed = obs_y,
    expected_deaths = year_mean,
    expected_deaths_se = year_se,
    expected_lower95 = exp_lower,
    expected_upper95 = exp_upper,
    excess = excess_point,
    excess_lower95 = excess_lower,
    excess_upper95 = excess_upper
  )
}

year2020 <- economist_year_excess(year=2020)
year2021 <- economist_year_excess(year=2021)
year2022 <- economist_year_excess(year=2022)
year2023 <- economist_year_excess(year=2023)
year2024 <- economist_year_excess(year=2024)

economist_yearly <- bind_rows(year2020,year2021,year2022,year2023,year2024)

write.csv(economist_yearly,file=paste0(outdir,"economist_yearly.csv"), row.names = FALSE)
write.csv(economist_monthly,file=paste0(outdir,"economist_monthly.csv"), row.names = FALSE)


# ------------------------ Acosta and Irizarry  ------------------------ 
library(excessmort)

dat_mo_f3<-dat_monthly%>%
  mutate(date = as.Date(format(as.Date(date_yr_mnth, format="%m/%d/%Y"),"%Y-%m-%d")))
  
### Compute expected counts
## Main results: using reference period of 2014-2019 and using extrapolate = FALSE
# Remove COVID years and periods with bad flu season
exclude_dates_v1 <- c(
  make_date(2015, 1, 1),
  make_date(2017, 1, 1),
  make_date(2018, 1, 1),
  seq(make_date(2020, 1, 1), make_date(2024, 12, 1), by = "month"))

# Flexible model to account for tappering off trend
expected_v1 <- excessmort::compute_expected(dat_mo_f3, exclude = exclude_dates_v1, harmonics = 2, trend.knots.per.year = 0.25, frequency = 12, weekday.effect = FALSE, extrapolate=FALSE)

# Estimate excess for each year and also the curve
em_v1 <- excess_model(expected_v1, 
                      start = make_date(2014, 1, 1), 
                      end = make_date(2024, 12, 1), 
                      exclude = c(make_date(2015, 1, 1),
                                  make_date(2017, 1, 1),
                                  make_date(2018, 1, 1)), # exclude months with abnormal flu
                      knots.per.year = 3, 
                      intervals = list(seq(make_date(2020, 1, 1), make_date(2020, 12, 1), by = "month"),
                                       seq(make_date(2021, 1, 1), make_date(2021, 12, 1), by = "month"),
                                       seq(make_date(2022, 1, 1), make_date(2022, 12, 1), by = "month"),
                                       seq(make_date(2023, 1, 1), make_date(2023, 12, 1), by = "month"),
                                       seq(make_date(2024, 1, 1), make_date(2024, 12, 1), by = "month")))

## Sensitivity results: using reference period of 2017-2019 and using extrapolate = TRUE for 2017-2019
dat_monthly_frm2017 <- dat_mo_f3 %>%
  filter(date >= as_date("2017-01-01"))

# Exclude COVID and months with bad flu season
exclude_dates_v2 <- c(
  make_date(2017, 1, 1),
  make_date(2018, 1, 1),
  seq(make_date(2020, 1, 1), make_date(2024, 12, 1), by = "month"))

# Extending linear trend from 2017-2019 for extrapolation
expected_v2 <- excessmort::compute_expected(dat_monthly_frm2017, exclude = exclude_dates_v2, harmonics = 2, trend.knots.per.year = 0.25, frequency = 12, weekday.effect = FALSE, extrapolate=TRUE)

# Estimate excess for each year and also the curve
em_v2 <- excess_model(expected_v2, start = make_date(2017, 1, 1), end = make_date(2024, 12, 1), 
                      exclude = c(make_date(2017, 1, 1),
                                  make_date(2018, 1, 1)), # exclude months with abnormal flu
                      knots.per.year = 3, 
                      intervals = list(seq(make_date(2020, 1, 1), make_date(2020, 12, 1), by = "month"),
                                       seq(make_date(2021, 1, 1), make_date(2021, 12, 1), by = "month"),
                                       seq(make_date(2022, 1, 1), make_date(2022, 12, 1), by = "month"),
                                       seq(make_date(2023, 1, 1), make_date(2023, 12, 1), by = "month"),
                                       seq(make_date(2024, 1, 1), make_date(2024, 12, 1), by = "month")))

## Additional diagnostics
diagnostics_v1 <- expected_diagnostic(expected_v1)
diagnostics_v2 <- expected_diagnostic(expected_v2)
# Fitted model
diagnostics_v1$expected
diagnostics_v2$expected
# Examine the trend
diagnostics_v1$trend
diagnostics_v2$trend

tibble(date = em_v1$date, obs = em_v1$observed, excess = em_v1$expected*(exp(em_v1$fitted) - 1), 
       upper = em_v1$expected*(exp(em_v1$fitted + 1.96*em_v1$log_expected_se) - 1),
       lower = em_v1$expected*(exp(em_v1$fitted - 1.96*em_v1$log_expected_se) - 1)) |>
  filter(date >= make_date(2022,1,1)) |>
  ggplot(aes(date, excess)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5, fill = "purple") +
  geom_line(color = "purple") +
  geom_hline(yintercept = 0, lty = 2) +
  theme_bw()

tibble(date = em_v2$date, obs = em_v2$observed, excess = em_v2$expected*(exp(em_v2$fitted) - 1), 
       upper = em_v2$expected*(exp(em_v2$fitted + 1.96*em_v2$log_expected_se) - 1),
       lower = em_v2$expected*(exp(em_v2$fitted - 1.96*em_v2$log_expected_se) - 1)) |>
  filter(date >= make_date(2022,1,1)) |>
  ggplot(aes(date, excess)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5, fill = "purple") +
  geom_line(color = "purple") +
  geom_hline(yintercept = 0, lty = 2) +
  theme_bw()

## Expected and excess deaths for each year
yearly_est_v1 <- em_v1$excess |> mutate(year = year(start)) |> 
  mutate(expected_lower95 = expected - 1.96*sd, expected_upper95 = expected + 1.96*sd,
         excess_lower95 = excess - 1.96*sd, excess_upper95 = excess + 1.96*sd) |>
  select(year,  observed, expected, sd, expected_lower95, expected_upper95, 
         excess, excess_lower95, excess_upper95) |>
  mutate(across(-year, function(x) round(x,.1) )) |>
  arrange(desc(year)) 

yearly_est_v2 <- em_v2$excess |> mutate(year = year(start)) |> 
  mutate(expected_lower95 = expected - 1.96*sd, expected_upper95 = expected + 1.96*sd,
         excess_lower95 = excess - 1.96*sd, excess_upper95 = excess + 1.96*sd) |>
  select(year,  observed, expected, sd, expected_lower95, expected_upper95, 
         excess, excess_lower95, excess_upper95) |>
  mutate(across(-year, function(x) round(x,.1) )) |>
  arrange(desc(year)) 

## Expected and excess deaths for each month
mo_tb_v1 <- tibble(date = em_v1$date, observed = em_v1$observed, 
                   expected_deaths = em_v1$expected,
                   expected_lower95 = (exp(log(em_v1$expected) - 1.96*em_v1$log_expected_se) - 1),
                   expected_upper95 = (exp(log(em_v1$expected) + 1.96*em_v1$log_expected_se) - 1),
                   excess = em_v1$expected*(exp(em_v1$fitted) - 1), 
                   excess_lower95 = em_v1$expected*(exp(em_v1$fitted - 1.96*em_v1$log_expected_se) - 1),
                   excess_upper95 = em_v1$expected*(exp(em_v1$fitted + 1.96*em_v1$log_expected_se) - 1)) |> 
  mutate(across(-date, function(x) round(x,.1)))
monthly_est_v1<-as.data.frame(mo_tb_v1)

mo_tb_v2 <- tibble(date = em_v2$date, observed = em_v2$observed, 
                   expected_deaths = em_v2$expected,
                   expected_lower95 = (exp(log(em_v2$expected) - 1.96*em_v2$log_expected_se) - 1),
                   expected_upper95 = (exp(log(em_v2$expected) + 1.96*em_v2$log_expected_se) - 1),
                   excess = em_v2$expected*(exp(em_v2$fitted) - 1), 
                   excess_lower95 = em_v2$expected*(exp(em_v2$fitted - 1.96*em_v2$log_expected_se) - 1),
                   excess_upper95 = em_v2$expected*(exp(em_v2$fitted + 1.96*em_v2$log_expected_se) - 1)) |> 
  mutate(across(-date, function(x) round(x,.1)))
monthly_est_v2<-as.data.frame(mo_tb_v2)

write.csv(yearly_est_v1,file=paste0(outdir,"harvard_yearly.csv"), row.names = FALSE)
write.csv(yearly_est_v2,file=paste0(outdir,"harvard_sensitivity_yearly.csv"), row.names = FALSE)
write.csv(monthly_est_v1,file=paste0(outdir,"harvard_monthly.csv"), row.names = FALSE)
write.csv(monthly_est_v2,file=paste0(outdir,"harvard_sensitivity_monthly.csv"), row.names = FALSE)

# ------------------------ IHME  ------------------------ 
## We referenced the source R and Python scripts from 
##  01_data_prep, 02_data_processing, and 
##  03_ensemble_excess_model from the following GitHub repository: 
##  https://github.com/ihmeuw-demographics/publication_covid_em

## We translated the Python code to R for the four base models, 
##  which use regmod in Python in the source script: 02b_fit_regmod_model.py.

set.seed(123)
library(splines)

dat_mo_f4 <- dat_monthly %>%
  rename(
    total_deaths = outcome,
    year_month   = date_yr_mnth
  ) %>%
  mutate(
    year_month = as.Date(year_month),
    month      = as.numeric(format(year_month, "%m")),
    year       = as.numeric(format(year_month, "%Y")),
    log_pop    = log(population)
  )

# Create index starting in 2014-01
dat_mo_f4$chron_index <- (dat_mo_f4$year - 2014) * 12 + dat_mo_f4$month


### --- Weight models ---


# --- Set up knot placements for the weight models (evaluation window March–Dec 2019) ---

t0_wt <- dat_mo_f4 %>%
  filter(year_month == as.Date("2019-03-01")) %>%
  pull(chron_index) %>%
  unique()
stopifnot(length(t0_wt) == 1)

knot_6_wt  <- t0_wt - 6
knot_12_wt <- t0_wt - 12
knot_18_wt <- t0_wt - 18
knot_24_wt <- t0_wt - 24

# --- First layer of estimation: seasonal model ---

seasonal_mod <- glm(
  formula = total_deaths ~ bs(month, degree = 3) + log_pop,
  data    = dat_mo_f4 %>% filter(year_month < as.Date("2019-03-01")),
  family  = poisson(link = "log")
)

dat_mo_f4$preds <- predict(seasonal_mod, newdata = dat_mo_f4)  # log-scale linear predictor

# --- Alternative knot placements for *weight* models ---

weight_mod_1 <- glm(
  total_deaths ~ preds + bs(chron_index, degree = 1, knots = knot_6_wt),
  data   = dat_mo_f4 %>% filter(year_month < as.Date("2019-03-01")),
  family = poisson(link = "log")
)

weight_mod_2 <- glm(
  total_deaths ~ preds + bs(chron_index, degree = 1, knots = knot_12_wt),
  data   = dat_mo_f4 %>% filter(year_month < as.Date("2019-03-01")),
  family = poisson(link = "log")
)

weight_mod_3 <- glm(
  total_deaths ~ preds + bs(chron_index, degree = 1, knots = knot_18_wt),
  data   = dat_mo_f4 %>% filter(year_month < as.Date("2019-03-01")),
  family = poisson(link = "log")
)

weight_mod_4 <- glm(
  total_deaths ~ preds + bs(chron_index, degree = 1, knots = knot_24_wt),
  data   = dat_mo_f4 %>% filter(year_month < as.Date("2019-03-01")),
  family = poisson(link = "log")
)

weight_mod_5 <- glm(
  total_deaths ~ factor(month) + factor(year),
  data   = dat_mo_f4 %>% filter(year_month < as.Date("2019-03-01")),
  family = poisson(link = "log")
)

# --- RMSE-based weights using March–Dec 2019 ---

eval_range <- dat_mo_f4 %>%
  filter(year_month > as.Date("2019-02-01"),
         year_month < as.Date("2020-01-01")) %>%
  arrange(year_month)

true_2019 <- eval_range$total_deaths

pred_1 <- predict(weight_mod_1, newdata = eval_range)  # log-scale
pred_2 <- predict(weight_mod_2, newdata = eval_range)
pred_3 <- predict(weight_mod_3, newdata = eval_range)
pred_4 <- predict(weight_mod_4, newdata = eval_range)
pred_5 <- predict(weight_mod_5, newdata = eval_range)

# Last-year-carried-forward
ly_range <- dat_mo_f4 %>%
  filter(year_month > as.Date("2018-02-01"),
         year_month < as.Date("2019-01-01")) %>%
  arrange(year_month)

pred_ly_eval <- ly_range$total_deaths  # length 10 (Mar–Dec)

stopifnot(length(pred_ly_eval) == length(true_2019))

# Fixed population used in RMSE scaling 
US_pop <- 328239523

get_weight <- function(pred_log) {
  mu <- exp(pred_log)
  errors <- mu - true_2019
  rmse <- sqrt(sum((errors / US_pop)^2) / length(mu))
  1 / (rmse^2)
}

# Last-year RMSE
errors_ly <- pred_ly_eval - true_2019
rmse_ly   <- sqrt(sum((errors_ly / US_pop)^2) / length(pred_ly_eval))
w_ly      <- 1 / (rmse_ly^2)

weights <- c(
  w_ly,
  get_weight(pred_1),
  get_weight(pred_2),
  get_weight(pred_3),
  get_weight(pred_4),
  get_weight(pred_5)
)

weight_vec <- weights / sum(weights)


### --- Final models for 2020+ predictions ---


# Knot placements for final baseline models (relative to first COVID month)
t0_final <- dat_mo_f4 %>%
  filter(year_month == as.Date("2020-01-01")) %>%
  pull(chron_index) %>%
  unique()
stopifnot(length(t0_final) == 1)

knot_6_final  <- t0_final - 6
knot_12_final <- t0_final - 12
knot_18_final <- t0_final - 18
knot_24_final <- t0_final - 24

# Train on all pre-2020 months
train_final <- dat_mo_f4 %>%
  filter(year_month < as.Date("2020-01-01"))

final_mod_1 <- glm(
  total_deaths ~ preds + bs(chron_index, degree = 1, knots = knot_6_final),
  data   = train_final,
  family = poisson(link = "log")
)

final_mod_2 <- glm(
  total_deaths ~ preds + bs(chron_index, degree = 1, knots = knot_12_final),
  data   = train_final,
  family = poisson(link = "log")
)

final_mod_3 <- glm(
  total_deaths ~ preds + bs(chron_index, degree = 1, knots = knot_18_final),
  data   = train_final,
  family = poisson(link = "log")
)

final_mod_4 <- glm(
  total_deaths ~ preds + bs(chron_index, degree = 1, knots = knot_24_final),
  data   = train_final,
  family = poisson(link = "log")
)

final_mod_5 <- glm(
  total_deaths ~ factor(month) + year,
  data   = dat_mo_f4 %>% filter(year_month < as.Date("2020-01-01")),
  family = poisson(link = "log")
)


### --- Prediction dataset (2020+) ---


pred_data <- dat_mo_f4 %>%
  filter(year_month >= as.Date("2020-01-01")) %>%
  arrange(year_month)

# Ensure year, month are present as integers
pred_data <- pred_data %>%
  mutate(
    year  = year(year_month),
    month = month(year_month)
  )

# Build last-year-carried-forward: for each (year, month), take total_deaths from previous year same month
ly_lookup <- dat_mo_f4 %>%
  mutate(
    year  = year(year_month),
    month = month(year_month)
  ) %>%
  dplyr::select(year, month, total_deaths)

pred_data <- pred_data %>%
  left_join(
    ly_lookup %>%
      mutate(year = year + 1) %>%    # deaths in 2019 -> baseline for 2020, etc.
      rename(ly_deaths = total_deaths),
    by = c("year", "month")
  )

pred_ly <- pred_data$ly_deaths
n_months <- nrow(pred_data)


### --- Draw-based prediction from GLMs ---


n_draws <- 1000

predict_draws_counts <- function(fit, newdata, n_draws = 1000) {
  # model matrix based on the fitted model's terms
  X <- model.matrix(terms(fit), newdata)
  beta_hat <- coef(fit)
  V <- vcov(fit)
  
  beta_draws <- MASS::mvrnorm(n = n_draws, mu = beta_hat, Sigma = V)  
  eta_draws  <- X %*% t(beta_draws)                                   
  mu_draws   <- exp(eta_draws)                                        
  mu_draws
}

mu1_draws <- predict_draws_counts(final_mod_1, pred_data, n_draws)
mu2_draws <- predict_draws_counts(final_mod_2, pred_data, n_draws)
mu3_draws <- predict_draws_counts(final_mod_3, pred_data, n_draws)
mu4_draws <- predict_draws_counts(final_mod_4, pred_data, n_draws)
mu5_draws <- predict_draws_counts(final_mod_5, pred_data, n_draws)

# Last-year baseline treated as fixed (no parameter uncertainty)
mu_ly_draws <- matrix(pred_ly, nrow = n_months, ncol = n_draws, byrow = FALSE)


### --- Ensemble draws, according to the developer's approach ---


stopifnot(length(weight_vec) == 6)

mu_ens_draws <- (
  weight_vec[1] * mu_ly_draws +
    weight_vec[2] * mu1_draws +
    weight_vec[3] * mu2_draws +
    weight_vec[4] * mu3_draws +
    weight_vec[5] * mu4_draws +
    weight_vec[6] * mu5_draws
)

obs_monthly <- pred_data$total_deaths

# Excess = observed - expected
excess_draws <- matrix(obs_monthly, nrow = n_months, ncol = n_draws, byrow = FALSE) - mu_ens_draws


### --- Monthly summaries (expected + excess) ---


month_summary <- t(apply(mu_ens_draws, 1, function(v) {
  c(
    expected_deaths  = mean(v),
    expected_lower = quantile(v, 0.025),
    expected_upper = quantile(v, 0.975)
  )
}))

excess_month_summary <- t(apply(excess_draws, 1, function(v) {
  c(
    excess  = mean(v),
    excess_lower = quantile(v, 0.025),
    excess_upper = quantile(v, 0.975)
  )
}))

ihme_monthly <- cbind(
  pred_data %>% 
    dplyr::select(year_month, year, month, total_deaths) %>% rename(observed = total_deaths),
  as.data.frame(month_summary),
  as.data.frame(excess_month_summary)
) %>%
  arrange(year_month)


### --- Yearly summaries (expected + excess) ---


years <- sort(unique(pred_data$year))

year_summ_list_expected <- lapply(years, function(yy) {
  idx <- which(pred_data$year == yy)
  yearly_draws <- colSums(mu_ens_draws[idx, , drop = FALSE])
  c(
    year           = yy,
    expected_deaths  = mean(yearly_draws),
    expected_lower = quantile(yearly_draws, 0.025),
    expected_upper = quantile(yearly_draws, 0.975)
  )
})

year_summ_list_excess <- lapply(years, function(yy) {
  idx <- which(pred_data$year == yy)
  yearly_excess_draws <- colSums(excess_draws[idx, , drop = FALSE])
  c(
    year         = yy,
    excess  = mean(yearly_excess_draws),
    excess_lower = quantile(yearly_excess_draws, 0.025),
    excess_upper = quantile(yearly_excess_draws, 0.975)
  )
})

year_expected <- do.call(rbind, year_summ_list_expected) %>% as.data.frame()
year_excess   <- do.call(rbind, year_summ_list_excess)   %>% as.data.frame()

ihme_yearly <- year_expected %>%
  left_join(year_excess, by = "year") %>%
  mutate(year = as.integer(year)) %>%
  arrange(year)

write.csv(ihme_yearly,file=paste0(outdir,"ihme_yearly.csv"), row.names = FALSE)
write.csv(ihme_monthly,file=paste0(outdir,"ihme_monthly.csv"), row.names = FALSE)
