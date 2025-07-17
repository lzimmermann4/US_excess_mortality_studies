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
dat_monthly <- read.csv(paste0(datadir,"monthly data.csv"), sep=',', header=TRUE)
dat_weekly <- read.csv(paste0(datadir,"weekly data.csv"), sep=',', header=TRUE)

# Read in monthly estimates that were generated in Python for IHME
monthly_regmod6 <- read.csv(paste0(indir,"regmod_ts_6_months.csv"), sep=',', header=TRUE)
monthly_regmod12 <- read.csv(paste0(indir,"regmod_ts_12_months.csv"), sep=',', header=TRUE)
monthly_regmod18 <- read.csv(paste0(indir,"regmod_ts_18_months.csv"), sep=',', header=TRUE)
monthly_regmod24 <- read.csv(paste0(indir,"regmod_ts_24_months.csv"), sep=',', header=TRUE)


# Read in estimates from the CDC
cdc_dat <- read.csv(paste0(datadir,"Provisional_COVID-19_Death_Counts_by_Week_Ending_Date_and_State_20250418.csv"), sep=',', header=TRUE) %>%
  dplyr::rename(observed = Total.Deaths)%>%
  dplyr::mutate(year = as.numeric(Year),
                month = Month,
                pct_of_expected_deaths = Percent.of.Expected.Deaths/100,
                expected_deaths = observed / pct_of_expected_deaths,
                excess = observed - expected_deaths
  )%>%
  dplyr::filter(year <= 2024)

# Read in estimates from Karlinsky & Kobak 
wmd_excess <- read.csv(paste0(indir,"WMD/excess-mortality-timeseries.csv"), sep=',', header=TRUE) %>%
  filter(country_name == "United States")

# ------------------------ CDC ------------------------
# Subset to monthly nationwide estimates from the CDC
cdc_monthly <- cdc_dat %>%
  dplyr::filter((State == "United States") & (Group=="By Month")) %>%
  dplyr::mutate( month_year=format(as.Date(paste0(year,"-",month,"-","01")), "%Y-%m"))%>%
  dplyr::select(month_year,year,month,expected_deaths,excess)

# Subset to yearly nationwide estimates from the CDC
cdc_yearly <- cdc_dat %>%
  dplyr::filter((State == "United States") & (Group=="By Year")) %>%
  dplyr::select(year,expected_deaths,excess)

write.csv(cdc_yearly,file=paste0(outdir,"cdc_yearly.csv"), row.names = FALSE)
write.csv(cdc_monthly,file=paste0(outdir,"cdc_monthly.csv"), row.names = FALSE)

# ------------------------ World Mortality Dataset ------------------------
dat_wk_f0 <- dat_weekly %>% group_by(year) %>% mutate(time = ifelse(year == 2020, row_number() + 1, row_number())) %>% ungroup() %>%  rename(observed = outcome)

# Merge observed deaths onto excess deaths and back-out expected deaths
wmd_weekly <- merge(wmd_excess, dat_wk_f0, by=c("year","time"),all.x=T) %>%
  dplyr::mutate(expected = observed - excess.deaths) %>%
  dplyr::rename(excess = excess.deaths) %>%
  dplyr::select(date, observed, excess, expected)

wmd_yearly <- wmd_excess %>%
  group_by(year) %>%
  summarise(excess = sum(excess.deaths)) %>%
  ungroup()

write.csv(wmd_yearly,file=paste0(outdir,"wmd_yearly.csv"), row.names = FALSE)
write.csv(wmd_weekly,file=paste0(outdir,"wmd_weekly.csv"), row.names = FALSE)

# ------------------------ Msemburi et. al  ------------------------ 
# Link to developer's repository: https://github.com/WHOexcessc19/Codebase
# The below code is adapted from spline_model_monthly.R and Final_Sampling.R, for the US case study.

dat_mo_f1 <- dat_monthly %>% rename(observed = outcome)
dat_wk_f1 <- dat_weekly %>% rename(observed = outcome)
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
      #theta.i=theta.i,
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
dat_wk_f2 <- dat_weekly %>% rename(observed = outcome)
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
results_monthly <- results_monthly %>%
  mutate(expected_lower95 = expected_deaths - (1.96*expected_deaths_se),
         expected_upper95 = expected_deaths + (1.96*expected_deaths_se),
         excess_lower95 = excess - (1.96*expected_deaths_se),
         excess_upper95 = excess + (1.96*expected_deaths_se))%>%
  rename(observed = total_deaths) %>%
  select(year,month,observed,expected_deaths, expected_deaths_se, expected_lower95, expected_upper95, excess, excess_lower95, excess_upper95)

## Obtain yearly counts & UIs 
dat20 <- subset(results_monthly, year==2020)
dat21 <- subset(results_monthly, year==2021)
dat22 <- subset(results_monthly, year==2022)
dat23 <- subset(results_monthly, year==2023)
dat24 <- subset(results_monthly, year==2024)
monthly_variances20 <- dat20$expected_deaths_se^2
monthly_variances21 <- dat21$expected_deaths_se^2
monthly_variances22 <- dat22$expected_deaths_se^2
monthly_variances23 <- dat23$expected_deaths_se^2
monthly_variances24 <- dat24$expected_deaths_se^2

# Obtain covariance between each pair of months
library(tidyr)
wide20 <- dat_wk_f2 %>%
  select(year, month, week, observed) %>%
  filter(year==2020) %>%
  pivot_wider(names_from = month, values_from = observed) %>%
  select(-c(year, week)) 
cov_matrix20 <- cov(data.matrix(wide20), use ="complete.obs")
total_variance20 <- sum(monthly_variances20) + 2*sum(cov_matrix20[lower.tri(cov_matrix20)])
se20 <- sqrt(total_variance20)

wide21 <- dat_wk_f2 %>%
  select(year, month, week, observed) %>%
  filter(year==2021) %>%
  pivot_wider(names_from = month, values_from = observed) %>%
  select(-c(year, week)) 
cov_matrix21 <- cov(data.matrix(wide21), use ="complete.obs")
total_variance21 <- sum(monthly_variances21) + 2*sum(cov_matrix21[lower.tri(cov_matrix21)])
se21<-sqrt(total_variance21)

wide22 <- dat_wk_f2 %>%
  select(year, month, week, observed) %>%
  filter(year==2022) %>%
  pivot_wider(names_from = month, values_from = observed) %>%
  select(-c(year, week)) 
cov_matrix22 <- cov(data.matrix(wide22), use ="complete.obs")
total_variance22 <- sum(monthly_variances22) + 2*sum(cov_matrix22[lower.tri(cov_matrix22)])
se22<-sqrt(total_variance22)

wide23 <- dat_wk_f2 %>%
  select(year, month, week, observed) %>%
  filter(year==2023) %>%
  pivot_wider(names_from = month, values_from = observed) %>%
  select(-c(year, week)) 
cov_matrix23 <- cov(data.matrix(wide23), use ="complete.obs")
total_variance23 <- sum(monthly_variances23) + 2*sum(cov_matrix23[lower.tri(cov_matrix23)])
se23<-sqrt(total_variance23)

wide24 <- dat_wk_f2 %>%
  select(year, month, week, observed) %>%
  filter(year==2024) %>%
  pivot_wider(names_from = month, values_from = observed) %>%
  select(-c(year, week)) 
cov_matrix24 <- cov(data.matrix(wide24), use ="complete.obs")
total_variance24 <- sum(monthly_variances24) + 2*sum(cov_matrix24[lower.tri(cov_matrix24)])
se24<-sqrt(total_variance24)

economist_yearly <- results_monthly%>%
  group_by(year)%>%
  summarise(expected_deaths = sum(expected_deaths),
            excess = sum(excess))%>%
  ungroup()%>%
  mutate(expected_deaths_se = case_when(year == 2020 ~ se20,
                                        year == 2021 ~ se21,
                                        year == 2022 ~ se22,
                                        year == 2023 ~ se23,
                                        year == 2024 ~ se24),
         expected_lower95 = expected_deaths - (1.96*expected_deaths_se),
         expected_upper95 = expected_deaths + (1.96*expected_deaths_se),
         excess_lower95 = excess - (1.96*expected_deaths_se),
         excess_upper95 = excess + (1.96*expected_deaths_se))%>%
  relocate(excess,.before = excess_lower95)

write.csv(economist_yearly,file=paste0(outdir,"economist_yearly.csv"), row.names = FALSE)
write.csv(economist_monthly,file=paste0(outdir,"economist_monthly.csv"), row.names = FALSE)


# ------------------------ Acosta and Irizarry  ------------------------ 
library(excessmort)

dat_wk_f3<-dat_weekly %>%
  mutate(date = as.Date(date, format="%Y-%m-%d"))
  
### Compute expected counts
## Main results: using reference period of 2014-2019 and using extrapolate = FALSE
# Remove COVID years and periods with bad flu season
exclude_dates_v1 <- c(
  seq(make_date(2015, 1, 1), make_date(2015, 2, 15), by = "day"),
  seq(make_date(2017, 1, 1), make_date(2017, 2, 15), by = "day"),
  seq(make_date(2018, 1, 1), make_date(2018, 2, 15), by = "day"),
  seq(make_date(2020, 1, 1), make_date(2024, 12, 31), by = "day"))

# Flexible model to account for tappering off trend
expected_v1 <- excessmort::compute_expected(dat_wk_f3, exclude = exclude_dates_v1, harmonics = 2, trend.knots.per.year = 1/4, extrapolate=FALSE)

# Estimate excess for each year and also the curve
em_v1 <- excess_model(expected_v1, start = make_date(2014, 1, 1), end = make_date(2024, 12, 31), exclude = exclude_dates, 
                      knots.per.year = 18, 
                      intervals = list(seq(make_date(2020, 1, 1), make_date(2020, 12, 31), by = "day"),
                                       seq(make_date(2021, 1, 1), make_date(2021, 12, 31), by = "day"),
                                       seq(make_date(2022, 1, 1), make_date(2022, 12, 31), by = "day"),
                                       seq(make_date(2023, 1, 1), make_date(2023, 12, 31), by = "day"),
                                       seq(make_date(2024, 1, 1), make_date(2024, 12, 31), by = "day")))

## Sensitivity results: using reference period of 2017-2019 and using extrapolate = TRUE for 2017-2019
dat_weekly_frm2017 <- dat_wk_f3 %>%
  filter(date >= as_date("2017-01-07"))

# Exclude COVID and years with bad flu season
exclude_dates_v2 <- c(
  seq(make_date(2017, 1, 1), make_date(2017, 2, 15), by = "day"),
  seq(make_date(2018, 1, 1), make_date(2018, 2, 15), by = "day"),
  seq(make_date(2020, 1, 1), make_date(2024, 12, 31), by = "day"))

# Extending linear trend from 2017-2019 for extrapolation
expected_v2 <- excessmort::compute_expected(dat_weekly_frm2017, exclude = exclude_dates_v2, harmonics = 2, trend.knots.per.year = 1/4, extrapolate=TRUE)

# Estimate excess for each year and also the curve
em_v2 <- excess_model(expected_v2, start = make_date(2017, 1, 1), end = make_date(2024, 12, 31), exclude = exclude_dates, 
                      knots.per.year = 18, 
                      intervals = list(seq(make_date(2020, 1, 1), make_date(2020, 12, 31), by = "day"),
                                       seq(make_date(2021, 1, 1), make_date(2021, 12, 31), by = "day"),
                                       seq(make_date(2022, 1, 1), make_date(2022, 12, 31), by = "day"),
                                       seq(make_date(2023, 1, 1), make_date(2023, 12, 31), by = "day"),
                                       seq(make_date(2024, 1, 1), make_date(2024, 12, 31), by = "day")))

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

## Expected and excess deaths for each week
wk_tb_v1 <- tibble(date = em_v1$date, observed = em_v1$observed, 
                   expected = em_v1$expected,
                   expected_lower95 = (exp(log(em_v1$expected) - 1.96*em_v1$log_expected_se) - 1),
                   expected_upper95 = (exp(log(em_v1$expected) + 1.96*em_v1$log_expected_se) - 1),
                   excess = em_v1$expected*(exp(em_v1$fitted) - 1), 
                   excess_lower95 = em_v1$expected*(exp(em_v1$fitted - 1.96*em_v1$log_expected_se) - 1),
                   excess_upper95 = em_v1$expected*(exp(em_v1$fitted + 1.96*em_v1$log_expected_se) - 1)) |> 
  mutate(across(-date, function(x) round(x,.1)))
weekly_est_v1<-as.data.frame(wk_tb_v1)

wk_tb_v2 <- tibble(date = em_v2$date, observed = em_v2$observed, 
                   expected = em_v2$expected,
                   expected_lower95 = (exp(log(em_v2$expected) - 1.96*em_v2$log_expected_se) - 1),
                   expected_upper95 = (exp(log(em_v2$expected) + 1.96*em_v2$log_expected_se) - 1),
                   excess = em_v2$expected*(exp(em_v2$fitted) - 1), 
                   excess_lower95 = em_v2$expected*(exp(em_v2$fitted - 1.96*em_v2$log_expected_se) - 1),
                   excess_upper95 = em_v2$expected*(exp(em_v2$fitted + 1.96*em_v2$log_expected_se) - 1)) |> 
  mutate(across(-date, function(x) round(x,.1)))
weekly_est_v2<-as.data.frame(wk_tb_v2)

write.csv(yearly_est_v1,file=paste0(outdir,"harvard_yearly.csv"), row.names = FALSE)
write.csv(yearly_est_v2,file=paste0(outdir,"harvard_sensitivity_yearly.csv"), row.names = FALSE)
write.csv(weekly_est_v1,file=paste0(outdir,"harvard_weekly.csv"), row.names = FALSE)
write.csv(weekly_est_v2,file=paste0(outdir,"harvard_sensitivity_weekly.csv"), row.names = FALSE)

# ------------------------ IHME  ------------------------ 
# Link to developer's repository: https://github.com/ihmeuw-demographics/publication_covid_em
# The below code is adapted from 01_data_prep.R, 02_data_processing.R, and 03_ensemble_excess_model.R, for the US case study.

## Read in estimates from four base models, which use regmod in Python. 
# We adapted the code in 02b_fit_regmod_model.py to obtain these estimates.

## Obtain estimates for two of the base models
library(argparse)
library(arrow)
library(assertable)
library(data.table)
library(demUtils)
library(fs)

# Set seed
set.seed(123)
n_draws = 1
if (n_draws == 1) n_draws <- 1000

ihme_loc <- "USA"

# For weekly estimates
dat_wk_f4 <- dat_weekly %>%
  rename(deaths = outcome) %>%
  mutate(ihme_loc = ihme_loc,
         age_start = 0, age_end = 125,
         age_name = "0 to 125",
         sex = "all",
         time_unit = "week",
         date =  format(as.Date(date, format="%Y-%m-%d"),"%m/%d/%Y"),
         month = format(as.Date(date, format="%m/%d/%Y"),"%m"),
         month = as.numeric(month),
         week = as.Date(date, format="%m/%d/%Y"),
         year = format(as.Date(date, format="%m/%d/%Y"),"%Y"),
         year = as.numeric(year)) %>%
  group_by(year) %>%
  mutate(week = row_number())%>%
  ungroup()%>%
  mutate(year_start = as.integer(year), week_start = week, time_start = week_start)

# For monthly estimates
dat_mo_f4 <- dat_monthly %>%
  rename(deaths = outcome,
         date = date_yr_mnth) %>%
  mutate(ihme_loc = ihme_loc,
         age_start = 0, age_end = 125,
         age_name = "0 to 125",
         sex = "all",
         time_unit = "month",
         date =  format(as.Date(date, format="%Y-%m-%d"),"%m/%d/%Y"),
         month = format(as.Date(date, format="%m/%d/%Y"),"%m"),
         month = as.numeric(month),
         year = format(as.Date(date, format="%m/%d/%Y"),"%Y"),
         year = as.numeric(year),
         offset = log(population)) %>%
  group_by(year) %>%
  mutate(month = row_number())%>%
  ungroup()%>%
  mutate(year_start = as.integer(year), month_start = month, time_start = month_start)


run_model <- function(model_data, model_type, ss, predict_data, all_age) {
  
  message(paste0(Sys.time(), " | Fit ", model_type, " model for sex ", ss))
  
  # subset on age and sex
  model_data <- copy(model_data)
  predict_data <- copy(predict_data)
  model_data <- model_data %>% filter(sex == ss)
  predict_data <- predict_data %>% filter(sex == ss)
  if (all_age) {
    model_data <- model_data %>% filter(age_name == "0 to 125")
    predict_data <- predict_data %>% filter(age_name == "0 to 125")
    one_age_group <- TRUE
  } else {
    model_data <- model_data %>% filter(age_name != "0 to 125")
    predict_data <- predict_data %>% filter(age_name != "0 to 125")
    # check if more than 1 detailed age group
    one_age_group <- length(unique(model_data$age_name)) == 1
  }
  
  # create formula. only use age effect if age-specific
  if (model_type == "logistic") {
    model_data[, death_rate := deaths / population]
    form <- paste0("death_rate ~ factor(time_start) + factor(year_start) + ",
                   ifelse(one_age_group, "", " factor(age_name)"))
  } else {
    form <- paste0("deaths ~ factor(time_start) + factor(year_start) + ",
                   ifelse(one_age_group, "", " factor(age_name) + "),
                   "offset(log(population))")
  }
  message(paste0("Using formula: ", form))
  form <- as.formula(form)
  
  # fit Poisson model
  if (model_type == "poisson") {
    fit <- stats::glm(
      formula = form,
      family = poisson(link = log),
      data = model_data
    )
  } else if (model_type == "neg_binom") {
    fit <- MASS::glm.nb(
      formula = form,
      link = log,
      data = model_data
    )
  } else if (model_type == "logistic") {
    fit <- stats::glm(
      formula = form,
      data = model_data,
      family = "binomial"
    )
  } else {
    stop(paste0("'", model_type, "' is an unsupported model type."))
  }
  
  # get fixed effects draws from variance covariance matrix
  fe <- stats::coef(fit)
  fe_vcov <- stats::vcov(fit)
  
  fe <- round(fe, 10)
  fe_vcov <- round(fe_vcov, 10)
  
  fe_draws <- MASS::mvrnorm(n = n_draws, mu = fe, Sigma = fe_vcov)
  fe_draws <- as.data.table(fe_draws)
  fe_draws$draw <- c(1:n_draws)
  
  # reshape fixed effects draws
  fe_draws <- melt(fe_draws, id.vars = "draw")
  fe_draws[, type := tstrsplit(variable, "\\(|\\)", keep = 2)]
  fe_draws[type != "Intercept", group := tstrsplit(variable, "\\(|\\)", keep = 3)]
  fe_draws_intercept <- fe_draws[
    type == "Intercept",
    list(draw, fe_intercept = value)
  ]
  fe_draws_time <- fe_draws[
    type == "time_start",
    list(draw, time_start = as.numeric(group), fe_time = value)
  ]
  fe_draws_year <- fe_draws[
    type == "year_start",
    list(draw, year_start = as.numeric(group), fe_year = value)
  ]
  fe_draws_age <- fe_draws[
    type == "age_name",
    list(draw, age_name = group, fe_age = value)
  ]
  
  # fill in reference groups w/ fixed effects zeros
  ref_time <- setdiff(unique(model_data$time_start), unique(fe_draws_time$time_start))
  ref_yr <- setdiff(unique(model_data$year_start), unique(fe_draws_year$year_start))
  ref_age <- setdiff(unique(model_data$age_name), unique(fe_draws_age$age_name))
  fe_draws_time <- rbind(
    fe_draws_time,
    data.table::CJ(draw = 1:n_draws, time_start = ref_time, fe_time = 0)
  )
  fe_draws_year <- rbind(
    fe_draws_year,
    data.table::CJ(draw = 1:n_draws, year_start = ref_yr, fe_year = 0)
  )
  fe_draws_age <- rbind(
    fe_draws_age,
    data.table::CJ(draw = 1:n_draws, age_name = ref_age, fe_age = 0)
  )
  
  # extrapolate year effect to any prediction years not in model fit dataset
  fill_yrs <- setdiff(
    unique(predict_data$year_start),
    unique(fe_draws_year$year_start)
  )
  
  if (length(fill_yrs) > 0) {
    if (max(fill_yrs) > max(fe_draws_year$year_start)) {
      fe_draws_year <- demUtils::extrapolate(
        fe_draws_year,
        id_cols = c("draw", "year_start"),
        extrapolate_col = "year_start",
        value_col = "fe_year",
        extrapolate_vals = unique(predict_data$year_start),
        method = "uniform",
        n_groups_fit = 2
      )
    } else {
      fe_draws_year <- demUtils::interpolate(
        fe_draws_year,
        id_cols = c("draw", "year_start"),
        interpolate_col = "year_start",
        value_col = "fe_year",
        interpolate_vals = unique(predict_data$year_start)
      )
    }
  }
  
  # merge fixed effects onto predict data
  d <- merge(predict_data, fe_draws_time, by = "time_start", allow.cartesian = TRUE)
  d <- merge(d, fe_draws_year, by = c("year_start", "draw"))
  d <- merge(d, fe_draws_age, by = c("age_name", "draw"))
  d <- merge(d, fe_draws_intercept, by = c("draw"))
  
  # predict expected deaths from regression formula
  if (model_type == "logistic") {
    d[, logit_death_rate := fe_intercept + fe_time + fe_year + fe_age]
    d[, death_rate := demUtils::invlogit(logit_death_rate)]
    d[, deaths_expected := death_rate * population]
    d[, c("death_rate", "logit_death_rate") := NULL]
  } else {
    d <- d %>%
      mutate(deaths_expected = exp(fe_intercept + fe_time + fe_year + fe_age + log(population)))
  }
  
  # value = excess mortality rate
  d <- d %>% 
    mutate(value = (deaths - deaths_expected) / population,
           model_type = model_type)
  
  return(d)
}

# Loop over sex and model type (excluding 'regmod' model type)
run_model_types <- c("poisson")

draws <- rbindlist(lapply(run_model_types, function(model_type) {
  draws_mt <- rbindlist(lapply(unique(dat_mo_f4$sex), function(ss) {  
    model_data <- dat_mo_f4 %>% filter((year <= 2019) & (sex == ss)) # fit model on 2014-2019
    predict_data <- dat_mo_f4 %>% filter((sex == ss))
    
    only_all_age <- length(unique(model_data$age_name)) == 1
    
    # all-age model
    d_all_age <- run_model(
      model_data = model_data,
      model_type = model_type,
      ss = ss,
      predict_data = predict_data,
      all_age = TRUE
    )
    
    # age-specific model
    if (!only_all_age) {
      d_age_specific <- run_model(
        model_data = model_data,
        model_type = model_type,
        ss = ss,
        predict_data = predict_data,
        all_age = FALSE
      )
      d <- rbind(d_all_age, d_age_specific)
    } else {
      d <- d_all_age
    }
    return(d)
  }))
  return(draws_mt)
}))

draws_poisson_yearmonth <- draws

draws_poisson_yearmonth <- draws_poisson_yearmonth %>% 
  mutate(death_rate_expected = deaths_expected / population,
         death_rate_observed = deaths / population) %>%
  filter(year>=2020)
draws <- draws_poisson_yearmonth

write.csv(draws_poisson_yearmonth,file=paste0(outdir, "draws/","draws_poisson_yearmonth.csv"), row.names = FALSE)

library(hierarchyUtils)
library(lubridate)
library(magrittr)

ihme_loc <- "USA"
draw_id_cols <- unique(c("date", "year", "month", "ihme_loc", "draw"))

## Obtain estimates for base model: "last year model"
# use previous year's observed mortality as the "expected mortality" for
# the current year
covid_years <- c(2020, 2021, 2022, 2023, 2024)

dt_last_yr <- setDT(dat_mo_f4) 
dt_last_yr[, model_type := "previous_year"]
dt_last_yr[, location_id := "USA"]
dt_last_yr[, death_rate_observed := deaths / population] # We replaced "person_years" with "population"
setkeyv(dt_last_yr, c("location_id", "year_start", "time_start"))
dt_last_yr[
  !(year_start %in% covid_years),
  death_rate_observed_non_covid := death_rate_observed
]
dt_last_yr[,
           death_rate_expected := shift(death_rate_observed_non_covid, 1),
           by = setdiff(c("location_id", "year_start", "time_start"), "year_start")
]
dt_last_yr <- dt_last_yr %>%
  dplyr::group_by(location_id, time_start, age_name, sex) %>%
  tidyr::fill(death_rate_expected, .direction = "down") %>%
  tidyr::fill(death_rate_expected, .direction = "up") %>%
  dplyr::ungroup() %>% setDT

# Compute excess mortality rate
dt_last_yr[, value := death_rate_observed - death_rate_expected]

# Clean up
dt_last_yr <- dt_last_yr[!is.na(value)]
dt_last_yr[, `:=` (draw = 1, process_stage = "final", groupings = "estimated")]

# Expand to draws
draw_frame <- data.table(draw = 1, new_draw = 1:n_draws)
dt_last_yr <- merge(dt_last_yr, draw_frame, by = "draw", allow.cartesian = T)
dt_last_yr[, `:=` (draw = new_draw, new_draw = NULL)]
dt_last_yr <- dt_last_yr %>% filter(year>=2020)

write.csv(dt_last_yr,file=paste0(outdir, "draws/","dt_last_yr (monthly).csv"), row.names = FALSE)

# Combine with other models
draws <- bind_rows(draws_poisson_yearmonth, dt_last_yr) 

monthly_poisson_prevyear <- draws %>%
  select(model_type, year_start, death_rate_expected, death_rate_observed)

draws_regmod<-rbind(monthly_regmod6,monthly_regmod12,monthly_regmod18,monthly_regmod24)%>% 
  rename(deaths_expected = expected_deaths) %>%
  mutate(death_rate_expected = deaths_expected / population,
         death_rate_observed = deaths / population,
         value = excess_deaths / population,
         year_start = year)

# Combine with other models
draws <- bind_rows(draws, draws_regmod)

# Check yearly estimates of each model
modspecific_yearly <- draws%>%
  dplyr::group_by(model_type, year) %>%
  dplyr::summarise(observed = mean(deaths),
                   expected_deaths = mean(deaths_expected),
                   excess = mean(excess_deaths))%>%
  ungroup()

# Out-of-sample, only keep 2019 post time cutoff, note: we instead keep 2020 onwards
dt_oos_rmse <- draws[
  year_start >= 2020 &
    !is.na(death_rate_expected) & !is.na(death_rate_observed)
]

# Compute rmse
dt_oos_rmse[, sq_error := (death_rate_expected - death_rate_observed) ^2]

# Get one rmse for all location-time
rmse_oos <- dt_oos_rmse[, list(rmse = sqrt(mean(sq_error))), by = "model_type"]

# Repeat with separate rmse values by location
rmse_oos_2 <- dt_oos_rmse[,
                          list(rmse = sqrt(mean(sq_error))),
                          by = c("model_type")
]

# Ensemble ----------------------------------------------------------------

# Compute weights
weights <- copy(rmse_oos)
weights[, weight := (1/(rmse^2))]
weights[, weight := weight / (sum(weight))]

# Save weights
readr::write_csv(
  weights,
  fs::path(indir, "stage1_ensemble_weights.csv")
)

## Ensemble across base models
draws <- merge(draws, weights, by = "model_type")
table(draws$model_type)

# Rescale weights in case a model type is missing
draws[,
      weight := weight / sum(weight),
      by = setdiff(draw_id_cols, "model_type")
]

# Compute model-level means and variances
model_stats <- draws %>%
  rename(observed = deaths) %>%
  group_by(year, month, model_type, weight) %>%
  summarise(
    model_mean = mean(value),
    model_var = var(value),
    .groups = "drop"
  )

# Apply Rubin's Rule
rubins_weighted <- model_stats %>%
  group_by(year, month) %>%
  summarise(
    pooled_mean = sum(weight * model_mean),
    within_var = sum(weight * model_var),
    between_var = sum(weight * (model_mean - pooled_mean)^2),
    total_var = within_var + between_var,
    .groups = "drop"
  )

# Convert to count-scale using population
dat_mo_f4 <- dat_mo_f4 %>% rename(observed = deaths)
rubins_weighted <- merge(rubins_weighted, dat_mo_f4[,c("year","month","population","observed")], by=c("year","month"))
results_monthly <- rubins_weighted %>%
  mutate(
    excess = pooled_mean * population,
    excess_death_var = total_var * population^2,
    excess_death_sd = sqrt(excess_death_var),
    excess_lower95 = excess - 1.96 * excess_death_sd,
    excess_upper95 = excess + 1.96 * excess_death_sd,
    expected_deaths = observed - excess,
    expected_lower95 = observed - excess_upper95,
    expected_upper95 = observed - excess_lower95
  ) 

ihme_monthly <- results_monthly %>%
  mutate(month_year=format(as.Date(paste0(year,"-",month,"-","01")), "%Y-%m")) %>%
  select(month_year,year, month, observed, expected_deaths, expected_lower95, expected_upper95, excess, excess_lower95, excess_upper95)%>%
  arrange(year,month)

## Obtain yearly 
dat20 <- subset(results_monthly, year==2020)
dat21 <- subset(results_monthly, year==2021)
dat22 <- subset(results_monthly, year==2022)
dat23 <- subset(results_monthly, year==2023)
dat24 <- subset(results_monthly, year==2024)
monthly_variances20 <- dat20$excess_death_sd^2
monthly_variances21 <- dat21$excess_death_sd^2
monthly_variances22 <- dat22$excess_death_sd^2
monthly_variances23 <- dat23$excess_death_sd^2
monthly_variances24 <- dat24$excess_death_sd^2

# Obtain covariance between each pair of months
library(tidyr)
dat_wk_f2 <- dat_weekly %>% rename(observed = outcome)
wide20 <- dat_wk_f2 %>%
  select(year, month, week, observed) %>%
  filter(year==2020) %>%
  pivot_wider(names_from = month, values_from = observed) %>%
  select(-c(year, week)) 
cov_matrix20 <- cov(data.matrix(wide20), use ="complete.obs")
total_variance20 <- sum(monthly_variances20) + 2*sum(cov_matrix20[lower.tri(cov_matrix20)])
se20 <- sqrt(total_variance20)

wide21 <- dat_wk_f2 %>%
  select(year, month, week, observed) %>%
  filter(year==2021) %>%
  pivot_wider(names_from = month, values_from = observed) %>%
  select(-c(year, week)) 
cov_matrix21 <- cov(data.matrix(wide21), use ="complete.obs")
total_variance21 <- sum(monthly_variances21) + 2*sum(cov_matrix21[lower.tri(cov_matrix21)])
se21<-sqrt(total_variance21)

wide22 <- dat_wk_f2 %>%
  select(year, month, week, observed) %>%
  filter(year==2022) %>%
  pivot_wider(names_from = month, values_from = observed) %>%
  select(-c(year, week)) 
cov_matrix22 <- cov(data.matrix(wide22), use ="complete.obs")
total_variance22 <- sum(monthly_variances22) + 2*sum(cov_matrix22[lower.tri(cov_matrix22)])
se22<-sqrt(total_variance22)

wide23 <- dat_wk_f2 %>%
  select(year, month, week, observed) %>%
  filter(year==2023) %>%
  pivot_wider(names_from = month, values_from = observed) %>%
  select(-c(year, week)) 
cov_matrix23 <- cov(data.matrix(wide23), use ="complete.obs")
total_variance23 <- sum(monthly_variances23) + 2*sum(cov_matrix23[lower.tri(cov_matrix23)])
se23<-sqrt(total_variance23)

wide24 <- dat_wk_f2 %>%
  select(year, month, week, observed) %>%
  filter(year==2024) %>%
  pivot_wider(names_from = month, values_from = observed) %>%
  select(-c(year, week)) 
cov_matrix24 <- cov(data.matrix(wide24), use ="complete.obs")
total_variance24 <- sum(monthly_variances24) + 2*sum(cov_matrix24[lower.tri(cov_matrix24)])
se24<-sqrt(total_variance24)

ihme_yearly <- results_monthly%>%
  group_by(year)%>%
  summarise(expected_deaths = sum(expected_deaths),
            excess = sum(excess))%>%
  ungroup()%>%
  mutate(excess_deaths_se = case_when(year == 2020 ~ se20,
                                      year == 2021 ~ se21,
                                      year == 2022 ~ se22,
                                      year == 2023 ~ se23,
                                      year == 2024 ~ se24),
         expected_lower95 = expected_deaths - (1.96*excess_deaths_se),
         expected_upper95 = expected_deaths + (1.96*excess_deaths_se),
         excess_lower95 = excess - (1.96*excess_deaths_se),
         excess_upper95 = excess + (1.96*excess_deaths_se))%>%
  relocate(excess,.before = excess_lower95) %>%
  select(-c(excess_deaths_se))

write.csv(ihme_yearly,file=paste0(outdir,"ihme_yearly.csv"), row.names = FALSE)
write.csv(ihme_monthly,file=paste0(outdir,"ihme_monthly.csv"), row.names = FALSE)
