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


# ------------------------ Msemburi et. al  ------------------------ 
# Link to .R script: https://github.com/WHOexcessc19/Codebase/blob/main/Expected/spline_model_monthly.R

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
countries <- unique(exp_obs_bymonthyear$iso3)
acm_predictions <- data.frame(iso3 = rep(countries, each = 48),
                              year = rep(c(rep(2020, 12), 
                                           rep(2021, 12),
                                           rep(2022, 12),
                                           rep(2023, 12)), 
                                         times = length(countries)),
                              month = rep(c(1:12, 1:12, 1:12, 1:12), 
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

### Predict expected monthly mortality in 2020-2021
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
                                                rep(2023, 12)),
                                       month = c(1:12, 1:12, 1:12, 1:12)))
  whichs_pred <- which(acm_predictions$country_num == i)
  acm_predictions[whichs_pred, "expected_acm"] <- pred$fit
  acm_predictions[whichs_pred, "expected_acm_se"] <- pred$se.fit
  pred_log <- predict(annual_model,
                      se.fit = TRUE,
                      newdata = data.frame(year = c(rep(2020, 12),
                                                    rep(2021, 12),
                                                    rep(2022, 12),
                                                    rep(2023, 12)),
                                           month = c(1:12, 1:12, 1:12, 1:12)))
  acm_predictions[whichs_pred, "expected_log_acm"] <- pred_log$fit
  acm_predictions[whichs_pred, "expected_log_acm_se"] <- pred_log$se.fit
  
  # Get gamma parameters
  gamma_E <- rep(0, 48)
  gamma_delta <- rep(0, 48)
  gamma_E_nb <- rep(0, 48)
  gamma_delta_nb <- rep(0, 48)
  for(j in 1:48){
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

# Get monthly counts & CIs
acm_predictions$expected_upper95 <- acm_predictions$expected_acm+1.96*acm_predictions$expected_acm_se
acm_predictions$expected_lower95 <- acm_predictions$expected_acm-1.96*acm_predictions$expected_acm_se
acm_predictions$excess <- observed2020$observed-acm_predictions$expected_acm
acm_predictions$excess_upper95 <- acm_predictions$excess+1.96*acm_predictions$expected_acm_se
acm_predictions$excess_lower95 <- acm_predictions$excess-1.96*acm_predictions$expected_acm_se

# Get cumulative counts & UIs 
dat20 <- subset(acm_predictions, year==2020)
dat21 <- subset(acm_predictions, year==2021)
dat22 <- subset(acm_predictions, year==2022)
dat23 <- subset(acm_predictions, year==2023)
monthly_variances20 <- dat20$expected_acm_se^2
monthly_variances21 <- dat21$expected_acm_se^2
monthly_variances22 <- dat22$expected_acm_se^2
monthly_variances23 <- dat23$expected_acm_se^2

# Obtain covariance between each pair of months
library(tidyr)
wide20 <- dat_wk_f1 %>%
  select(year, month, week, observed) %>%
  filter(year==2020) %>%
  pivot_wider(names_from = month, values_from = observed) %>%
  select(-c(year, week)) 
cov_matrix20 <- cov(data.matrix(wide20), use ="complete.obs")
total_variance20 <- sum(monthly_variances20) + 2*sum(cov_matrix20[lower.tri(cov_matrix20)])
se20 <- sqrt(total_variance20)

wide21 <- dat_wk_f1 %>%
  select(year, month, week, observed) %>%
  filter(year==2021) %>%
  pivot_wider(names_from = month, values_from = observed) %>%
  select(-c(year, week)) 
cov_matrix21 <- cov(data.matrix(wide21), use ="complete.obs")
total_variance21 <- sum(monthly_variances21) + 2*sum(cov_matrix21[lower.tri(cov_matrix21)])
se21<-sqrt(total_variance21)

wide22 <- dat_wk_f1 %>%
  select(year, month, week, observed) %>%
  filter(year==2022) %>%
  pivot_wider(names_from = month, values_from = observed) %>%
  select(-c(year, week)) 
cov_matrix22 <- cov(data.matrix(wide22), use ="complete.obs")
total_variance22 <- sum(monthly_variances22) + 2*sum(cov_matrix22[lower.tri(cov_matrix22)])
se22<-sqrt(total_variance22)

wide23 <- dat_wk_f1 %>%
  select(year, month, week, observed) %>%
  filter(year==2023) %>%
  pivot_wider(names_from = month, values_from = observed) %>%
  select(-c(year, week)) 
cov_matrix23 <- cov(data.matrix(wide23), use ="complete.obs")
total_variance23 <- sum(monthly_variances23) + 2*sum(cov_matrix23[lower.tri(cov_matrix23)])
se23<-sqrt(total_variance23)

who_monthly <- acm_predictions %>%
  rename(expected_deaths = expected_acm) %>%
  mutate( month_year=format(as.Date(paste0(year,"-",month,"-","01")), "%Y-%m"))%>%
  select(month_year, year, month, expected_deaths, expected_lower95, expected_upper95, 
         excess, excess_lower95, excess_upper95)

who_yearly <- acm_predictions%>%
  group_by(year)%>%
  summarise(expected_acm = sum(expected_acm),
            excess = sum(excess))%>%
  ungroup()%>%
  mutate(expected_acm_se = case_when(year == 2020 ~ se20,
                                     year == 2021 ~ se21,
                                     year == 2022 ~ se22,
                                     year == 2023 ~ se23),
         expected_lower95 = expected_acm - (1.96*expected_acm_se),
         expected_upper95 = expected_acm + (1.96*expected_acm_se),
         excess_lower95 = excess - (1.96*expected_acm_se),
         excess_upper95 = excess + (1.96*expected_acm_se))%>%
  relocate(excess,.before = excess_lower95)

write.csv(who_yearly,file=paste0(outdir,"/who_yearly.csv"), row.names = FALSE)
write.csv(who_monthly,file=paste0(outdir,"/who_monthly.csv"), row.names = FALSE)


# ------------------------ The Economist  ------------------------ 

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
      filter(end_date <= as.Date("2019-12-31")) %>%
      mutate(total_deaths_per_day = total_deaths / days)
    expected_deaths_model <- lm(expected_deaths_formula,train_df)
    expected_deaths <- df %>% filter(year >= 2020) %>%
      mutate(expected_deaths = predict(expected_deaths_model,newdata=.) * days)
    
  } else if(frequency %in% c("monthly","quarterly")) {
    
    # Train an Economist monthly or quarterly model
    train_df <- df %>% 
      filter(end_date <= as.Date("2019-12-31")) %>%
      mutate(total_deaths_per_day = total_deaths / days)
    expected_deaths_model <- lm(expected_deaths_formula,train_df)
    expected_deaths <- df %>% filter(year >= 2020) %>%
      mutate(expected_deaths = predict(expected_deaths_model,newdata=.) * days)
    
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
  mutate(start_date = c(seq(make_date(2013, 12, 30), make_date(2023, 11, 30), by = "month")),
         end_date = c(seq(make_date(2014, 01, 30), make_date(2023, 12, 30), by = "month")),
         days = lead(as.numeric(interval(ymd(start_date),ymd(end_date)) %/% days(1))),
         days = ifelse(is.na(days),31,days))

res_mo <- get_excess_deaths(dat_mo_f2,
                            expected_deaths_model = NULL,
                            frequency = "monthly")
results_monthly <- res_mo[[2]]
results_monthly$excess <- results_monthly$total_deaths-results_monthly$expected_deaths
results_monthly <- results_monthly

results_yearly <- results_monthly %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(total_deaths = sum(total_deaths),
                   expected_deaths = sum(expected_deaths))%>%
  dplyr::ungroup()
results_yearly$excess <- results_yearly$total_deaths-results_yearly$expected_deaths

## Obtain monthly uncertainty intervals for relevant years
# obtain mean and population standard deviation over-time for each month
dat_wk_f2 <- dat_weekly %>%
  mutate(country = "United States",
         region = "United States",
         region_code = 0,
         days = 7,
         expected_deaths = "TBC") %>%
  group_by(year) %>%
  mutate(week = row_number()) %>% # week <- c(1:52,1:52,1:53,1:52,1:52,1:52,1:52,1:52,1:53,1:52)
  ungroup() %>%
  rename(total_deaths = outcome) %>%
  relocate(population, .after=total_deaths) %>%
  mutate(start_date = c(seq(make_date(2013, 12, 29), make_date(2023, 12, 24), by = "week")),
         end_date = date)

monthly_param <- dat_wk_f2%>%
  dplyr::mutate(month_year=format(as.Date(end_date), "%Y-%m"),
                year = year(as.Date(end_date)))%>%
  dplyr::group_by(month_year)%>%
  dplyr::mutate(wk_index = row_number(), # to identify the first week in the month
                num_wks = n(), # total number of weeks in the month for each year
                mean_param_t = mean(total_deaths),
                sqdiff = (total_deaths-mean_param_t)^2,
                sum_sqdiff = sum(sqdiff),
                # divide by number of weeks in month
                sd_param_t = sqrt((1/(n())) * ( sum_sqdiff )) 
  )%>%
  dplyr::ungroup()%>%
  dplyr::filter(wk_index == 1)%>%
  dplyr::mutate(m_index = row_number())

# Plot mean and standard deviation over-time to examine whether appear to be time-varying
require(ggplot)
require(cowplot)

plot1 <- ggplot(data=monthly_param,aes(x=month_year, group = 1)) +
  geom_line(aes(y = mean_param_t),linetype="solid",linewidth=1.5,color="grey70") + 
  theme_bw() +
  scale_y_continuous(limits = c(0, 100000), breaks = seq(0, 100000, by = 25000)) +
  theme(text = element_text(family = "Helvetica",size = 14),
        axis.text.x = element_text(angle = 320),
        plot.title        = element_text(size = 16, face = "bold")) + 
  ggtitle(expression(paste("Mean parameter (",mu,") over time"))) + 
  labs(x = "Month-Year", y = "mu")

plot2 <- ggplot(data=monthly_param,aes(x=month_year, group = 1)) +
  geom_line(aes(y = sd_param_t),linetype="solid",linewidth=1.5,color="grey70") + 
  theme_bw() +
  theme(text = element_text(family = "Helvetica",size = 14),
        axis.text.x = element_text(angle = 320),
        plot.title        = element_text(size = 16, face = "bold")) + 
  ggtitle(expression(paste("Standard deviation paramater (",sigma,") over time"))) + 
  labs(x = "Month-Year", y = "sigma")

png(file=paste0(outdir,"Parameters_Overtime_Plot.png"), height = 4000, width = 9000, 
    res = 300)
cowplot::plot_grid( 
  plot1, plot2,
  nrow=2,
  align="v")
dev.off()

# Create look-up of all variables except observed deaths
lookup_df <- dat_mo_f2 %>%
  dplyr::select(-c(total_deaths))

# Initiate dataframe to simulate n samples from specified
#  distribution, for each month, and estimate 1,000 times
B <- 1000
M <- 120 # number of months in sample

# initialize output
boot_excess_out <- data.frame()

for(i in 1:B){
  # Set seed for reproducibility and set w/in loop for random draws to differ
  set.seed(i) 
  
  # Initialize output
  dat_boot <- data.frame()
  boot_excess_temp <- data.frame()
  
  # Loop through months to draw from normal distribution 
  #  based on mu and sigma for specific month-year, to align with modeling in The Economist
  for(m in 1:M){
    dat_boot_m <- data.frame()
    mean_param <- subset(monthly_param, m_index==m)$mean_param_t
    sd_param <- subset(monthly_param, m_index==m)$sd_param_t
    num_wks <- subset(monthly_param, m_index==m)$num_wks
    dat_boot_m <- rnorm(n = num_wks, 
                        mean = mean_param, 
                        sd = sd_param)  
    dat_boot_m <- as.data.frame(dat_boot_m)
    
    # Iteratively append results for each iteration
    dat_boot <- bind_rows(dat_boot, dat_boot_m)
  }  
  
  dat_boot$date <- c(seq(make_date(2014, 1, 4),
                         make_date(2023, 12, 31), by =
                           "week"))
  dat_boot <- dat_boot %>%
    mutate(date = as.character(ymd(date)),
          year = as.numeric(format(as.Date(date, format="%Y-%m-%d"),"%Y")),
          month = as.numeric(format(as.Date(date, format="%Y-%m-%d"),"%m")),
          date_yr_mnth = make_date(year = year, month = month)) %>%
    group_by(date_yr_mnth,year,month)%>%
    summarise(total_deaths = sum(dat_boot_m,na.rm=T)) %>%
    ungroup

  # Merge on population counts
  dat_boot <- merge(x=dat_boot, y=lookup_df, by = c("date_yr_mnth","year","month"),all.x=T)

  # Re-fit the excess deaths for the sample of counts drawn from distribution above
  # Note also obtain estimates for 2020 and 2021 for figure
  res_temp <- get_excess_deaths(dat_boot,
                                expected_deaths_model = NULL,
                                frequency = "monthly")
  results_temp <- res_temp[[2]]
  
  boot_excess_temp <- results_temp %>% 
    dplyr::rename(excess_boot = excess_deaths,
           expected_boot = expected_deaths) %>%
    dplyr::mutate(iteration_b = i)
  
  # Iteratively append results for each iteration
  boot_excess_out <- bind_rows(boot_excess_out, boot_excess_temp)
  
}

# Compute monthly estimates for each bootstrap iteration (1,000 iter.)
boot_monthly <- boot_excess_out%>%
  dplyr::mutate(month_year=format(as.Date(date_yr_mnth), "%Y-%m"))%>%
  dplyr:: rename( monthly_expected = expected_boot,
                  monthly_excess = excess_boot)

# Obtain yearly 
boot_yearly <- boot_monthly%>%
  dplyr::group_by(iteration_b, year) %>%
  dplyr::summarise(yearly_expected = sum(monthly_expected),
                   yearly_excess = sum(monthly_excess)) %>%
  dplyr::ungroup()

# For plots, sort the point estimates across bootstrap iterations in ascending order within each month 
boot_monthly_expected <- boot_monthly %>%
  dplyr::arrange(month_year,monthly_expected) 

boot_monthly_excess <- boot_monthly %>%
  dplyr::arrange(month_year,monthly_excess) 

# Plot the sampling distribution for first two months and last two months 
#  of years 2022 and 2023, respectively
options(scipen = 999)
color1 <- "#dfcac3"
color2 <- "#ab715e"
color3 <- "#c2c7b9"
color4 <- "#788067"
color_dash <- "grey50"

hist_a <- function(dat, monthyr_no, monthyr_char, color_hex){
  ggplot(data=filter(dat,month_year==monthyr_no),aes(x= monthly_expected )) +
    geom_histogram(alpha=.6,color=color_hex, fill=color_hex) + 
    geom_vline(aes(xintercept=quantile(monthly_expected,0.025)), 
               color=color_dash, linetype="dashed") +
    geom_vline(aes(xintercept=quantile(monthly_expected,0.975)), 
               color=color_dash, linetype="dashed") +
    theme_bw() +
    scale_x_continuous() +
    theme(text = element_text(family = "Helvetica",size = 14),
          axis.text.x = element_text(angle = 320),
          plot.title        = element_text(size = 14, face = "bold")) + 
    ggtitle(paste0("Sampling distribution of monthly estimated expected deaths: ",
                   monthyr_char," (1000 iterations)")) + 
    labs(x = "Estimated Expected Deaths", y = "Frequency")
}

hist_b <- function(dat, monthyr_no, monthyr_char, color_hex){
  ggplot(data=filter(dat,month_year==monthyr_no),aes(x= monthly_excess )) +
    geom_histogram(alpha=.6,color=color_hex, fill=color_hex) + 
    geom_vline(aes(xintercept=quantile(monthly_excess,0.025)), 
               color=color_dash, linetype="dashed") +
    geom_vline(aes(xintercept=quantile(monthly_excess,0.975)), 
               color=color_dash, linetype="dashed") +
    theme_bw() +
    scale_x_continuous() +
    theme(text = element_text(family = "Helvetica",size = 14),
          axis.text.x = element_text(angle = 320),
          plot.title        = element_text(size = 14, face = "bold")) + 
    ggtitle(paste0("Sampling distribution of monthly estimated excess deaths: ",
                   monthyr_char," (1000 iterations)")) + 
    labs(x = "Estimated Excess Deaths", y = "Frequency")
}

# Expected deaths plots
hist_plot1a = hist_a(dat=boot_monthly_expected,monthyr_no="2022-01",monthyr_char="Jan 2022",color_hex=color1)
hist_plot2a = hist_a(dat=boot_monthly_expected,monthyr_no="2022-02",monthyr_char="Feb 2022",color_hex=color1)
hist_plot3a = hist_a(dat=boot_monthly_expected,monthyr_no="2022-11",monthyr_char="Nov 2022",color_hex=color1)
hist_plot4a = hist_a(dat=boot_monthly_expected,monthyr_no="2022-12",monthyr_char="Dec 2022",color_hex=color1)
hist_plot5a = hist_a(dat=boot_monthly_expected,monthyr_no="2023-01",monthyr_char="Jan 2023",color_hex=color2)
hist_plot6a = hist_a(dat=boot_monthly_expected,monthyr_no="2023-02",monthyr_char="Feb 2023",color_hex=color2)
hist_plot7a = hist_a(dat=boot_monthly_expected,monthyr_no="2023-11",monthyr_char="Nov 2023",color_hex=color2)
hist_plot8a = hist_a(dat=boot_monthly_expected,monthyr_no="2023-12",monthyr_char="Dec 2023",color_hex=color2)

# Excess deaths plots
hist_plot1b = hist_b(dat=boot_monthly_excess,monthyr_no="2022-01",monthyr_char="Jan 2022",color_hex=color3)
hist_plot2b = hist_b(dat=boot_monthly_excess,monthyr_no="2022-02",monthyr_char="Feb 2022",color_hex=color3)
hist_plot3b = hist_b(dat=boot_monthly_excess,monthyr_no="2022-11",monthyr_char="Nov 2022",color_hex=color3)
hist_plot4b = hist_b(dat=boot_monthly_excess,monthyr_no="2022-12",monthyr_char="Dec 2022",color_hex=color3)
hist_plot5b = hist_b(dat=boot_monthly_excess,monthyr_no="2023-01",monthyr_char="Jan 2023",color_hex=color4)
hist_plot6b = hist_b(dat=boot_monthly_excess,monthyr_no="2023-02",monthyr_char="Feb 2023",color_hex=color4)
hist_plot7b = hist_b(dat=boot_monthly_excess,monthyr_no="2023-11",monthyr_char="Nov 2023",color_hex=color4)
hist_plot8b = hist_b(dat=boot_monthly_excess,monthyr_no="2023-12",monthyr_char="Dec 2023",color_hex=color4)

png(file=paste0(outdir,"Expected_Histogram_Panel_Plot - monthly est.png"), height = 3400, width = 5550, res = 300)
cowplot::plot_grid( 
  hist_plot1a,
  hist_plot2a,
  hist_plot3a,
  hist_plot4a,
  hist_plot5a,
  hist_plot6a,
  hist_plot7a,
  hist_plot8a,
  nrow=4
)
dev.off()

png(file=paste0(outdir,"Excess_Histogram_Panel_Plot - monthly est.png"), height = 3400, width = 5550, res = 300)
cowplot::plot_grid( 
  hist_plot1b,
  hist_plot2b,
  hist_plot3b,
  hist_plot4b,
  hist_plot5b,
  hist_plot6b,
  hist_plot7b,
  hist_plot8b,
  nrow=4
)
dev.off()

# Calculate empirical 95% confidence intervals for excess death estimates
#  as lower bound is 2.5% percentile and upper bound is 97.5% percentile 
#  of the sampling distributions of the estimates from get_excess_deaths() 
#  for the 1,000 bootstrap repetitions for each month respectively
monthly_UI <- boot_monthly %>%
  group_by(month_year, year, month) %>%
  summarise(
            expected_lower95 = quantile(monthly_expected,0.025),
            expected_upper95 = quantile(monthly_expected,0.975),
            excess_lower95 = quantile(monthly_excess,0.025),
            excess_upper95 = quantile(monthly_excess,0.975))%>%
  ungroup()
yearly_UI <- boot_yearly %>%
  group_by(year) %>%
  summarise(
            expected_lower95 = quantile(yearly_expected,0.025),
            expected_upper95 = quantile(yearly_expected,0.975),
            excess_lower95 = quantile(yearly_excess,0.025),
            excess_upper95 = quantile(yearly_excess,0.975))%>%
  ungroup()

# Merge to point estimates
economist_monthly <- merge(results_monthly,monthly_UI,by=c("year","month")) %>%
  relocate(c(year, month), .after=month_year)%>%
  relocate(excess,.before=excess_lower95) %>% 
  select(year,month,expected_deaths, expected_lower95, expected_upper95, excess, excess_lower95, excess_upper95)
economist_yearly <- merge(results_yearly,yearly_UI,by=c("year"))%>%
  relocate(excess,.before=excess_lower95)

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
  seq(make_date(2020, 1, 1), make_date(2023, 12, 31), by = "day"))

# Flexible model to account for tappering off trend
expected_v1 <- excessmort::compute_expected(dat_wk_f3, exclude = exclude_dates_v1, harmonics = 2, trend.knots.per.year = 1/4, extrapolate=FALSE)

# Estimate excess for each year and also the curve
em_v1 <- excess_model(expected_v1, start = make_date(2014, 1, 1), end = make_date(2023, 12, 31), exclude = exclude_dates, 
                      knots.per.year = 18, 
                      intervals = list(seq(make_date(2020, 1, 1), make_date(2020, 12, 31), by = "day"),
                                       seq(make_date(2021, 1, 1), make_date(2021, 12, 31), by = "day"),
                                       seq(make_date(2022, 1, 1), make_date(2022, 12, 31), by = "day"),
                                       seq(make_date(2023, 1, 1), make_date(2023, 12, 31), by = "day")))


## Sensitivity results: using reference period of 2017-2019 and using extrapolate = TRUE for 2017-2019
dat_weekly_frm2017 <- dat_wk_f3 %>%
  filter(date >= as_date("2017-01-07"))

# Exclude COVID and years with bad flu season
exclude_dates_v2 <- c(
  seq(make_date(2017, 1, 1), make_date(2017, 2, 15), by = "day"),
  seq(make_date(2018, 1, 1), make_date(2018, 2, 15), by = "day"),
  seq(make_date(2020, 1, 1), make_date(2023, 12, 31), by = "day"))

# Extending linear trend from 2017-2019 for extrapolation
expected_v2 <- excessmort::compute_expected(dat_weekly_frm2017, exclude = exclude_dates_v2, harmonics = 2, trend.knots.per.year = 1/4, extrapolate=TRUE)

# Estimate excess for each year and also the curve
em_v2 <- excess_model(expected_v2, start = make_date(2017, 1, 1), end = make_date(2023, 12, 31), exclude = exclude_dates, 
                   knots.per.year = 18, 
                   intervals = list(seq(make_date(2020, 1, 1), make_date(2020, 12, 31), by = "day"),
                                    seq(make_date(2021, 1, 1), make_date(2021, 12, 31), by = "day"),
                                    seq(make_date(2022, 1, 1), make_date(2022, 12, 31), by = "day"),
                                    seq(make_date(2023, 1, 1), make_date(2023, 12, 31), by = "day")))

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
write.csv(yearly_est_v2,file=paste0(outdir,"harvard_yearly_sensitivity.csv"), row.names = FALSE)
write.csv(weekly_est_v1,file=paste0(outdir,"harvard_weekly.csv"), row.names = FALSE)
write.csv(weekly_est_v2,file=paste0(outdir,"harvard_weekly_sensitivity.csv"), row.names = FALSE)


# ------------------------ IHME  ------------------------ 

## We utilized the source R and Python scripts from 
##  01_data_prep, 02_data_processing, and 
##  03_ensemble_excess_model from the following GitHub repository: 
##  https://github.com/ihmeuw-demographics/publication_covid_em

### Read in estimates from four base models, which use regmod in Python. 
## We used the 02b_fit_regmod_model.py source Python script to obtain these estimates.

### The below code is adapted from the source R script: 02_fit_model.R
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


### The below code is adapted from the source R script: 03_complete_stage1.R
library(hierarchyUtils)
library(lubridate)
library(magrittr)

ihme_loc <- "USA"
draw_id_cols <- unique(c("date", "year", "month", "ihme_loc", "draw"))

## Obtain estimates for base model: "last year model"
# use previous year's observed mortality as the "expected mortality" for
# the current year
covid_years <- c(2020, 2021, 2022, 2023)

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

# Check regmod's yearly estimates
regmod_yearly <- draws_regmod%>%
  dplyr::group_by(model_type, year) %>%
  dplyr::summarise(observed = mean(deaths),
                   expected_deaths = mean(deaths_expected),
                   excess = mean(excess_deaths))%>%
  ungroup()

dt_oos <- bind_rows(draws, draws_regmod) 

# Out-of-sample, only keep 2019 post time cutoff, note: we instead keep 2020 onwards
dt_oos_rmse <- dt_oos[
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

# Combine with other models
draws <- bind_rows(draws, draws_regmod)

## Ensemble across base models
draws <- merge(draws, weights, by = "model_type")
table(draws$model_type)

# Rescale weights in case a model type is missing
draws[,
      weight := weight / sum(weight),
      by = setdiff(draw_id_cols, "model_type")
]

# Get weighted average
ensemble_draws <- draws %>%
  group_by(draw,year,month) %>%
  summarise(value = sum(value * weight),
            wt_check = sum(weight)) %>%
  ungroup()

# Merge on population and deaths to obtain absolute counts
ensemble_full <- merge(ensemble_draws,dat_monthly,by=c("year","month"), all.x=T) %>%
  filter(year>=2020) %>% 
  rename(deaths = outcome) %>%
  mutate(excess = value * population,
         expected_deaths = deaths - excess)

# Obtain empirical UIs 
ihme_monthly <- ensemble_full %>%
  mutate( month_year=format(as.Date(paste0(year,"-",month,"-","01")), "%Y-%m"))%>%
  dplyr::group_by(month_year, year, month) %>%
  dplyr::summarise(
    expected_lower95 = quantile(expected_deaths,0.025),
    expected_upper95 = quantile(expected_deaths,0.975),
    expected_deaths = mean(expected_deaths),
    excess_lower95 = quantile(excess,0.025),
    excess_upper95 = quantile(excess,0.975),
    excess = mean(excess))%>%
  ungroup()

# Obtain yearly 
ensemble_yearly <- ensemble_full%>%
  dplyr::group_by(draw, year) %>%
  dplyr::summarise(observed = sum(deaths),
                   expected_deaths = sum(expected_deaths),
                   excess = sum(excess)) %>%
  dplyr::ungroup()

ihme_yearly <- ensemble_yearly %>%
  group_by(year, observed) %>%
  summarise(
    expected_lower95 = quantile(expected_deaths,0.025),
    expected_upper95 = quantile(expected_deaths,0.975),
    excess_lower95 = quantile(excess,0.025),
    excess_upper95 = quantile(excess,0.975),
    expected_deaths = mean(expected_deaths),
    excess = mean(excess))%>%
  ungroup()

write.csv(ihme_yearly,file=paste0(outdir,"ihme_yearly.csv"), row.names = FALSE)
write.csv(ihme_monthly,file=paste0(outdir,"ihme_monthly.csv"), row.names = FALSE)

