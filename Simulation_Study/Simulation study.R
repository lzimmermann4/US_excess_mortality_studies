library(tidyverse)
library(lubridate)
library(mgcv)
library(splines)
library(excessmort)
library(reshape2)
library(ggpubr)
library(xtable)

#### For visualizations: a colorblind-friendly palette ####

okabe_ito_palette <- c(
  "#009E73", # bluish green
  "#0072B2", # blue
  "#D55E00", # vermillion
  "#CC79A7",# reddish purple
  "#E69F00", # orange
  "#56B4E9", # sky blue
  "#F0E442", # yellow
  "#0072B2", # blue
  "#D55E00", # vermillion
  "#000000"  # black
)

options(ggplot2.discrete.colour= okabe_ito_palette) ; options(ggplot2.discrete.fill= okabe_ito_palette) 


#### SIMULATION FUNCTIONS ####
# simulate_mortality is the main simulation function.
# simulate_mortality_hypo and _hypo2 change elements of the function:
# hypo1 has a younger age structure, perturbed cause-specific mortality, no summer peak, and a 5% inccrease in COVID impact.
# hypo2 has an older age structure, perturbed cause-specific mortality, and a 2% decrease in COVID impact.
# aggregate_to_monthly_totals does exactly what it says on the tin.

simulate_mortality <- function(start_year = 2014, end_year = 2024, divergence_year = 2020, seed_set = 1) {
  base_death_rate <- 7.5 # per 1000
  
  # Pre-calculate all date information
  all_years <- start_year:end_year
  n_years <- length(all_years)
  
  # Create complete date grid efficiently
  date_grid <- data.frame()
  for (year in all_years) {
    year_start <- as.Date(paste0(year, "-01-01"))
    weeks_in_year <- min(53, ceiling(as.numeric(difftime(
      as.Date(paste0(year, "-12-31")), year_start, units = "weeks"))))
    
    year_dates <- data.frame(
      year = year,
      week = 1:weeks_in_year,
      week_start_date = year_start + weeks(0:(weeks_in_year-1)),
      week_index = pmin(1:weeks_in_year, 52)
    )
    date_grid <- rbind(date_grid, year_dates)
  }
  
  # Pre-calculate population by year (from Census Bureau)
  years_elapsed <- all_years - start_year
  populations <- c(318301008, 320635163, 322941311, 324985539, 326687501, 328239523,
                   331577720, 332099760, 334017321, 336806231, 340110988)
  #populations <- base_population * (1 + annual_population_growth)^years_elapsed
  names(populations) <- all_years
  
  # Age groups in proportion
  #age_groups <- c("0-17" = 0.21, "18-34" = 0.26, "35-54" = 0.32, 
  #                "55-74" = 0.11, "75+" = 0.10)
  
  calculate_age_groups <- function(year) {
    # Base year 2014 proportions
    base_proportions <- c("0-17" = 0.21, "18-34" = 0.26, "35-54" = 0.32, 
                          "55-74" = 0.11, "75+" = 0.10)
    # Years since 2014
    years_since <- year - 2014
    
    # Youth declining due to lower birth rates
    prop_0_17 <- base_proportions["0-17"] * (0.997^years_since)  # -0.3% per year
    
    # Young adults declining due to lower birth rates in 1990s-2000s
    prop_18_34 <- base_proportions["18-34"] * (0.998^years_since)  # -0.2% per year
    
    # Middle-aged relatively stable but slight decline
    prop_35_54 <- base_proportions["35-54"] * (0.9995^years_since)  # -0.05% per year
    
    # Older adults increasing (baby boomers aging)
    prop_55_74 <- base_proportions["55-74"] * (1.015^years_since)  # +1.5% per year
    
    # Elderly increasing (life expectancy + baby boomers)
    prop_75_plus <- base_proportions["75+"] * (1.025^years_since)  # +2.5% per year
    
    # Normalize to ensure they sum to 1.0
    raw_props <- c(prop_0_17, prop_18_34, prop_35_54, prop_55_74, prop_75_plus)
    normalized_props <- raw_props / sum(raw_props)
    
    names(normalized_props) <- c("0-17", "18-34", "35-54", "55-74", "75+")
    return(normalized_props)
  }
  
  # Age-specific death rates per 1000 people
  age_death_rates <- c("0-17" = 1/1.1, "18-34" = 2.3/1.3, "35-54" = 7/1.3, 
                       "55-74" = 28.7/1.1, "75+" = 35.0/1.0)
  
  # Cause distributions 
  pre_2020_causes <- c(
    "Heart disease" = 0.235, "Cancer" = 0.225, "Accidents" = 0.052,
    "Stroke" = 0.052, "Chronic respiratory" = 0.052, "Alzheimer's" = 0.030,
    "Diabetes" = 0.022, "Influenza/Pneumonia" = 0.021, "Kidney disease" = 0.020,
    "Suicide" = 0.012, "Septicemia" = 0.015, "Chronic liver disease" = 0.014,
    "Hypertension" = 0.008, "Parkinson's disease" = 0.004, "Other" = 0.238
  )
  
  post_2020_causes_covid <- c(
    "Heart disease" = 0.23, "Cancer" = 0.21, "COVID-19" = 0.0,
    "Accidents" = 0.06, "Stroke" = 0.05, "Chronic respiratory" = 0.05,
    "Alzheimer's" = 0.04, "Diabetes" = 0.03, "Influenza/Pneumonia" = 0.02,
    "Kidney disease" = 0.02, "Suicide" = 0.015, "Other" = 0.275
  )
  
  # Seasonal variation
  seasonal_base <- 1 + 0.03 * cos(2 * pi * (1:52 - 4) / 52)  
  
  # Pre-calculate cause-specific seasonal pattern
  create_seasonal_fast <- function(winter_mult, summer_mult, peak_week = 4) {
    weeks <- 1:52
    
    winter_mult_adj <- 1 + (winter_mult - 1) * 0.6
    summer_mult_adj <- 1 + (summer_mult - 1) * 0.6
    
    amplitude <- (winter_mult_adj - summer_mult_adj) / 2
    baseline <- (winter_mult_adj + summer_mult_adj) / 2
    pmax(baseline + amplitude * cos(2 * pi * (weeks - peak_week) / 52), 0.1)
  }
  
  seasonal_patterns <- list(
    "Heart disease" = create_seasonal_fast(1.2, 0.9, 4),
    "Cancer" = rep(1.0, 52),
    "Accidents" = create_seasonal_fast(0.9, 1.2, 30),
    "Stroke" = create_seasonal_fast(1.2, 0.9, 4),
    "Chronic respiratory" = create_seasonal_fast(1.4, 0.7, 6),
    "Alzheimer's" = rep(1.0, 52),
    "Diabetes" = rep(1.0, 52),
    "Influenza/Pneumonia" = create_seasonal_fast(2.5, 0.3, 8),
    "Kidney disease" = rep(1.0, 52),
    "Suicide" = create_seasonal_fast(0.9, 1.1, 20),
    "Septicemia" = create_seasonal_fast(1.1, 0.9, 6),
    "Chronic liver disease" = rep(1.0, 52),
    "Hypertension" = create_seasonal_fast(1.1, 0.9, 4),
    "Parkinson's disease" = rep(1.0, 52),
    "COVID-19" = rep(1.0, 52),
    "Other" = rep(1.0, 52)
  )
  
  # Pre-calculate secular trends for all years
  calculate_trends <- function(years, base_year = 2014) {
    elapsed <- (years - base_year)
    trends <- data.frame(
      year = years,
      Heart_disease = 0.98^elapsed,
      Cancer = 0.985^elapsed,
      Accidents = 1.01^elapsed,
      Stroke = 0.98^elapsed,
      Chronic_respiratory = 0.99^elapsed,
      Alzheimers = 1.03^elapsed,
      Diabetes = 1.015^elapsed,
      Influenza_Pneumonia = 0.995^elapsed,
      Kidney_disease = 1.005^elapsed,
      Suicide = 1.02^elapsed,
      Septicemia = 1.01^elapsed,
      Chronic_liver_disease = 1.025^elapsed,
      Hypertension = 1.01^elapsed,
      Parkinson_disease = 1.015^elapsed,
      Other = 1.0
    )
    return(trends)
  }
  
  pre_trends <- calculate_trends(all_years)
  
  # COVID impact functions
  covid_impacts <- c("2020" = 1.05, "2021" = 1.12, "2022" = 1.05, 
                     "2023" = 0.88, "2024" = 0.85)
  covid_props <- c("2020" = 0.11, "2021" = 0.2, "2022" = 0.05, 
                   "2023" = 0.02, "2024" = 0.01)
  
  # STEP 1: Generate identical pre-2020 data for BOTH scenarios
  #cat("Generating identical pre-2020 data...\n")
  
  pre_2020_data <- data.frame()
  
  for (year in start_year:(divergence_year - 1)) {
    # Calculate age groups for this specific year
    age_groups <- calculate_age_groups(year)
    
    year_data <- date_grid[date_grid$year == year, ]
    current_pop <- populations[as.character(year)]
    base_trend <- pre_trends[pre_trends$year == year, ]
    
    # Calculate base weekly deaths
    weeks_per_year <- 365.25 / 7
    
    # Create all combinations
    year_results <- expand.grid(
      week = year_data$week,
      age_group = names(age_groups),
      cause = names(pre_2020_causes),
      stringsAsFactors = FALSE
    )
    
    # Add calculated values
    year_results$year <- year
    year_results$week_start_date <- year_data$week_start_date[match(year_results$week, year_data$week)]
    year_results$week_index <- year_data$week_index[match(year_results$week, year_data$week)]
    year_results$age_prop <- age_groups[year_results$age_group]
    year_results$age_death_rate <- age_death_rates[year_results$age_group]
    year_results$cause_prop <- pre_2020_causes[year_results$cause]
    year_results$population <- round(current_pop * year_results$age_prop)
    
    year_results$seasonal_mult <- seasonal_base[year_results$week_index]
    year_results$holiday_mult <- ifelse(year_results$week %in% c(1, 47, 52), 
                                        c(1.04, 1.02, 1.06)[match(year_results$week, c(1, 47, 52))], 1.0)  
    year_results$seasonal_mult <- year_results$seasonal_mult * year_results$holiday_mult
    
    # Cause-specific seasonal patterns
    year_results$cause_seasonal <- 1.0
    for (cause in names(pre_2020_causes)) {
      mask <- year_results$cause == cause
      if (any(mask)) {
        year_results$cause_seasonal[mask] <- seasonal_patterns[[cause]][year_results$week_index[mask]]
      }
    }
    
    # Secular trends
    cause_to_trend <- list(
      "Heart disease" = "Heart_disease", "Cancer" = "Cancer", "Accidents" = "Accidents",
      "Stroke" = "Stroke", "Chronic respiratory" = "Chronic_respiratory", 
      "Alzheimer's" = "Alzheimers", "Diabetes" = "Diabetes",
      "Influenza/Pneumonia" = "Influenza_Pneumonia", "Kidney disease" = "Kidney_disease",
      "Suicide" = "Suicide", "Septicemia" = "Septicemia", 
      "Chronic liver disease" = "Chronic_liver_disease", "Hypertension" = "Hypertension",
      "Parkinson's disease" = "Parkinson_disease", "Other" = "Other"
    )
    
    year_results$secular_trend <- 1.0
    for (cause in names(pre_2020_causes)) {
      mask <- year_results$cause == cause
      if (any(mask)) {
        trend_col <- cause_to_trend[[cause]]
        if (!is.null(trend_col)) {
          year_results$secular_trend[mask] <- base_trend[[trend_col]]
        }
      }
    }
    
    # Calculate base deaths
    year_results$base_deaths <- (year_results$population * year_results$age_death_rate / 1000 / weeks_per_year)
    year_results$adjusted_deaths <- year_results$base_deaths * year_results$cause_prop * 
      year_results$seasonal_mult * year_results$cause_seasonal * 
      year_results$secular_trend
    
    # CRITICAL: Use SAME seed for pre-2020 data to ensure identical results
    set.seed(year * 1000 + seed_set)  # Same seed regardless of scenario
    variations <- rnorm(nrow(year_results), mean = 0, sd = 0.1)  
    year_results$final_deaths <- pmax(0, round(year_results$adjusted_deaths * (1 + variations)))
    
    # Adjust deaths and variation for specific years and peak weeks
    if (year == 2015 || year == 2018) {
      peak_weeks <- 44:52
      reduction_factor <- 0.92
      variation_reduction_factor <- 0.2
      year_results$final_deaths[year_results$week %in% peak_weeks] <- year_results$final_deaths[year_results$week %in% peak_weeks] * reduction_factor
      variations[year_results$week %in% peak_weeks] <- variations[year_results$week %in% peak_weeks] * variation_reduction_factor
    }
    
    if (year == 2016 || year == 2019) {
      peak_weeks <- 1:15
      reduction_factor <- 0.92
      variation_reduction_factor <- 0.2
      year_results$final_deaths[year_results$week %in% peak_weeks] <- year_results$final_deaths[year_results$week %in% peak_weeks] * reduction_factor
      variations[year_results$week %in% peak_weeks] <- variations[year_results$week %in% peak_weeks] * variation_reduction_factor
    }
    # Keep only non-zero deaths
    year_final <- year_results[year_results$final_deaths > 0, 
                               c("year", "week", "age_group", "cause", "final_deaths", "population", "week_start_date")]
    names(year_final)[5] <- "deaths"
    
    pre_2020_data <- rbind(pre_2020_data, year_final)
  }
  
  #cat("Pre-2020 data generated. \n")
  
  # STEP 2: Generate post-2020 data for each scenario separately
  generate_post_2020 <- function(scenario_name) {
    # cat("Generating post-2020 data for", scenario_name, "scenario...\n")
    
    post_data <- data.frame()
    
    for (year in divergence_year:end_year) {
      # Calculate age groups for this specific year
      age_groups <- calculate_age_groups(year)
      
      year_data <- date_grid[date_grid$year == year, ]
      current_pop <- populations[as.character(year)]
      
      # Choose parameters based on scenario
      
      if (scenario_name == "COVID") {
        causes <- post_2020_causes_covid
        base_trend <- pre_trends[pre_trends$year == year, ]
        # Add COVID disruption to trends
        years_since <- year - 2020
        base_trend$Heart_disease <- base_trend$Heart_disease * 1.01^years_since
        base_trend$Cancer <- base_trend$Cancer * 1.005^years_since
        base_trend$Suicide <- base_trend$Suicide * 1.03^years_since
      } else {
        causes <- pre_2020_causes  # Continue pre-pandemic patterns
        base_trend <- pre_trends[pre_trends$year == year, ]
      }
      
      # Create all combinations
      year_results <- expand.grid(
        week = year_data$week,
        age_group = names(age_groups),
        cause = names(causes),
        stringsAsFactors = FALSE
      )
      
      # Add calculated values
      year_results$year <- year
      year_results$week_start_date <- year_data$week_start_date[match(year_results$week, year_data$week)]
      year_results$week_index <- year_data$week_index[match(year_results$week, year_data$week)]
      year_results$age_prop <- age_groups[year_results$age_group]
      year_results$age_death_rate <- age_death_rates[year_results$age_group]
      year_results$cause_prop <- causes[year_results$cause]
      year_results$population <- round(current_pop * year_results$age_prop)
      
      # Seasonal multipliers
      year_results$seasonal_mult <- seasonal_base[year_results$week_index] 
      year_results$holiday_mult <- ifelse(year_results$week %in% c(1, 47, 52), 
                                          c(1.04, 1.02, 1.06)[match(year_results$week, c(1, 47, 52))], 1.0)
      year_results$seasonal_mult <- year_results$seasonal_mult * year_results$holiday_mult
      
      # Cause-specific seasonal patterns
      year_results$cause_seasonal <- 1.0
      for (cause in names(causes)) {
        mask <- year_results$cause == cause
        if (any(mask) && cause %in% names(seasonal_patterns)) {
          year_results$cause_seasonal[mask] <- seasonal_patterns[[cause]][year_results$week_index[mask]]
        }
      }
      
      # COVID-specific adjustments for COVID scenario
      
      if (scenario_name == "COVID") {
        year_str <- as.character(year)
        covid_mult <- ifelse(year_str %in% names(covid_impacts), covid_impacts[year_str], 1.0)
        covid_prop_base <- ifelse(year_str %in% names(covid_props), covid_props[year_str], 0.0)
        
        if (covid_prop_base > 0) {
          
          if (year == 2020) {
            # Sharp, dramatic peaks with very low baseline (reduced amplitudes)
            spring_wave <- pmax(0, 3 * exp(-0.20 * ((year_results$week - 17)^2) / (1^2)))  
            summer_wave <- pmax(0, 2.0 * exp(-0.5 * ((year_results$week - 32)^2) / (3^2)))     
            winter_wave <- pmax(0, 5.0 * exp(-0.6 * ((year_results$week - 50)^2) / (3^2))) 
            wave_pattern <- 0.03 + spring_wave + summer_wave + winter_wave  
            
          } else if (year == 2021) {
            winter_surge <- pmax(0, 4.5 * exp(-0.8 * ((year_results$week - 4)^2) / (3^2)))   
            delta_wave <- pmax(0, 4.5 * exp(-0.75 * ((year_results$week - 35)^2) / (3^2)))     
            late_winter <- pmax(0, 2 * exp(-0.25 * ((year_results$week - 50)^2) / (3^2))) 
            wave_pattern <- 0 + winter_surge + delta_wave + late_winter
            
          } else if (year == 2022) {
            omicron_wave <- pmax(0, 3.5 * exp(-0.5 * ((year_results$week - 8)^2) / (2^2)))   
            summer_wave <- pmax(0, 1.8 * exp(-0.5 * ((year_results$week - 28)^2) / (2^2)))   
            wave_pattern <- 0.2 + omicron_wave + summer_wave
            
          } else if (year == 2023) {
            # SMALLER TROUGHS: Higher baseline, smaller amplitude
            seasonal_covid <- 1 + 0.7 * cos(2 * pi * (year_results$week - 4) / 52)  # 1.2 → 0.8 amplitude
            wave_pattern <- pmax(1, seasonal_covid)  # 0.3 → 0.7 baseline (much higher floor)
            
          } else if (year == 2024) {
            # EVEN SMALLER TROUGHS: Higher baseline, minimal amplitude
            seasonal_covid <- 1 + 0.4 * cos(2 * pi * (year_results$week - 4) / 52)  # Even smaller amplitude
            wave_pattern <- pmax(1, seasonal_covid)  # 0.3 → 0.8 baseline (very high floor) 
          }
          else {
            # Still dramatic but smaller
            seasonal_covid <- 1 + 0.7 * cos(2 * pi * (year_results$week - 4) / 52)  
            wave_pattern <- 0.8*pmax(0.3, seasonal_covid)
          }
          
          # Reduced maximum COVID proportion
          covid_prop_weekly <- pmin(covid_prop_base * wave_pattern, 0.30)  
          
          # Reduced overall mortality surge during COVID peaks
          surge_intensity <- covid_prop_weekly / max(covid_prop_weekly, 0.01)
          surge_multiplier <- 1 + 0.2 * surge_intensity  
          
          year_results$covid_prop <- covid_prop_weekly
          year_results$cause_prop <- ifelse(year_results$cause == "COVID-19", 
                                            year_results$covid_prop,
                                            year_results$cause_prop * (1 - year_results$covid_prop))
          year_results$covid_mult <- covid_mult * surge_multiplier
          
          covid_activity <- mean(covid_prop_weekly)
          variation_scale <- 0.1 * (1 + 1.0 * covid_activity)  
          
        } else {
          year_results$covid_prop <- 0
          year_results$covid_mult <- covid_mult
          variation_scale <- 0.15
        }
      } else {
        year_results$covid_mult <- 1.0
        variation_scale <- 0.15
      }
      
      # Secular trends
      cause_to_trend <- list(
        "Heart disease" = "Heart_disease", "Cancer" = "Cancer", "Accidents" = "Accidents",
        "Stroke" = "Stroke", "Chronic respiratory" = "Chronic_respiratory", 
        "Alzheimer's" = "Alzheimers", "Diabetes" = "Diabetes",
        "Influenza/Pneumonia" = "Influenza_Pneumonia", "Kidney disease" = "Kidney_disease",
        "Suicide" = "Suicide", "Septicemia" = "Septicemia", 
        "Chronic liver disease" = "Chronic_liver_disease", "Hypertension" = "Hypertension",
        "Parkinson's disease" = "Parkinson_disease", "COVID-19" = "Other", "Other" = "Other"
      )
      
      year_results$secular_trend <- 1.0
      for (cause in names(causes)) {
        mask <- year_results$cause == cause
        if (any(mask)) {
          trend_col <- cause_to_trend[[cause]]
          if (!is.null(trend_col)) {
            year_results$secular_trend[mask] <- base_trend[[trend_col]]
          }
        }
      }
      
      # Calculate deaths
      weeks_per_year <- 365.25 / 7
      year_results$base_deaths <- (year_results$population * year_results$age_death_rate / 1000 / weeks_per_year)
      year_results$adjusted_deaths <- year_results$base_deaths * year_results$cause_prop * 
        year_results$seasonal_mult * year_results$cause_seasonal * 
        year_results$secular_trend * year_results$covid_mult
      
      # Different seeds for different scenarios in post-2020 period
      scenario_seed_offset <- ifelse(scenario_name == "COVID", 100, 200)
      set.seed(year * 1000 + scenario_seed_offset + seed_set)
      variations <- rnorm(nrow(year_results), mean = 0, sd = variation_scale)
      year_results$final_deaths <- pmax(0, round(year_results$adjusted_deaths * (1 + variations)))
      
      # Keep only non-zero deaths
      year_final <- year_results[year_results$final_deaths > 0, 
                                 c("year", "week", "age_group", "cause", "final_deaths", "population", "week_start_date")]
      names(year_final)[5] <- "deaths"
      
      post_data <- rbind(post_data, year_final)
    }
    
    return(post_data)
  }
  
  # Generate post-2020 data for both scenarios
  covid_post_2020 <- generate_post_2020("COVID")
  no_covid_post_2020 <- generate_post_2020("No COVID")
  
  # STEP 3: Combine pre-2020 (identical) with post-2020 (different) for each scenario
  covid_scenario <- rbind(
    pre_2020_data %>% mutate(scenario = "COVID"),
    covid_post_2020 %>% mutate(scenario = "COVID")
  ) %>% mutate(deaths = as.integer(deaths)) # ensure counts
  
  no_covid_scenario <- rbind(
    pre_2020_data %>% mutate(scenario = "No COVID"),  # Same pre-2020 data!
    no_covid_post_2020 %>% mutate(scenario = "No COVID")
  ) %>% mutate(deaths = as.integer(deaths)) # ensure counts
  
  #cat("Final datasets created.")
  
  return(list(
    covid_scenario = covid_scenario,
    no_covid_scenario = no_covid_scenario,
    divergence_year = divergence_year
  ))
}
simulate_mortality_hypo1 <- function(start_year = 2014, end_year = 2024, divergence_year = 2020, seed_set = 1) {
  base_death_rate <- 7.5 # per 1000
  
  # Pre-calculate all date information
  all_years <- start_year:end_year
  n_years <- length(all_years)
  
  # Create complete date grid efficiently
  date_grid <- data.frame()
  for (year in all_years) {
    year_start <- as.Date(paste0(year, "-01-01"))
    weeks_in_year <- min(53, ceiling(as.numeric(difftime(
      as.Date(paste0(year, "-12-31")), year_start, units = "weeks"))))
    
    year_dates <- data.frame(
      year = year,
      week = 1:weeks_in_year,
      week_start_date = year_start + weeks(0:(weeks_in_year-1)),
      week_index = pmin(1:weeks_in_year, 52)
    )
    date_grid <- rbind(date_grid, year_dates)
  }
  
  # Pre-calculate population by year (from Census Bureau)
  years_elapsed <- all_years - start_year
  populations <- c(318301008, 320635163, 322941311, 324985539, 326687501, 328239523,
                   331577720, 332099760, 334017321, 336806231, 340110988)
  #populations <- base_population * (1 + annual_population_growth)^years_elapsed
  names(populations) <- all_years
  
  # Age groups in proportion
  
  calculate_age_groups <- function(year) {
    # Hypothetical proportions for a much younger population
    base_proportions <- c("0-17" = 0.31, "18-34" = 0.35, "35-54" = 0.25, 
                          "55-74" = 0.05, "75+" = 0.04)
    # Years since 2014
    years_since <- year - 2014
    
    # Youth declining due to lower birth rates
    prop_0_17 <- base_proportions["0-17"] * (1.01^years_since)  # +1% per year
    
    # Young adults declining due to lower birth rates in 1990s-2000s
    prop_18_34 <- base_proportions["18-34"] * (1.005^years_since)  # +0.5% per year
    
    # Middle-aged relatively stable but slight decline
    prop_35_54 <- base_proportions["35-54"] * (0.9995^years_since)  # +0.05% per year
    
    # Older adults increasing (baby boomers aging)
    prop_55_74 <- base_proportions["55-74"] * (0.985^years_since)  # -1.5% per year
    
    # Elderly increasing (life expectancy + baby boomers)
    prop_75_plus <- base_proportions["75+"] * (.975^years_since)  # -2.5% per year
    
    # Normalize to ensure they sum to 1.0
    raw_props <- c(prop_0_17, prop_18_34, prop_35_54, prop_55_74, prop_75_plus)
    normalized_props <- raw_props / sum(raw_props)
    
    names(normalized_props) <- c("0-17", "18-34", "35-54", "55-74", "75+")
    return(normalized_props)
  }
  
  # Age-specific death rates per 1000 people
  age_death_rates <- c("0-17" = 1/1.1, "18-34" = 2.3/1.3, "35-54" = 7/1.3, 
                       "55-74" = 28.7/1.1, "75+" = 35.0/1.0)
  
  # Cause distributions 
  pre_2020_causes <- c(
    "Heart disease" = 0.229, "Cancer" = 0.222, "Accidents" = 0.052,
    "Stroke" = 0.051, "Chronic respiratory" = 0.049, "Alzheimer's" = 0.032,
    "Diabetes" = 0.020, "Influenza/Pneumonia" = 0.028, "Kidney disease" = 0.028,
    "Suicide" = 0.009, "Septicemia" = 0.020, "Chronic liver disease" = 0.014,
    "Hypertension" = 0.009, "Parkinson's disease" = 0.012, "Other" = 0.225
  )
  
  post_2020_causes_covid <- c(
    "Heart disease" = 0.220, "Cancer" = 0.198, "COVID-19" = 0.007,
    "Accidents" = 0.065, "Stroke" = 0.057, "Chronic respiratory" = 0.059,
    "Alzheimer's" = 0.046, "Diabetes" = 0.029, "Influenza/Pneumonia" = 0.025,
    "Kidney disease" = 0.012, "Suicide" = 0.018, "Other" = 0.264
  )
  
  # Seasonal variation
  seasonal_base <- 1 + 0.03 * cos(2 * pi * (1:52 - 4) / 52)  
  
  # Pre-calculate cause-specific seasonal pattern
  create_seasonal_fast <- function(winter_mult, summer_mult, peak_week = 4) {
    weeks <- 1:52
    
    winter_mult_adj <- 1 + (winter_mult - 1) * 0.6
    summer_mult_adj <- 1 + (summer_mult - 1) * 0.6
    
    amplitude <- (winter_mult_adj - summer_mult_adj) / 2
    baseline <- (winter_mult_adj + summer_mult_adj) / 2
    pmax(baseline + amplitude * cos(2 * pi * (weeks - peak_week) / 52), 0.1)
  }
  
  seasonal_patterns <- list(
    "Heart disease" = create_seasonal_fast(1.2, 0.9, 4),
    "Cancer" = rep(1.0, 52),
    "Accidents" = create_seasonal_fast(0.9, 1.2, 30),
    "Stroke" = create_seasonal_fast(1.2, 0.9, 4),
    "Chronic respiratory" = create_seasonal_fast(1.4, 0.7, 6),
    "Alzheimer's" = rep(1.0, 52),
    "Diabetes" = rep(1.0, 52),
    "Influenza/Pneumonia" = create_seasonal_fast(2.5, 0.3, 8),
    "Kidney disease" = rep(1.0, 52),
    "Suicide" = create_seasonal_fast(0.9, 1.1, 20),
    "Septicemia" = create_seasonal_fast(1.1, 0.9, 6),
    "Chronic liver disease" = rep(1.0, 52),
    "Hypertension" = create_seasonal_fast(1.1, 0.9, 4),
    "Parkinson's disease" = rep(1.0, 52),
    "COVID-19" = rep(1.0, 52),
    "Other" = rep(1.0, 52)
  )
  
  # Pre-calculate secular trends for all years
  calculate_trends <- function(years, base_year = 2014) {
    elapsed <- (years - base_year)
    trends <- data.frame(
      year = years,
      Heart_disease = 0.98^elapsed,
      Cancer = 0.985^elapsed,
      Accidents = 1.01^elapsed,
      Stroke = 0.98^elapsed,
      Chronic_respiratory = 0.99^elapsed,
      Alzheimers = 1.03^elapsed,
      Diabetes = 1.015^elapsed,
      Influenza_Pneumonia = 0.995^elapsed,
      Kidney_disease = 1.005^elapsed,
      Suicide = 1.02^elapsed,
      Septicemia = 1.01^elapsed,
      Chronic_liver_disease = 1.025^elapsed,
      Hypertension = 1.01^elapsed,
      Parkinson_disease = 1.015^elapsed,
      Other = 1.0
    )
    return(trends)
  }
  
  pre_trends <- calculate_trends(all_years)
  
  # COVID impact functions
  covid_impacts <- c("2020" = 1.1, "2021" = 1.12, "2022" = 1.05, 
                     "2023" = 0.88, "2024" = 0.85)
  covid_props <- c("2020" = 0.11, "2021" = 0.2, "2022" = 0.05, 
                   "2023" = 0.02, "2024" = 0.01)
  
  # STEP 1: Generate identical pre-2020 data for BOTH scenarios
  #cat("Generating identical pre-2020 data...\n")
  
  pre_2020_data <- data.frame()
  
  for (year in start_year:(divergence_year - 1)) {
    # Calculate age groups for this specific year
    age_groups <- calculate_age_groups(year)
    
    year_data <- date_grid[date_grid$year == year, ]
    current_pop <- populations[as.character(year)]
    base_trend <- pre_trends[pre_trends$year == year, ]
    
    # Calculate base weekly deaths
    weeks_per_year <- 365.25 / 7
    
    # Create all combinations
    year_results <- expand.grid(
      week = year_data$week,
      age_group = names(age_groups),
      cause = names(pre_2020_causes),
      stringsAsFactors = FALSE
    )
    
    # Add calculated values
    year_results$year <- year
    year_results$week_start_date <- year_data$week_start_date[match(year_results$week, year_data$week)]
    year_results$week_index <- year_data$week_index[match(year_results$week, year_data$week)]
    year_results$age_prop <- age_groups[year_results$age_group]
    year_results$age_death_rate <- age_death_rates[year_results$age_group]
    year_results$cause_prop <- pre_2020_causes[year_results$cause]
    year_results$population <- round(current_pop * year_results$age_prop)
    
    year_results$seasonal_mult <- seasonal_base[year_results$week_index]
    year_results$holiday_mult <- ifelse(year_results$week %in% c(1, 47, 52), 
                                        c(1.04, 1.02, 1.06)[match(year_results$week, c(1, 47, 52))], 1.0)  
    year_results$seasonal_mult <- year_results$seasonal_mult * year_results$holiday_mult
    
    # Cause-specific seasonal patterns
    year_results$cause_seasonal <- 1.0
    for (cause in names(pre_2020_causes)) {
      mask <- year_results$cause == cause
      if (any(mask)) {
        year_results$cause_seasonal[mask] <- seasonal_patterns[[cause]][year_results$week_index[mask]]
      }
    }
    
    # Secular trends
    cause_to_trend <- list(
      "Heart disease" = "Heart_disease", "Cancer" = "Cancer", "Accidents" = "Accidents",
      "Stroke" = "Stroke", "Chronic respiratory" = "Chronic_respiratory", 
      "Alzheimer's" = "Alzheimers", "Diabetes" = "Diabetes",
      "Influenza/Pneumonia" = "Influenza_Pneumonia", "Kidney disease" = "Kidney_disease",
      "Suicide" = "Suicide", "Septicemia" = "Septicemia", 
      "Chronic liver disease" = "Chronic_liver_disease", "Hypertension" = "Hypertension",
      "Parkinson's disease" = "Parkinson_disease", "Other" = "Other"
    )
    
    year_results$secular_trend <- 1.0
    for (cause in names(pre_2020_causes)) {
      mask <- year_results$cause == cause
      if (any(mask)) {
        trend_col <- cause_to_trend[[cause]]
        if (!is.null(trend_col)) {
          year_results$secular_trend[mask] <- base_trend[[trend_col]]
        }
      }
    }
    
    # Calculate base deaths
    year_results$base_deaths <- (year_results$population * year_results$age_death_rate / 1000 / weeks_per_year)
    year_results$adjusted_deaths <- year_results$base_deaths * year_results$cause_prop * 
      year_results$seasonal_mult * year_results$cause_seasonal * 
      year_results$secular_trend
    
    # CRITICAL: Use SAME seed for pre-2020 data to ensure identical results
    set.seed(year * 1000 + seed_set)  # Same seed regardless of scenario
    variations <- rnorm(nrow(year_results), mean = 0, sd = 0.1)  
    year_results$final_deaths <- pmax(0, round(year_results$adjusted_deaths * (1 + variations)))
    
    # Adjust deaths and variation for specific years and peak weeks
    if (year == 2015 || year == 2018) {
      peak_weeks <- 44:52
      reduction_factor <- 0.92
      variation_reduction_factor <- 0.2
      year_results$final_deaths[year_results$week %in% peak_weeks] <- year_results$final_deaths[year_results$week %in% peak_weeks] * reduction_factor
      variations[year_results$week %in% peak_weeks] <- variations[year_results$week %in% peak_weeks] * variation_reduction_factor
    }
    
    if (year == 2016 || year == 2019) {
      peak_weeks <- 1:15
      reduction_factor <- 0.92
      variation_reduction_factor <- 0.2
      year_results$final_deaths[year_results$week %in% peak_weeks] <- year_results$final_deaths[year_results$week %in% peak_weeks] * reduction_factor
      variations[year_results$week %in% peak_weeks] <- variations[year_results$week %in% peak_weeks] * variation_reduction_factor
    }
    # Keep only non-zero deaths
    year_final <- year_results[year_results$final_deaths > 0, 
                               c("year", "week", "age_group", "cause", "final_deaths", "population", "week_start_date")]
    names(year_final)[5] <- "deaths"
    
    pre_2020_data <- rbind(pre_2020_data, year_final)
  }
  
  #cat("Pre-2020 data generated. \n")
  
  # STEP 2: Generate post-2020 data for each scenario separately
  generate_post_2020 <- function(scenario_name) {
    # cat("Generating post-2020 data for", scenario_name, "scenario...\n")
    
    post_data <- data.frame()
    
    for (year in divergence_year:end_year) {
      # Calculate age groups for this specific year
      age_groups <- calculate_age_groups(year)
      
      year_data <- date_grid[date_grid$year == year, ]
      current_pop <- populations[as.character(year)]
      
      # Choose parameters based on scenario
      
      if (scenario_name == "COVID") {
        causes <- post_2020_causes_covid
        base_trend <- pre_trends[pre_trends$year == year, ]
        # Add COVID disruption to trends
        years_since <- year - divergence_year
        base_trend$Heart_disease <- base_trend$Heart_disease * 1.01^years_since
        base_trend$Cancer <- base_trend$Cancer * 1.005^years_since
        base_trend$Suicide <- base_trend$Suicide * 1.03^years_since
      } else {
        causes <- pre_2020_causes  # Continue pre-pandemic patterns
        base_trend <- pre_trends[pre_trends$year == year, ]
      }
      
      # Create all combinations
      year_results <- expand.grid(
        week = year_data$week,
        age_group = names(age_groups),
        cause = names(causes),
        stringsAsFactors = FALSE
      )
      
      # Add calculated values
      year_results$year <- year
      year_results$week_start_date <- year_data$week_start_date[match(year_results$week, year_data$week)]
      year_results$week_index <- year_data$week_index[match(year_results$week, year_data$week)]
      year_results$age_prop <- age_groups[year_results$age_group]
      year_results$age_death_rate <- age_death_rates[year_results$age_group]
      year_results$cause_prop <- causes[year_results$cause]
      year_results$population <- round(current_pop * year_results$age_prop)
      
      # Seasonal multipliers
      year_results$seasonal_mult <- seasonal_base[year_results$week_index] 
      year_results$holiday_mult <- ifelse(year_results$week %in% c(1, 47, 52), 
                                          c(1.04, 1.02, 1.06)[match(year_results$week, c(1, 47, 52))], 1.0)
      year_results$seasonal_mult <- year_results$seasonal_mult * year_results$holiday_mult
      
      # Cause-specific seasonal patterns
      year_results$cause_seasonal <- 1.0
      for (cause in names(causes)) {
        mask <- year_results$cause == cause
        if (any(mask) && cause %in% names(seasonal_patterns)) {
          year_results$cause_seasonal[mask] <- seasonal_patterns[[cause]][year_results$week_index[mask]]
        }
      }
      
      # COVID-specific adjustments for COVID scenario
      
      if (scenario_name == "COVID") {
        year_str <- as.character(year)
        covid_mult <- ifelse(year_str %in% names(covid_impacts), covid_impacts[year_str], 1.0)
        covid_prop_base <- ifelse(year_str %in% names(covid_props), covid_props[year_str], 0.0)
        
        if (covid_prop_base > 0) {
          
          if (year == 2020) {
            winter_surge <- pmax(0, 9 * exp(-1 * ((year_results$week - 7)^2) / (2^2)))   
            delta_wave <- pmax(0, 4.8 * exp(-0.8 * ((year_results$week - 35)^2) / (3^2)))     
            late_winter <- pmax(0, 5 * exp(-0.4 * ((year_results$week - 52)^2) / (2^2))) 
            wave_pattern <- 0 + winter_surge + late_winter
            
          } else if (year == 2021) {
            # Sharp, dramatic peaks with very low baseline (reduced amplitudes)
            spring_wave <- pmax(0, 3.3 * exp(-0.25 * ((year_results$week - 15)^2) / (2^2)))  
            winter_wave <- pmax(0, 5.0 * exp(-0.25 * ((year_results$week - 47)^2) / (2.5^2))) 
            wave_pattern <- 0.03 + spring_wave + winter_wave  
            
          } else if (year == 2022) {
            omicron_wave <- pmax(0, 3.5 * exp(-0.5 * ((year_results$week - 8)^2) / (2^2)))   
            summer_wave <- pmax(0, 1.8 * exp(-0.5 * ((year_results$week - 28)^2) / (2^2)))   
            wave_pattern <- 0.2 + omicron_wave + summer_wave
            
          } else if (year == 2023) {
            # SMALLER TROUGHS: Higher baseline, smaller amplitude
            seasonal_covid <- 1 + 0.7 * cos(2 * pi * (year_results$week - 4) / 52)  # 1.2 → 0.8 amplitude
            wave_pattern <- pmax(1, seasonal_covid)  # 0.3 → 0.7 baseline (much higher floor)
            
          } else if (year == 2024) {
            # EVEN SMALLER TROUGHS: Higher baseline, minimal amplitude
            seasonal_covid <- 1 + 0.4 * cos(2 * pi * (year_results$week - 4) / 52)  # Even smaller amplitude
            wave_pattern <- pmax(1, seasonal_covid)  # 0.3 → 0.8 baseline (very high floor) 
          }
          else {
            # Still dramatic but smaller
            seasonal_covid <- 1 + 0.7 * cos(2 * pi * (year_results$week - 4) / 52)  
            wave_pattern <- 0.8*pmax(0.3, seasonal_covid)
          }
          
          # Reduced maximum COVID proportion
          covid_prop_weekly <- pmin(covid_prop_base * wave_pattern, 0.30)  
          
          # Reduced overall mortality surge during COVID peaks
          surge_intensity <- covid_prop_weekly / max(covid_prop_weekly, 0.01)
          surge_multiplier <- 1 + 0.2 * surge_intensity  
          
          year_results$covid_prop <- covid_prop_weekly
          year_results$cause_prop <- ifelse(year_results$cause == "COVID-19", 
                                            year_results$covid_prop,
                                            year_results$cause_prop * (1 - year_results$covid_prop))
          year_results$covid_mult <- covid_mult * surge_multiplier
          
          covid_activity <- mean(covid_prop_weekly)
          variation_scale <- 0.1 * (1 + 1.0 * covid_activity)  
          
        } else {
          year_results$covid_prop <- 0
          year_results$covid_mult <- covid_mult
          variation_scale <- 0.15
        }
      } else {
        year_results$covid_mult <- 1.0
        variation_scale <- 0.15
      }
      
      # Secular trends
      cause_to_trend <- list(
        "Heart disease" = "Heart_disease", "Cancer" = "Cancer", "Accidents" = "Accidents",
        "Stroke" = "Stroke", "Chronic respiratory" = "Chronic_respiratory", 
        "Alzheimer's" = "Alzheimers", "Diabetes" = "Diabetes",
        "Influenza/Pneumonia" = "Influenza_Pneumonia", "Kidney disease" = "Kidney_disease",
        "Suicide" = "Suicide", "Septicemia" = "Septicemia", 
        "Chronic liver disease" = "Chronic_liver_disease", "Hypertension" = "Hypertension",
        "Parkinson's disease" = "Parkinson_disease", "COVID-19" = "Other", "Other" = "Other"
      )
      
      year_results$secular_trend <- 1.0
      for (cause in names(causes)) {
        mask <- year_results$cause == cause
        if (any(mask)) {
          trend_col <- cause_to_trend[[cause]]
          if (!is.null(trend_col)) {
            year_results$secular_trend[mask] <- base_trend[[trend_col]]
          }
        }
      }
      
      # Calculate deaths
      weeks_per_year <- 365.25 / 7
      year_results$base_deaths <- (year_results$population * year_results$age_death_rate / 1000 / weeks_per_year)
      year_results$adjusted_deaths <- year_results$base_deaths * year_results$cause_prop * 
        year_results$seasonal_mult * year_results$cause_seasonal * 
        year_results$secular_trend * year_results$covid_mult
      
      # Different seeds for different scenarios in post-2020 period
      scenario_seed_offset <- ifelse(scenario_name == "COVID", 100, 200)
      set.seed(year * 1000 + scenario_seed_offset + seed_set)
      variations <- rnorm(nrow(year_results), mean = 0, sd = variation_scale)
      year_results$final_deaths <- pmax(0, round(year_results$adjusted_deaths * (1 + variations)))
      
      # Keep only non-zero deaths
      year_final <- year_results[year_results$final_deaths > 0, 
                                 c("year", "week", "age_group", "cause", "final_deaths", "population", "week_start_date")]
      names(year_final)[5] <- "deaths"
      
      post_data <- rbind(post_data, year_final)
    }
    
    return(post_data)
  }
  
  # Generate post-2020 data for both scenarios
  covid_post_2020 <- generate_post_2020("COVID")
  no_covid_post_2020 <- generate_post_2020("No COVID")
  
  # STEP 3: Combine pre-2020 (identical) with post-2020 (different) for each scenario
  covid_scenario <- rbind(
    pre_2020_data %>% mutate(scenario = "COVID"),
    covid_post_2020 %>% mutate(scenario = "COVID")
  ) %>% mutate(deaths = as.integer(deaths)) # ensure counts
  
  no_covid_scenario <- rbind(
    pre_2020_data %>% mutate(scenario = "No COVID"),  # Same pre-2020 data!
    no_covid_post_2020 %>% mutate(scenario = "No COVID")
  ) %>% mutate(deaths = as.integer(deaths)) # ensure counts
  
  #cat("Final datasets created.")
  
  return(list(
    covid_scenario = covid_scenario,
    no_covid_scenario = no_covid_scenario,
    divergence_year = divergence_year
  ))
}
simulate_mortality_hypo2 <- function(start_year = 2014, end_year = 2024, divergence_year = 2020, seed_set = 1) {
  base_death_rate <- 7.5 # per 1000
  
  # Pre-calculate all date information
  all_years <- start_year:end_year
  n_years <- length(all_years)
  
  # Create complete date grid efficiently
  date_grid <- data.frame()
  for (year in all_years) {
    year_start <- as.Date(paste0(year, "-01-01"))
    weeks_in_year <- min(53, ceiling(as.numeric(difftime(
      as.Date(paste0(year, "-12-31")), year_start, units = "weeks"))))
    
    year_dates <- data.frame(
      year = year,
      week = 1:weeks_in_year,
      week_start_date = year_start + weeks(0:(weeks_in_year-1)),
      week_index = pmin(1:weeks_in_year, 52)
    )
    date_grid <- rbind(date_grid, year_dates)
  }
  
  # Pre-calculate population by year (from Census Bureau)
  years_elapsed <- all_years - start_year
  populations <- c(318301008, 320635163, 322941311, 324985539, 326687501, 328239523,
                   331577720, 332099760, 334017321, 336806231, 340110988)
  #populations <- base_population * (1 + annual_population_growth)^years_elapsed
  names(populations) <- all_years
  
  # Age groups in proportion
  calculate_age_groups <- function(year) {
    # Hypothetical proportions for a much older population
    base_proportions <- c("0-17" = 0.15, "18-34" = 0.20, "35-54" = 0.25, 
                          "55-74" = 0.25, "75+" = 0.15)
    # Years since 2014
    years_since <- year - 2014
    
    # Youth declining due to lower birth rates
    prop_0_17 <- base_proportions["0-17"] * (0.99^years_since)  # -1% per year
    
    # Young adults stable but slight decline
    prop_18_34 <- base_proportions["18-34"] * (0.995^years_since)  # -0.5% per year
    
    # Middle-aged relatively stable
    prop_35_54 <- base_proportions["35-54"] * (1.00^years_since)  # 0% per year
    
    # Older adults increasing
    prop_55_74 <- base_proportions["55-74"] * (1.01^years_since)  # +1% per year
    
    # Elderly significantly increasing (life expectancy + baby boomers)
    prop_75_plus <- base_proportions["75+"] * (1.015^years_since)  # +1.5% per year
    
    # Normalize to ensure they sum to 1.0
    raw_props <- c(prop_0_17, prop_18_34, prop_35_54, prop_55_74, prop_75_plus)
    normalized_props <- raw_props / sum(raw_props)
    
    names(normalized_props) <- c("0-17", "18-34", "35-54", "55-74", "75+")
    return(normalized_props)
  }
  
  # Age-specific death rates per 1000 people
  age_death_rates <- c("0-17" = 1/1.1, "18-34" = 2.3/1.3, "35-54" = 7/1.3, 
                       "55-74" = 28.7/1.1, "75+" = 35.0/1.0)
  
  # Cause distributions 
  pre_2020_causes <- c(
    "Heart disease" = 0.231, "Cancer" = 0.221, "Accidents" = 0.055,
    "Stroke" = 0.059, "Chronic respiratory" = 0.048, "Alzheimer's" = 0.031,
    "Diabetes" = 0.028, "Influenza/Pneumonia" = 0.032, "Kidney disease" = 0.026,
    "Suicide" = 0.012, "Septicemia" = 0.018, "Chronic liver disease" = 0.012,
    "Hypertension" = 0.010, "Parkinson's disease" = 0.002, "Other" = 0.216
  )
  
  post_2020_causes_covid <- c(
    "Heart disease" = 0.219, "Cancer" = 0.194, "COVID-19" = 0.010,
    "Accidents" = 0.071, "Stroke" = 0.055, "Chronic respiratory" = 0.057,
    "Alzheimer's" = 0.053, "Diabetes" = 0.032, "Influenza/Pneumonia" = 0.022,
    "Kidney disease" = 0.014, "Suicide" = 0.015, "Other" = 0.257
  )
  
  # Seasonal variation
  seasonal_base <- 1 + 0.03 * cos(2 * pi * (1:52 - 4) / 52)  
  
  # Pre-calculate cause-specific seasonal pattern
  create_seasonal_fast <- function(winter_mult, summer_mult, peak_week = 4) {
    weeks <- 1:52
    
    winter_mult_adj <- 1 + (winter_mult - 1) * 0.6
    summer_mult_adj <- 1 + (summer_mult - 1) * 0.6
    
    amplitude <- (winter_mult_adj - summer_mult_adj) / 2
    baseline <- (winter_mult_adj + summer_mult_adj) / 2
    pmax(baseline + amplitude * cos(2 * pi * (weeks - peak_week) / 52), 0.1)
  }
  
  seasonal_patterns <- list(
    "Heart disease" = create_seasonal_fast(1.2, 0.9, 4),
    "Cancer" = rep(1.0, 52),
    "Accidents" = create_seasonal_fast(0.9, 1.2, 30),
    "Stroke" = create_seasonal_fast(1.2, 0.9, 4),
    "Chronic respiratory" = create_seasonal_fast(1.4, 0.7, 6),
    "Alzheimer's" = rep(1.0, 52),
    "Diabetes" = rep(1.0, 52),
    "Influenza/Pneumonia" = create_seasonal_fast(2.5, 0.3, 8),
    "Kidney disease" = rep(1.0, 52),
    "Suicide" = create_seasonal_fast(0.9, 1.1, 20),
    "Septicemia" = create_seasonal_fast(1.1, 0.9, 6),
    "Chronic liver disease" = rep(1.0, 52),
    "Hypertension" = create_seasonal_fast(1.1, 0.9, 4),
    "Parkinson's disease" = rep(1.0, 52),
    "COVID-19" = rep(1.0, 52),
    "Other" = rep(1.0, 52)
  )
  
  # Pre-calculate secular trends for all years
  calculate_trends <- function(years, base_year = 2014) {
    elapsed <- (years - base_year)
    trends <- data.frame(
      year = years,
      Heart_disease = 0.98^elapsed,
      Cancer = 0.985^elapsed,
      Accidents = 1.01^elapsed,
      Stroke = 0.98^elapsed,
      Chronic_respiratory = 0.99^elapsed,
      Alzheimers = 1.03^elapsed,
      Diabetes = 1.015^elapsed,
      Influenza_Pneumonia = 0.995^elapsed,
      Kidney_disease = 1.005^elapsed,
      Suicide = 1.02^elapsed,
      Septicemia = 1.01^elapsed,
      Chronic_liver_disease = 1.025^elapsed,
      Hypertension = 1.01^elapsed,
      Parkinson_disease = 1.015^elapsed,
      Other = 1.0
    )
    return(trends)
  }
  
  pre_trends <- calculate_trends(all_years)
  
  # COVID impact functions
  covid_impacts <- c("2020" = 1.03, "2021" = 1.12, "2022" = 1.05, 
                     "2023" = 0.88, "2024" = 0.85)
  covid_props <- c("2020" = 0.11, "2021" = 0.2, "2022" = 0.05, 
                   "2023" = 0.02, "2024" = 0.01)
  
  # STEP 1: Generate identical pre-2020 data for BOTH scenarios
  #cat("Generating identical pre-2020 data...\n")
  
  pre_2020_data <- data.frame()
  
  for (year in start_year:(divergence_year - 1)) {
    # Calculate age groups for this specific year
    age_groups <- calculate_age_groups(year)
    
    year_data <- date_grid[date_grid$year == year, ]
    current_pop <- populations[as.character(year)]
    base_trend <- pre_trends[pre_trends$year == year, ]
    
    # Calculate base weekly deaths
    weeks_per_year <- 365.25 / 7
    
    # Create all combinations
    year_results <- expand.grid(
      week = year_data$week,
      age_group = names(age_groups),
      cause = names(pre_2020_causes),
      stringsAsFactors = FALSE
    )
    
    # Add calculated values
    year_results$year <- year
    year_results$week_start_date <- year_data$week_start_date[match(year_results$week, year_data$week)]
    year_results$week_index <- year_data$week_index[match(year_results$week, year_data$week)]
    year_results$age_prop <- age_groups[year_results$age_group]
    year_results$age_death_rate <- age_death_rates[year_results$age_group]
    year_results$cause_prop <- pre_2020_causes[year_results$cause]
    year_results$population <- round(current_pop * year_results$age_prop)
    
    year_results$seasonal_mult <- seasonal_base[year_results$week_index]
    year_results$holiday_mult <- ifelse(year_results$week %in% c(1, 47, 52), 
                                        c(1.04, 1.02, 1.06)[match(year_results$week, c(1, 47, 52))], 1.0)  
    year_results$seasonal_mult <- year_results$seasonal_mult * year_results$holiday_mult
    
    # Cause-specific seasonal patterns
    year_results$cause_seasonal <- 1.0
    for (cause in names(pre_2020_causes)) {
      mask <- year_results$cause == cause
      if (any(mask)) {
        year_results$cause_seasonal[mask] <- seasonal_patterns[[cause]][year_results$week_index[mask]]
      }
    }
    
    # Secular trends
    cause_to_trend <- list(
      "Heart disease" = "Heart_disease", "Cancer" = "Cancer", "Accidents" = "Accidents",
      "Stroke" = "Stroke", "Chronic respiratory" = "Chronic_respiratory", 
      "Alzheimer's" = "Alzheimers", "Diabetes" = "Diabetes",
      "Influenza/Pneumonia" = "Influenza_Pneumonia", "Kidney disease" = "Kidney_disease",
      "Suicide" = "Suicide", "Septicemia" = "Septicemia", 
      "Chronic liver disease" = "Chronic_liver_disease", "Hypertension" = "Hypertension",
      "Parkinson's disease" = "Parkinson_disease", "Other" = "Other"
    )
    
    year_results$secular_trend <- 1.0
    for (cause in names(pre_2020_causes)) {
      mask <- year_results$cause == cause
      if (any(mask)) {
        trend_col <- cause_to_trend[[cause]]
        if (!is.null(trend_col)) {
          year_results$secular_trend[mask] <- base_trend[[trend_col]]
        }
      }
    }
    
    # Calculate base deaths
    year_results$base_deaths <- (year_results$population * year_results$age_death_rate / 1000 / weeks_per_year)
    year_results$adjusted_deaths <- year_results$base_deaths * year_results$cause_prop * 
      year_results$seasonal_mult * year_results$cause_seasonal * 
      year_results$secular_trend
    
    # CRITICAL: Use SAME seed for pre-2020 data to ensure identical results
    set.seed(year * 1000 + seed_set)  # Same seed regardless of scenario
    variations <- rnorm(nrow(year_results), mean = 0, sd = 0.1)  
    year_results$final_deaths <- pmax(0, round(year_results$adjusted_deaths * (1 + variations)))
    
    # Adjust deaths and variation for specific years and peak weeks
    if (year == 2015 || year == 2018) {
      peak_weeks <- 44:52
      reduction_factor <- 0.92
      variation_reduction_factor <- 0.2
      year_results$final_deaths[year_results$week %in% peak_weeks] <- year_results$final_deaths[year_results$week %in% peak_weeks] * reduction_factor
      variations[year_results$week %in% peak_weeks] <- variations[year_results$week %in% peak_weeks] * variation_reduction_factor
    }
    
    if (year == 2016 || year == 2019) {
      peak_weeks <- 1:15
      reduction_factor <- 0.92
      variation_reduction_factor <- 0.2
      year_results$final_deaths[year_results$week %in% peak_weeks] <- year_results$final_deaths[year_results$week %in% peak_weeks] * reduction_factor
      variations[year_results$week %in% peak_weeks] <- variations[year_results$week %in% peak_weeks] * variation_reduction_factor
    }
    # Keep only non-zero deaths
    year_final <- year_results[year_results$final_deaths > 0, 
                               c("year", "week", "age_group", "cause", "final_deaths", "population", "week_start_date")]
    names(year_final)[5] <- "deaths"
    
    pre_2020_data <- rbind(pre_2020_data, year_final)
  }
  
  #cat("Pre-2020 data generated. \n")
  
  # STEP 2: Generate post-2020 data for each scenario separately
  generate_post_2020 <- function(scenario_name) {
    # cat("Generating post-2020 data for", scenario_name, "scenario...\n")
    
    post_data <- data.frame()
    
    for (year in divergence_year:end_year) {
      # Calculate age groups for this specific year
      age_groups <- calculate_age_groups(year)
      
      year_data <- date_grid[date_grid$year == year, ]
      current_pop <- populations[as.character(year)]
      
      # Choose parameters based on scenario
      
      if (scenario_name == "COVID") {
        causes <- post_2020_causes_covid
        base_trend <- pre_trends[pre_trends$year == year, ]
        # Add COVID disruption to trends
        years_since <- year - divergence_year
        base_trend$Heart_disease <- base_trend$Heart_disease * 1.01^years_since
        base_trend$Cancer <- base_trend$Cancer * 1.005^years_since
        base_trend$Suicide <- base_trend$Suicide * 1.03^years_since
      } else {
        causes <- pre_2020_causes  # Continue pre-pandemic patterns
        base_trend <- pre_trends[pre_trends$year == year, ]
      }
      
      # Create all combinations
      year_results <- expand.grid(
        week = year_data$week,
        age_group = names(age_groups),
        cause = names(causes),
        stringsAsFactors = FALSE
      )
      
      # Add calculated values
      year_results$year <- year
      year_results$week_start_date <- year_data$week_start_date[match(year_results$week, year_data$week)]
      year_results$week_index <- year_data$week_index[match(year_results$week, year_data$week)]
      year_results$age_prop <- age_groups[year_results$age_group]
      year_results$age_death_rate <- age_death_rates[year_results$age_group]
      year_results$cause_prop <- causes[year_results$cause]
      year_results$population <- round(current_pop * year_results$age_prop)
      
      # Seasonal multipliers
      year_results$seasonal_mult <- seasonal_base[year_results$week_index] 
      year_results$holiday_mult <- ifelse(year_results$week %in% c(1, 47, 52), 
                                          c(1.04, 1.02, 1.06)[match(year_results$week, c(1, 47, 52))], 1.0)
      year_results$seasonal_mult <- year_results$seasonal_mult * year_results$holiday_mult
      
      # Cause-specific seasonal patterns
      year_results$cause_seasonal <- 1.0
      for (cause in names(causes)) {
        mask <- year_results$cause == cause
        if (any(mask) && cause %in% names(seasonal_patterns)) {
          year_results$cause_seasonal[mask] <- seasonal_patterns[[cause]][year_results$week_index[mask]]
        }
      }
      
      # COVID-specific adjustments for COVID scenario
      
      if (scenario_name == "COVID") {
        year_str <- as.character(year)
        covid_mult <- ifelse(year_str %in% names(covid_impacts), covid_impacts[year_str], 1.0)
        covid_prop_base <- ifelse(year_str %in% names(covid_props), covid_props[year_str], 0.0)
        
        if (covid_prop_base > 0) {
          
          if (year == 2020) {
            winter_surge <- pmax(0, 9 * exp(-1 * ((year_results$week - 7)^2) / (2^2)))   
            delta_wave <- pmax(0, 4.8 * exp(-0.8 * ((year_results$week - 35)^2) / (3^2)))     
            late_winter <- pmax(0, 5 * exp(-0.4 * ((year_results$week - 52)^2) / (2^2))) 
            wave_pattern <- 0 + winter_surge + delta_wave + late_winter
            
          } else if (year == 2021) {
            # Sharp, dramatic peaks with very low baseline (reduced amplitudes)
            spring_wave <- pmax(0, 3.3 * exp(-0.25 * ((year_results$week - 15)^2) / (2^2)))  
            winter_wave <- pmax(0, 5.0 * exp(-0.25 * ((year_results$week - 47)^2) / (2.5^2))) 
            wave_pattern <- 0.03 + spring_wave + winter_wave  
            
          } else if (year == 2022) {
            omicron_wave <- pmax(0, 3.5 * exp(-0.5 * ((year_results$week - 8)^2) / (2^2)))   
            summer_wave <- pmax(0, 1.8 * exp(-0.5 * ((year_results$week - 28)^2) / (2^2)))   
            wave_pattern <- 0.2 + omicron_wave + summer_wave
            
          } else if (year == 2023) {
            # SMALLER TROUGHS: Higher baseline, smaller amplitude
            seasonal_covid <- 1 + 0.7 * cos(2 * pi * (year_results$week - 4) / 52)  # 1.2 → 0.8 amplitude
            wave_pattern <- pmax(1, seasonal_covid)  # 0.3 → 0.7 baseline (much higher floor)
            
          } else if (year == 2024) {
            # EVEN SMALLER TROUGHS: Higher baseline, minimal amplitude
            seasonal_covid <- 1 + 0.4 * cos(2 * pi * (year_results$week - 4) / 52)  # Even smaller amplitude
            wave_pattern <- pmax(1, seasonal_covid)  # 0.3 → 0.8 baseline (very high floor) 
          }
          else {
            # Still dramatic but smaller
            seasonal_covid <- 1 + 0.7 * cos(2 * pi * (year_results$week - 4) / 52)  
            wave_pattern <- 0.8*pmax(0.3, seasonal_covid)
          }
          
          # Reduced maximum COVID proportion
          covid_prop_weekly <- pmin(covid_prop_base * wave_pattern, 0.30)  
          
          # Reduced overall mortality surge during COVID peaks
          surge_intensity <- covid_prop_weekly / max(covid_prop_weekly, 0.01)
          surge_multiplier <- 1 + 0.2 * surge_intensity  
          
          year_results$covid_prop <- covid_prop_weekly
          year_results$cause_prop <- ifelse(year_results$cause == "COVID-19", 
                                            year_results$covid_prop,
                                            year_results$cause_prop * (1 - year_results$covid_prop))
          year_results$covid_mult <- covid_mult * surge_multiplier
          
          covid_activity <- mean(covid_prop_weekly)
          variation_scale <- 0.1 * (1 + 1.0 * covid_activity)  
          
        } else {
          year_results$covid_prop <- 0
          year_results$covid_mult <- covid_mult
          variation_scale <- 0.15
        }
      } else {
        year_results$covid_mult <- 1.0
        variation_scale <- 0.15
      }
      
      # Secular trends
      cause_to_trend <- list(
        "Heart disease" = "Heart_disease", "Cancer" = "Cancer", "Accidents" = "Accidents",
        "Stroke" = "Stroke", "Chronic respiratory" = "Chronic_respiratory", 
        "Alzheimer's" = "Alzheimers", "Diabetes" = "Diabetes",
        "Influenza/Pneumonia" = "Influenza_Pneumonia", "Kidney disease" = "Kidney_disease",
        "Suicide" = "Suicide", "Septicemia" = "Septicemia", 
        "Chronic liver disease" = "Chronic_liver_disease", "Hypertension" = "Hypertension",
        "Parkinson's disease" = "Parkinson_disease", "COVID-19" = "Other", "Other" = "Other"
      )
      
      year_results$secular_trend <- 1.0
      for (cause in names(causes)) {
        mask <- year_results$cause == cause
        if (any(mask)) {
          trend_col <- cause_to_trend[[cause]]
          if (!is.null(trend_col)) {
            year_results$secular_trend[mask] <- base_trend[[trend_col]]
          }
        }
      }
      
      # Calculate deaths
      weeks_per_year <- 365.25 / 7
      year_results$base_deaths <- (year_results$population * year_results$age_death_rate / 1000 / weeks_per_year)
      year_results$adjusted_deaths <- year_results$base_deaths * year_results$cause_prop * 
        year_results$seasonal_mult * year_results$cause_seasonal * 
        year_results$secular_trend * year_results$covid_mult
      
      # Different seeds for different scenarios in post-2020 period
      scenario_seed_offset <- ifelse(scenario_name == "COVID", 100, 200)
      set.seed(year * 1000 + scenario_seed_offset + seed_set)
      variations <- rnorm(nrow(year_results), mean = 0, sd = variation_scale)
      year_results$final_deaths <- pmax(0, round(year_results$adjusted_deaths * (1 + variations)))
      
      # Keep only non-zero deaths
      year_final <- year_results[year_results$final_deaths > 0, 
                                 c("year", "week", "age_group", "cause", "final_deaths", "population", "week_start_date")]
      names(year_final)[5] <- "deaths"
      
      post_data <- rbind(post_data, year_final)
    }
    
    return(post_data)
  }
  
  # Generate post-2020 data for both scenarios
  covid_post_2020 <- generate_post_2020("COVID")
  no_covid_post_2020 <- generate_post_2020("No COVID")
  
  # STEP 3: Combine pre-2020 (identical) with post-2020 (different) for each scenario
  covid_scenario <- rbind(
    pre_2020_data %>% mutate(scenario = "COVID"),
    covid_post_2020 %>% mutate(scenario = "COVID")
  ) %>% mutate(deaths = as.integer(deaths)) # ensure counts
  
  no_covid_scenario <- rbind(
    pre_2020_data %>% mutate(scenario = "No COVID"),  # Same pre-2020 data!
    no_covid_post_2020 %>% mutate(scenario = "No COVID")
  ) %>% mutate(deaths = as.integer(deaths)) # ensure counts
  
  #cat("Final datasets created.")
  
  return(list(
    covid_scenario = covid_scenario,
    no_covid_scenario = no_covid_scenario,
    divergence_year = divergence_year
  ))
}

aggregate_to_monthly_totals <- function(weekly_data, date_col = "week_start_date", deaths_col = "deaths") {
  weekly_data %>%
    mutate(
      month = month(.data[[date_col]]),
      year_month = as.Date(paste(.data[["year"]], sprintf("%02d", month(.data[[date_col]])), "01", sep = "-"))
    ) %>%
    group_by(scenario, year, month, year_month) %>%
    summarise(
      total_deaths = sum(.data[[deaths_col]], na.rm = TRUE),
      weeks_in_month = n_distinct(week),
      .groups = 'drop'
    ) %>%
    mutate(
      month_name = month.name[month],
    ) %>%
    select(scenario, year, month, month_name, year_month, 
           total_deaths, weeks_in_month) %>%
    arrange(scenario, year, month)
}


population <- data.frame(year=c(rep(2014, 12), rep(2015, 12), rep(2016, 12), rep(2017, 12), rep(2018, 12), rep(2019, 12),
                                rep(2020, 12), rep(2021, 12), rep(2022, 12), rep(2023, 12), rep(2024, 12)),
                         pop=c(rep(318301008, 12), rep(320635163, 12), rep(322941311, 12),
                               rep(324985539, 12), rep(326687501, 12), rep(328239523, 12),
                               rep(331577720, 12), rep(332099760, 12), rep(334017321, 12),
                               rep(336806231, 12), rep(340110988, 12)))
#### MORTALITY CURVE VISUALIZATIONS ####

results <- simulate_mortality(2014, 2023, 2020, seed_set = 1)
combined2 <- rbind(results$covid_scenario %>% mutate(scenario = "COVID"), 
                   results$no_covid_scenario %>% mutate(scenario = "No COVID"))
weekly_totals  <- combined2 %>% group_by(scenario, year, week, week_start_date) %>%
  summarise(total_deaths = sum(deaths), .groups = 'drop')

results_hypo1 <- simulate_mortality_hypo1(2014, 2023, 2020, seed_set = 1)
combined2 <- rbind(results_hypo1$covid_scenario %>% mutate(scenario = "COVID"), 
                   results_hypo1$no_covid_scenario %>% mutate(scenario = "No COVID"))
weekly_totals_hypo1 <- combined2 %>% group_by(scenario, year, week, week_start_date) %>%
  summarise(total_deaths = sum(deaths), .groups = 'drop')

results_hypo2 <- simulate_mortality_hypo2(2014, 2023, 2020, seed_set = 1)
combined2 <- rbind(results_hypo2$covid_scenario %>% mutate(scenario = "COVID"), 
                   results_hypo2$no_covid_scenario %>% mutate(scenario = "No COVID"))
weekly_totals_hypo2 <- combined2 %>% group_by(scenario, year, week, week_start_date) %>%
  summarise(total_deaths = sum(deaths), .groups = 'drop')


ggplot(rbind.data.frame(weekly_totals %>% filter(scenario=="Without COVID") %>% mutate(Scenario="Scenario 1"),
                        weekly_totals_hypo1 %>% filter(scenario=="No COVID") %>% mutate(Scenario="Scenario 2"),
                        weekly_totals_hypo2 %>% filter(scenario=="No COVID") %>% mutate(Scenario="Scenario 3"),
                        weekly_totals %>% filter(scenario=="With COVID") %>% mutate(Scenario="Scenario 1 (COVID)"),
                        weekly_totals_hypo1 %>% filter(scenario=="COVID") %>% mutate(Scenario="Scenario 2 (COVID)"),
                        weekly_totals_hypo2 %>% filter(scenario=="COVID") %>% mutate(Scenario="Scenario 3 (COVID)")), 
       aes(x=week_start_date, y=total_deaths, color = Scenario)) + geom_line(alpha=0.7, size=1)+
  geom_vline(xintercept = as.Date("2020-01-01"), linetype = "dashed", alpha = 0.5, size=0.75) +
  geom_label(data = data.frame(x = as.Date("2014-01-01"), y = 104000, label = "Scenario 3: Older population, milder pandemic"), 
             aes(x = x, y = y, label = label), 
             color = "black", fill = "white", hjust = 0,
             linewidth=0.5, size = 4.5) +
  geom_label(data = data.frame(x = as.Date("2014-01-01"), y = 39000, label = "Scenario 2: Younger population, deadlier pandemic"), 
             aes(x = x, y = y, label = label), 
             color = "black", fill = "white", hjust = 0,
             linewidth=0.5, size = 4.5) +
  geom_label(data = data.frame(x = as.Date("2014-01-01"), y = 67000, label = "Scenario 1: Realistic for United States"), 
             aes(x = x, y = y, label = label), 
             color = "black", fill = "white", hjust = 0,
             linewidth=0.5, size = 4.5) +
  geom_label(data = data.frame(x = as.Date("2019-03-01"), y = 130000, label = "COVID-19 pandemic begins"), 
             aes(x = x, y = y, label = label), 
             color = "black", fill = "white", hjust = 0,
             linewidth=0.5, size = 4.5) +
  labs(  
    x = "Date", y = "Total deaths", color = "Scenario") +# scale_color_manual(values=c("coral1", "palegreen3", "black")) +
  theme_minimal() +  theme(plot.title = element_text(size = 14),
                           plot.subtitle = element_text(size=12),
                           axis.title.x = element_text(size = 13), axis.title.y = element_text(size = 13), 
                           axis.text = element_text(size = 10),
                           legend.text=element_text(size=11), legend.title=element_text(size=12)) + guides(color="none")


#### Functions to calculate performance metrics for simulations ####

# yearly functions:

get_prmse_y <- function(x, y){
  prmse_list <- list()
  for (n in 1:N){
    meths <- x[[n]]$y_est %>% filter(year == y) %>% select(all_of(methods)) 
    truth <- x[[n]]$y_est %>% filter(year == y) %>% select(truth) %>% unlist %>% as.vector
    prmse_list[[n]] <- apply(meths, 2, function(k) (k - truth)^2)
  }
  prmse_y <- do.call(rbind.data.frame, prmse_list) 
  prmse_y_df <- apply(prmse_y, 2, median) %>% sqrt
  names(prmse_y_df) <- methods
  return(prmse_y_df)
}
get_emp_sd_y <- function(x, y){
  y_est_list <- list()
  for (n in 1:N){
    y_est_list[[n]] <- x[[n]]$y_est %>% filter(year == y) %>% select(all_of(methods)) 
  }
  res <- do.call(rbind.data.frame, y_est_list)
  return(apply(res, 2, sd))
}
get_mape_y <- function(x, y){
  mape_list <- list()
  for (n in 1:N){
    meths <- x[[n]]$y_est %>% filter(year == y) %>% select(all_of(methods)) 
    truth <- x[[n]]$y_est %>% filter(year == y) %>% select(truth) %>% unlist %>% as.vector
    mape_list[[n]] <- abs((meths - truth)/truth)*100
  }
  return(apply(do.call(rbind.data.frame, mape_list),2,median))
  #return(do.call(rbind.data.frame, mape_list))
}

## monthly functions:

get_prmse_m <- function(x, y){
  avg_m_list <- c()
  for (m in 1:12){
    m_est_list <- list()
    for (n in 1:N){
      meths <- x[[n]]$m_est %>% filter(year == y, month == m) %>% select(all_of(methods))
      truth <- x[[n]]$m_est %>% filter(year == y, month == m) %>% select(truth) %>% unlist %>% as.vector
      m_est_list[[n]] <- apply(meths, 2, function(x) (x-truth)^2)
    }
    # m_est_list is the list of sq. pred. errors for one month across all runs
    hold_df <- do.call(rbind.data.frame, m_est_list)
    avg_m_list[[m]] <- apply(hold_df, 2, median) # get median across runs for month t
  }
  hold2 <- do.call(rbind.data.frame, avg_m_list) ; colnames(hold2) <- methods
  return(sqrt(colMeans(hold2)))
}
get_emp_sd_m <- function(x, y){
  avg_m_list <- c()
  for (m in 1:12){
    m_est_list <- list()
    for (n in 1:N){
      m_est_list[[n]] <- x[[n]]$m_est %>% filter(year == y, month == m) %>% select(all_of(methods))
    }
    
    hold_df <- do.call(rbind.data.frame, m_est_list)
    avg_m_list[[m]] <- apply(hold_df, 2, var) 
    # each entry is 1/(k-1)*sum(delta_kmt - bar delta_kmt)^2 for one t
  }
  hold2 <- do.call(rbind.data.frame, avg_m_list) ; colnames(hold2) <- methods
  return(sqrt(colMeans(hold2)))
}
get_mape_m <- function(x, y){
  mape_list <- c()
  for (m in 1:12){
    m_est_list <- list()
    for (n in 1:N){
      meths <- x[[n]]$m_est %>% filter(year == y, month == m) %>% select(all_of(methods))
      truth <- x[[n]]$m_est %>% filter(year == y, month == m) %>% select(truth) %>% unlist %>% as.vector
      m_est_list[[n]] <- apply(meths, 2, function(x) abs((x - truth)/truth)*100)
    }
    hold_df <- do.call(rbind.data.frame, m_est_list)
    mape_list[[m]] <- apply(hold_df, 2, median) 
    # each entry is the median abs. pred. error for month t across runs. 
  }
  hold2 <- do.call(rbind.data.frame, mape_list) ; colnames(hold2) <- methods
  return(colMeans(hold2))}

# For formatting data tables with metrics relative to WMD:
format_data <- function(df, df_n){
  out <- df_n
  for (i in 1:3){
    out[1, i] <- paste0(round(df[1,i],2), " (1.00)")
  }
  return(out)
}

#### SCENARIO 1: 2014 to 2017, 2020 and 2021 ####
N <- 1000
sim_loop_sc1 <- function(n, start_year, end_year){
  
  # Simulate data
  results <- simulate_mortality(start_year, end_year, 2020, seed_set = n) # Documented seeds
  obs_data_m <- aggregate_to_monthly_totals(results$covid_scenario)
  cf_data_m <- aggregate_to_monthly_totals(results$no_covid_scenario)
  
  # Set-up for method estimates
  obs_death_m <- obs_data_m %>% filter(year > 2019) %>% select(total_deaths) %>% unname %>% unlist
  cf_death_m <- cf_data_m %>% filter(year > 2019) %>% select(total_deaths) %>% unname %>% unlist
  obs_death_y <- obs_data_m %>% filter(year > 2019) %>% aggregate(total_deaths ~ year, sum) %>% select(total_deaths) %>% unname %>% unlist
  cf_death_y <- cf_data_m %>% filter(year > 2019) %>% aggregate(total_deaths ~ year, sum) %>% select(total_deaths) %>% unname %>% unlist
  get_pop <-  population %>% filter(year >= start_year & year <= end_year)
  obs_data_m$population <- get_pop[,"pop"]
  
  if (end_year==2020){
  year_shell <- data.frame(year = rep(2020, 12), month = 1:12)}
  
  if (end_year==2021){
    year_shell <- data.frame(year = c(rep(2020, 12), rep(2021, 12)), month = rep(1:12, 2))
  }
  
  if (end_year==2023){
    year_shell <- data.frame(year = c(rep(2020, 12), rep(2021, 12),
                                      rep(2022, 12), rep(2023, 12)), month = rep(1:12, 4))
  }
  
  # WMD method: Baseline method against which methods are compared
  
  wmd_mod <- lm(total_deaths ~ factor(month) + year, data=obs_data_m %>% filter(year < 2020))
  wmd_est_m <- cbind.data.frame(year_shell,
                                predict(wmd_mod, newdata = obs_data_m %>% filter(year > 2019)))
  colnames(wmd_est_m) <- c("year", "month", "estimate")
  wmd_est_y <- wmd_est_m %>% aggregate(estimate ~ year, sum)
  
  # WHO method

  model <- gam(total_deaths ~ s(year, k=3, bs="tp") + 
                   s(month, bs = "cc"),
                 data=obs_data_m %>% filter(year <= 2019),
                 family=nb(theta = NULL, link= "log"))
  
  who_est_m <- cbind.data.frame(year_shell, 
                                predict(model, se.fit=FALSE, type="response", newdata = obs_data_m %>% filter(year > 2019)))
  colnames(who_est_m) <- c("year", "month", "estimate")
  who_est_y <- who_est_m %>% aggregate(estimate ~ year, sum)
  
  # The Economist method
  
  obs_data_m$days_in_month <- days_in_month(obs_data_m$year_month)
  obs_data_m$deaths_per_day <- obs_data_m$total_deaths / obs_data_m$days_in_month
  expected_mod <- lm(deaths_per_day ~  year + factor(month), data=obs_data_m %>% filter(year <= 2019))
  
  months <- c(obs_data_m %>% filter(year == 2020) %>% select(days_in_month) %>% as.matrix %>% as.vector)
  econ_est_m <- data.frame(year_shell, predict(expected_mod, newdata = obs_data_m %>% filter(year >= 2020))*months)
  colnames(econ_est_m) <- c("year", "month", "estimate") 
  econ_est_y <- econ_est_m %>% aggregate(estimate ~ year, sum)
  
  # IHME method (most intensive)
  
  # Step 1
  
  # Data preparation:
  obs_data_m$log_pop <- log(obs_data_m$population)
  obs_data_m$chron_index <- (obs_data_m$year - 2014)*12 + obs_data_m$month # for time trend
  
  # First layer of estimation
  seasonal_mod <- glm(formula = total_deaths ~ bs(month, degree = 3) + log_pop, data = obs_data_m %>% 
                        filter(year_month < "2019-03-01"), family = poisson(link = "log"))
  obs_data_m$preds <- predict(seasonal_mod, newdata=obs_data_m) # only pre-2020 data will be used in training model
  
  # Alternative knot placements
  weight_mod_1 <- glm(formula = total_deaths ~ preds + bs(chron_index, degree = 1, knots = (63 - 6)), # Last knot at 6 months prior
                      data = obs_data_m %>% filter(year_month < "2019-03-01"),
                      family = poisson(link = "log"))
  
  weight_mod_2 <- glm(formula = total_deaths ~ preds + bs(chron_index, degree = 1, knots = (63 - 12)), # Last knot at 12 months prior
                      data = obs_data_m %>% filter(year_month < "2019-03-01"),
                      family = poisson(link = "log"))
  
  weight_mod_3 <- glm(formula = total_deaths ~ preds + bs(chron_index, degree = 1, knots = (63 - 18)), # Last knot at 18 months prior
                      data = obs_data_m %>% filter(year_month < "2019-03-01"),
                      family = poisson(link = "log"))
  
  weight_mod_4 <- glm(formula = total_deaths ~ preds + bs(chron_index, degree = 1, knots = (63 - 24)), # Last knot at 24 months prior
                      data = obs_data_m %>% filter(year_month < "2019-03-01"),
                      family = poisson(link = "log"))
  weight_mod_5 <- glm(formula = total_deaths ~ factor(month) + factor(year), 
                      data =  obs_data_m %>% filter(year_month < "2019-03-01"),
                      family = poisson(link = "log"))
  
  # As the data is new every time, we need to recalculate weights every time. Weights calculated based on 
  # predicting Match to December 2019.
  
  # Knot-model predictions:
  pred_1 <- predict(weight_mod_1, newdata = obs_data_m %>%  filter(year_month > "2019-02-01" & year_month < "2020-01-01"))
  pred_2 <- predict(weight_mod_2, newdata = obs_data_m %>%  filter(year_month > "2019-02-01" & year_month < "2020-01-01"))
  pred_3 <- predict(weight_mod_3, newdata = obs_data_m %>%  filter(year_month > "2019-02-01" & year_month < "2020-01-01"))
  pred_4 <- predict(weight_mod_4, newdata = obs_data_m %>%  filter(year_month > "2019-02-01" & year_month < "2020-01-01"))
  
  # Dummy variable prediction:
  pred_5 <- predict(weight_mod_5, newdata = obs_data_m %>%  filter(year_month > "2019-02-01" & year_month < "2020-01-01"))
  
  # Last year carried forward:
  pred_ly <- obs_data_m %>% filter(year_month > "2018-02-01" & year_month < "2019-01-01") %>% select(total_deaths)
  
  # See Wang et al. for details on the calculation of these weights.
  
  get_weight <- function(pred){
    errors <- (pred %>% exp) - (obs_data_m %>%  filter(year_month > "2019-02-01" & year_month < "2020-01-01") %>% select(total_deaths))
    rmse <- sqrt( sum( (errors / 328239523) ^2 )/length(pred) )
    return( 1 / (rmse^2))}
  
  errors_ly <- pred_ly - obs_data_m %>% filter(year_month > "2019-02-01" & year_month < "2020-01-01") %>% select(total_deaths)
  rmse_ly <- sqrt( sum( (errors_ly / 328239523) ^2 )/9 )
  weights <- c(1/(rmse_ly^2), get_weight(pred_1), get_weight(pred_2), get_weight(pred_3), get_weight(pred_4),
               get_weight(pred_5))
  weight_vec <- weights / sum(weights) # ensemble weights for averages.
  
  # Now to the actual prediction:
  
  final_mod_1 <- glm(formula = total_deaths ~ preds + bs(chron_index, degree = 1, knots = (73 - 6)), # Last knot at 6 months prior
                     data = obs_data_m %>% filter(year_month < "2020-01-01"),
                     family = poisson(link = "log"))
  
  final_mod_2 <- glm(formula = total_deaths ~ preds + bs(chron_index, degree = 1, knots = (73 - 12)), # Last knot at 12 months prior
                     data = obs_data_m %>% filter(year_month < "2020-01-01"),
                     family = poisson(link = "log"))
  
  final_mod_3 <- glm(formula = total_deaths ~ preds + bs(chron_index, degree = 1, knots = (73 - 18)), # Last knot at 18 months prior
                     data = obs_data_m %>% filter(year_month < "2020-01-01"),
                     family = poisson(link = "log"))
  
  final_mod_4 <- glm(formula = total_deaths ~ preds + bs(chron_index, degree = 1, knots = (73 - 24)), # Last knot at 24 months prior
                     data = obs_data_m %>% filter(year_month < "2020-01-01"),
                     family = poisson(link = "log"))
  final_mod_5 <- glm(formula = total_deaths ~ factor(month) + factor(year), 
                     data =  obs_data_m %>% filter(year_month < "2020-03-01"),
                     family = poisson(link = "log"))
  
  
  pred_ly <-obs_data_m %>%  filter(year_month >= "2019-01-01" & year_month <= "2019-12-01") %>% select(total_deaths) %>% unlist %>% as.vector
  pred_1 <- predict(final_mod_1, newdata = obs_data_m %>%  filter(year_month >= "2020-01-01"))
  pred_2 <- predict(final_mod_2, newdata = obs_data_m %>%  filter(year_month >= "2020-01-01"))
  pred_3 <- predict(final_mod_3, newdata = obs_data_m %>%  filter(year_month >= "2020-01-01")) 
  pred_4 <- predict(final_mod_4, newdata = obs_data_m %>%  filter(year_month >= "2020-01-01"))
  pred_5 <- predict(final_mod_5, newdata = obs_data_m %>%  filter(year==2020))
  # year effect for 2021 is equal to 2020, so the same predictions are made for the 5th model for both years.
  
  # Predictions at the monthly level:
  preds_df <- data.frame(ly = pred_ly, knot_1 = pred_1 %>% exp, knot_2 = pred_2 %>% exp, 
                         knot_3 = pred_3 %>% exp, knot_4 = pred_4 %>% exp,
                         dummy = pred_5 %>% exp)
  
  ihme_est_m <- data.frame(year_shell, estimate=apply(preds_df, 1, function(x) sum(x*weight_vec)))
  ihme_est_y <- ihme_est_m %>% aggregate(estimate ~ year, sum)
  
  # Acosta and Irizarry method
  
  obs_data_m <- obs_data_m %>% rename(outcome=total_deaths, date=year_month)
  obs_data_m$outcome <- as.integer(obs_data_m$outcome)
  
  exclude_dates_v1 <- c(
    make_date(2015, 1, 1),
    make_date(2017, 1, 1),
    make_date(2018, 1, 1),
    seq(make_date(2020, 1, 1), make_date(2024, 12, 1), by = "month"))
  
  # Flexible model to account for tapering off trend
  
  if (2020 - start_year < 5){
    ai_ed <- compute_expected(obs_data_m, exclude = exclude_dates_v1, harmonics = 2, 
                              trend.knots.per.year = 1/4, frequency = 12, weekday.effect = FALSE, 
                              extrapolate=FALSE, verbose=FALSE, include.trend=FALSE)   
  }else{
    ai_ed <- compute_expected(obs_data_m, exclude = exclude_dates_v1, harmonics = 2, 
                              trend.knots.per.year = 1/4, frequency = 12, weekday.effect = FALSE, 
                              extrapolate=FALSE, verbose=FALSE)   
  }
   
  if (end_year == 2020){
  model <- excess_model(ai_ed, start = make_date(2014, 1, 1), end = make_date(2020, 12, 31),
                        exclude = c(make_date(2015, 1, 1),
                                    make_date(2017, 1, 1), 
                                    make_date(2018, 1, 1)), # exclude months with abnormal flu
                        knots.per.year = 3,
                        intervals = seq(make_date(2020, 1, 1), make_date(2020, 12, 31), by = "month"))}
  if (end_year > 2020){
    model <- excess_model(ai_ed, start = make_date(2014, 1, 1), end = make_date(2021, 12, 31),
                          exclude = c(make_date(2015, 1, 1),
                                      make_date(2017, 1, 1), 
                                      make_date(2018, 1, 1)), # exclude months with abnormal flu
                          knots.per.year = 3,
                          intervals = seq(make_date(2020, 1, 1), make_date(2023, 12, 31), by = "month"))
  }
  
  ai_est_m <- model$excess %>% select(start, excess) %>% mutate(year=year(start), month=month(start)) %>% filter(year > 2019) %>% 
    filter(excess != 0) %>% select(year, month, excess)
  colnames(ai_est_m) <- c("year", "month", "estimate")
  
  ai_est_y <- ai_est_m %>% aggregate(estimate ~ year, sum)
  
  ## Storing results from each method for each different scale.
  
  m_est <- cbind.data.frame(year_shell, 
                            truth = obs_death_m - cf_death_m, # true excess deaths
                            wmd_est = obs_death_m - wmd_est_m$estimate, # baseline
                            who_est = obs_death_m - who_est_m$estimate, 
                            econ_est = obs_death_m - econ_est_m$estimate,
                            ai_est = ai_est_m$estimate,
                            ihme_est = obs_death_m - ihme_est_m$estimate)
  
  y_est <- cbind.data.frame(year = 2020:end_year, 
                            truth = obs_death_y - cf_death_y, # true excess deaths
                            wmd_est = obs_death_y - wmd_est_y$estimate, # baseline
                            who_est = obs_death_y - who_est_y$estimate, 
                            econ_est = obs_death_y - econ_est_y$estimate,
                            ai_est = ai_est_y$estimate,
                            ihme_est = obs_death_y - ihme_est_y$estimate)
  return(list(m_est = m_est, y_est = y_est))
}

methods <- c("wmd_est", "who_est", "econ_est", "ai_est", "ihme_est")

### Scenario 1: 2014 to 2021
sc1_2014 <- list()
for (n in 1:N){
  sc1_2014[[n]] <- suppressWarnings(sim_loop_sc1(n, start_year = 2014, end_year = 2023))
  print(paste0("Simulation run ", n, " completed."))
}

sc1_2014_2020_m <- cbind(get_prmse_m(sc1_2014, 2020), get_mape_m(sc1_2014, 2020), get_emp_sd_m(sc1_2014, 2020)) 
sc1_2014_2020_y <- cbind(get_prmse_y(sc1_2014, 2020), get_mape_y(sc1_2014, 2020), get_emp_sd_y(sc1_2014, 2020)) 

sc1_2014_2020_m_rel <- apply(sc1_2014_2020_m, 1, function(x) sc1_2014_2020_m[1,]/x) %>% t %>% round(2)
sc1_2014_2020_y_rel <- apply(sc1_2014_2020_y, 1, function(x) sc1_2014_2020_y[1,]/x) %>% t %>% round(2)

sc1_2014_2021_m <- cbind(get_prmse_m(sc1_2014, 2021), get_mape_m(sc1_2014, 2021), get_emp_sd_m(sc1_2014, 2021)) 
sc1_2014_2021_y <- cbind(get_prmse_y(sc1_2014, 2021), get_mape_y(sc1_2014, 2021), get_emp_sd_y(sc1_2014, 2021)) 

sc1_2014_2021_m_rel <- apply(sc1_2014_2021_m, 1, function(x) sc1_2014_2021_m[1,]/x) %>% t %>% round(2)
sc1_2014_2021_y_rel <- apply(sc1_2014_2021_y, 1, function(x) sc1_2014_2021_y[1,]/x) %>% t %>% round(2)

sc1_2014_2022_m <- cbind(get_prmse_m(sc1_2014, 2022), get_mape_m(sc1_2014, 2022), get_emp_sd_m(sc1_2014, 2022)) 
sc1_2014_2022_y <- cbind(get_prmse_y(sc1_2014, 2022), get_mape_y(sc1_2014, 2022), get_emp_sd_y(sc1_2014, 2022)) 

sc1_2014_2022_m_rel <- apply(sc1_2014_2022_m, 1, function(x) sc1_2014_2022_m[1,]/x) %>% t %>% round(2)
sc1_2014_2022_y_rel <- apply(sc1_2014_2022_y, 1, function(x) sc1_2014_2022_y[1,]/x) %>% t %>% round(2)

sc1_2014_2023_m <- cbind(get_prmse_m(sc1_2014, 2023), get_mape_m(sc1_2014, 2023), get_emp_sd_m(sc1_2014, 2023)) 
sc1_2014_2023_y <- cbind(get_prmse_y(sc1_2014, 2023), get_mape_y(sc1_2014, 2023), get_emp_sd_y(sc1_2014, 2023)) 

sc1_2014_2023_m_rel <- apply(sc1_2014_2023_m, 1, function(x) sc1_2014_2023_m[1,]/x) %>% t %>% round(2)
sc1_2014_2023_y_rel <- apply(sc1_2014_2023_y, 1, function(x) sc1_2014_2023_y[1,]/x) %>% t %>% round(2)

#### Scenario 1: 2017 to 2021 
sc1_2017 <- list()
for (n in 1:N){
  sc1_2017[[n]] <- suppressWarnings(sim_loop_sc1(n, start_year = 2017, end_year = 2023))
  print(paste0("Simulation run ", n, " completed."))
}

sc1_2017_2020_m <- cbind(get_prmse_m(sc1_2017, 2020), get_mape_m(sc1_2017, 2020), get_emp_sd_m(sc1_2017, 2020)) 
sc1_2017_2020_y <- cbind(get_prmse_y(sc1_2017, 2020), get_mape_y(sc1_2017, 2020), get_emp_sd_y(sc1_2017, 2020)) 

sc1_2017_2020_m_rel <- apply(sc1_2017_2020_m, 1, function(x) sc1_2017_2020_m[1,]/x) %>% t %>% round(2)
sc1_2017_2020_y_rel <- apply(sc1_2017_2020_y, 1, function(x) sc1_2017_2020_y[1,]/x) %>% t %>% round(2)

sc1_2017_2021_m <- cbind(get_prmse_m(sc1_2017, 2021), get_mape_m(sc1_2017, 2021), get_emp_sd_m(sc1_2017, 2021)) 
sc1_2017_2021_y <- cbind(get_prmse_y(sc1_2017, 2021), get_mape_y(sc1_2017, 2021), get_emp_sd_y(sc1_2017, 2021)) 

sc1_2017_2021_m_rel <- apply(sc1_2017_2021_m, 1, function(x) sc1_2017_2021_m[1,]/x) %>% t %>% round(2)
sc1_2017_2021_y_rel <- apply(sc1_2017_2021_y, 1, function(x) sc1_2017_2021_y[1,]/x) %>% t %>% round(2)

sc1_2017_2022_m <- cbind(get_prmse_m(sc1_2017, 2022), get_mape_m(sc1_2017, 2022), get_emp_sd_m(sc1_2017, 2022)) 
sc1_2017_2022_y <- cbind(get_prmse_y(sc1_2017, 2022), get_mape_y(sc1_2017, 2022), get_emp_sd_y(sc1_2017, 2022)) 

sc1_2017_2022_m_rel <- apply(sc1_2017_2022_m, 1, function(x) sc1_2017_2022_m[1,]/x) %>% t %>% round(2)
sc1_2017_2022_y_rel <- apply(sc1_2017_2022_y, 1, function(x) sc1_2017_2022_y[1,]/x) %>% t %>% round(2)

sc1_2017_2023_m <- cbind(get_prmse_m(sc1_2017, 2023), get_mape_m(sc1_2017, 2023), get_emp_sd_m(sc1_2017, 2023)) 
sc1_2017_2023_y <- cbind(get_prmse_y(sc1_2017, 2023), get_mape_y(sc1_2017, 2023), get_emp_sd_y(sc1_2017, 2023)) 

sc1_2017_2023_m_rel <- apply(sc1_2017_2023_m, 1, function(x) sc1_2017_2023_m[1,]/x) %>% t %>% round(2)
sc1_2017_2023_y_rel <- apply(sc1_2017_2023_y, 1, function(x) sc1_2017_2023_y[1,]/x) %>% t %>% round(2)

#### SCENARIO 2: 2014 to 2017, 2020 and 2021 ####

sim_loop_sc2 <- function(n, start_year, end_year){
  
  # Simulate data
  results <- simulate_mortality_hypo1(start_year, end_year, 2020, seed_set = n) # Documented seeds
  obs_data_m <- aggregate_to_monthly_totals(results$covid_scenario)
  cf_data_m <- aggregate_to_monthly_totals(results$no_covid_scenario)
  
  # Set-up for method estimates
  obs_death_m <- obs_data_m %>% filter(year > 2019) %>% select(total_deaths) %>% unname %>% unlist
  cf_death_m <- cf_data_m %>% filter(year > 2019) %>% select(total_deaths) %>% unname %>% unlist
  obs_death_y <- obs_data_m %>% filter(year > 2019) %>% aggregate(total_deaths ~ year, sum) %>% select(total_deaths) %>% unname %>% unlist
  cf_death_y <- cf_data_m %>% filter(year > 2019) %>% aggregate(total_deaths ~ year, sum) %>% select(total_deaths) %>% unname %>% unlist
  get_pop <-  population %>% filter(year >= start_year & year <= end_year)
  obs_data_m$population <- get_pop[,"pop"]
  
  if (end_year==2020){
    year_shell <- data.frame(year = rep(2020, 12), month = 1:12)}
  
  if (end_year==2021){
    year_shell <- data.frame(year = c(rep(2020, 12), rep(2021, 12)), month = rep(1:12, 2))
  }
  
  if (end_year==2023){
    year_shell <- data.frame(year = c(rep(2020, 12), rep(2021, 12),
                                      rep(2022, 12), rep(2023, 12)), month = rep(1:12, 4))
  }
  
  # WMD method: Baseline method against which methods are compared
  
  wmd_mod <- lm(total_deaths ~ factor(month) + year, data=obs_data_m %>% filter(year < 2020))
  wmd_est_m <- cbind.data.frame(year_shell,
                                predict(wmd_mod, newdata = obs_data_m %>% filter(year > 2019)))
  colnames(wmd_est_m) <- c("year", "month", "estimate")
  wmd_est_y <- wmd_est_m %>% aggregate(estimate ~ year, sum)
  
  # WHO method
  
  model <- gam(total_deaths ~ s(year, k=3, bs="tp") + 
                 s(month, bs = "cc"),
               data=obs_data_m %>% filter(year <= 2019),
               family=nb(theta = NULL, link= "log"))
  
  who_est_m <- cbind.data.frame(year_shell, 
                                predict(model, se.fit=FALSE, type="response", newdata = obs_data_m %>% filter(year > 2019)))
  colnames(who_est_m) <- c("year", "month", "estimate")
  who_est_y <- who_est_m %>% aggregate(estimate ~ year, sum)
  
  # The Economist method
  
  obs_data_m$days_in_month <- days_in_month(obs_data_m$year_month)
  obs_data_m$deaths_per_day <- obs_data_m$total_deaths / obs_data_m$days_in_month
  expected_mod <- lm(deaths_per_day ~  year + factor(month), data=obs_data_m %>% filter(year <= 2019))
  
  months <- c(obs_data_m %>% filter(year == 2020) %>% select(days_in_month) %>% as.matrix %>% as.vector)
  econ_est_m <- data.frame(year_shell, predict(expected_mod, newdata = obs_data_m %>% filter(year >= 2020))*months)
  colnames(econ_est_m) <- c("year", "month", "estimate") 
  econ_est_y <- econ_est_m %>% aggregate(estimate ~ year, sum)
  
  # IHME method (most intensive)
  
  # Step 1
  
  # Data preparation:
  obs_data_m$log_pop <- log(obs_data_m$population)
  obs_data_m$chron_index <- (obs_data_m$year - 2014)*12 + obs_data_m$month # for time trend
  
  # First layer of estimation
  seasonal_mod <- glm(formula = total_deaths ~ bs(month, degree = 3) + log_pop, data = obs_data_m %>% 
                        filter(year_month < "2019-03-01"), family = poisson(link = "log"))
  obs_data_m$preds <- predict(seasonal_mod, newdata=obs_data_m) # only pre-2020 data will be used in training model
  
  # Alternative knot placements
  weight_mod_1 <- glm(formula = total_deaths ~ preds + bs(chron_index, degree = 1, knots = (63 - 6)), # Last knot at 6 months prior
                      data = obs_data_m %>% filter(year_month < "2019-03-01"),
                      family = poisson(link = "log"))
  
  weight_mod_2 <- glm(formula = total_deaths ~ preds + bs(chron_index, degree = 1, knots = (63 - 12)), # Last knot at 12 months prior
                      data = obs_data_m %>% filter(year_month < "2019-03-01"),
                      family = poisson(link = "log"))
  
  weight_mod_3 <- glm(formula = total_deaths ~ preds + bs(chron_index, degree = 1, knots = (63 - 18)), # Last knot at 18 months prior
                      data = obs_data_m %>% filter(year_month < "2019-03-01"),
                      family = poisson(link = "log"))
  
  weight_mod_4 <- glm(formula = total_deaths ~ preds + bs(chron_index, degree = 1, knots = (63 - 24)), # Last knot at 24 months prior
                      data = obs_data_m %>% filter(year_month < "2019-03-01"),
                      family = poisson(link = "log"))
  weight_mod_5 <- glm(formula = total_deaths ~ factor(month) + factor(year), 
                      data =  obs_data_m %>% filter(year_month < "2019-03-01"),
                      family = poisson(link = "log"))
  
  # As the data is new every time, we need to recalculate weights every time. Weights calculated based on 
  # predicting Match to December 2019.
  
  # Knot-model predictions:
  pred_1 <- predict(weight_mod_1, newdata = obs_data_m %>%  filter(year_month > "2019-02-01" & year_month < "2020-01-01"))
  pred_2 <- predict(weight_mod_2, newdata = obs_data_m %>%  filter(year_month > "2019-02-01" & year_month < "2020-01-01"))
  pred_3 <- predict(weight_mod_3, newdata = obs_data_m %>%  filter(year_month > "2019-02-01" & year_month < "2020-01-01"))
  pred_4 <- predict(weight_mod_4, newdata = obs_data_m %>%  filter(year_month > "2019-02-01" & year_month < "2020-01-01"))
  
  # Dummy variable prediction:
  pred_5 <- predict(weight_mod_5, newdata = obs_data_m %>%  filter(year_month > "2019-02-01" & year_month < "2020-01-01"))
  
  # Last year carried forward:
  pred_ly <- obs_data_m %>% filter(year_month > "2018-02-01" & year_month < "2019-01-01") %>% select(total_deaths)
  
  # See Wang et al. for details on the calculation of these weights.
  
  get_weight <- function(pred){
    errors <- (pred %>% exp) - (obs_data_m %>%  filter(year_month > "2019-02-01" & year_month < "2020-01-01") %>% select(total_deaths))
    rmse <- sqrt( sum( (errors / 328239523) ^2 )/length(pred) )
    return( 1 / (rmse^2))}
  
  errors_ly <- pred_ly - obs_data_m %>% filter(year_month > "2019-02-01" & year_month < "2020-01-01") %>% select(total_deaths)
  rmse_ly <- sqrt( sum( (errors_ly / 328239523) ^2 )/9 )
  weights <- c(1/(rmse_ly^2), get_weight(pred_1), get_weight(pred_2), get_weight(pred_3), get_weight(pred_4),
               get_weight(pred_5))
  weight_vec <- weights / sum(weights) # ensemble weights for averages.
  
  # Now to the actual prediction:
  
  final_mod_1 <- glm(formula = total_deaths ~ preds + bs(chron_index, degree = 1, knots = (73 - 6)), # Last knot at 6 months prior
                     data = obs_data_m %>% filter(year_month < "2020-01-01"),
                     family = poisson(link = "log"))
  
  final_mod_2 <- glm(formula = total_deaths ~ preds + bs(chron_index, degree = 1, knots = (73 - 12)), # Last knot at 12 months prior
                     data = obs_data_m %>% filter(year_month < "2020-01-01"),
                     family = poisson(link = "log"))
  
  final_mod_3 <- glm(formula = total_deaths ~ preds + bs(chron_index, degree = 1, knots = (73 - 18)), # Last knot at 18 months prior
                     data = obs_data_m %>% filter(year_month < "2020-01-01"),
                     family = poisson(link = "log"))
  
  final_mod_4 <- glm(formula = total_deaths ~ preds + bs(chron_index, degree = 1, knots = (73 - 24)), # Last knot at 24 months prior
                     data = obs_data_m %>% filter(year_month < "2020-01-01"),
                     family = poisson(link = "log"))
  final_mod_5 <- glm(formula = total_deaths ~ factor(month) + factor(year), 
                     data =  obs_data_m %>% filter(year_month < "2020-03-01"),
                     family = poisson(link = "log"))
  
  
  pred_ly <-obs_data_m %>%  filter(year_month >= "2019-01-01" & year_month <= "2019-12-01") %>% select(total_deaths) %>% unlist %>% as.vector
  pred_1 <- predict(final_mod_1, newdata = obs_data_m %>%  filter(year_month >= "2020-01-01"))
  pred_2 <- predict(final_mod_2, newdata = obs_data_m %>%  filter(year_month >= "2020-01-01"))
  pred_3 <- predict(final_mod_3, newdata = obs_data_m %>%  filter(year_month >= "2020-01-01")) 
  pred_4 <- predict(final_mod_4, newdata = obs_data_m %>%  filter(year_month >= "2020-01-01"))
  pred_5 <- predict(final_mod_5, newdata = obs_data_m %>%  filter(year==2020))
  # year effect for 2021 is equal to 2020, so the same predictions are made for the 5th model for both years.
  
  # Predictions at the monthly level:
  preds_df <- data.frame(ly = pred_ly, knot_1 = pred_1 %>% exp, knot_2 = pred_2 %>% exp, 
                         knot_3 = pred_3 %>% exp, knot_4 = pred_4 %>% exp,
                         dummy = pred_5 %>% exp)
  
  ihme_est_m <- data.frame(year_shell, estimate=apply(preds_df, 1, function(x) sum(x*weight_vec)))
  ihme_est_y <- ihme_est_m %>% aggregate(estimate ~ year, sum)
  
  # Acosta and Irizarry method
  
  obs_data_m <- obs_data_m %>% rename(outcome=total_deaths, date=year_month)
  obs_data_m$outcome <- as.integer(obs_data_m$outcome)
  
  exclude_dates_v1 <- c(
    make_date(2015, 1, 1),
    make_date(2017, 1, 1),
    make_date(2018, 1, 1),
    seq(make_date(2020, 1, 1), make_date(2024, 12, 1), by = "month"))
  
  # Flexible model to account for tapering off trend
  
  if (2020 - start_year < 5){
    ai_ed <- compute_expected(obs_data_m, exclude = exclude_dates_v1, harmonics = 2, 
                              trend.knots.per.year = 1/4, frequency = 12, weekday.effect = FALSE, 
                              extrapolate=FALSE, verbose=FALSE, include.trend=FALSE)   
  }else{
    ai_ed <- compute_expected(obs_data_m, exclude = exclude_dates_v1, harmonics = 2, 
                              trend.knots.per.year = 1/4, frequency = 12, weekday.effect = FALSE, 
                              extrapolate=FALSE, verbose=FALSE)   
  }
  
  if (end_year == 2020){
    model <- excess_model(ai_ed, start = make_date(2014, 1, 1), end = make_date(2020, 12, 31),
                          exclude = c(make_date(2015, 1, 1),
                                      make_date(2017, 1, 1), 
                                      make_date(2018, 1, 1)), # exclude months with abnormal flu
                          knots.per.year = 3,
                          intervals = seq(make_date(2020, 1, 1), make_date(2020, 12, 31), by = "month"))}
  if (end_year > 2020){
    model <- excess_model(ai_ed, start = make_date(2014, 1, 1), end = make_date(2021, 12, 31),
                          exclude = c(make_date(2015, 1, 1),
                                      make_date(2017, 1, 1), 
                                      make_date(2018, 1, 1)), # exclude months with abnormal flu
                          knots.per.year = 3,
                          intervals = seq(make_date(2020, 1, 1), make_date(2023, 12, 31), by = "month"))
  }
  
  ai_est_m <- model$excess %>% select(start, excess) %>% mutate(year=year(start), month=month(start)) %>% filter(year > 2019) %>% 
    filter(excess != 0) %>% select(year, month, excess)
  colnames(ai_est_m) <- c("year", "month", "estimate")
  
  ai_est_y <- ai_est_m %>% aggregate(estimate ~ year, sum)
  
  ## Storing results from each method for each different scale.
  
  m_est <- cbind.data.frame(year_shell, 
                            truth = obs_death_m - cf_death_m, # true excess deaths
                            wmd_est = obs_death_m - wmd_est_m$estimate, # baseline
                            who_est = obs_death_m - who_est_m$estimate, 
                            econ_est = obs_death_m - econ_est_m$estimate,
                            ai_est = ai_est_m$estimate,
                            ihme_est = obs_death_m - ihme_est_m$estimate)
  
  y_est <- cbind.data.frame(year = 2020:end_year, 
                            truth = obs_death_y - cf_death_y, # true excess deaths
                            wmd_est = obs_death_y - wmd_est_y$estimate, # baseline
                            who_est = obs_death_y - who_est_y$estimate, 
                            econ_est = obs_death_y - econ_est_y$estimate,
                            ai_est = ai_est_y$estimate,
                            ihme_est = obs_death_y - ihme_est_y$estimate)
  return(list(m_est = m_est, y_est = y_est))
}

methods <- c("wmd_est", "who_est", "econ_est", "ai_est", "ihme_est")

#### Scenario 2: 2014 to 2021
sc2_2014 <- list()
for (n in 1:N){
  sc2_2014[[n]] <- suppressWarnings(sim_loop_sc2(n, start_year = 2014, end_year = 2023))
  print(paste0("Simulation run ", n, " completed."))
}

sc2_2014_2020_m <- cbind(get_prmse_m(sc2_2014, 2020), get_mape_m(sc2_2014, 2020), get_emp_sd_m(sc2_2014, 2020)) 
sc2_2014_2020_y <- cbind(get_prmse_y(sc2_2014, 2020), get_mape_y(sc2_2014, 2020), get_emp_sd_y(sc2_2014, 2020)) 

sc2_2014_2020_m_rel <- apply(sc2_2014_2020_m, 1, function(x) sc2_2014_2020_m[1,]/x) %>% t %>% round(2)
sc2_2014_2020_y_rel <- apply(sc2_2014_2020_y, 1, function(x) sc2_2014_2020_y[1,]/x) %>% t %>% round(2)

sc2_2014_2021_m <- cbind(get_prmse_m(sc2_2014, 2021), get_mape_m(sc2_2014, 2021), get_emp_sd_m(sc2_2014, 2021)) 
sc2_2014_2021_y <- cbind(get_prmse_y(sc2_2014, 2021), get_mape_y(sc2_2014, 2021), get_emp_sd_y(sc2_2014, 2021)) 

sc2_2014_2021_m_rel <- apply(sc2_2014_2021_m, 1, function(x) sc2_2014_2021_m[1,]/x) %>% t %>% round(2)
sc2_2014_2021_y_rel <- apply(sc2_2014_2021_y, 1, function(x) sc2_2014_2021_y[1,]/x) %>% t %>% round(2)

sc2_2014_2022_m <- cbind(get_prmse_m(sc2_2014, 2022), get_mape_m(sc2_2014, 2022), get_emp_sd_m(sc2_2014, 2022)) 
sc2_2014_2022_y <- cbind(get_prmse_y(sc2_2014, 2022), get_mape_y(sc2_2014, 2022), get_emp_sd_y(sc2_2014, 2022)) 

sc2_2014_2022_m_rel <- apply(sc2_2014_2022_m, 1, function(x) sc2_2014_2022_m[1,]/x) %>% t %>% round(2)
sc2_2014_2022_y_rel <- apply(sc2_2014_2022_y, 1, function(x) sc2_2014_2022_y[1,]/x) %>% t %>% round(2)

sc2_2014_2023_m <- cbind(get_prmse_m(sc2_2014, 2023), get_mape_m(sc2_2014, 2023), get_emp_sd_m(sc2_2014, 2023)) 
sc2_2014_2023_y <- cbind(get_prmse_y(sc2_2014, 2023), get_mape_y(sc2_2014, 2023), get_emp_sd_y(sc2_2014, 2023)) 

sc2_2014_2023_m_rel <- apply(sc2_2014_2023_m, 1, function(x) sc2_2014_2023_m[1,]/x) %>% t %>% round(2)
sc2_2014_2023_y_rel <- apply(sc2_2014_2023_y, 1, function(x) sc2_2014_2023_y[1,]/x) %>% t %>% round(2)

#### Scenario 2: 2017 to 2021 
sc2_2017 <- list()
for (n in 1:N){
  sc2_2017[[n]] <- suppressWarnings(sim_loop_sc2(n, start_year = 2017, end_year = 2023))
  print(paste0("Simulation run ", n, " completed."))
}

sc2_2017_2020_m <- cbind(get_prmse_m(sc2_2017, 2020), get_mape_m(sc2_2017, 2020), get_emp_sd_m(sc2_2017, 2020)) 
sc2_2017_2020_y <- cbind(get_prmse_y(sc2_2017, 2020), get_mape_y(sc2_2017, 2020), get_emp_sd_y(sc2_2017, 2020)) 

sc2_2017_2020_m_rel <- apply(sc2_2017_2020_m, 1, function(x) sc2_2017_2020_m[1,]/x) %>% t %>% round(2)
sc2_2017_2020_y_rel <- apply(sc2_2017_2020_y, 1, function(x) sc2_2017_2020_y[1,]/x) %>% t %>% round(2)

sc2_2017_2021_m <- cbind(get_prmse_m(sc2_2017, 2021), get_mape_m(sc2_2017, 2021), get_emp_sd_m(sc2_2017, 2021)) 
sc2_2017_2021_y <- cbind(get_prmse_y(sc2_2017, 2021), get_mape_y(sc2_2017, 2021), get_emp_sd_y(sc2_2017, 2021)) 

sc2_2017_2021_m_rel <- apply(sc2_2017_2021_m, 1, function(x) sc2_2017_2021_m[1,]/x) %>% t %>% round(2)
sc2_2017_2021_y_rel <- apply(sc2_2017_2021_y, 1, function(x) sc2_2017_2021_y[1,]/x) %>% t %>% round(2)

sc2_2017_2022_m <- cbind(get_prmse_m(sc2_2017, 2022), get_mape_m(sc2_2017, 2022), get_emp_sd_m(sc2_2017, 2022)) 
sc2_2017_2022_y <- cbind(get_prmse_y(sc2_2017, 2022), get_mape_y(sc2_2017, 2022), get_emp_sd_y(sc2_2017, 2022)) 

sc2_2017_2022_m_rel <- apply(sc2_2017_2022_m, 1, function(x) sc2_2017_2022_m[1,]/x) %>% t %>% round(2)
sc2_2017_2022_y_rel <- apply(sc2_2017_2022_y, 1, function(x) sc2_2017_2022_y[1,]/x) %>% t %>% round(2)

sc2_2017_2023_m <- cbind(get_prmse_m(sc2_2017, 2023), get_mape_m(sc2_2017, 2023), get_emp_sd_m(sc2_2017, 2023)) 
sc2_2017_2023_y <- cbind(get_prmse_y(sc2_2017, 2023), get_mape_y(sc2_2017, 2023), get_emp_sd_y(sc2_2017, 2023)) 

sc2_2017_2023_m_rel <- apply(sc2_2017_2023_m, 1, function(x) sc2_2017_2023_m[1,]/x) %>% t %>% round(2)
sc2_2017_2023_y_rel <- apply(sc2_2017_2023_y, 1, function(x) sc2_2017_2023_y[1,]/x) %>% t %>% round(2)

#### SCENARIO 3: 2014 to 2017, 2020 and 2021 ####

sim_loop_sc3 <- function(n, start_year, end_year){
  
  # Simulate data
  results <- simulate_mortality_hypo2(start_year, end_year, 2020, seed_set = n) # Documented seeds
  obs_data_m <- aggregate_to_monthly_totals(results$covid_scenario)
  cf_data_m <- aggregate_to_monthly_totals(results$no_covid_scenario)
  
  # Set-up for method estimates
  obs_death_m <- obs_data_m %>% filter(year > 2019) %>% select(total_deaths) %>% unname %>% unlist
  cf_death_m <- cf_data_m %>% filter(year > 2019) %>% select(total_deaths) %>% unname %>% unlist
  obs_death_y <- obs_data_m %>% filter(year > 2019) %>% aggregate(total_deaths ~ year, sum) %>% select(total_deaths) %>% unname %>% unlist
  cf_death_y <- cf_data_m %>% filter(year > 2019) %>% aggregate(total_deaths ~ year, sum) %>% select(total_deaths) %>% unname %>% unlist
  get_pop <-  population %>% filter(year >= start_year & year <= end_year)
  obs_data_m$population <- get_pop[,"pop"]
  
  if (end_year==2020){
    year_shell <- data.frame(year = rep(2020, 12), month = 1:12)}
  
  if (end_year==2021){
    year_shell <- data.frame(year = c(rep(2020, 12), rep(2021, 12)), month = rep(1:12, 2))
  }
  
  if (end_year==2023){
    year_shell <- data.frame(year = c(rep(2020, 12), rep(2021, 12),
                                      rep(2022, 12), rep(2023, 12)), month = rep(1:12, 4))
  }
  
  # WMD method: Baseline method against which methods are compared
  
  wmd_mod <- lm(total_deaths ~ factor(month) + year, data=obs_data_m %>% filter(year < 2020))
  wmd_est_m <- cbind.data.frame(year_shell,
                                predict(wmd_mod, newdata = obs_data_m %>% filter(year > 2019)))
  colnames(wmd_est_m) <- c("year", "month", "estimate")
  wmd_est_y <- wmd_est_m %>% aggregate(estimate ~ year, sum)
  
  # WHO method
  
  model <- gam(total_deaths ~ s(year, k=3, bs="tp") + 
                 s(month, bs = "cc"),
               data=obs_data_m %>% filter(year <= 2019),
               family=nb(theta = NULL, link= "log"))
  
  who_est_m <- cbind.data.frame(year_shell, 
                                predict(model, se.fit=FALSE, type="response", newdata = obs_data_m %>% filter(year > 2019)))
  colnames(who_est_m) <- c("year", "month", "estimate")
  who_est_y <- who_est_m %>% aggregate(estimate ~ year, sum)
  
  # The Economist method
  
  obs_data_m$days_in_month <- days_in_month(obs_data_m$year_month)
  obs_data_m$deaths_per_day <- obs_data_m$total_deaths / obs_data_m$days_in_month
  expected_mod <- lm(deaths_per_day ~  year + factor(month), data=obs_data_m %>% filter(year <= 2019))
  
  months <- c(obs_data_m %>% filter(year == 2020) %>% select(days_in_month) %>% as.matrix %>% as.vector)
  econ_est_m <- data.frame(year_shell, predict(expected_mod, newdata = obs_data_m %>% filter(year >= 2020))*months)
  colnames(econ_est_m) <- c("year", "month", "estimate") 
  econ_est_y <- econ_est_m %>% aggregate(estimate ~ year, sum)
  
  # IHME method (most intensive)
  
  # Step 1
  
  # Data preparation:
  obs_data_m$log_pop <- log(obs_data_m$population)
  obs_data_m$chron_index <- (obs_data_m$year - 2014)*12 + obs_data_m$month # for time trend
  
  # First layer of estimation
  seasonal_mod <- glm(formula = total_deaths ~ bs(month, degree = 3) + log_pop, data = obs_data_m %>% 
                        filter(year_month < "2019-03-01"), family = poisson(link = "log"))
  obs_data_m$preds <- predict(seasonal_mod, newdata=obs_data_m) # only pre-2020 data will be used in training model
  
  # Alternative knot placements
  weight_mod_1 <- glm(formula = total_deaths ~ preds + bs(chron_index, degree = 1, knots = (63 - 6)), # Last knot at 6 months prior
                      data = obs_data_m %>% filter(year_month < "2019-03-01"),
                      family = poisson(link = "log"))
  
  weight_mod_2 <- glm(formula = total_deaths ~ preds + bs(chron_index, degree = 1, knots = (63 - 12)), # Last knot at 12 months prior
                      data = obs_data_m %>% filter(year_month < "2019-03-01"),
                      family = poisson(link = "log"))
  
  weight_mod_3 <- glm(formula = total_deaths ~ preds + bs(chron_index, degree = 1, knots = (63 - 18)), # Last knot at 18 months prior
                      data = obs_data_m %>% filter(year_month < "2019-03-01"),
                      family = poisson(link = "log"))
  
  weight_mod_4 <- glm(formula = total_deaths ~ preds + bs(chron_index, degree = 1, knots = (63 - 24)), # Last knot at 24 months prior
                      data = obs_data_m %>% filter(year_month < "2019-03-01"),
                      family = poisson(link = "log"))
  weight_mod_5 <- glm(formula = total_deaths ~ factor(month) + factor(year), 
                      data =  obs_data_m %>% filter(year_month < "2019-03-01"),
                      family = poisson(link = "log"))
  
  # As the data is new every time, we need to recalculate weights every time. Weights calculated based on 
  # predicting Match to December 2019.
  
  # Knot-model predictions:
  pred_1 <- predict(weight_mod_1, newdata = obs_data_m %>%  filter(year_month > "2019-02-01" & year_month < "2020-01-01"))
  pred_2 <- predict(weight_mod_2, newdata = obs_data_m %>%  filter(year_month > "2019-02-01" & year_month < "2020-01-01"))
  pred_3 <- predict(weight_mod_3, newdata = obs_data_m %>%  filter(year_month > "2019-02-01" & year_month < "2020-01-01"))
  pred_4 <- predict(weight_mod_4, newdata = obs_data_m %>%  filter(year_month > "2019-02-01" & year_month < "2020-01-01"))
  
  # Dummy variable prediction:
  pred_5 <- predict(weight_mod_5, newdata = obs_data_m %>%  filter(year_month > "2019-02-01" & year_month < "2020-01-01"))
  
  # Last year carried forward:
  pred_ly <- obs_data_m %>% filter(year_month > "2018-02-01" & year_month < "2019-01-01") %>% select(total_deaths)
  
  # See Wang et al. for details on the calculation of these weights.
  
  get_weight <- function(pred){
    errors <- (pred %>% exp) - (obs_data_m %>%  filter(year_month > "2019-02-01" & year_month < "2020-01-01") %>% select(total_deaths))
    rmse <- sqrt( sum( (errors / 328239523) ^2 )/length(pred) )
    return( 1 / (rmse^2))}
  
  errors_ly <- pred_ly - obs_data_m %>% filter(year_month > "2019-02-01" & year_month < "2020-01-01") %>% select(total_deaths)
  rmse_ly <- sqrt( sum( (errors_ly / 328239523) ^2 )/9 )
  weights <- c(1/(rmse_ly^2), get_weight(pred_1), get_weight(pred_2), get_weight(pred_3), get_weight(pred_4),
               get_weight(pred_5))
  weight_vec <- weights / sum(weights) # ensemble weights for averages.
  
  # Now to the actual prediction:
  
  final_mod_1 <- glm(formula = total_deaths ~ preds + bs(chron_index, degree = 1, knots = (73 - 6)), # Last knot at 6 months prior
                     data = obs_data_m %>% filter(year_month < "2020-01-01"),
                     family = poisson(link = "log"))
  
  final_mod_2 <- glm(formula = total_deaths ~ preds + bs(chron_index, degree = 1, knots = (73 - 12)), # Last knot at 12 months prior
                     data = obs_data_m %>% filter(year_month < "2020-01-01"),
                     family = poisson(link = "log"))
  
  final_mod_3 <- glm(formula = total_deaths ~ preds + bs(chron_index, degree = 1, knots = (73 - 18)), # Last knot at 18 months prior
                     data = obs_data_m %>% filter(year_month < "2020-01-01"),
                     family = poisson(link = "log"))
  
  final_mod_4 <- glm(formula = total_deaths ~ preds + bs(chron_index, degree = 1, knots = (73 - 24)), # Last knot at 24 months prior
                     data = obs_data_m %>% filter(year_month < "2020-01-01"),
                     family = poisson(link = "log"))
  final_mod_5 <- glm(formula = total_deaths ~ factor(month) + factor(year), 
                     data =  obs_data_m %>% filter(year_month < "2020-03-01"),
                     family = poisson(link = "log"))
  
  
  pred_ly <-obs_data_m %>%  filter(year_month >= "2019-01-01" & year_month <= "2019-12-01") %>% select(total_deaths) %>% unlist %>% as.vector
  pred_1 <- predict(final_mod_1, newdata = obs_data_m %>%  filter(year_month >= "2020-01-01"))
  pred_2 <- predict(final_mod_2, newdata = obs_data_m %>%  filter(year_month >= "2020-01-01"))
  pred_3 <- predict(final_mod_3, newdata = obs_data_m %>%  filter(year_month >= "2020-01-01")) 
  pred_4 <- predict(final_mod_4, newdata = obs_data_m %>%  filter(year_month >= "2020-01-01"))
  pred_5 <- predict(final_mod_5, newdata = obs_data_m %>%  filter(year==2020))
  # year effect for 2021 is equal to 2020, so the same predictions are made for the 5th model for both years.
  
  # Predictions at the monthly level:
  preds_df <- data.frame(ly = pred_ly, knot_1 = pred_1 %>% exp, knot_2 = pred_2 %>% exp, 
                         knot_3 = pred_3 %>% exp, knot_4 = pred_4 %>% exp,
                         dummy = pred_5 %>% exp)
  
  ihme_est_m <- data.frame(year_shell, estimate=apply(preds_df, 1, function(x) sum(x*weight_vec)))
  ihme_est_y <- ihme_est_m %>% aggregate(estimate ~ year, sum)
  
  # Acosta and Irizarry method
  
  obs_data_m <- obs_data_m %>% rename(outcome=total_deaths, date=year_month)
  obs_data_m$outcome <- as.integer(obs_data_m$outcome)
  
  exclude_dates_v1 <- c(
    make_date(2015, 1, 1),
    make_date(2017, 1, 1),
    make_date(2018, 1, 1),
    seq(make_date(2020, 1, 1), make_date(2024, 12, 1), by = "month"))
  
  # Flexible model to account for tapering off trend
  
  if (2020 - start_year < 5){
    ai_ed <- compute_expected(obs_data_m, exclude = exclude_dates_v1, harmonics = 2, 
                              trend.knots.per.year = 1/4, frequency = 12, weekday.effect = FALSE, 
                              extrapolate=FALSE, verbose=FALSE, include.trend=FALSE)   
  }else{
    ai_ed <- compute_expected(obs_data_m, exclude = exclude_dates_v1, harmonics = 2, 
                              trend.knots.per.year = 1/4, frequency = 12, weekday.effect = FALSE, 
                              extrapolate=FALSE, verbose=FALSE)   
  }
  
  if (end_year == 2020){
    model <- excess_model(ai_ed, start = make_date(2014, 1, 1), end = make_date(2020, 12, 31),
                          exclude = c(make_date(2015, 1, 1),
                                      make_date(2017, 1, 1), 
                                      make_date(2018, 1, 1)), # exclude months with abnormal flu
                          knots.per.year = 3,
                          intervals = seq(make_date(2020, 1, 1), make_date(2020, 12, 31), by = "month"))}
  if (end_year > 2020){
    model <- excess_model(ai_ed, start = make_date(2014, 1, 1), end = make_date(2021, 12, 31),
                          exclude = c(make_date(2015, 1, 1),
                                      make_date(2017, 1, 1), 
                                      make_date(2018, 1, 1)), # exclude months with abnormal flu
                          knots.per.year = 3,
                          intervals = seq(make_date(2020, 1, 1), make_date(2023, 12, 31), by = "month"))
  }
  
  ai_est_m <- model$excess %>% select(start, excess) %>% mutate(year=year(start), month=month(start)) %>% filter(year > 2019) %>% 
    filter(excess != 0) %>% select(year, month, excess)
  colnames(ai_est_m) <- c("year", "month", "estimate")
  
  ai_est_y <- ai_est_m %>% aggregate(estimate ~ year, sum)
  
  ## Storing results from each method for each different scale.
  
  m_est <- cbind.data.frame(year_shell, 
                            truth = obs_death_m - cf_death_m, # true excess deaths
                            wmd_est = obs_death_m - wmd_est_m$estimate, # baseline
                            who_est = obs_death_m - who_est_m$estimate, 
                            econ_est = obs_death_m - econ_est_m$estimate,
                            ai_est = ai_est_m$estimate,
                            ihme_est = obs_death_m - ihme_est_m$estimate)
  
  y_est <- cbind.data.frame(year = 2020:end_year, 
                            truth = obs_death_y - cf_death_y, # true excess deaths
                            wmd_est = obs_death_y - wmd_est_y$estimate, # baseline
                            who_est = obs_death_y - who_est_y$estimate, 
                            econ_est = obs_death_y - econ_est_y$estimate,
                            ai_est = ai_est_y$estimate,
                            ihme_est = obs_death_y - ihme_est_y$estimate)
  return(list(m_est = m_est, y_est = y_est))
}

methods <- c("wmd_est", "who_est", "econ_est", "ai_est", "ihme_est")

#### Scenario 3: 2014 to 2021
sc3_2014 <- list()
for (n in 1:N){
  sc3_2014[[n]] <- suppressWarnings(sim_loop_sc3(n, start_year = 2014, end_year = 2023))
  print(paste0("Simulation run ", n, " completed."))
}

sc3_2014_2020_m <- cbind(get_prmse_m(sc3_2014, 2020), get_mape_m(sc3_2014, 2020), get_emp_sd_m(sc3_2014, 2020)) 
sc3_2014_2020_y <- cbind(get_prmse_y(sc3_2014, 2020), get_mape_y(sc3_2014, 2020), get_emp_sd_y(sc3_2014, 2020)) 

sc3_2014_2020_m_rel <- apply(sc3_2014_2020_m, 1, function(x) sc3_2014_2020_m[1,]/x) %>% t %>% round(2)
sc3_2014_2020_y_rel <- apply(sc3_2014_2020_y, 1, function(x) sc3_2014_2020_y[1,]/x) %>% t %>% round(2)

sc3_2014_2021_m <- cbind(get_prmse_m(sc3_2014, 2021), get_mape_m(sc3_2014, 2021), get_emp_sd_m(sc3_2014, 2021)) 
sc3_2014_2021_y <- cbind(get_prmse_y(sc3_2014, 2021), get_mape_y(sc3_2014, 2021), get_emp_sd_y(sc3_2014, 2021)) 

sc3_2014_2021_m_rel <- apply(sc3_2014_2021_m, 1, function(x) sc3_2014_2021_m[1,]/x) %>% t %>% round(2)
sc3_2014_2021_y_rel <- apply(sc3_2014_2021_y, 1, function(x) sc3_2014_2021_y[1,]/x) %>% t %>% round(2)

sc3_2014_2022_m <- cbind(get_prmse_m(sc3_2014, 2022), get_mape_m(sc3_2014, 2022), get_emp_sd_m(sc3_2014, 2022)) 
sc3_2014_2022_y <- cbind(get_prmse_y(sc3_2014, 2022), get_mape_y(sc3_2014, 2022), get_emp_sd_y(sc3_2014, 2022)) 

sc3_2014_2022_m_rel <- apply(sc3_2014_2022_m, 1, function(x) sc3_2014_2022_m[1,]/x) %>% t %>% round(2)
sc3_2014_2022_y_rel <- apply(sc3_2014_2022_y, 1, function(x) sc3_2014_2022_y[1,]/x) %>% t %>% round(2)

sc3_2014_2023_m <- cbind(get_prmse_m(sc3_2014, 2023), get_mape_m(sc3_2014, 2023), get_emp_sd_m(sc3_2014, 2023)) 
sc3_2014_2023_y <- cbind(get_prmse_y(sc3_2014, 2023), get_mape_y(sc3_2014, 2023), get_emp_sd_y(sc3_2014, 2023)) 

sc3_2014_2023_m_rel <- apply(sc3_2014_2023_m, 1, function(x) sc3_2014_2023_m[1,]/x) %>% t %>% round(2)
sc3_2014_2023_y_rel <- apply(sc3_2014_2023_y, 1, function(x) sc3_2014_2023_y[1,]/x) %>% t %>% round(2)

#### Scenario 3: 2017 to 2021 
sc3_2017 <- list()
for (n in 1:N){
  sc3_2017[[n]] <- suppressWarnings(sim_loop_sc3(n, start_year = 2017, end_year = 2023))
  print(paste0("Simulation run ", n, " completed."))
}

sc3_2017_2020_m <- cbind(get_prmse_m(sc3_2017, 2020), get_mape_m(sc3_2017, 2020), get_emp_sd_m(sc3_2017, 2020)) 
sc3_2017_2020_y <- cbind(get_prmse_y(sc3_2017, 2020), get_mape_y(sc3_2017, 2020), get_emp_sd_y(sc3_2017, 2020)) 

sc3_2017_2020_m_rel <- apply(sc3_2017_2020_m, 1, function(x) sc3_2017_2020_m[1,]/x) %>% t %>% round(2)
sc3_2017_2020_y_rel <- apply(sc3_2017_2020_y, 1, function(x) sc3_2017_2020_y[1,]/x) %>% t %>% round(2)

sc3_2017_2021_m <- cbind(get_prmse_m(sc3_2017, 2021), get_mape_m(sc3_2017, 2021), get_emp_sd_m(sc3_2017, 2021)) 
sc3_2017_2021_y <- cbind(get_prmse_y(sc3_2017, 2021), get_mape_y(sc3_2017, 2021), get_emp_sd_y(sc3_2017, 2021)) 

sc3_2017_2021_m_rel <- apply(sc3_2017_2021_m, 1, function(x) sc3_2017_2021_m[1,]/x) %>% t %>% round(2)
sc3_2017_2021_y_rel <- apply(sc3_2017_2021_y, 1, function(x) sc3_2017_2021_y[1,]/x) %>% t %>% round(2)

sc3_2017_2022_m <- cbind(get_prmse_m(sc3_2017, 2022), get_mape_m(sc3_2017, 2022), get_emp_sd_m(sc3_2017, 2022)) 
sc3_2017_2022_y <- cbind(get_prmse_y(sc3_2017, 2022), get_mape_y(sc3_2017, 2022), get_emp_sd_y(sc3_2017, 2022)) 

sc3_2017_2022_m_rel <- apply(sc3_2017_2022_m, 1, function(x) sc3_2017_2022_m[1,]/x) %>% t %>% round(2)
sc3_2017_2022_y_rel <- apply(sc3_2017_2022_y, 1, function(x) sc3_2017_2022_y[1,]/x) %>% t %>% round(2)

sc3_2017_2023_m <- cbind(get_prmse_m(sc3_2017, 2023), get_mape_m(sc3_2017, 2023), get_emp_sd_m(sc3_2017, 2023)) 
sc3_2017_2023_y <- cbind(get_prmse_y(sc3_2017, 2023), get_mape_y(sc3_2017, 2023), get_emp_sd_y(sc3_2017, 2023)) 

sc3_2017_2023_m_rel <- apply(sc3_2017_2023_m, 1, function(x) sc3_2017_2023_m[1,]/x) %>% t %>% round(2)
sc3_2017_2023_y_rel <- apply(sc3_2017_2023_y, 1, function(x) sc3_2017_2023_y[1,]/x) %>% t %>% round(2)

#### VISUALIZATION: Method performance at yearly scale ####
### Visualizing scenario 1

sc1_2014_res_y <- rbind.data.frame(cbind.data.frame(melt(sc1_2014_2020_y_rel), start="2014–2019", est=2020),
                                 cbind.data.frame(melt(sc1_2014_2021_y_rel), start="2014–2019", est=2021),
                                 cbind.data.frame(melt(sc1_2014_2022_y_rel), start="2014–2019", est=2022),
                                 cbind.data.frame(melt(sc1_2014_2023_y_rel), start="2014–2019", est=2023))

colnames(sc1_2014_res_y) <- c("Method", "Metric", "Value", "Starting year", "Year estimated")
sc1_2014_res_y$Metric <- factor(sc1_2014_res_y$Metric)
sc1_2014_res_y$Method <- factor(sc1_2014_res_y$Method)
levels(sc1_2014_res_y$Method) <- c("WMD", "WHO", "Economist", "A-I", "IHME")
levels(sc1_2014_res_y$Metric) <- c("RMSPE", "MAPE", "MCSD")


sc1_2017_res_y <- rbind.data.frame(cbind.data.frame(melt(sc1_2017_2020_y_rel), start="2017–2019", est=2020),
                                 cbind.data.frame(melt(sc1_2017_2021_y_rel), start="2017–2019", est=2021),
                                 cbind.data.frame(melt(sc1_2017_2022_y_rel), start="2017–2019", est=2022),
                                 cbind.data.frame(melt(sc1_2017_2023_y_rel), start="2017–2019", est=2023))

colnames(sc1_2017_res_y) <- c("Method", "Metric", "Value", "Starting year", "Year estimated")
sc1_2017_res_y$Metric <- factor(sc1_2017_res_y$Metric)
sc1_2017_res_y$Method <- factor(sc1_2017_res_y$Method)
levels(sc1_2017_res_y$Method) <- c("WMD", "WHO", "Economist", "A-I", "IHME")
levels(sc1_2017_res_y$Metric) <- c("RMSPE", "MAPE", "MCSD")

#### Visualizing scenario 2

sc2_2014_res_y <- rbind.data.frame(cbind.data.frame(melt(sc2_2014_2020_y_rel), start="2014–2019", est=2020),
                                   cbind.data.frame(melt(sc2_2014_2021_y_rel), start="2014–2019", est=2021),
                                   cbind.data.frame(melt(sc2_2014_2022_y_rel), start="2014–2019", est=2022),
                                   cbind.data.frame(melt(sc2_2014_2023_y_rel), start="2014–2019", est=2023))

colnames(sc2_2014_res_y) <- c("Method", "Metric", "Value", "Starting year", "Year estimated")
sc2_2014_res_y$Metric <- factor(sc2_2014_res_y$Metric)
sc2_2014_res_y$Method <- factor(sc2_2014_res_y$Method)
levels(sc2_2014_res_y$Method) <- c("WMD", "WHO", "Economist", "A-I", "IHME")
levels(sc2_2014_res_y$Metric) <- c("RMSPE", "MAPE", "MCSD")

sc2_2017_res_y <- rbind.data.frame(cbind.data.frame(melt(sc2_2017_2020_y_rel), start="2017–2019", est=2020),
                                   cbind.data.frame(melt(sc2_2017_2021_y_rel), start="2017–2019", est=2021),
                                   cbind.data.frame(melt(sc2_2017_2022_y_rel), start="2017–2019", est=2022),
                                   cbind.data.frame(melt(sc2_2017_2023_y_rel), start="2017–2019", est=2023))

colnames(sc2_2017_res_y) <- c("Method", "Metric", "Value", "Starting year", "Year estimated")
sc2_2017_res_y$Metric <- factor(sc2_2017_res_y$Metric)
sc2_2017_res_y$Method <- factor(sc2_2017_res_y$Method)
levels(sc2_2017_res_y$Method) <- c("WMD", "WHO", "Economist", "A-I", "IHME")
levels(sc2_2017_res_y$Metric) <- c("RMSPE", "MAPE", "MCSD")

#### Visualizing scenario 3

sc3_2014_res_y <- rbind.data.frame(cbind.data.frame(melt(sc3_2014_2020_y_rel), start="2014–2019", est=2020),
                                   cbind.data.frame(melt(sc3_2014_2021_y_rel), start="2014–2019", est=2021),
                                   cbind.data.frame(melt(sc3_2014_2022_y_rel), start="2014–2019", est=2022),
                                   cbind.data.frame(melt(sc3_2014_2023_y_rel), start="2014–2019", est=2023))

colnames(sc3_2014_res_y) <- c("Method", "Metric", "Value", "Starting year", "Year estimated")
sc3_2014_res_y$Metric <- factor(sc3_2014_res_y$Metric)
sc3_2014_res_y$Method <- factor(sc3_2014_res_y$Method)
levels(sc3_2014_res_y$Method) <- c("WMD", "WHO", "Economist", "A-I", "IHME")
levels(sc3_2014_res_y$Metric) <- c("RMSPE", "MAPE", "MCSD")

sc3_2017_res_y <- rbind.data.frame(cbind.data.frame(melt(sc3_2017_2020_y_rel), start="2017–2019", est=2020),
                                   cbind.data.frame(melt(sc3_2017_2021_y_rel), start="2017–2019", est=2021),
                                   cbind.data.frame(melt(sc3_2017_2022_y_rel), start="2017–2019", est=2022),
                                   cbind.data.frame(melt(sc3_2017_2023_y_rel), start="2017–2019", est=2023))

colnames(sc3_2017_res_y) <- c("Method", "Metric", "Value", "Starting year", "Year estimated")
sc3_2017_res_y$Metric <- factor(sc3_2017_res_y$Metric)
sc3_2017_res_y$Method <- factor(sc3_2017_res_y$Method)
levels(sc3_2017_res_y$Method) <- c("WMD", "WHO", "Economist", "A-I", "IHME")
levels(sc3_2017_res_y$Metric) <- c("RMSPE", "MAPE", "MCSD")

sc1_y_plot <- 
  ggplot(rbind(sc1_2014_res_y, sc1_2017_res_y), 
       aes(y=Value, x=`Year estimated`, group=Metric, color=Method))+
  geom_hline(yintercept=1, linetype="dashed", size=0.8) + geom_point(show.legend=TRUE, size=2.5) + 
  geom_line(show.legend=TRUE, aes(group=Method), size=1.5 ) + 
  facet_grid(rows=vars(Metric), cols=vars(`Starting year`), axes="margins") +
  theme(legend.text = element_text(size = 18), 
        legend.title = element_text(size = 18), 
        legend.key.size = unit(1, 'cm'), 
        legend.position="bottom",
        axis.text = element_text(size=16),  # Adjust axis text size
        axis.title = element_text(size=18), # Adjust axis labels text size
        strip.text = element_text(size=18), # Adjust facet group text size
        plot.title = element_text(size=18))  +  theme(strip.background = element_rect(fill="lightblue", size=1, color="lightblue4")) +
  labs(title="Scenario 1: Realistic United States", y="Relative performance to WMD baseline")

sc2_y_plot <- ggplot(rbind(sc2_2014_res_y, sc2_2017_res_y), 
                     aes(y=Value, x=`Year estimated`, group=Metric, color=Method))+
  geom_hline(yintercept=1, linetype="dashed", size=0.8) + geom_point(show.legend=TRUE, size=2.5) + 
  geom_line(show.legend=TRUE, aes(group=Method), size=1.5 ) + 
  facet_grid(rows=vars(Metric), cols=vars(`Starting year`), axes="margins") +
  theme(legend.text = element_text(size = 18), 
        legend.title = element_text(size = 18), 
        legend.key.size = unit(1, 'cm'), 
        legend.position="bottom",
        axis.text = element_text(size=16),  # Adjust axis text size
        axis.title = element_text(size=18), # Adjust axis labels text size
        strip.text = element_text(size=18), # Adjust facet group text size
        plot.title = element_text(size=18))  +  theme(strip.background = element_rect(fill="lightblue", size=1, color="lightblue4")) +
  labs(title="Scenario 2: Younger population, deadlier pandemic", y="")

sc3_y_plot <- ggplot(rbind(sc3_2014_res_y, sc3_2017_res_y), 
       aes(y=Value, x=`Year estimated`, group=Metric, color=Method))+
  geom_hline(yintercept=1, linetype="dashed", size=0.8) + geom_point(show.legend=TRUE, size=2.5) + 
  geom_line(show.legend=TRUE, aes(group=Method), size=1.5 ) + 
  facet_grid(rows=vars(Metric), cols=vars(`Starting year`), axes="margins") +
  theme(legend.text = element_text(size = 18), 
        legend.title = element_text(size = 18), 
        legend.key.size = unit(1, 'cm'), 
        legend.position="bottom",
        axis.text = element_text(size=16),  # Adjust axis text size
        axis.title = element_text(size=18), # Adjust axis labels text size
        strip.text = element_text(size=18), # Adjust facet group text size
        plot.title = element_text(size=18))  +  theme(strip.background = element_rect(fill="lightblue", size=1, color="lightblue4")) +
  labs(title="Scenario 3: Older population, milder pandemic", y="")

ggarrange(sc1_y_plot, sc2_y_plot, sc3_y_plot, ncol=3, common.legend=TRUE)

#### VISUALIZATION: Method performance at monthly scale ####
sc1_2014_res_m <- rbind.data.frame(cbind.data.frame(melt(sc1_2014_2020_m_rel), start="2014–2019", est=2020),
                                   cbind.data.frame(melt(sc1_2014_2021_m_rel), start="2014–2019", est=2021),
                                   cbind.data.frame(melt(sc1_2014_2022_m_rel), start="2014–2019", est=2022),
                                   cbind.data.frame(melt(sc1_2014_2023_m_rel), start="2014–2019", est=2023))

colnames(sc1_2014_res_m) <- c("Method", "Metric", "Value", "Starting year", "Year estimated")
sc1_2014_res_m$Metric <- factor(sc1_2014_res_m$Metric)
sc1_2014_res_m$Method <- factor(sc1_2014_res_m$Method)
levels(sc1_2014_res_m$Method) <- c("WMD", "WHO", "Economist", "A-I", "IHME")
levels(sc1_2014_res_m$Metric) <- c("RMSPE", "MAPE", "MCSD")


sc1_2017_res_m <- rbind.data.frame(cbind.data.frame(melt(sc1_2017_2020_m_rel), start="2017–2019", est=2020),
                                   cbind.data.frame(melt(sc1_2017_2021_m_rel), start="2017–2019", est=2021),
                                   cbind.data.frame(melt(sc1_2017_2022_m_rel), start="2017–2019", est=2022),
                                   cbind.data.frame(melt(sc1_2017_2023_m_rel), start="2017–2019", est=2023))

colnames(sc1_2017_res_m) <- c("Method", "Metric", "Value", "Starting year", "Year estimated")
sc1_2017_res_m$Metric <- factor(sc1_2017_res_m$Metric)
sc1_2017_res_m$Method <- factor(sc1_2017_res_m$Method)
levels(sc1_2017_res_m$Method) <- c("WMD", "WHO", "Economist", "A-I", "IHME")
levels(sc1_2017_res_m$Metric) <- c("RMSPE", "MAPE", "MCSD")

sc1_2014_plot_m <- ggplot(sc1_2014_res_m, 
                          aes(y=Value, x=`Year estimated`, group=Metric, color=Method)) + 
  geom_hline(yintercept=1, linetype="dashed", size=0.8) + geom_point(show.legend=TRUE, size=2.5) + 
  geom_line(show.legend=TRUE, aes(group=Method), size=1.5) + facet_wrap(~Metric, ncol=5)  + 
  labs(title="Scenario 1: Reference period 2014–2019", y="Relative performance",x="Year estimated") + 
  theme(legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        legend.key.size = unit(1, 'cm'), 
        axis.text = element_text(size=12),  # Adjust axis text size
        axis.title = element_text(size=14), # Adjust axis labels text size
        strip.text = element_text(size=12), # Adjust facet group text size
        plot.title = element_text(size=16)) 

sc1_2017_plot_m <- ggplot(sc1_2017_res_m, 
                          aes(y=Value, x=`Year estimated`, group=Metric, color=Method)) + 
  geom_hline(yintercept=1, linetype="dashed", size=0.8) + geom_point(show.legend=TRUE, size=2.5) + 
  geom_line(show.legend=TRUE, aes(group=Method), size=1.5) + facet_wrap(~Metric, ncol=5)  + 
  labs(title="Scenario 1: Reference period 2017–2019", y="Relative performance",x="Year estimated") + 
  theme(legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        legend.key.size = unit(1, 'cm'), 
        axis.text = element_text(size=12),  # Adjust axis text size
        axis.title = element_text(size=14), # Adjust axis labels text size
        strip.text = element_text(size=12), # Adjust facet group text size
        plot.title = element_text(size=16)) 


sc1_plots <- ggarrange(sc1_2014_plot_m, sc1_2017_plot_m, nrow=2, common.legend=TRUE)
sc1_plots

#### Visualizing scenario 2

sc2_2014_res_m <- rbind.data.frame(cbind.data.frame(melt(sc2_2014_2020_m_rel), start="2014–2019", est=2020),
                                   cbind.data.frame(melt(sc2_2014_2021_m_rel), start="2014–2019", est=2021),
                                   cbind.data.frame(melt(sc2_2014_2022_m_rel), start="2014–2019", est=2022),
                                   cbind.data.frame(melt(sc2_2014_2023_m_rel), start="2014–2019", est=2023))

colnames(sc2_2014_res_m) <- c("Method", "Metric", "Value", "Starting year", "Year estimated")
sc2_2014_res_m$Metric <- factor(sc2_2014_res_m$Metric)
sc2_2014_res_m$Method <- factor(sc2_2014_res_m$Method)
levels(sc2_2014_res_m$Method) <- c("WMD", "WHO", "Economist", "A-I", "IHME")
levels(sc2_2014_res_m$Metric) <- c("RMSPE", "MAPE", "MCSD")


sc2_2017_res_m <- rbind.data.frame(cbind.data.frame(melt(sc2_2017_2020_m_rel), start="2017–2019", est=2020),
                                   cbind.data.frame(melt(sc2_2017_2021_m_rel), start="2017–2019", est=2021),
                                   cbind.data.frame(melt(sc2_2017_2022_m_rel), start="2017–2019", est=2022),
                                   cbind.data.frame(melt(sc2_2017_2023_m_rel), start="2017–2019", est=2023))

colnames(sc2_2017_res_m) <- c("Method", "Metric", "Value", "Starting year", "Year estimated")
sc2_2017_res_m$Metric <- factor(sc2_2017_res_m$Metric)
sc2_2017_res_m$Method <- factor(sc2_2017_res_m$Method)
levels(sc2_2017_res_m$Method) <- c("WMD", "WHO", "Economist", "A-I", "IHME")
levels(sc2_2017_res_m$Metric) <- c("RMSPE", "MAPE", "MCSD")

sc2_2014_plot_m <- ggplot(sc2_2014_res_m, 
                          aes(y=Value, x=`Year estimated`, group=Metric, color=Method)) + 
  geom_hline(yintercept=1, linetype="dashed", size=0.8) + geom_point(show.legend=TRUE, size=2.5) + 
  geom_line(show.legend=TRUE, aes(group=Method), size=1.5) + facet_wrap(~Metric, ncol=5)  + 
  labs(title="Scenario 2: Reference period 2014–2019", y="Relative performance",x="Year estimated") + 
  theme(legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        legend.key.size = unit(1, 'cm'), 
        axis.text = element_text(size=12),  # Adjust axis text size
        axis.title = element_text(size=14), # Adjust axis labels text size
        strip.text = element_text(size=12), # Adjust facet group text size
        plot.title = element_text(size=16)) 

sc2_2017_plot_m <- ggplot(sc2_2017_res_m %>% filter(Method %in% c("WHO", "Economist", "A-I", "IHME")), 
                          aes(y=Value, x=`Year estimated`, group=Metric, color=Method)) + 
  geom_hline(yintercept=1, linetype="dashed", size=0.8) + geom_point(show.legend=TRUE, size=2.5) + 
  geom_line(show.legend=TRUE, aes(group=Method), size=1.5) + facet_wrap(~Metric, ncol=5)  + 
  labs(title="Scenario 2: Reference period 2017–2019", y="Relative performance",x="Year estimated") + 
  theme(legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        legend.key.size = unit(1, 'cm'), 
        axis.text = element_text(size=12),  # Adjust axis text size
        axis.title = element_text(size=14), # Adjust axis labels text size
        strip.text = element_text(size=12), # Adjust facet group text size
        plot.title = element_text(size=16)) 


sc2_plots <- ggarrange(sc2_2014_plot_m, sc2_2017_plot_m, nrow=2, common.legend=TRUE)
sc2_plots

#### Visualizing scenario 3

sc3_2014_res_m <- rbind.data.frame(cbind.data.frame(melt(sc3_2014_2020_m_rel), start="2014–2019", est=2020),
                                   cbind.data.frame(melt(sc3_2014_2021_m_rel), start="2014–2019", est=2021),
                                   cbind.data.frame(melt(sc3_2014_2022_m_rel), start="2014–2019", est=2022),
                                   cbind.data.frame(melt(sc3_2014_2023_m_rel), start="2014–2019", est=2023))

colnames(sc3_2014_res_m) <- c("Method", "Metric", "Value", "Starting year", "Year estimated")
sc3_2014_res_m$Metric <- factor(sc3_2014_res_m$Metric)
sc3_2014_res_m$Method <- factor(sc3_2014_res_m$Method)
levels(sc3_2014_res_m$Method) <- c("WMD", "WHO", "Economist", "A-I", "IHME")
levels(sc3_2014_res_m$Metric) <- c("RMSPE", "MAPE", "MCSD")

sc3_2017_res_m <- rbind.data.frame(cbind.data.frame(melt(sc3_2017_2020_m_rel), start="2017–2019", est=2020),
                                   cbind.data.frame(melt(sc3_2017_2021_m_rel), start="2017–2019", est=2021),
                                   cbind.data.frame(melt(sc3_2017_2022_m_rel), start="2017–2019", est=2022),
                                   cbind.data.frame(melt(sc3_2017_2023_m_rel), start="2017–2019", est=2023))

colnames(sc3_2017_res_m) <- c("Method", "Metric", "Value", "Starting year", "Year estimated")
sc3_2017_res_m$Metric <- factor(sc3_2017_res_m$Metric)
sc3_2017_res_m$Method <- factor(sc3_2017_res_m$Method)
levels(sc3_2017_res_m$Method) <- c("WMD", "WHO", "Economist", "A-I", "IHME")
levels(sc3_2017_res_m$Metric) <- c("RMSPE", "MAPE", "MCSD")

sc3_2014_plot_m <- ggplot(sc3_2014_res_m, 
                          aes(y=Value, x=`Year estimated`, group=Metric, color=Method)) + 
  geom_hline(yintercept=1, linetype="dashed", size=0.8) + geom_point(show.legend=TRUE, size=2.5) + 
  geom_line(show.legend=TRUE, aes(group=Method), size=1.5) + facet_wrap(~Metric, ncol=5)  + 
  labs(title="Scenario 3: Reference period 2014–2019", y="Relative performance",x="Year estimated") + 
  theme(legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        legend.key.size = unit(1, 'cm'), 
        axis.text = element_text(size=12),  # Adjust axis text size
        axis.title = element_text(size=14), # Adjust axis labels text size
        strip.text = element_text(size=12), # Adjust facet group text size
        plot.title = element_text(size=16)) 

sc3_2017_plot_m <- ggplot(sc3_2017_res_m, 
                          aes(y=Value, x=`Year estimated`, group=Metric, color=Method)) + 
  geom_hline(yintercept=1, linetype="dashed", size=0.8) + geom_point(show.legend=TRUE, size=2.5) + 
  geom_line(show.legend=TRUE, aes(group=Method), size=1.5) + facet_wrap(~Metric, ncol=5)  + 
  labs(title="Scenario 3: Reference period 2017–2019", y="Relative performance",x="Year estimated") + 
  theme(legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        legend.key.size = unit(1, 'cm'), 
        axis.text = element_text(size=12),  # Adjust axis text size
        axis.title = element_text(size=14), # Adjust axis labels text size
        strip.text = element_text(size=12), # Adjust facet group text size
        plot.title = element_text(size=16)) 


sc1_m_plot <- ggplot(rbind(sc1_2014_res_m, sc1_2017_res_m), 
                     aes(y=Value, x=`Year estimated`, group=Metric, color=Method))+
  geom_hline(yintercept=1, linetype="dashed", size=0.8) + geom_point(show.legend=TRUE, size=2.5) + 
  geom_line(show.legend=TRUE, aes(group=Method), size=1.5 ) + 
  facet_grid(rows=vars(Metric), cols=vars(`Starting year`), axes="margins") +
  theme(legend.text = element_text(size = 18), 
        legend.title = element_text(size = 18), 
        legend.key.size = unit(1, 'cm'), 
        legend.position="bottom",
        axis.text = element_text(size=16),  # Adjust axis text size
        axis.title = element_text(size=18), # Adjust axis labels text size
        strip.text = element_text(size=18), # Adjust facet group text size
        plot.title = element_text(size=18))  +  theme(strip.background = element_rect(fill="lightblue", size=1, color="lightblue4")) +
  labs(title="Scenario 1: Realistic United States", y="Relative performance to WMD baseline")

sc2_m_plot <- ggplot(rbind(sc2_2014_res_m, sc2_2017_res_m), 
                     aes(y=Value, x=`Year estimated`, group=Metric, color=Method))+
  geom_hline(yintercept=1, linetype="dashed", size=0.8) + geom_point(show.legend=TRUE, size=2.5) + 
  geom_line(show.legend=TRUE, aes(group=Method), size=1.5 ) + 
  facet_grid(rows=vars(Metric), cols=vars(`Starting year`), axes="margins") +
  theme(legend.text = element_text(size = 18), 
        legend.title = element_text(size = 18), 
        legend.key.size = unit(1, 'cm'), 
        legend.position="bottom",
        axis.text = element_text(size=16),  # Adjust axis text size
        axis.title = element_text(size=18), # Adjust axis labels text size
        strip.text = element_text(size=18), # Adjust facet group text size
        plot.title = element_text(size=18))  +  theme(strip.background = element_rect(fill="lightblue", size=1, color="lightblue4")) +
  labs(title="Scenario 2: Younger population, deadlier pandemic", y="")

sc3_m_plot <- ggplot(rbind(sc3_2014_res_m, sc3_2017_res_m), 
                     aes(y=Value, x=`Year estimated`, group=Metric, color=Method))+
  geom_hline(yintercept=1, linetype="dashed", size=0.8) + geom_point(show.legend=TRUE, size=2.5) + 
  geom_line(show.legend=TRUE, aes(group=Method), size=1.5 ) + 
  facet_grid(rows=vars(Metric), cols=vars(`Starting year`), axes="margins") +
  theme(legend.text = element_text(size = 18), 
        legend.title = element_text(size = 18), 
        legend.key.size = unit(1, 'cm'), 
        legend.position="bottom",
        axis.text = element_text(size=16),  # Adjust axis text size
        axis.title = element_text(size=18), # Adjust axis labels text size
        strip.text = element_text(size=18), # Adjust facet group text size
        plot.title = element_text(size=18))  +  theme(strip.background = element_rect(fill="lightblue", size=1, color="lightblue4")) +
  labs(title="Scenario 3: Older population, milder pandemic", y="")

ggarrange(sc1_m_plot, sc2_m_plot, sc3_m_plot, ncol=3, common.legend=TRUE)


#### Data table: 2014 reference period ####

### Scenario 1 table

sc1_2014_2020_df <- rbind(2020, format_data(sc1_2014_2020_y, sc1_2014_2020_y_rel))
colnames(sc1_2014_2020_df) <- c("RMSPE", "MAPE", "MCSD")
rownames(sc1_2014_2020_df) <- c("Year", "WMD", "WHO", "Economist", "A-I", "IHME")

sc1_2014_2021_df <- rbind(2021, format_data(sc1_2014_2021_y, sc1_2014_2021_y_rel))
colnames(sc1_2014_2021_df) <- c("RMSPE", "MAPE", "MCSD")
rownames(sc1_2014_2021_df) <- c("Year", "WMD", "WHO", "Economist", "A-I", "IHME")

sc1_2014_2022_df <- rbind(2022, format_data(sc1_2014_2022_y, sc1_2014_2022_y_rel))
colnames(sc1_2014_2022_df) <- c("RMSPE", "MAPE", "MCSD")
rownames(sc1_2014_2022_df) <- c("Year", "WMD", "WHO", "Economist", "A-I", "IHME")

sc1_2014_2023_df <- rbind(2023, format_data(sc1_2014_2023_y, sc1_2014_2023_y_rel))
colnames(sc1_2014_2023_df) <- c("RMSPE", "MAPE", "MCSD")
rownames(sc1_2014_2023_df) <- c("Year", "WMD", "WHO", "Economist", "A-I", "IHME")

sc1_y_df <- rbind(sc1_2014_2020_df, sc1_2014_2021_df, sc1_2014_2022_df, sc1_2014_2023_df)

xtable(sc1_y_df)

### Scenario 2 table

sc2_2014_2020_df <- rbind(2020, format_data(sc2_2014_2020_y, sc2_2014_2020_y_rel))
colnames(sc2_2014_2020_df) <- c("RMSPE", "MAPE", "MCSD")
rownames(sc2_2014_2020_df) <- c("Year", "WMD", "WHO", "Economist", "A-I", "IHME")

sc2_2014_2021_df <- rbind(2021, format_data(sc2_2014_2021_y, sc2_2014_2021_y_rel))
colnames(sc2_2014_2021_df) <- c("RMSPE", "MAPE", "MCSD")
rownames(sc2_2014_2021_df) <- c("Year", "WMD", "WHO", "Economist", "A-I", "IHME")

sc2_2014_2022_df <- rbind(2022, format_data(sc2_2014_2022_y, sc2_2014_2022_y_rel))
colnames(sc2_2014_2022_df) <- c("RMSPE", "MAPE", "MCSD")
rownames(sc2_2014_2022_df) <- c("Year", "WMD", "WHO", "Economist", "A-I", "IHME")

sc2_2014_2023_df <- rbind(2023, format_data(sc2_2014_2023_y, sc2_2014_2023_y_rel))
colnames(sc2_2014_2023_df) <- c("RMSPE", "MAPE", "MCSD")
rownames(sc2_2014_2023_df) <- c("Year", "WMD", "WHO", "Economist", "A-I", "IHME")

sc2_y_df <- rbind(sc2_2014_2020_df, sc2_2014_2021_df, sc2_2014_2022_df, sc2_2014_2023_df)

xtable(sc2_y_df)

### Scenario 3 table

sc3_2014_2020_df <- rbind(2020, format_data(sc3_2014_2020_y, sc3_2014_2020_y_rel))
colnames(sc3_2014_2020_df) <- c("RMSPE", "MAPE", "MCSD")
rownames(sc3_2014_2020_df) <- c("Year", "WMD", "WHO", "Economist", "A-I", "IHME")

sc3_2014_2021_df <- rbind(2021, format_data(sc3_2014_2021_y, sc3_2014_2021_y_rel))
colnames(sc3_2014_2021_df) <- c("RMSPE", "MAPE", "MCSD")
rownames(sc3_2014_2021_df) <- c("Year", "WMD", "WHO", "Economist", "A-I", "IHME")

sc3_2014_2022_df <- rbind(2022, format_data(sc3_2014_2022_y, sc3_2014_2022_y_rel))
colnames(sc3_2014_2022_df) <- c("RMSPE", "MAPE", "MCSD")
rownames(sc3_2014_2022_df) <- c("Year", "WMD", "WHO", "Economist", "A-I", "IHME")

sc3_2014_2023_df <- rbind(2023, format_data(sc3_2014_2023_y, sc3_2014_2023_y_rel))
colnames(sc3_2014_2023_df) <- c("RMSPE", "MAPE", "MCSD")
rownames(sc3_2014_2023_df) <- c("Year", "WMD", "WHO", "Economist", "A-I", "IHME")

sc3_y_df <- rbind(sc3_2014_2020_df, sc3_2014_2021_df, sc3_2014_2022_df, sc3_2014_2023_df)

xtable(sc3_y_df)

### Monthly scenario 1

sc1_2014_2020_df <- rbind(2020, format_data(sc1_2014_2020_m, sc1_2014_2020_m_rel))
colnames(sc1_2014_2020_df) <- c("RMSPE", "MAPE", "MCSD")
rownames(sc1_2014_2020_df) <- c("Year", "WMD", "WHO", "Economist", "A-I", "IHME")

sc1_2014_2021_df <- rbind(2021, format_data(sc1_2014_2021_m, sc1_2014_2021_m_rel))
colnames(sc1_2014_2021_df) <- c("RMSPE", "MAPE", "MCSD")
rownames(sc1_2014_2021_df) <- c("Year", "WMD", "WHO", "Economist", "A-I", "IHME")

sc1_2014_2022_df <- rbind(2022, format_data(sc1_2014_2022_m, sc1_2014_2022_m_rel))
colnames(sc1_2014_2022_df) <- c("RMSPE", "MAPE", "MCSD")
rownames(sc1_2014_2022_df) <- c("Year", "WMD", "WHO", "Economist", "A-I", "IHME")

sc1_2014_2023_df <- rbind(2023, format_data(sc1_2014_2023_m, sc1_2014_2023_m_rel))
colnames(sc1_2014_2023_df) <- c("RMSPE", "MAPE", "MCSD")
rownames(sc1_2014_2023_df) <- c("Year", "WMD", "WHO", "Economist", "A-I", "IHME")

sc1_m_df <- rbind(sc1_2014_2020_df, sc1_2014_2021_df, sc1_2014_2022_df, sc1_2014_2023_df)

xtable(sc1_m_df)

### Monthly scenario 2

sc2_2014_2020_df <- rbind(2020, format_data(sc2_2014_2020_m, sc2_2014_2020_m_rel))
colnames(sc2_2014_2020_df) <- c("RMSPE", "MAPE", "MCSD")
rownames(sc2_2014_2020_df) <- c("Year", "WMD", "WHO", "Economist", "A-I", "IHME")

sc2_2014_2021_df <- rbind(2021, format_data(sc2_2014_2021_m, sc2_2014_2021_m_rel))
colnames(sc2_2014_2021_df) <- c("RMSPE", "MAPE", "MCSD")
rownames(sc2_2014_2021_df) <- c("Year", "WMD", "WHO", "Economist", "A-I", "IHME")

sc2_2014_2022_df <- rbind(2022, format_data(sc2_2014_2022_m, sc2_2014_2022_m_rel))
colnames(sc2_2014_2022_df) <- c("RMSPE", "MAPE", "MCSD")
rownames(sc2_2014_2022_df) <- c("Year", "WMD", "WHO", "Economist", "A-I", "IHME")

sc2_2014_2023_df <- rbind(2023, format_data(sc2_2014_2023_m, sc2_2014_2023_m_rel))
colnames(sc2_2014_2023_df) <- c("RMSPE", "MAPE", "MCSD")
rownames(sc2_2014_2023_df) <- c("Year", "WMD", "WHO", "Economist", "A-I", "IHME")

sc2_m_df <- rbind(sc2_2014_2020_df, sc2_2014_2021_df, sc2_2014_2022_df, sc2_2014_2023_df)

xtable(sc2_m_df)

### Monthly scenario 3

sc3_2014_2020_df <- rbind(2020, format_data(sc3_2014_2020_m, sc3_2014_2020_m_rel))
colnames(sc3_2014_2020_df) <- c("RMSPE", "MAPE", "MCSD")
rownames(sc3_2014_2020_df) <- c("Year", "WMD", "WHO", "Economist", "A-I", "IHME")

sc3_2014_2021_df <- rbind(2021, format_data(sc3_2014_2021_m, sc3_2014_2021_m_rel))
colnames(sc3_2014_2021_df) <- c("RMSPE", "MAPE", "MCSD")
rownames(sc3_2014_2021_df) <- c("Year", "WMD", "WHO", "Economist", "A-I", "IHME")

sc3_2014_2022_df <- rbind(2022, format_data(sc3_2014_2022_m, sc3_2014_2022_m_rel))
colnames(sc3_2014_2022_df) <- c("RMSPE", "MAPE", "MCSD")
rownames(sc3_2014_2022_df) <- c("Year", "WMD", "WHO", "Economist", "A-I", "IHME")

sc3_2014_2023_df <- rbind(2023, format_data(sc3_2014_2023_m, sc3_2014_2023_m_rel))
colnames(sc3_2014_2023_df) <- c("RMSPE", "MAPE", "MCSD")
rownames(sc3_2014_2023_df) <- c("Year", "WMD", "WHO", "Economist", "A-I", "IHME")

sc3_m_df <- rbind(sc3_2014_2020_df, sc3_2014_2021_df, sc3_2014_2022_df, sc3_2014_2023_df)

xtable(sc3_m_df)

#### Data table: 2017 reference period ####

### Scenario 1 table

sc1_2017_2020_df <- rbind(2020, format_data(sc1_2017_2020_y, sc1_2017_2020_y_rel))
colnames(sc1_2017_2020_df) <- c("RMSPE", "MAPE", "MCSD")
rownames(sc1_2017_2020_df) <- c("Year", "WMD", "WHO", "Economist", "A-I", "IHME")

sc1_2017_2021_df <- rbind(2021, format_data(sc1_2017_2021_y, sc1_2017_2021_y_rel))
colnames(sc1_2017_2021_df) <- c("RMSPE", "MAPE", "MCSD")
rownames(sc1_2017_2021_df) <- c("Year", "WMD", "WHO", "Economist", "A-I", "IHME")

sc1_2017_2022_df <- rbind(2022, format_data(sc1_2017_2022_y, sc1_2017_2022_y_rel))
colnames(sc1_2017_2022_df) <- c("RMSPE", "MAPE", "MCSD")
rownames(sc1_2017_2022_df) <- c("Year", "WMD", "WHO", "Economist", "A-I", "IHME")

sc1_2017_2023_df <- rbind(2023, format_data(sc1_2017_2023_y, sc1_2017_2023_y_rel))
colnames(sc1_2017_2023_df) <- c("RMSPE", "MAPE", "MCSD")
rownames(sc1_2017_2023_df) <- c("Year", "WMD", "WHO", "Economist", "A-I", "IHME")

sc1_y_df <- rbind(sc1_2017_2020_df, sc1_2017_2021_df, sc1_2017_2022_df, sc1_2017_2023_df)

xtable(sc1_y_df)

### Scenario 2 table

sc2_2017_2020_df <- rbind(2020, format_data(sc2_2017_2020_y, sc2_2017_2020_y_rel))
colnames(sc2_2017_2020_df) <- c("RMSPE", "MAPE", "MCSD")
rownames(sc2_2017_2020_df) <- c("Year", "WMD", "WHO", "Economist", "A-I", "IHME")

sc2_2017_2021_df <- rbind(2021, format_data(sc2_2017_2021_y, sc2_2017_2021_y_rel))
colnames(sc2_2017_2021_df) <- c("RMSPE", "MAPE", "MCSD")
rownames(sc2_2017_2021_df) <- c("Year", "WMD", "WHO", "Economist", "A-I", "IHME")

sc2_2017_2022_df <- rbind(2022, format_data(sc2_2017_2022_y, sc2_2017_2022_y_rel))
colnames(sc2_2017_2022_df) <- c("RMSPE", "MAPE", "MCSD")
rownames(sc2_2017_2022_df) <- c("Year", "WMD", "WHO", "Economist", "A-I", "IHME")

sc2_2017_2023_df <- rbind(2023, format_data(sc2_2017_2023_y, sc2_2017_2023_y_rel))
colnames(sc2_2017_2023_df) <- c("RMSPE", "MAPE", "MCSD")
rownames(sc2_2017_2023_df) <- c("Year", "WMD", "WHO", "Economist", "A-I", "IHME")

sc2_y_df <- rbind(sc2_2017_2020_df, sc2_2017_2021_df, sc2_2017_2022_df, sc2_2017_2023_df)

xtable(sc2_y_df)

### Scenario 3 table

sc3_2017_2020_df <- rbind(2020, format_data(sc3_2017_2020_y, sc3_2017_2020_y_rel))
colnames(sc3_2017_2020_df) <- c("RMSPE", "MAPE", "MCSD")
rownames(sc3_2017_2020_df) <- c("Year", "WMD", "WHO", "Economist", "A-I", "IHME")

sc3_2017_2021_df <- rbind(2021, format_data(sc3_2017_2021_y, sc3_2017_2021_y_rel))
colnames(sc3_2017_2021_df) <- c("RMSPE", "MAPE", "MCSD")
rownames(sc3_2017_2021_df) <- c("Year", "WMD", "WHO", "Economist", "A-I", "IHME")

sc3_2017_2022_df <- rbind(2022, format_data(sc3_2017_2022_y, sc3_2017_2022_y_rel))
colnames(sc3_2017_2022_df) <- c("RMSPE", "MAPE", "MCSD")
rownames(sc3_2017_2022_df) <- c("Year", "WMD", "WHO", "Economist", "A-I", "IHME")

sc3_2017_2023_df <- rbind(2023, format_data(sc3_2017_2023_y, sc3_2017_2023_y_rel))
colnames(sc3_2017_2023_df) <- c("RMSPE", "MAPE", "MCSD")
rownames(sc3_2017_2023_df) <- c("Year", "WMD", "WHO", "Economist", "A-I", "IHME")

sc3_y_df <- rbind(sc3_2017_2020_df, sc3_2017_2021_df, sc3_2017_2022_df, sc3_2017_2023_df)

xtable(sc3_y_df)

### Monthly scenario 1

sc1_2017_2020_df <- rbind(2020, format_data(sc1_2017_2020_m, sc1_2017_2020_m_rel))
colnames(sc1_2017_2020_df) <- c("RMSPE", "MAPE", "MCSD")
rownames(sc1_2017_2020_df) <- c("Year", "WMD", "WHO", "Economist", "A-I", "IHME")

sc1_2017_2021_df <- rbind(2021, format_data(sc1_2017_2021_m, sc1_2017_2021_m_rel))
colnames(sc1_2017_2021_df) <- c("RMSPE", "MAPE", "MCSD")
rownames(sc1_2017_2021_df) <- c("Year", "WMD", "WHO", "Economist", "A-I", "IHME")

sc1_2017_2022_df <- rbind(2022, format_data(sc1_2017_2022_m, sc1_2017_2022_m_rel))
colnames(sc1_2017_2022_df) <- c("RMSPE", "MAPE", "MCSD")
rownames(sc1_2017_2022_df) <- c("Year", "WMD", "WHO", "Economist", "A-I", "IHME")

sc1_2017_2023_df <- rbind(2023, format_data(sc1_2017_2023_m, sc1_2017_2023_m_rel))
colnames(sc1_2017_2023_df) <- c("RMSPE", "MAPE", "MCSD")
rownames(sc1_2017_2023_df) <- c("Year", "WMD", "WHO", "Economist", "A-I", "IHME")

sc1_m_df <- rbind(sc1_2017_2020_df, sc1_2017_2021_df, sc1_2017_2022_df, sc1_2017_2023_df)

xtable(sc1_m_df)

### Monthly scenario 2

sc2_2017_2020_df <- rbind(2020, format_data(sc2_2017_2020_m, sc2_2017_2020_m_rel))
colnames(sc2_2017_2020_df) <- c("RMSPE", "MAPE", "MCSD")
rownames(sc2_2017_2020_df) <- c("Year", "WMD", "WHO", "Economist", "A-I", "IHME")

sc2_2017_2021_df <- rbind(2021, format_data(sc2_2017_2021_m, sc2_2017_2021_m_rel))
colnames(sc2_2017_2021_df) <- c("RMSPE", "MAPE", "MCSD")
rownames(sc2_2017_2021_df) <- c("Year", "WMD", "WHO", "Economist", "A-I", "IHME")

sc2_2017_2022_df <- rbind(2022, format_data(sc2_2017_2022_m, sc2_2017_2022_m_rel))
colnames(sc2_2017_2022_df) <- c("RMSPE", "MAPE", "MCSD")
rownames(sc2_2017_2022_df) <- c("Year", "WMD", "WHO", "Economist", "A-I", "IHME")

sc2_2017_2023_df <- rbind(2023, format_data(sc2_2017_2023_m, sc2_2017_2023_m_rel))
colnames(sc2_2017_2023_df) <- c("RMSPE", "MAPE", "MCSD")
rownames(sc2_2017_2023_df) <- c("Year", "WMD", "WHO", "Economist", "A-I", "IHME")

sc2_m_df <- rbind(sc2_2017_2020_df, sc2_2017_2021_df, sc2_2017_2022_df, sc2_2017_2023_df)

xtable(sc2_m_df)

### Monthly scenario 3

sc3_2017_2020_df <- rbind(2020, format_data(sc3_2017_2020_m, sc3_2017_2020_m_rel))
colnames(sc3_2017_2020_df) <- c("RMSPE", "MAPE", "MCSD")
rownames(sc3_2017_2020_df) <- c("Year", "WMD", "WHO", "Economist", "A-I", "IHME")

sc3_2017_2021_df <- rbind(2021, format_data(sc3_2017_2021_m, sc3_2017_2021_m_rel))
colnames(sc3_2017_2021_df) <- c("RMSPE", "MAPE", "MCSD")
rownames(sc3_2017_2021_df) <- c("Year", "WMD", "WHO", "Economist", "A-I", "IHME")

sc3_2017_2022_df <- rbind(2022, format_data(sc3_2017_2022_m, sc3_2017_2022_m_rel))
colnames(sc3_2017_2022_df) <- c("RMSPE", "MAPE", "MCSD")
rownames(sc3_2017_2022_df) <- c("Year", "WMD", "WHO", "Economist", "A-I", "IHME")

sc3_2017_2023_df <- rbind(2023, format_data(sc3_2017_2023_m, sc3_2017_2023_m_rel))
colnames(sc3_2017_2023_df) <- c("RMSPE", "MAPE", "MCSD")
rownames(sc3_2017_2023_df) <- c("Year", "WMD", "WHO", "Economist", "A-I", "IHME")

sc3_m_df <- rbind(sc3_2017_2020_df, sc3_2017_2021_df, sc3_2017_2022_df, sc3_2017_2023_df)

xtable(sc3_m_df)
