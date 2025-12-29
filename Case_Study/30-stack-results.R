library(tidyverse)
library(lubridate)
library(ggplot2)
library(ggtext)
setwd("[insert file path to working directory]")
datadir <- "./data/"
outdir <- "./output/"

dat_monthly <- read.csv(paste0(datadir,"monthly data.csv"), sep=',', header=TRUE)%>%
  rename(observed = outcome)

list_mo_methods <- c("who","economist","ihme","harvard")
list_bench_methods <- c("wmd","cdc")
list_all_methods <- c(list_mo_methods,list_bench_methods)

## Read in and append monthly results from each method
res_monthly <- data.frame()
read_fn_mo <- function(i) {
  res_monthly_i <- read.csv(paste0(outdir,i,"_monthly.csv"), sep=',', header=TRUE) %>%
    mutate(method = i)
  res_monthly <- bind_rows(res_monthly,res_monthly_i)
  
  return(res_monthly)
} 

for(m in list_mo_methods){
  res_monthly <- read_fn_mo(m)
}
# Merge on observed deaths for plots
res_m_wobs <- merge(dat_monthly,res_monthly,by=c("year","month"),all.x=T)
res_m_wobs$date <- make_date(year = res_m_wobs$year, month = res_m_wobs$month)

## Prepare data for plots
dat_plot_m <- res_m_wobs %>%
  arrange(method, date) %>%
  group_by(method) %>%
  mutate(
    excess           = zoo::rollmean(excess, k = 5, fill = NA, align = "center"),
    excess_lower95   = zoo::rollmean(excess_lower95, k = 5, fill = NA, align = "center"),
    excess_upper95   = zoo::rollmean(excess_upper95, k = 5, fill = NA, align = "center"),
    
    expected_deaths  = zoo::rollmean(expected_deaths, k = 5, fill = NA, align = "center"),
    expected_lower95 = zoo::rollmean(expected_lower95, k = 5, fill = NA, align = "center"),
    expected_upper95 = zoo::rollmean(expected_upper95, k = 5, fill = NA, align = "center")
  ) %>%
  fill(
    excess, excess_lower95, excess_upper95,
    expected_deaths, expected_lower95, expected_upper95,
    .direction = "downup"
  ) %>%
  ungroup()%>%
  mutate(Study = case_when(method == "who" ~ "World Health Organization",
                           method == "economist" ~ "The Economist",
                           method == "harvard" ~ "Acosta and Irizarry",
                           method == "harvard sensitivity" ~ "Acosta and Irizarry Sensitivity Results",
                           method == "ihme" ~ "Institute for Health Metrics and Evaluation",
                           method == "wmd" ~ "World Mortality Dataset",
                           method == "cdc" ~ "Centers for Disease Control and Prevention",
                           TRUE ~ "Missing"),
         Observed = "Observed",
         method = ifelse(is.na(method),"observed",method),
         Study = ifelse(method == "observed","Missing",Study),
         # recode bounds of uncertainty interval to be the point estimate, for plotting methods where there is no uncertainty interval
         expected_lower95 = ifelse(method %in% c("wmd", "cdc"), expected_deaths, expected_lower95),
         expected_upper95 = ifelse(method %in% c("wmd", "cdc"), expected_deaths, expected_upper95),
         excess_lower95 = ifelse(method %in% c("wmd", "cdc"), excess, excess_lower95),
         excess_upper95 = ifelse(method %in% c("wmd", "cdc"), excess, excess_upper95),
         excess_CI_width = excess_upper95 - excess_lower95,
         expected_CI_width = expected_upper95 - expected_lower95,
         across(starts_with("excess"), function(x) round(x, 0)),
         across(starts_with("expected"), function(x) round(x, 0)),
         deaths_per_100k = (observed / population) * 100000
  ) %>%
  arrange(year,month,Study)

## Prepare plot parameters
theme_expected <- theme(
  axis.text.x=element_text(angle = 45, vjust = 1, hjust=1),
  panel.background = element_rect(fill = "white"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major.y = element_line(color = "lightgrey"),
  panel.grid.major.x = element_line(color = "lightgrey"),
  legend.position = c(0.30, 0.92),
  legend.justification = c("right", "top"),
  legend.box.just = "right",
  legend.margin = margin(0.5,0.5,0.5,0.5),
  legend.box.background = element_rect(color="lightgrey",linewidth=1.2),
  legend.title=element_blank(),
  legend.text=element_markdown()
) 
theme_excess <-   theme(
  axis.text.x=element_text(angle = 45, vjust = 1, hjust=1),
  panel.background = element_rect(fill = "white"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major.y = element_line(color = "lightgrey"),
  panel.grid.major.x = element_line(color = "lightgrey"),
  legend.position = c(0.30, 0.97),
  legend.justification = c("right", "top"),
  legend.box.just = "right",
  legend.margin = margin(0.5,0.5,0.5,0.5),
  legend.box.background = element_rect(color="lightgrey",linewidth=1.2),
  legend.title=element_blank(),
  legend.text=element_markdown()
)

list_study = c(
  "World Health Organization",  
  "The Economist",
  "Acosta and Irizarry",
  "Institute for Health Metrics and Evaluation",
  "Centers for Disease Control and Prevention"
)

dat_plot_m$Study <- factor(dat_plot_m$Study, levels = list_study)

list_observed = c("Observed")
dat_plot_m$Observed <- factor(dat_plot_m$Observed, levels = list_observed)

okabe_ito_palette <- c(
  "#009E73",
  "#0072B2", 
  "#D55E00", 
  "#CC79A7",
  "#E69F00", 
  "#56B4E9", 
  "#F0E442", 
  "#0072B2", 
  "#D55E00", 
  "#000000" 
)


color_list <- c("World Health Organization"="#009E73",
                "The Economist"="#0072B2",
                "Acosta and Irizarry"="#D55E00",
                "Institute for Health Metrics and Evaluation"="#CC79A7",
                "Centers for Disease Control and Prevention"="grey60",
                "Observed"="black")

## Monthly plots for IHME, WHO, the Economist, and A-I
excess_p <- ggplot() +
  geom_ribbon(
    data = filter(dat_plot_m, dat_plot_m$year >= 2022 & !(dat_plot_m$method %in% c("wmd", "cdc","observed"))),
    aes(x = date, ymin = excess_lower95, ymax = excess_upper95, fill = Study),
    alpha = 0.3, show.legend = FALSE
  ) +
  geom_line(
    data = filter(dat_plot_m, dat_plot_m$year >= 2022),
    aes(x = date, y = excess, colour = Study),
    linewidth = 1.3
  ) +
  scale_x_date(breaks = scales::breaks_pretty(n = 24), date_labels = "%b %Y", expand = c(0, 0), limits=c(as.Date("2022-01-01"),as.Date("2024-12-01"))) +
  scale_y_continuous(breaks=scales::breaks_pretty(n = 10), limits=c(-15000, 100000)) +
  labs(title = "3B. Monthly Excess Deaths, January 2022-December 2024") +
  theme_excess +
  labs(x = "",
       y = "Excess Deaths",
       color="") +
  scale_color_manual(values = color_list) +
  scale_fill_manual(values = color_list) +
  geom_line(data=filter(dat_plot_m, dat_plot_m$year>=2022), 
            aes(x=date, y=0))


expected_p <- ggplot() +
  geom_vline(xintercept = as.numeric(as.Date("2020-01-01")), 
             linetype = "twodash", color = "#292929", size=0.7)+
  geom_vline(xintercept = as.numeric(as.Date("2022-01-01")), 
             linetype = "twodash", color = "#292929", size=0.7)+
  geom_vline(xintercept = as.numeric(as.Date("2024-12-01")), 
             linetype = "twodash", color = "#292929", size=0.7)+
  geom_ribbon(
    data = filter(dat_plot_m, !(dat_plot_m$method %in% c("cdc", "wmd","missing"))),
    aes(x = date, ymin = expected_lower95, ymax = expected_upper95, fill = Study),
    alpha = 0.3, show.legend = FALSE
  ) +
  geom_line(
    data = filter(dat_plot_m, !(dat_plot_m$method %in% c("cdc", "wmd","observed"))),
    aes(x = date, y = expected_deaths, colour = Study),
    linewidth = 1.1
  ) +
  scale_x_date(breaks = scales::breaks_pretty(n=20), date_labels = "%b %Y", expand = c(0,20), limits=c(as.Date("2014-01-01"),as.Date("2024-12-01"))) +
  scale_y_continuous(breaks=scales::breaks_pretty(n = 5),limits=c(150000, 430000)) +
  labs(title = "3A. Monthly Observed and Expected Deaths, January 2014-December 2024") +
  theme_expected +
  labs(x = "",y = "Deaths") +
  geom_point(data=dat_plot_m,
             aes(x=date, y=observed,color=Observed)) + 
  scale_color_manual(values=color_list) + 
  scale_fill_manual(values=color_list) +
  geom_text(aes(x=as.Date("2021-01-01"), 
                label=stringr::str_wrap("\nAcute Pandemic Period: Expected Deaths as Predicted by the Baseline Period",25), y=375000), 
            colour="black", size=16/.pt) + 
  geom_text(aes(x=as.Date("2023-06-01"), 
                label=stringr::str_wrap("\nPost Pandemic Period: Expected Deaths Estimation Period",25), y=375000), 
            colour="black", size=16/.pt)

## Create plot of monthly trend in observed deaths per capita
#  by first calculating deaths per 100,000, through dividing the total number of deaths by the total population, then multiplying by 100,000. 
#  this represents the mortality rate, which indicates the likelihood of death within a population
observed_p <- ggplot() +
  geom_vline(xintercept = as.numeric(as.Date("2020-01-01")), 
             linetype = "twodash", color = "#292929", size=0.7)+
  geom_vline(xintercept = as.numeric(as.Date("2022-01-01")), 
             linetype = "twodash", color = "#292929", size=0.7)+
  geom_vline(xintercept = as.numeric(as.Date("2024-12-01")), 
             linetype = "twodash", color = "#292929", size=0.7)+
  scale_x_date(breaks = scales::breaks_pretty(n=20), date_labels = "%b %Y", expand = c(0,20), limits=c(as.Date("2014-01-01"),as.Date("2024-12-01"))) +
  scale_y_continuous(breaks=scales::breaks_pretty(n = 5),limits=c(0, 150)) +
  labs(title = "") +
  theme(
    axis.text.x=element_text(angle = 45, vjust = 1, hjust=1),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "lightgrey"),
    panel.grid.major.x = element_line(color = "lightgrey"),
    legend.position="none") +
  labs(x = "",y = "Observed Mortality Rate (per 100,000 Population)") +
  geom_point(data=dat_plot_m, size=2,
             aes(x=date, y=deaths_per_100k, color=Observed)) + 
  scale_color_manual(values=color_list) + 
  scale_fill_manual(values=color_list) +
  geom_text(aes(x=as.Date("2021-01-01"), 
                label=stringr::str_wrap("\nAcute Pandemic Period: Expected Deaths as Predicted by the Baseline Period",22), y=143), 
            colour="black", size=13/.pt) + 
  geom_text(aes(x=as.Date("2023-06-01"), 
                label=stringr::str_wrap("\nPost Pandemic Period: Expected Deaths Estimation Period",29), y=143), 
            colour="black", size=13/.pt)
print(observed_p)

## Save plots 
library(cowplot)
theme_set(theme_cowplot(font_size = 16))
png(file=paste0(outdir,"Figure_Monthly_Excess.png"), height = 3800, width = 4560,  res = 300)
cowplot::plot_grid( 
  expected_p,
  excess_p, ncol=1
)
dev.off()

png(file=paste0(outdir,"Figure_Mortality_Rate.png"), height = 2200, width = 3550, res = 300)
observed_p
dev.off()

### Table with yearly estimates

## Read in and append monthly results from each method
res_yearly <- data.frame()
read_fn_yr <- function(i) {
  res_yearly_i <- read.csv(paste0(outdir,i,"_yearly.csv"), sep=',', header=TRUE) %>%
    mutate(method = gsub("_", " ",i))
  res_yearly <- bind_rows(res_yearly,res_yearly_i)
  
  return(res_yearly)
} 

for(m in list_all_methods){
  res_yearly <- read_fn_yr(m)
}

tbl_yrly <- res_yearly %>%
  dplyr::mutate(excess_estimate = paste0(round(excess/1000,0),"K (",round(excess_lower95/1000,0), "K, ",round(excess_upper95/1000,0),"K)")) %>%
  dplyr::select(method, year, excess_estimate) %>%
  dplyr::arrange(method, desc(year))

write.csv(tbl_yrly,file=paste0(outdir,"Table_2.csv"), row.names = FALSE)

tbl_monthly <- res_monthly %>%
  dplyr::mutate(excess_estimate = paste0(round(excess/1000,0),"K (",round(excess_lower95/1000,0), "K, ",round(excess_upper95/1000,0),"K)")) %>%
  dplyr::select(method, year, month, excess_estimate)

write.csv(tbl_monthly,file=paste0(outdir,"Table_S1.csv"), row.names = FALSE)

### Figure with Yearly Estimates
dat_plot_y <- res_yearly %>%
  mutate(Study = case_when(method == "who" ~ "World Health Organization",
                           method == "economist" ~ "The Economist",
                           method == "harvard" ~ "Acosta and Irizarry",
                           method == "harvard sensitivity" ~ "Acosta and Irizarry Sensitivity Results",
                           method == "ihme" ~ "Institute for Health Metrics and Evaluation",
                           method == "wmd" ~ "World Mortality Dataset",
                           method == "cdc" ~ "Centers for Disease Control and Prevention",
                           TRUE ~ "Missing")) %>%
  # remove partial estimates 
  filter(!(method == "cdc" & year == 2023))

list_study = c(
  "World Health Organization",  
  "The Economist",
  "Acosta and Irizarry",
  "Institute for Health Metrics and Evaluation",
  "Centers for Disease Control and Prevention",
  "World Mortality Dataset"
)
dat_plot_y$Study <- factor(dat_plot_y$Study, levels = list_study)

color_list_y <- c("World Health Organization"="#009E73",
                  "The Economist"="#0072B2",
                  "Acosta and Irizarry"="#D55E00",
                  "Institute for Health Metrics and Evaluation"="#CC79A7",
                  "Centers for Disease Control and Prevention"="grey60",
                  "World Mortality Dataset"="#b80f0a")

library(scales)

forest_p <- ggplot(filter(dat_plot_y,dat_plot_y$year>=2022 & dat_plot_y$method != "harvard sensitivity"), aes(x = excess, y = method, color = Study)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = excess_lower95, xmax = excess_upper95), height = 0.3, show.legend = F) +
  facet_wrap(~year, ncol = 1, strip.position = "left") +
  theme_minimal(base_size = 13) +
  labs(
    title = "",
    x = "Excess Death Estimate (95% UI)",
    y = "",
    color = ""
  ) +
  scale_x_continuous(breaks=scales::breaks_pretty(n = 5),labels = label_number(suffix = "K", scale = 1e-3)) +
  scale_color_manual(values = color_list_y) +
  theme(
    strip.text.y.left = element_text(angle = 90, size = 13),
    strip.placement = "outside",
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = c(0.99, 0.03),
    legend.justification = c(1, 0),
    legend.background = element_rect(fill = "white", color = "gray80"),
    legend.box.background = element_rect(color = "gray80"),
    legend.title=element_blank(),
    legend.text=element_markdown()
  ) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50")

png(file=paste0(outdir,"Figure_Forest_Plot.png"), height = 2400, width = 2700,  res = 300)
forest_p
dev.off()