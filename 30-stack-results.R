library(tidyverse)
library(lubridate)
library(ggplot2)
library(ggtext)
setwd("[insert file path to working directory]")
datadir <- "./data/"
outdir <- "./output/"

dat_monthly <- read.csv(paste0(datadir,"monthly data.csv"), sep=',', header=TRUE)%>%
  rename(observed = outcome)

list_mo_methods <- c("who","economist","ihme")
list_wk_methods <- c("harvard")

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

## Read in and append weekly results
res_weekly <- data.frame()
read_fn_wk <- function(i) {
  res_weekly_i <- read.csv(paste0(outdir,i,"_weekly.csv"), sep=',', header=TRUE) %>%
    mutate(method = gsub("_"," ",i))
  res_weekly <- bind_rows(res_weekly,res_weekly_i)
  
  return(res_weekly)
} 

for(w in list_wk_methods){
  res_weekly <- read_fn_wk(w)
}
res_weekly <- res_weekly %>% filter(!is.na(date))

## Prepare data for plots
dat_plot_m <- res_m_wobs %>%
  select(-c(observed.y))%>%
  rename(observed = observed.x) %>%
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
         # Recode bounds of uncertainty interval to be the point estimate, for plotting methods where there is no uncertainty interval
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

dat_plot_w <- res_weekly %>%
  rename(expected_deaths = expected)%>%
  mutate(Study = case_when(method == "who" ~ "World Health Organization",
                           method == "economist" ~ "The Economist",
                           method == "harvard" ~ "Acosta and Irizarry",
                           method == "harvard sensitivity" ~ "Acosta and Irizarry Sensitivity Results",
                           method == "ihme" ~ "Institute for Health Metrics and Evaluation",
                           method == "wmd" ~ "World Mortality Dataset",
                           method == "cdc" ~ "Centers for Disease Control and Prevention",
                           TRUE ~ "Missing"),
         Observed = "Observed",
         # Recode bounds of uncertainty interval to be the point estimate, for plotting methods where there is no uncertainty interval
         expected_lower95 = ifelse(method %in% c("wmd", "cdc"), expected_deaths, expected_lower95),
         expected_upper95 = ifelse(method %in% c("wmd", "cdc"), expected_deaths, expected_upper95),
         excess_lower95 = ifelse(method %in% c("wmd", "cdc"), excess, excess_lower95),
         excess_upper95 = ifelse(method %in% c("wmd", "cdc"), excess, excess_upper95),
         no_UI = ifelse(method %in% c("wmd", "cdc"), 1, 0),
         date = as.Date(date),
         year = as.numeric(year(date)),
         across(starts_with("excess"), function(x) ifelse(year<=2019,NA,x)),
         across(starts_with("expected"), function(x) ifelse(year<=2019,NA,x)))

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
  "Institute for Health Metrics and Evaluation",
  "World Health Organization",
  "The Economist",
  "Centers for Disease Control and Prevention"
)
list_study_w = c(
  "Acosta and Irizarry",
  "Acosta and Irizarry Sensitivity Results",
  "World Mortality Dataset"
)
dat_plot_m$Study <- factor(dat_plot_m$Study, levels = list_study)
dat_plot_w$Study <- factor(dat_plot_w$Study, levels = list_study_w)

list_observed = c("Observed")
dat_plot_m$Observed <- factor(dat_plot_m$Observed, levels = list_observed)
dat_plot_w$Observed <- factor(dat_plot_w$Observed, levels = list_observed)

color_list <- c("Institute for Health Metrics and Evaluation"="#66c2a5",
                "World Health Organization"="#8da0cb",
                "The Economist"="#fc8d62",
                "Centers for Disease Control and Prevention"="grey60",
                "Observed"="black")
color_list_w <- c("Acosta and Irizarry"="navy",
                  "Acosta and Irizarry Sensitivity Results" = "#6baed6",
                  "World Mortality Dataset"="#b80f0a",
                  "Observed"="black")


## Monthly plots for IHME, WHO, and the Economist
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
  scale_y_continuous(breaks=scales::breaks_pretty(n = 10), limits=c(-100000, 250000)) +
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


## Weekly plots for Acosta and Irizarry
excess_p_w <- ggplot() +
  stat_smooth(data=filter(dat_plot_w, dat_plot_w$year>=2022 & !(dat_plot_w$method %in% c("wmd","cdc"))),
              geom="ribbon",aes(x=date,y = excess,ymin = excess_lower95,ymax = excess_upper95,fill=Study),
              alpha = 0.3, span=0.3, show.legend = F
  ) +
  geom_smooth(data=filter(dat_plot_w, dat_plot_w$year>=2022),
              aes(x=date,y = excess, colour = Study),span=0.3, se=F,linewidth=1) + 
  scale_x_date(breaks = scales::breaks_pretty(n = 36), date_labels = "%b %Y", expand = c(0, 4)
  ) +
  scale_y_continuous(breaks=scales::breaks_pretty(n = 10)) +
  labs(title = "Fig A2-B. Weekly Excess Deaths, January 2022-December 2024") +
  theme_excess +
  labs(x = "",
       y = "Excess Deaths",
       color="") +
  scale_color_manual(values = color_list_w) +
  scale_fill_manual(values = color_list_w) +
  geom_line(data=filter(dat_plot_w, dat_plot_w$year>=2022),
            aes(x=date,y=0))

expected_p_w <- ggplot() +
  geom_vline(xintercept = as.numeric(as.Date("2020-01-01")), 
             linetype = "twodash", color = "#292929", size=0.7)+
  geom_vline(xintercept = as.numeric(as.Date("2022-01-01")), 
             linetype = "twodash", color = "#292929", size=0.7)+
  geom_vline(xintercept = as.numeric(as.Date("2024-12-28")), 
             linetype = "twodash", color = "#292929", size=0.7)+
  stat_smooth(data=filter(dat_plot_w, !(dat_plot_w$method %in% c("wmd","cdc"))),
              geom="ribbon",
              aes(x = date,
                  y = expected_deaths,
                  ymax = expected_upper95,
                  ymin = expected_lower95, fill=Study),
              alpha = 0.3, span=0.3, show.legend = F
  ) +
  geom_smooth(data=dat_plot_w, 
              aes(x=date, y = expected_deaths, colour = Study),span=0.3, se=F,linewidth=0.9) + 
  scale_x_date(breaks = scales::breaks_pretty(n=20), date_labels = "%b %Y", expand = c(0,0), limits=c(as.Date("2014-01-01"),as.Date("2024-12-28"))) +
  scale_y_continuous(breaks=scales::breaks_pretty(n = 10),limits=c(40000, 100000)) +
  labs(title = "Fig A2-A. Weekly Observed and Expected Deaths, January 2014-December 2024") +
  theme_expected +
  labs(x = "",y = "Deaths") +
  geom_point(data=dat_plot_w,
             aes(x=date, y=observed,color=Observed)) + 
  scale_color_manual(values=color_list_w) + 
  scale_fill_manual(values=color_list_w) +
  geom_text(aes(x=as.Date("2021-01-01"), 
                label=stringr::str_wrap("\nAcute Pandemic Period: Expected Deaths as Predicted by the Baseline Period",25), y=95200), 
            colour="black", size=16/.pt) + 
  geom_text(aes(x=as.Date("2023-06-15"), 
                label=stringr::str_wrap("\nPost Pandemic Period: Expected Deaths Estimation Period",25), y=95000), 
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

png(file=paste0(outdir,"Figure_Weekly_Excess.png"), height = 3800, width = 4550, res = 300)
cowplot::plot_grid( 
  expected_p_w, excess_p_w, ncol = 1
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
  mutate(excess_estimate = paste0(round(excess/1000,1),"K (",round(excess_lower95/1000,1), "K, ",round(excess_upper95/1000,1),"K)")) %>%
  select(method, year, excess_estimate) %>%
  arrange(method, desc(year))

write.csv(tbl_yrly,file=paste0(outdir,"Table_2.csv"), row.names = FALSE)

tbl_monthly <- res_monthly %>%
  mutate(excess_estimate = paste0(round(excess/1000,1),"K (",round(excess_lower95/1000,1), "K, ",round(excess_upper95/1000,1),"K)")) %>%
  select(method, year, month, excess_estimate)

write.csv(tbl_monthly,file=paste0(outdir,"Table_S1.csv"), row.names = FALSE)