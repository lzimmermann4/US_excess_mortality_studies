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
    mutate(method = i)
  res_weekly <- bind_rows(res_weekly,res_weekly_i)
  
  return(res_weekly)
} 

for(w in list_wk_methods){
  res_weekly <- read_fn_wk(w)
}

## Read in and append sensitivity results 
res_weekly_sens <- read.csv(paste0(outdir,"harvard_weekly_sensitivity.csv"), sep=',', header=TRUE) %>%
  mutate(method = "harvard sensitivity")
res_weekly <- bind_rows(res_weekly,res_weekly_sens)

## Prepare data for plots
dat_plot_m <- res_m_wobs %>%
  mutate(Study = case_when(method == "who" ~ "Msemburi et al<sup>1</sup>",
                           method == "economist" ~ "The Economist and Solstad<sup>2</sup>",
                           method == "harvard" ~ "Acosta and Irizarry<sup>3</sup>",
                           method == "harvard sensitivity" ~ "Acosta and Irizarry<sup>3</sup> Sensitivity Results",
                           method == "ihme" ~ "Institute for Health Metrics and Evaluation<sup>4</sup>",
                           TRUE ~ "Missing"),
         Observed = "Observed",
         excess_CI_width = excess_upper95 - excess_lower95,
         expected_CI_width = expected_upper95 - expected_lower95,
         across(starts_with("excess"), function(x) round(x, 0)),
         across(starts_with("expected"), function(x) round(x, 0))
         ) %>%
  arrange(year,month,Study)
  
dat_plot_w <- res_weekly %>%
  rename(expected_deaths = expected)%>%
  mutate(Study = case_when(method == "who" ~ "Msemburi et al<sup>1</sup>",
                           method == "economist" ~ "The Economist and Solstad<sup>2</sup>",
                           method == "harvard" ~ "Acosta and Irizarry<sup>3</sup>",
                           method == "harvard sensitivity" ~ "Acosta and Irizarry<sup>3</sup> Sensitivity Results",
                           method == "ihme" ~ "Institute for Health Metrics and Evaluation<sup>4</sup>",
                           TRUE ~ "Missing"),
         Observed = "Observed",
         date = as.Date(date),
         year = as.numeric(year(date)),
         across(starts_with("excess"), function(x) ifelse(year<=2019,NA,x)),
         across(starts_with("expected"), function(x) ifelse(year<=2019,NA,x)))

## Prepare plot parameters
theme_expected <- theme(plot.title=element_text(size=10, family="Arial", color="black"), 
                        axis.title=element_text(size=10, family="Arial",color="black"),
                        axis.text.x=element_text(size=10, family="Arial",color="black", angle = 45, vjust = 1, hjust=1),
                        axis.text.y=element_text(size=10, family="Arial",color="black"),
                        panel.background = element_rect(fill = "white"),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.grid.major.y = element_line(color = "lightgrey"),
                        panel.grid.major.x = element_line(color = "lightgrey"),
                        legend.position = c(0.193, 0.98),
                        legend.justification = c("right", "top"),
                        legend.box.just = "right",
                        legend.margin = margin(0.5,0.5,0.5,0.5),
                        legend.box.background = element_rect(color="lightgrey",linewidth=1),
                        legend.title=element_blank(),
                        legend.text=element_markdown()
) 
theme_excess <-   theme(plot.title=element_text(size=10, family="Arial", color="black"), 
                        axis.title=element_text(size=10, family="Arial",color="black"),
                        axis.text.x=element_text(size=10, family="Arial",color="black", angle = 45, vjust = 1, hjust=1),
                        axis.text.y=element_text(size=10, family="Arial",color="black"),
                        panel.background = element_rect(fill = "white"),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.grid.major.y = element_line(color = "lightgrey"),
                        panel.grid.major.x = element_line(color = "lightgrey"),
                        legend.position = c(0.955, 0.98),
                        legend.justification = c("right", "top"),
                        legend.box.just = "right",
                        legend.margin = margin(0.5,0.5,0.5,0.5),
                        legend.box.background = element_rect(color="lightgrey",linewidth=1),
                        legend.title=element_blank(),
                        legend.text=element_markdown()
)

list_study = c(
  "Institute for Health Metrics and Evaluation<sup>4</sup>",
  "Msemburi et al<sup>1</sup>",
  "The Economist and Solstad<sup>2</sup>"
)
list_study_w = c(
  "Acosta and Irizarry<sup>3</sup>",
  "Acosta and Irizarry<sup>3</sup> Sensitivity Results"
)
dat_plot_m$Study <- factor(dat_plot_m$Study, levels = list_study)
dat_plot_w$Study <- factor(dat_plot_w$Study, levels = list_study_w)

list_observed = c("Observed")
dat_plot_m$Observed <- factor(dat_plot_m$Observed, levels = list_observed)
dat_plot_w$Observed <- factor(dat_plot_w$Observed, levels = list_observed)

color_list <- c("Institute for Health Metrics and Evaluation<sup>4</sup>"="#66c2a5",
                "Msemburi et al<sup>1</sup>"="#8da0cb",
                "The Economist and Solstad<sup>2</sup>"="#fc8d62",
                "Observed"="black")
color_list_w <- c("Acosta and Irizarry<sup>3</sup>"="#6baed6",
                  "Acosta and Irizarry<sup>3</sup> Sensitivity Results" = "grey60",
                  "Observed"="black")

## Monthly plots for IHME, WHO, and the Economist
excess_p <- ggplot(data=filter(dat_plot_m, dat_plot_m$year>=2022), aes(x=date)) +
  stat_smooth(
    geom="ribbon",aes(y = excess,ymin = excess_lower95,ymax = excess_upper95,fill=Study),
    alpha = 0.3, span=0.3, show.legend = F
  ) +
  geom_smooth(aes(y = excess, colour = Study),span=0.3, se=F,linewidth=1) + 
  scale_x_date(breaks = scales::breaks_pretty(n = 24), date_labels = "%b %Y", expand = c(0, 0), limits=c(as.Date("2022-01-01"),as.Date("2023-12-01"))) +
  scale_y_continuous(breaks=scales::breaks_pretty(n = 10), limits=c(-100000, 250000)) +
  labs(title = "1B. Monthly Excess Deaths, January 2022-December 2023") +
  theme_excess +
  labs(x = "",
       y = "Excess Deaths",
       color="") +
  scale_color_manual(values = color_list) +
  scale_fill_manual(values = color_list) +
  geom_line(aes(y=0))

expected_p <- ggplot(dat_plot_m, aes(x=date)) +
  geom_vline(xintercept = as.numeric(as.Date("2020-01-01")), 
             linetype = "twodash", color = "#292929", size=0.7)+
  geom_vline(xintercept = as.numeric(as.Date("2022-01-01")), 
             linetype = "twodash", color = "#292929", size=0.7)+
  geom_vline(xintercept = as.numeric(as.Date("2023-12-01")), 
             linetype = "twodash", color = "#292929", size=0.7)+
  stat_smooth(
    geom="ribbon",
    aes(y = expected_deaths,
        ymax = expected_upper95,
        ymin = expected_lower95, fill=Study),
    alpha = 0.3, span=0.3, show.legend = F
  ) +
  geom_smooth(aes(y = expected_deaths, colour = Study),span=0.3, se=F,linewidth=0.9) + 
  scale_x_date(breaks = scales::breaks_pretty(n=20), date_labels = "%b %Y", expand = c(0,20), limits=c(as.Date("2014-01-01"),as.Date("2023-12-01"))) +
  scale_y_continuous(breaks=scales::breaks_pretty(n = 5),limits=c(150000, 430000)) +
  labs(title = "1A. Monthly Observed and Expected Deaths, January 2014-December 2023") +
  theme_expected +
  labs(x = "",y = "Deaths") +
  geom_point(aes(y=observed,color=Observed)) + 
  scale_color_manual(values=color_list) + 
  scale_fill_manual(values=color_list) +
  geom_text(aes(x=as.Date("2021-01-01"), 
                label=stringr::str_wrap("\nAcute Pandemic Period: Expected Deaths as Predicted by the Baseline Period",25), y=375000), 
            colour="black", size=11/.pt) + 
  geom_text(aes(x=as.Date("2023-01-01"), 
                label=stringr::str_wrap("\nPost Pandemic Period: Expected Deaths Estimation Period",25), y=375000), 
            colour="black", size=11/.pt)


## Weekly plots for Acosta and Irizarry
excess_p_w <- ggplot(data=filter(dat_plot_w, dat_plot_w$year>=2022), aes(x=date)) +
  stat_smooth(
    geom="ribbon",aes(y = excess,ymin = excess_lower95,ymax = excess_upper95,fill=Study),
    alpha = 0.3, span=0.3, show.legend = F
  ) +
  geom_smooth(aes(y = excess, colour = Study),span=0.3, se=F,linewidth=1) + 
  scale_x_date(breaks = scales::breaks_pretty(n = 24), date_labels = "%b %Y", expand = c(0, 12), limits=c(as.Date("2022-01-01"),as.Date("2023-12-01"))) +
  scale_y_continuous(breaks=scales::breaks_pretty(n = 10)) +
  labs(title = "S2B. Weekly Excess Deaths, January 2022-December 2023") +
  theme_excess +
  labs(x = "",
       y = "Excess Deaths",
       color="") +
  scale_color_manual(values = color_list_w) +
  scale_fill_manual(values = color_list_w) +
  geom_line(aes(y=0))

expected_p_w <- ggplot(dat_plot_w, aes(x=date)) +
  geom_vline(xintercept = as.numeric(as.Date("2020-01-01")), 
             linetype = "twodash", color = "#292929", size=0.7)+
  geom_vline(xintercept = as.numeric(as.Date("2022-01-01")), 
             linetype = "twodash", color = "#292929", size=0.7)+
  geom_vline(xintercept = as.numeric(as.Date("2023-12-01")), 
             linetype = "twodash", color = "#292929", size=0.7)+
  stat_smooth(
    geom="ribbon",
    aes(y = expected_deaths,
        ymax = expected_upper95,
        ymin = expected_lower95, fill=Study),
    alpha = 0.3, span=0.3, show.legend = F
  ) +
  geom_smooth(aes(y = expected_deaths, colour = Study),span=0.3, se=F,linewidth=0.9) + 
  scale_x_date(breaks = scales::breaks_pretty(n=20), date_labels = "%b %Y", expand = c(0,0), limits=c(as.Date("2014-01-01"),as.Date("2023-12-01"))) +
  scale_y_continuous(breaks=scales::breaks_pretty(n = 10),limits=c(40000, 100000)) +
  labs(title = "S2A. Weekly Observed and Expected Deaths, January 2014-December 2023") +
  theme_expected +
  labs(x = "",y = "Deaths") +
  geom_point(aes(y=observed,color=Observed)) + 
  scale_color_manual(values=color_list_w) + 
  scale_fill_manual(values=color_list_w) +
  geom_text(aes(x=as.Date("2021-01-01"), 
                label=stringr::str_wrap("\nAcute Pandemic Period: Expected Deaths as Predicted by the Baseline Period",25), y=95000), 
            colour="black", size=11/.pt) + 
  geom_text(aes(x=as.Date("2023-01-01"), 
                label=stringr::str_wrap("\nPost Pandemic Period: Expected Deaths Estimation Period",25), y=95000), 
            colour="black", size=11/.pt)

## Save plots 
png(file=paste0(outdir,"Figure_1.png"), height = 4800, width = 5550,  res = 300)
cowplot::plot_grid( 
  expected_p,
  excess_p, ncol=1
)
dev.off()

png(file=paste0(outdir,"Figure_S2.png"), height = 4800, width = 5550, res = 300)
cowplot::plot_grid( 
  expected_p_w, excess_p_w, ncol = 1
)
dev.off()