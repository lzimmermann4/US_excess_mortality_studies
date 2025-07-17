library(httr2)
library(jsonlite)
library(data.table)
library(lubridate)
library(excessmort)
library(tidyverse)

setwd("[insert file path to working directory]")
codedir<-"./code/"
datadir<-"./data/"

end_date_ymd <- "2024-12-31"

# ------------------------ Data Preparation ------------------------
## Download data from CDC
mmwr_to_date <- function(mmwr_year, mmwr_week) {
  first_day <- floor_date(make_date(mmwr_year, 1, 4) , unit = "week")
  date <- first_day + weeks(mmwr_week - 1) + days(6) 
  return(date)
}

api <- "https://data.cdc.gov/resource/3yf8-kanr.json"
dt2014_2019 <- request(api) |> 
  req_url_query("$limit" = 10000000) |>
  req_perform() |> resp_body_string() |> 
  fromJSON(flatten = TRUE) |>
  setDT()

dat1 <- dt2014_2019[jurisdiction_of_occurrence == "United States"]

if (!identical(with(dat1, mmwr_to_date(as.numeric(mmwryear), as.numeric(mmwrweek))), as_date(dat1$weekendingdate))) stop("something wrong with dates.")
dat1[, date := as_date(weekendingdate)]
dat1 <- dat1[, c("date", "allcause")]
setnames(dat1, "allcause", "outcome")

api <- "https://data.cdc.gov/resource/r8kw-7aab.json"
dt2020_present <- request(api) |> 
  req_url_query("$limit" = 10000000) |>
  req_perform() |> resp_body_string() |> 
  fromJSON(flatten = TRUE) |>
  setDT()

dat2 <- dt2020_present[state == "United States" & group == "By Week"]
dat2 <- dat2[, c("end_date", "total_deaths","covid_19_deaths")]
setnames(dat2, c("end_date", "total_deaths","covid_19_deaths"), c("date", "outcome","covid_deaths"))
dat2[, date := as_date(date)]

dat <- bind_rows(list(dat1, dat2))
dat[, outcome := as.numeric(outcome)]
dat[, covid_deaths := as.numeric(covid_deaths)]
dat <- dat[date <= as_date(end_date_ymd)]
dat <- dat[order(date),]

## Get population from United States Census
source(paste0(codedir,"census-key.R"))  

api <- "https://api.census.gov/data/2019/pep/population"
pop1 <- request(api) |>  
  req_url_query(get = I("POP,DATE_CODE,DATE_DESC"), 
                `for` = I("us:*"),
                key = census_key) |>
  req_perform() |>
  resp_body_string() |> 
  fromJSON(flatten = TRUE) |>
  as.data.table()

pop1 <- pop1 |> janitor::row_to_names(1) 
pop1 <- pop1[!grepl("4/1", DATE_DESC)]
pop1 <- data.table(year = 2010 + as.numeric(pop1$DATE_CODE) - 3, population = as.numeric(pop1$POP))

pop2 <- fread("https://www2.census.gov/programs-surveys/popest/datasets/2020-2024/state/totals/NST-EST2024-ALLDATA.csv") 
years <- 2020:2024
pop2 <- pop2[NAME == "United States", ]
cols <- paste0("POPESTIMATE", years)
pop2 <- data.table(year = years, population = unlist(pop2[,..cols]))

pop <- rbindlist(list(pop1,pop2))[order(year)]
pop <- approx_demographics(pop, first_day = min(dat$date), last_day = max(dat$date), )

dat <- merge(dat, pop, by = "date", all.x = TRUE)

## Create additional variables for analysis
dat_weekly <- dat %>%
  mutate(year = as.numeric(format(as.Date(date, format="%m/%d/%Y"),"%Y")),
         month = as.numeric(format(as.Date(date, format="%m/%d/%Y"),"%m")),
         date_yr_mnth = make_date(year = year, month = month)
        )%>%
  group_by(year, month) %>%
  mutate(week = row_number())%>%
  ungroup()%>%
  relocate(covid_deaths,.after = week)

# Obtain monthly observed deaths
dat_monthly <- dat_weekly %>% 
  group_by(date_yr_mnth, year, month) %>%
  summarise(outcome = sum(outcome),
            covid_deaths = sum(covid_deaths),
            population = max(population)) %>%
  ungroup()

write.csv(dat_weekly,file=paste0(datadir,"/weekly data.csv"), row.names = FALSE)
write.csv(dat_monthly,file=paste0(datadir,"/monthly data.csv"), row.names = FALSE)
