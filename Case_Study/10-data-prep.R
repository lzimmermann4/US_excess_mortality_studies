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
## Download data from CDC Wonder
library(cdcwonder)

df_final <- cdcwonder::wonder(
  database   = "d76_mortality",       # final 1999–2020
  group_by   = c("Year", "Month"),    # monthly aggregation
  start_year = 2014,
  end_year   = 2020,
  measure    = "Deaths",              # all-cause mortality
  api_key    = NULL                 
)


df_prov <- cdcwonder::wonder(
  database = "mcd_provisional",   # provisional
  group_by = c("Month", "Year"),
  start_year = 2021,
  end_year = 2024,
  measure = "Deaths", 
  api_key = NULL  
)

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
dat <- merge(dat, pop, by = "date", all.x = TRUE)

## Create additional variables for analysis
dat_monthly <- dat %>%
  mutate(year = as.numeric(format(as.Date(date, format="%m/%d/%Y"),"%Y")),
         month = as.numeric(format(as.Date(date, format="%m/%d/%Y"),"%m")),
         date_yr_mnth = make_date(year = year, month = month)
  )

write.csv(dat_monthly,file=paste0(datadir,"/monthly data.csv"), row.names = FALSE)
