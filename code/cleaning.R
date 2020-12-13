library(dplyr)
library(readr)

default_covs <-
  c(
    "e(0)",
    "LA1and10",
    "PovertyRate",
    "MedianFamilyIncome",
    "TractSNAPper",
    "Urban",
    "TractHUNVper",
    "TractBlackper",
    "TractHispanicper"
  )


clean_and_load <- function(myvars=default_covs) {
  #USDA Data on food accessability by census tract
  USDA <- read_csv("../data/USDA_data.zip")
  #CDC Data on average life expectancy by census tract
  life_exp <- read_csv("../data/US_A.zip")
  #removes some track data do to missing observations in life expectancy data, 74,000 rows to 65000 rows
  USDA_Data <-
    dplyr::inner_join(life_exp, USDA, by = c("Tract ID" = "CensusTract"))
  #normalize SNAP as a % of hoursing units in the tract, % of houses with no vehicle in the tract and % of each race in tract
  USDA_Data$TractSNAPper = USDA_Data$TractSNAP / USDA_Data$OHU2010
  USDA_Data$TractHUNVper = USDA_Data$TractHUNV / USDA_Data$OHU2010
  
  USDA_Data$TractBlackper = USDA_Data$TractBlack/USDA_Data$POP2010
  USDA_Data$TractWhiteper = USDA_Data$TractWhite/USDA_Data$POP2010
  USDA_Data$TractHispanicper = USDA_Data$TractHispanic/USDA_Data$POP2010
  
  #create new data and filter NA,NAN and Inf
  USDA_new <- USDA_Data[myvars]
  #make sure if its not systematically missing data
  USDA_new <- USDA_new[Reduce(`&`, lapply(USDA_new, function(x) !is.na(x)  & is.finite(x))),]
  USDA_new <- USDA_new[USDA_new$MedianFamilyIncome !=0,]
  colnames(USDA_new) <-
    c(
      "lifeexp",
      "LA1and10",
      "PovertyRate"   ,
      "MedianFamilyIncome",
      "TractSNAPper",
      "Urban",
      "TractHUNVper",
      "TractBlackper",
      "TractHispanicper"
    )
  USDA_new <<- USDA_new[USDA_new$TractSNAPper <= 1 & USDA_new$TractHUNVper <= 1,]
}

