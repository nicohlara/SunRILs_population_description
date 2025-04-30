###Script to take a large dataframe with raw phenotypic values and turn them into BLUE phenotypes
## Nicolas A. H. Lara

## set location
dir = 'G:/My Drive/Nico_PhD.lnk/data/'
setwd(dir)

## load in packages
library("dplyr")
library("tidyr")
library("asreml")


## load in data
cross_info <- read.delim('cross_info/cross_info.csv', sep=",") %>% filter(Cross_ID != 'UX1999')
phenotype <- read.delim("phenotype/phenotype.csv", sep=",")


## filter and subset data
cross_names <- c(cross_info$Cross_ID, 'Parent')
locations <- c('Kinston', 'Raleigh')
years <- c('2022', '2023')
effect_variables <- c("Location", "Year", "Cross_ID", "Entry",  "row", "column")
traits <- data.frame(trait = c("WDR", "flowering", "Powdery_mildew", "Height"),
                     type = c("qualitative", "quantitative", "qualitative", "quantitative"))
phenotype <- phenotype %>%
  mutate(Location = as.factor(Location), Cross_ID = as.factor(Cross_ID), Entry=as.factor(Entry), Year = as.factor(Year),
         row = as.factor(row), column = as.factor(column)) %>%
  # mutate(WDR = round(WDR, 0), Powdery_mildew = round(Powdery_mildew, 0)) %>%
  dplyr::filter(Year %in% years & Cross_ID %in% cross_names, Location %in% locations) %>%
  dplyr::select(all_of(c(effect_variables, traits$trait)))


for (row in 1:nrow(traits)) {
  trait <- traits[row, "trait"]
  pf <- dplyr::select(phenotype, all_of(c(effect_variables, trait))) %>%
    filter(!is.na(!!as.symbol(trait))) %>%
    droplevels()
  fix_formula <- as.formula(paste0(trait, "~", "Entry"))
  random_formula <- as.formula("~Location:Year")
  if (traits[row, "type"] == "qualitative") {
    pf[trait] <- ordered(round(pf[[trait]], 0))
    model <- asreml(fixed=fix_formula,
                    random=random_formula,
                    # residual = as.formula("~Location:Year:ar1(row):ar1(column)"),
                    na.action=na.method(x='include', y='include'),
                    data=pf, workspace="8gb",
                    family=asr_multinomial())
  } else {
    pf <- dplyr::mutate(pf, Env = as.factor(paste0(Location, "_", Year))) %>%
      complete(Env, row, column)
    level_num <<- length(unique(pf$Env))
    model <- asreml(fixed=fix_formula,
                    random=random_formula,
                    residual = ~dsum(~ar1(row):ar1(column)| Env, levels=1:level_num),
                    na.action=na.method(x='include', y='include'),
                    data=pf, workspace="8gb")
  }
  BLUE <- summary(model, coef=TRUE)$coef.fixed %>%
    data.frame()
  if (traits[row, "type"] == "quantitative") {
    BLUE <- mutate(BLUE, solution = solution + BLUE['(Intercept)', 'solution'])
  }
  BLUE <- BLUE %>%
    mutate(Entry = sapply(strsplit(rownames(BLUE), "_"), '[', 2)) %>%
    filter(!is.na(Entry)) %>%
    select(c(Entry, solution))
  colnames(BLUE)[colnames(BLUE)=='solution'] <- trait
  if (!exists("blues")) {blues <- BLUE} else {blues <- merge(blues, BLUE, by='Entry')}
}

cross_ID <- unique(select(phenotype, c(Entry, Cross_ID))) %>%
  mutate(Entry = as.character(Entry), Cross_ID = as.character(Cross_ID)) %>%
  unique()
blues <- merge(blues, cross_ID, by="Entry")
blues <- blues[c(ncol(blues), 1:(ncol(blues)-1))]

blues <- blues %>% rename(HD = flowering, PM = Powdery_mildew)

write.table(blues, "C:/Users/nalara/Documents/GitHub/SunRILs_population/data/blues.csv", quote=F, sep=",", row.names=F)
