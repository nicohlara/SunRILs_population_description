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
  dplyr::filter(Year %in% years & Cross_ID %in% cross_names, Location %in% locations) %>%
  dplyr::select(all_of(c(effect_variables, traits$trait)))


for (row in 1:nrow(traits)) {
  trait <- traits[row, "trait"]
  pf <- dplyr::select(phenotype, all_of(c(effect_variables, trait))) %>%
    filter(!is.na(!!as.symbol(trait))) %>%
    droplevels() #%>%
  #   mutate(Env = as.factor(paste0(Location, "_", Year)),
  #          Powdery_mildew = ordered(!!as.symbol(trait))) %>%
  #   # complete(Location, Year, row, column)
  #   complete(Env)
  fix_formula <- as.formula(paste0(trait, "~", "Entry"))
  random_formula <- as.formula("~Location:Year")
  # residual_formula <- as.formula("~Location:Year:ar1(row):ar1(column)")
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
write.table(blues, "C:/Users/nalara/Documents/GitHub/SunRILs_population/data/blues.csv", quote=F, sep=",", row.names=F)


# traits = c('Powdery_mildew', "WDR")
# for (trait in traits) {
#   pheno_file <- select(phenotype, c(Location, Year, Cross_ID, Entry, row, column, !!as.symbol(trait))) %>%
#     filter(!is.na(!!as.symbol(trait))) %>%
#     droplevels() %>%
#     mutate(Env = as.factor(paste0(Location, "_", Year)), row = as.factor(row), column = as.factor(column),
#            Powdery_mildew = ordered(!!as.symbol(trait))) %>%
#     #complete(Env, row, column)
#     complete(Env)
#   level_num <<- length(unique(pheno_file$Env))
#   print(level_num)
#   fix.formula <- as.formula(paste0(trait, '~', 'Entry'))
#   row_col <- asreml(fixed = fix.formula,
#                     random = ~Location:Year,
#                     residual = ~Location:Year:ar1(row):ar1(column),
#                     #residual=~dsum(~ar1(row):ar1(column)| Env, levels=1:level_num),
#                     data=pheno_file,
#                     na.action = na.method(x='include', y='include'),
#                     workspace="8gb",
#                     family=asr_multinomial())
#   BLUE <- summary(row_col, coef=TRUE)$coef.fixed %>%
#     data.frame()
#   BLUE <- BLUE %>%
#     mutate(#solution = solution + BLUE['(Intercept)', 'solution'],
#       Entry = sapply(strsplit(rownames(BLUE), "_"), '[', 2)) %>%
#     filter(!is.na(Entry)) %>%
#     select(c(Entry, solution))
#   cross_ID <- unique(select(pheno_file, c(Entry, Cross_ID))) %>%
#     mutate(Entry = as.character(Entry), Cross_ID = as.character(Cross_ID))
#   BLUEs <- merge(BLUE, cross_ID, by="Entry")
# }
# 
# make_blue <- function(pheno_file, trait) {
#   pheno_file <- select(pheno_file, c(Location, Year, Cross_ID, Entry, row, column, !!as.symbol(trait))) %>%
#     filter(!is.na(!!as.symbol(trait))) %>%
#     droplevels() %>%
#     mutate(Env = as.factor(paste0(Location, "_", Year))) %>%
#     complete(Env, row, column)
#   level_num <<- length(unique(pheno_file$Env))
#   print(level_num)
#   fix.formula <- as.formula(paste0(trait, '~', 'Entry'))
#   row_col <- asreml(fixed = fix.formula,
#                     random = ~Location:Year,
#                     #residual = ~Location:Year:ar1(row):ar1(column),
#                     residual=~dsum(~ar1(row):ar1(column)| Env, levels=1:level_num),
#                     data=pheno_file,
#                     na.action = na.method(x='include', y='include'),
#                     workspace="8gb")
#  
# }
# 
# BLUE_pheno <- phenotype %>%
#   filter(Cross_ID %in% cross_names & Location %in% locations & Year %in% years) %>%
#   mutate(row = as.factor(row), column = as.factor(column)) %>%
#   select(-c(Awns, Timepoint_height, Zadok.s.Growth.Scale, Frozen, Sample_date, Weight, Plot_number, Moisture_protein)) %>%
#   droplevels
