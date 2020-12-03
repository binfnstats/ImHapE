# This R script does:
# -------------------
# 1. Organizes geographical locations for each of the SARS-CoV-2 sequences obtained from GISAID and aligned with EMBOSS stretcher
# 2. Get hamming distances for each sequence (relative to Wuhan reference)
# 3. Gets VAFs for SARS-CoV-2 sequences from latest timepoints available given location and 
# 4. Saves aligned empirical SARS-CoV-2 sequences as binary-encoded text files

library(data.table)
library(dplyr)

# ==================================
# 1. ORGANIZE GEOGRAPHICAL LOCATIONS
# ==================================
names <- list.files('./alignments')
data <- data.table(filename = names, country = sapply(strsplit(names, '+', 3), '[', 1), id = sapply(strsplit(names, '+', 3), '[', 2), date = sapply(strsplit(sapply(strsplit(names, '+', 3), '[', 3), '.', 2), '[', 1))

library(maps)
library(countrycode)

# City and countries
data(world.cities)
colnames(world.cities) <- c('city', 'country', 'pop', 'lat', 'long', 'capital')
data$country <- ifelse(data$country == "ITALY", 'Italy', data$country)
data$country <- ifelse(data$country == "Bahrein", 'Bahrain', data$country)
data$country <- ifelse(data$country == "NanChang", 'Nanchang', data$country)
data$country <- ifelse(data$country == "Bucuresti", 'Bucharest', data$country)
data$city <- ifelse(data$country %in% world.cities$city & !data$country %in% world.cities$country, data$country, NA)
data$country <- ifelse(data$country %in% world.cities$city & !data$country %in% world.cities$country, world.cities$country, data$country)
data$city <- ifelse(data$country %in% c('Fujian', 'Guangdong', 'Henan', 'Jiangsu', 'Jiangxi', 'Liaoning', 'Shandong', 'Sichuan', 'Yunnan', 'Zhejiang', 'Yingtan'), data$country, data$city)
data$country <- ifelse(data$country %in% c('Fujian', 'Guangdong', 'Henan', 'Jiangsu', 'Jiangxi', 'Liaoning', 'Shandong', 'Sichuan', 'Yunnan', 'Zhejiang', 'Yingtan'), 'China', data$country)
data <- data[!data$country %in% c('tiger', 'mink', 'pangolin', 'env')]

# Continents
data$continent <- countrycode(sourcevar = data$country, origin = 'country.name', destination = 'continent')
data$continent <- ifelse(data$country == 'England', 'Europe', data$continent)
data$continent <- ifelse(data$country == 'Scotland', 'Europe', data$continent)
data$continent <- ifelse(data$country == 'Serbia and Montenegro', 'Europe', data$continent)
data$continent <- ifelse(data$country == 'Wales', 'Europe', data$continent)

# Regions
data$region <- countrycode(sourcevar = data$country, origin = 'country.name', destination = 'region')
data$region <- ifelse(data$country == 'England', 'Europe & Central Asia', data$region)
data$region <- ifelse(data$country == 'Scotland', 'Europe & Central Asia', data$region)
data$region <- ifelse(data$country == 'Wales', 'Europe & Central Asia', data$region)
data$month <- sapply(strsplit(data$date, '-'), '[', 2)

write.table(data, file = './data/covid_ids.csv', sep = ',', quote = FALSE, row.names = FALSE)

# ==========================================
# 2. GET HAMMING DISTANCES FOR EACH SEQUENCE
# ==========================================
vcf <- fread('./data/covid_vcf.csv')
haplotypes <- apply(vcf[,3:ncol(vcf)], 2, function(x) as.numeric(vcf$reference != x & x != 'N'))
hamming_distance <- apply(haplotypes, 2, function(x) sum(x))
hamming <- data.table(filename = colnames(vcf[,3:ncol(vcf)]), hamming = hamming_distance)
write.table(hamming, file = './processed/covid_hamming.csv', quote = FALSE, row.names = FALSE, sep =',')


# ================================================================
# 3. GET VAFS FOR LATEST TIMEPOINTS FOR SPECIFIC COUNTRIES/REGIONS
# ================================================================
ids <- c('England', 'Scotland', 'Wales', 'USA', 'India', 'Europe & Central Asia', 'North America', 'South Asia')
data$month <- as.numeric(data$month)
data <- filter(data, month != 12)
for (i in ids) {
  if (i %in% c('England', 'Scotland', 'Wales', 'USA', 'India')) {
    files <- filter(data, country == i & month == max(filter(data, country == i)$month))$filename  
  } else {
    files <- filter(data, region == i & month == max(filter(data, region == i)$month))$filename  
  }
  files <- files[files %in% colnames(vcf)]
  haps <- haplotypes[,files]
  vafs <- apply(haps, 1, function(x) sum(x))
  vafs <- vafs/ncol(haps)
  vafs[is.na(vafs)==TRUE] <- 0
  print(quantile(vafs))
  write.table(vafs, file = paste('./data/covid_vaf_', gsub(" ", "", i), '.tsv', sep = ''), quote = FALSE, row.names = FALSE, sep = '\t')
}

# ==============================================================================
# 4. STORE EMPIRICAL SARS-COV-2 SEQUENCES AS BINARY ENCODED HAPLOTYPE TEXT FILES
# ==============================================================================
# Haplotypes
for (i in colnames(haplotypes)) {
	seq <- paste(haplotypes[,i], collapse = '')
	save_seq <- file(paste('./haplotypes/', i, sep = ''))
	writeLines(seq, save_seq)
	close(save_seq)
}
