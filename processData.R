library(dplyr)
library(RColorBrewer)
data.files <- dir("./Data/")
### Processing error pretty much all of 2020 is missing
ignore_files <- c("2020_Tuberculosis.csv",
                  "fitting_data",
                  "fitting_data.csv",
                  "test_data.csv")
data.files <- setdiff(data.files, ignore_files)
all_dat <- sapply(paste0("./Data/", data.files),
                 read.csv)
clean_dat <- lapply(all_dat, 
                      function(dat, col_names, regions){
                        colnames(dat) <- toupper(colnames(dat))
                        col_names <- toupper(col_names)
                        
                        dat %>% 
                          dplyr::select(any_of(!!col_names) | 
                                   ends_with('CurrentQuarter') |
                                   ends_with('CurrentWeek')) %>%
                          dplyr::filter(REPORTINGAREA %in% !!regions) %>%
                          distinct()
                        
                      },
                      col_names = c('ReportingArea', 
                                    'MMWRYear',
                                    'MMWRQuarter',
                                    'MMWRWeek',
                                    'TuberculosisCum2021'),
                      regions = c('UNITED STATES',
                                  'TOTAL'))
X <- Reduce(rbind, clean_dat[-length(clean_dat)])
X <- X %>% arrange(MMWRYEAR, MMWRQUARTER)
X$MMWRYEAR <- factor(X$MMWRYEAR)

cols <- brewer.pal(length(levels(X$MMWRYEAR)),
                   'Pastel2')
jpeg(file="QuarterlyCases.jpeg", quality = 100)
plot(X$TUBERCULOSISCURRENTQUARTER, 
     type = 'l',
     ylim = c(1000, 2000),
     main = "Quarterly Cases Over Time",
     ylab = 'TB Counts',
     xlab = 'Time Index')
points(X$TUBERCULOSISCURRENTQUARTER,
       pch = 19,
       col = cols[as.integer(X$MMWRYEAR)])
legend("top",
       legend = levels(X$MMWRYEAR),
       pch = 19,
       col = cols,
       horiz = TRUE,
       cex = 0.7)
dev.off()
X <- rename(X, TBCOUNT = TUBERCULOSISCURRENTQUARTER)
X2 <- clean_dat$`./Data/2021_Tuberculosis_to_Tularemia.csv` %>% 
  arrange(MMWRWEEK) %>%
  dplyr::filter(MMWRWEEK %in% (13*(1:4)))
count2021 <- c(X2$TUBERCULOSISCUM2021[1], diff(X2$TUBERCULOSISCUM2021))
# %>% 
  # select(-REPORTINGAREA, -MMWRYEAR)
X2 <- X2 %>% mutate(quarter = 1 + floor((MMWRWEEK-1)/13))
X2 <- X2 %>% group_by(quarter) %>% 
  summarise(Count = sum(TUBERCULOSISCURRENTWEEK)) %>%
  ungroup()
# Reduce(union, sapply(clean_dat[-length(clean_dat)], colnames))
write.csv(X, file = "./Data/fitting_data.csv")
write.csv(data.frame(quarter = 1:4,
                     count = count2021), file = "./Data/test_data.csv")       
