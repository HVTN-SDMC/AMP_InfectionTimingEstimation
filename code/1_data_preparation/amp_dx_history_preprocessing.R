

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
#
# amp_dx_history_preprocessing.R
#
# -- Code originally by Philip Labuschagne
# -- Adapted for use and executed by C.A. Magaret (cmagaret@fredhutch.org)
#
# Originally run on 19 June 2020 on the study data set:
#   "AMP_diagnostic_testing_history_primary_2020_May_v1_all.csv"
#
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #


# refresh our workspace
rm (list=ls ())

# Initialization
library (tidyverse)

# CONFIG:  define the path to the directory of support files (including the
# prototype "tsic2" package)
path.support <- "/file/path/code/0_support"

# load tsic2 files and supplement
setwd (file.path (path.support, "tsic2/data"))
load ("all_assay_dynamics.rda")
setwd (file.path (path.support, "tsic2/R"))
source ("assay_dynamics.R")
source ("result_interpretation.R")
source ("interpretation_processing.R")
source ("DSHIVInfectionTiming.R")
source ("interpretation_presentation.R")
setwd (path.support)
source ("tsic2_supplement.Rlib")

# initialize our report object
report <- NULL


################################
# 2.a Load the data
################################


# CONFIG:  select which data set we need to use, and the destination of our output files
input_file <- '/file/path/AMP_diagnostic_testing_history_primary_2020_May_v1_all.csv'
output_folder <- '/output/directory'

# load input data
dat <- read.csv (input_file, stringsAsFactors = FALSE)


################################
# 2.b Basic summary metrics
################################


report <- c(report, paste0('Filename: ', basename(input_file), '\n'))
report <- c(report, paste0('Number of participants: ', length(unique(dat$Subject)), '\n'))
report <- c(report, paste0('Number of rows: ', nrow(dat), '\n'))
report <- c(report, paste0('Column names: ', paste(names(dat), sep = ', ', collapse = ', '), '\n'))

# Counts of test result
report <- c(report, 'Test Name, Number Negative, Number Positive')
tmp <- tapply(dat$Subject, list(dat$Test, dat$Result), length, default = 0)
for (i in 1:dim(tmp)[1]){
  report <- c(report, paste(attr(tmp, 'dimnames')[[1]][i],
                            tmp[i, 1],
                            tmp[i, 2],
                            sep = ', '))
}
report <- c(report, '\n')
write.table(report,
            paste0(output_folder, '/',
                   gsub('.csv$', '_preprocessing_report.txt', basename(input_file))), 
            row.names = FALSE,
            sep='',
            quote=FALSE,
            col.names = FALSE)

# Store the number of unique ptids for checks later in script
check_unique_ptids <- length(unique(dat$Subject))

# Very basic checks of the file format
stopifnot(all(names(dat) == c("Subject", "Date", "Test", "Result", "artstartyn", "artstartdtpartial", "artstartdt")))
stopifnot(all(sort(unique(dat$Result)) == c("Negative", "Positive")))


##########################################
# 2.c Data Cleaning - Stage 1: basic remapping
##########################################


# Extract ART initiation data
art_dat <- dat[, c('Subject', 'artstartyn', 'artstartdtpartial', 'artstartdt')]
art_dat <- unique(art_dat)
names(art_dat)[1] <- c('ptid')
art_dat$ptid <- paste0('p_', art_dat$ptid)
art_dat$artstartdt <- as.numeric(as.Date(art_dat$artstartdt))
art_dat$artstartdt[is.na(art_dat$artstartdt)] <- Inf
art_dat$artstartdt <- art_dat$artstartdt + 0.5 #put treatment start at noon
stopifnot(nrow(art_dat) == check_unique_ptids)

# Drop columns
dat <- dat[, c("Subject", "Date", "Test", "Result")]

# Rename columns
names(dat) <- c('ptid', 'sample_date', 'test', 'result')

# Place the tests at noon
dat$sample_date <- as.numeric(as.Date(dat$sample_date)) + 0.5

# Drop rows associated with tests that are not used for the Dx based infection timing.
valid_tests <- c(
  "Abbott ARCHITECT HIV Ag/Ab Combo", "Abbott Real Time HIV01 v1.0 m2000sp/m2000rt",
  "BioRad Geenius Fully Reactive", "BioRad Geenius Indeterminate",
  "BioRad GS HIV Combo Ag/Ab EIA", 'iSCAv2',
  "Roche Taqman v2.0"
  )
dat <- subset(dat, test %in% valid_tests)

valid_test_rows <- nrow(dat) # To check later that no additional rows were dropped

# Remap the test names
## Setup data.frame with the mapping
test_mapping <- matrix(
c(
"Abbott ARCHITECT HIV Ag/Ab Combo", 'architect_weib3_delaney',
"Abbott Real Time HIV01 v1.0 m2000sp/m2000rt", 'abbott_real_time_weib3_delaney_and_manufacturer',
"BioRad Geenius Fully Reactive", 'geenius_fr_weib3_delaney',
"BioRad Geenius Indeterminate", 'geenius_indet_weib3_delaney',
"BioRad GS HIV Combo Ag/Ab EIA", 'gs_combo_weib3_delaney',
'iSCAv2', 'iscav2_weib3_delaney_and_tosiano',
"Roche Taqman v2.0", 'taqman_weib3_delaney_and_manufacturer'
),
  ncol = 2,
  byrow = TRUE)
test_mapping <- data.frame(
  AMP_name = test_mapping[,1],
  tsic_name = test_mapping[,2],
  stringsAsFactors = FALSE)

## Join to map and drop extra columns
dat <- merge(dat, test_mapping,
  by.x = 'test', by.y = 'AMP_name')
dat <- dat[, c("ptid", "sample_date", "tsic_name", "result")]
names(dat) <- c('ptid', 'sample_date', 'test', 'result')
stopifnot(nrow(dat) == valid_test_rows) # Ensure that no rows were lost

# Remap the results
dat$result <- ifelse(dat$result == 'Negative', '-', '+')

# Remap the ptids
dat$ptid <- paste0('p_', dat$ptid)


######################################################
# 2.d Data Cleaning - Stage 2: Restricting to valid visits
######################################################


# 2.d.i Find important for each participant
# Find FP dates
ptid_important_dates <-
  subset(dat, result == '+') %>%
    group_by(ptid) %>% 
    summarize(FP = min(sample_date))

# Find LN dates
t_dat <- merge(dat, ptid_important_dates, by = 'ptid')
ptid_important_dates <- merge(ptid_important_dates, 
  subset(t_dat, sample_date <= FP & result == '-') %>%
    group_by(ptid) %>% 
    summarize(LN = max(sample_date))
  )
rm(t_dat)

# Find First positive Geenius Fully reactive dates
ptid_important_dates <- merge(ptid_important_dates, 
  subset(dat, result == '+' & test == 'geenius_fr_weib3_delaney') %>%
    group_by(ptid) %>% 
    summarize(FP_geenius_fr = min(sample_date)),
  all.x = TRUE) 

# Find last Geenius result
ptid_important_dates <- merge(ptid_important_dates, 
  subset(dat, test %in% c('geenius_indet_weib3_delaney', 'geenius_fr_weib3_delaney')) %>%
    group_by(ptid) %>% 
    summarize(last_geenius = max(sample_date)),
  all.x = TRUE) 

# Find the last useful Geenius test
ptid_important_dates$geenius_cutoff <- 
  with(ptid_important_dates, ifelse(is.na(FP_geenius_fr), last_geenius, FP_geenius_fr))

# Add the ART initiation date into the important date data structure
ptid_important_dates <-
merge(ptid_important_dates, art_dat[, c('ptid', 'artstartdt')], 
      by.x = 'ptid', by.y = 'ptid', all.x = TRUE)

# Dropping visits
## Make a backup
full_dat <- dat
## Merge in the dates that will be used to decide what gets dropped
dat <- merge(dat, ptid_important_dates)
stopifnot(nrow(full_dat) == nrow(dat))

# 2.d.ii Drop all visits before the LN date
dat <- subset(dat, sample_date >= LN)

# 2.d.iii Drop all visits after the ART initiation date
dat <- subset(dat, sample_date <= artstartdt)

# 2.d.iv Drop all visits after the visit that has the last useful Geenius result
dat <- subset(dat, sample_date <= geenius_cutoff)
stopifnot(length(unique(dat$ptid)) == check_unique_ptids)

# Drop columns that are no longer necessary
dat <- dat[, c('ptid', 'sample_date', 'test', 'result')]


#######################################
# 2.e Flag invalid diagnostic histories
#######################################


# 2.e.i.1 Once positive, always positive
violations <- NULL
for (c_ptid in sort(unique(dat$ptid))){
  c_dat <- subset(dat, ptid == c_ptid)
  for (c_test in sort(unique(c_dat$test))){
    cc_dat <- subset(c_dat, test == c_test)
    fp <- min(subset(cc_dat, result == '+')$sample_date)
    any_neg_after_fp <- '-' %in% subset(cc_dat, sample_date > fp)$result
    if (any_neg_after_fp){
      violations <- rbind(violations,
        data.frame(ptid = c_ptid,
                   test = c_test,
                   type = '- after +'))
    }
  }
}

# Build structure that tracks the order of the assays
## based on their window periods
assay_class_ordered <-
data.frame(
matrix(
       c('iscav2_weib3_delaney_and_tosiano', 1,
         'taqman_weib3_delaney_and_manufacturer', 1,
         'abbott_real_time_weib3_delaney_and_manufacturer', 1,
         'architect_weib3_delaney', 2,
         'gs_combo_weib3_delaney', 2,
         'geenius_indet_weib3_delaney', 3,
         'geenius_fr_weib3_delaney', 4),
       ncol = 2,
       byrow = TRUE))
names(assay_class_ordered) <- c('test', 'class')
assay_class_ordered$test <- as.character(assay_class_ordered$test)
assay_class_ordered$class <- as.numeric(assay_class_ordered$class)

# now loop down and perform the checks
c_ptid <- sort(unique(dat$ptid))
for (c_ptid in sort(unique(dat$ptid))){
  c_dat <- subset(dat, ptid == c_ptid)
  # 2.e.i.2 - If a test with a short window period is negative, then all tests from classes with longer window periods must also be negative. 
  for (c_date in sort(unique(c_dat$sample_date))){
    cc_dat <- subset(c_dat, sample_date == c_date)
    for (c_class in sort(unique(assay_class_ordered$class))){
      assays_of_c_class <- subset(assay_class_ordered, class == c_class)$test
      all_slower_must_be_neg <- ('-' %in% subset(cc_dat, test %in% assays_of_c_class)$result)
      violation <- FALSE
      if (all_slower_must_be_neg){
        all_slower_assays <- subset(assay_class_ordered, class > c_class)$test
        violation <- ('+' %in% subset(cc_dat, test %in% all_slower_assays)$result)
        if (violation){
          violations <- rbind(violations,
            data.frame(ptid = c_ptid,
                       test = c_test,
                       type = 'slower assay is +'))
        }
      }
  
    }
  }
  # 2.e.i.3 If a test with a long window period is positive, then all tests from classes with shorter window periods must also be positive. 
  for (c_date in sort(unique(c_dat$sample_date))){
    cc_dat <- subset(c_dat, sample_date == c_date)
    for (c_class in rev(sort(unique(assay_class_ordered$class)))){
      assays_of_c_class <- subset(assay_class_ordered, class == c_class)$test
      all_faster_must_be_pos <- ('+' %in% subset(cc_dat, test %in% assays_of_c_class)$result)
      violation <- FALSE
      if (all_faster_must_be_pos){
        all_faster_assays <- subset(assay_class_ordered, class < c_class)$test
        violation <- ('-' %in% subset(cc_dat, test %in% all_slower_assays)$result)
        if (violation){
          violations <- rbind(violations,
            data.frame(ptid = c_ptid,
                       test = c_test,
                       type = 'faster assay is -'))
        }
      }
  
    }
  }
}

# 2.e.i.4 At least 1 positive, at least 1 negative
## sometimes treatment discards the first visit (just 1 ptid)
for (c_ptid in sort(unique(dat$ptid))){
  c_dat <- subset(dat, ptid == c_ptid)
  at_least_1_neg <- '-' %in% c_dat$result
  at_least_1_pos <- '+' %in% c_dat$result
  if (!(at_least_1_neg & at_least_1_pos)){
    violations <- rbind(violations,
      data.frame(ptid = c_ptid,
                 test = 'All',
                 type = 'Need at least one - and +'))
  }
}

# 2.e.ii Output a list of the ptids that have invalid diagnostic histories
invalid_diagnostic_histories <- data.frame(
  ptid = as.character(unique(violations$ptid)),
  note = 'Invalid pattern of test results',
  stringsAsFactors = FALSE)
write.csv(invalid_diagnostic_histories, 
          paste0(output_folder, '/',
                 gsub('.csv$', '_invalid_histories.csv', basename(input_file))), 
          row.names = FALSE)


#########################################
# 2.f Output dataset for hBayes approach.
#########################################


write.csv(dat, 
          paste0(output_folder, '/',
                 gsub('.csv$', '_for_hBayes.csv', basename(input_file))), 
          row.names = FALSE)


#################################################################
# 2.g Restrict to results that the Weibull3 IDT approach may use.
#################################################################


fastest_to_slowest <- c("iscav2_weib3_delaney_and_tosiano",
            "taqman_weib3_delaney_and_manufacturer",
            "abbott_real_time_weib3_delaney_and_manufacturer",
            "architect_weib3_delaney", 
            "gs_combo_weib3_delaney",
            "geenius_indet_weib3_delaney",
            "geenius_fr_weib3_delaney")

# Reduce full dataset so that each visit only retains the most informative test result
## Makes the next step easier
mi_dat <- NULL
for (c_ptid in unique(dat$ptid)){
  for (c_date in unique(subset(dat, ptid == c_ptid)$sample_date)){

    mi_dat <- rbind(mi_dat,
      select_most_informative_results(subset(dat, ptid == c_ptid & sample_date == c_date),
                                      fastest_to_slowest)$kept_ihist,
      stringsAsFactors = FALSE)
  }
}

tmp_dat <- NULL
for (c_ptid in sort(unique(dat$ptid))){
  c_mi_dat <- subset(mi_dat, c_ptid == ptid)
  FP <- subset(ptid_important_dates, ptid == c_ptid)$FP
  LN <- subset(ptid_important_dates, ptid == c_ptid)$LN
  
  # 2.g.i
  FP_discordant <- '-' %in% subset(c_dat, sample_date == FP)$result
  if (FP_discordant){ # 2.g.i.1
    tmp_dat <- rbind(tmp_dat,
      subset(c_mi_dat, sample_date == FP),
      stringsAsFactors = FALSE
      )
  } else { # 2.g.i.2
    tmp_dat <- rbind(tmp_dat,
      subset(c_mi_dat, sample_date %in% c(LN, FP)),
      stringsAsFactors = FALSE
      )
  }
}
weib3IDT_dat <- tmp_dat
rm(tmp_dat)


##############################################
# 2.h Output dataset for Weibull3 IDT approach
##############################################


write.csv(weib3IDT_dat, 
          paste0(output_folder, '/',
                 gsub('.csv$', '_for_weib3IDT.csv', basename(input_file))), 
          row.names = FALSE)


