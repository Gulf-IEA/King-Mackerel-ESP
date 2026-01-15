
# example script
rm(list = ls())

load("indicator_processing/spec_file.RData")


# load data -------------------------------------

url <- ""

# process data ------------------



# save final object ----------------------

write.csv(dat, file = "data/formatted/formatted_csvs/<filename>.csv")  # csv file

save(ind, file = "data/formatted/final_objects/<filename>.RData")      # object for IEAanalyzer

#############################  END  ####################################

print("xxx -- SUCCESSFULLY RUN")


