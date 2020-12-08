# Load libraries
library(readxl)
library(dplyr)
library(reshape2)

# Data path
dh_path <- 'Documents/Massion_lab/Others/Dalton/data'

# Function to read and compile data into one dataframe
dh_comp <- function(path, sheet_ls, name_ls){

    UD <- function(x){ 
        res <- ifelse(x < 0, 'DN', 'UP')
        res 
    }

    comp <- matrix(nrow=0, ncol=5)

    for (i in 1:length(sheet_ls)){
        x <- read_excel(path = path, sheet = sheet_ls[i])
        # The log2(fold-change) is the log-ratio of a gene's or a transcript's expression values in two different conditions.
        x <- x %>% 
            tibble(.) %>%
            mutate(
                L2R = log2(OE / WT),
                state = UD(L2R),
                subtype = name_ls[i]
            ) %>%
            arrange(desc(abs(L2R))) %>%
            relocate(OE, .after=WT) %>%
            relocate(subtype, .after=metabolite) 

        comp <- rbind(comp, x)
    }
    comp
}

# Vector with excel sheet names 
sheet_ls <- excel_sheets(path = file.path(dh_path, 'MetabolonAllSubtypes.xlsx'))

# Vector with names to label rows (same names but replacing spaces with underscores)
name_ls <- chartr(" ", "_", sheet_ls) # alterative use for name_ls argument

data_DH <- dh_comp(path = file.path(dh_path, 'MetabolonAllSubtypes.xlsx'),
                   sheet_ls = sheet_ls,
                   name_ls = sheet_ls)

# Reshaped data
data_DH_rs <- data_DH
data_DH_rs$L2R <- NULL
data_DH_rs <- reshape2::melt(data_DH_rs, value.name = "value",  variable.name='condition')

# Save data
save(data_DH, data_DH_rs, file=file.path(dh_path, 'compdata.RData'))