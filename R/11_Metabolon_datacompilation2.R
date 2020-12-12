# Load libraries
library(readxl)
library(dplyr)
library(reshape2)

met_cells <- read_excel("Documents/Massion_lab/Others/Dalton/data/Metabolon_cells_MF121120.xlsx")

# Super pathway names
k <- which(!is.na(met_cells$Super_Pathway))
sp <- met_cells$Super_Pathway[k]
for(i in 1:length(k)){
    if(i<length(k)){
        met_cells$Super_Pathway[k[i]:(k[i+1]-1)] <- sp[i]
    } else {
        met_cells$Super_Pathway[k[i]:nrow(met_cells)] <- sp[i]
    }
}

# Sub pathway names
k2 <- which(!is.na(met_cells$Sub_Pathway))
subp <- met_cells$Sub_Pathway[k2]
for(i in 1:length(k2)){
    if(i<length(k2)){
        met_cells$Sub_Pathway[k2[i]:(k2[i+1]-1)] <- subp[i]
    } else {
        met_cells$Sub_Pathway[k2[i]:nrow(met_cells)] <- subp[i]
    }
}



UD <- function(x){ 
    res <- ifelse(x < 0, 'DN', 'UP')
    res 
}

# Compute normalization and state labels
met_cells <- met_cells %>%
    tibble(.) %>%
    mutate(
        A549_L2R = log2(A549_KO / A549_WT),
        A549_state = UD(A549_L2R),
        BEAS2B_L2R = log2(BEAS2B_OE / BEAS2B_WT),
        BEAS2B_state = UD(BEAS2B_L2R),
    ) %>%
    relocate(A549_KO, .after=A549_WT) %>%
    relocate(A549_L2R, .after=A549_KO) %>%
    relocate(A549_state, .after=A549_L2R) %>%
    relocate(BEAS2B_OE, .after=BEAS2B_WT) %>%
    relocate(BEAS2B_L2R, .after=BEAS2B_OE) %>%
    relocate(BEAS2B_state, .after=BEAS2B_L2R) 


# Create summary data frame (mean subpathway values)
met_pth<- met_cells%>%
    group_by(Sub_Pathway)%>% 
    summarise(A549_WT=mean(A549_WT), A549_KO=mean(A549_KO), 
        BEAS2B_WT=mean(BEAS2B_WT), BEAS2B_OE=mean(BEAS2B_OE)) %>%
    mutate(Super_Pathway = met_cells$Super_Pathway[match(x$Sub_Pathway, met_cells$Sub_Pathway)]) %>%
    relocate(Super_Pathway, .before=Sub_Pathway) %>%
    mutate(
        A549_L2R = log2(A549_KO / A549_WT),
        A549_state = UD(A549_L2R),
        BEAS2B_L2R = log2(BEAS2B_OE / BEAS2B_WT),
        BEAS2B_state = UD(BEAS2B_L2R),
    ) %>%
    relocate(A549_KO, .after=A549_WT) %>%
    relocate(A549_L2R, .after=A549_KO) %>%
    relocate(A549_state, .after=A549_L2R) %>%
    relocate(BEAS2B_OE, .after=BEAS2B_WT) %>%
    relocate(BEAS2B_L2R, .after=BEAS2B_OE) %>%
    relocate(BEAS2B_state, .after=BEAS2B_L2R) 

save(met_cells, met_pth, file= 'Documents/Massion_lab/Others/Dalton/data/compdata2.RData')
