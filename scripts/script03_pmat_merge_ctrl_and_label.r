rm(list= ls())

library(tidyverse)

WT_ctrl_lst_rep1 <- list.files(pattern = "WT_Ctrl.*rep1.*pmat")

WT_ctrl_df_rep1 <- NULL
for(lst in WT_ctrl_lst_rep1){
    data <- read.table(lst,header=T)
    WT_ctrl_df_rep1 <- rbind(WT_ctrl_df_rep1,data) 
}

WT_ctrl_lst_rep2 <- list.files(pattern = "WT_Ctrl.*rep2.*pmat")
WT_ctrl_df_rep2 <- NULL
for(lst in WT_ctrl_lst_rep2){
    data <- read.table(lst,header=T)
    WT_ctrl_df_rep2 <- rbind(WT_ctrl_df_rep2,data) 
}
WT_label_lst_rep1 <- list.files(pattern = "WT_Label.*rep1.*pmat")
WT_label_df_rep1 <- NULL
for(lst in WT_label_lst_rep1){
    data <- read.table(lst,header=T)
    WT_label_df_rep1 <- rbind(WT_label_df_rep1,data) 
}
WT_label_lst_rep2 <- list.files(pattern = "WT_Label.*rep2.*pmat")
WT_label_df_rep2 <- NULL
for(lst in WT_label_lst_rep2){
    data <- read.table(lst,header=T)
    WT_label_df_rep2 <- rbind(WT_label_df_rep2,data) 
}


WT_label_df_rep2 <- WT_label_df_rep2 %>% mutate(id = str_c(chrom,pos,sep="_"))
WT_label_df_rep1 <- WT_label_df_rep1 %>% mutate(id = str_c(chrom,pos,sep="_"))

WT_ctrl_df_rep1 <- WT_ctrl_df_rep1 %>% mutate(id = str_c(chrom,pos,sep="_"))
WT_ctrl_df_rep2 <- WT_ctrl_df_rep2 %>% mutate(id = str_c(chrom,pos,sep="_"))

WT_rep1 <- dplyr::inner_join(WT_ctrl_df_rep1,WT_label_df_rep1,by="id",suffix=c("_ctrl","_label"))

WT_rep2 <- dplyr::inner_join(WT_ctrl_df_rep2,WT_label_df_rep2,by="id",suffix=c("_ctrl","_label"))

WT <- dplyr::inner_join(WT_rep1,WT_rep2,by="id",suffix=c("_R1","_R2"))

write.csv(WT,"WT_ctrl_label_pmat_merge.csv",quote=F,row.names=F)