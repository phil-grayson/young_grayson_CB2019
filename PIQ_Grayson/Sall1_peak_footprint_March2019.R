library(tidyverse)
setwd("~/Desktop")

# NEVER try to load master.  It can't be done.  Use grep to grab what you want - even that takes 10 minutes to collect a scaffold.
#all.significant.calls <- read.csv("calls.Master.csv",sep=",",header=T) %>% tbl_df() %>% select(-X)
chick.chrom <- read_csv("chick_SALL1_chrom.csv",col_names = F) %>% tbl_df() %>% select(-X1) 
emu.chrom <-  read_csv("emu_SALL1_chrom.csv",col_names = F) %>% tbl_df() %>% select(-X1)

chick.chrom.fixed <- chick.chrom %>% separate(X7, c("purity","metadata"), sep=" results_") %>% 
  separate(metadata, c("library","pwmMeta","fileMeta"), sep="/") %>%
  separate(library, c("species","stage","tissue","Verts")) %>% select(-stage) %>%
  separate(fileMeta, c("pwmID","pwmTF.raw"), sep=" ") %>%
  separate(pwmTF.raw, c("pwmTF","toss"), sep="-calls") %>% select(-toss) %>% separate(pwmTF, c("pwmTF","RC"))

emu.chrom.fixed <- emu.chrom %>% separate(X7, c("purity","metadata"), sep=" results_") %>% 
  separate(metadata, c("library","pwmMeta","fileMeta"), sep="/") %>%
  separate(library, c("species","stage","tissue","Verts")) %>% select(-stage) %>%
  separate(fileMeta, c("pwmID","pwmTF.raw"), sep=" ") %>%
  separate(pwmTF.raw, c("pwmTF","toss"), sep="-calls") %>% select(-toss) %>% separate(pwmTF, c("pwmTF","RC"))


# first interested in significant footprints in the SALL1 peak region
# in chick: NC_006098.3:5332702-5333369
# in emu:  emu - scaffold_25:1990718-1991386

chick.SALL1.peak.footprints <- chick.chrom.fixed %>% filter(X3 >= 5332702, X3 <= 5333369) %>% filter(X6>0,purity>0.7)
write_tsv(chick.SALL1.peak.footprints,"Sall1_chick_footprint_scores.txt")
chick.SALL1.peak.footprints <- read_tsv("Sall1_chick_footprint_scores.txt")  %>% filter(X6>0,purity>0.7)
emu.SALL1.peak.footprints <- emu.chrom.fixed %>% filter(X3 >= 1990718, X3 <= 1991386) %>% filter(X6>0,purity>0.7)
emu.SALL1.peak.footprints <- read_tsv("Sall1_emu_footprint_scores.txt")  %>% filter(X6>0,purity>0.7)
write_tsv(emu.SALL1.peak.footprints,"Sall1_emu_footprint_scores.txt")
#emu.SALL1.peak.footprints %>% separate(pwmTF, c("pwmTF","RC"))

table(chick.SALL1.peak.footprints$tissue)
table(emu.SALL1.peak.footprints$tissue)

#look at chick Footprints
chick_SALL1.count.FL <- chick.SALL1.peak.footprints %>%filter(tissue == "FL") %>% dplyr::count(pwmTF) %>% dplyr::rename(chick.FL=n)
chick_SALL1.count.flank <- chick.SALL1.peak.footprints %>%filter(tissue == "flank") %>% dplyr::count(pwmTF)%>% dplyr::rename(chick.flank=n)
chick_SALL1.count.HL <- chick.SALL1.peak.footprints %>%filter(tissue == "HL") %>% dplyr::count(pwmTF)%>% dplyr::rename(chick.HL=n)

chick_SALL1.temp <- full_join(chick_SALL1.count.FL,chick_SALL1.count.flank)
chick_SALL1.count <- full_join(chick_SALL1.count.HL,chick_SALL1.temp)
#chick_SALL1.noflank <- chick_SALL1.count %>% filter(chick.FL == T, chick.HL == T, is.na(chick.flank))

# look at emu Footprints
emu_SALL1.count.FL <- emu.SALL1.peak.footprints %>%filter(tissue == "FL") %>% dplyr::count(pwmTF)%>% dplyr::rename(emu.FL=n)
emu_SALL1.count.flank <- emu.SALL1.peak.footprints %>%filter(tissue == "flank") %>% dplyr::count(pwmTF)%>% dplyr::rename(emu.flank=n)
emu_SALL1.count.HL <- emu.SALL1.peak.footprints %>%filter(tissue == "HL") %>% dplyr::count(pwmTF)%>% dplyr::rename(emu.HL=n)

emu_SALL1.temp <- full_join(emu_SALL1.count.FL,emu_SALL1.count.flank)
emu_SALL1.count <- full_join(emu_SALL1.count.HL,emu_SALL1.temp)
#emu_SALL1.onlyHL <- emu_SALL1.count %>% filter(emu.HL == T, is.na(emu.flank), is.na(emu.FL))

# compare!
SALL1_peak_together <- full_join(chick_SALL1.count,emu_SALL1.count)%>% replace(is.na(.),0)
SALL1_peak_cand <- SALL1_peak_together %>% filter(emu.HL >0, emu.FL == 0, emu.flank == 0, chick.FL > 0, chick.HL >0, chick.flank ==0)

setwd("~/Dropbox/Doctorate/atac/hh18/enhancer_screen/footPrint_SALL1peak/")
write_tsv(SALL1_peak_cand,"Sall1_peak_candidates.txt")
write_tsv(SALL1_peak_together,"Sall1_peak_all_footprint_count.txt")


