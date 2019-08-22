# This script contains information on the final run of candidate generation for the HH18 RNA/ATAC-seq

library(tidyverse)
library(DESeq2)

setwd("~/Desktop/chapter2/")
rosetta <- read.table(file="emu_chick_species_GeneID_rosetta_12315.txt",sep=" ",header = T) %>% tbl_df() %>% mutate(galGal = as.character(galGal))

# final rosetta is by geneID, but John's featurecounts were by geneNumber, so:
# grep "mRNA" /Volumes/SeaBlue/rna-seq-2018/chicken_files/GCF_000002315.3_Gallus_gallus-4.0_genomic.gff | grep -v "exon" - > /Volumes/SeaBlue/rna-seq-2018/chicken_files/GCF_000002315.3_Gallus_gallus-4.0_genomic_GREPmRNA.gff 

# gff <- read_tsv(file='/Volumes/SeaBlue/rna-seq-2018/chicken_files/galGal4_topLevel_mRNAgrep_exonVgrep.gff',col_names = F) %>% separate(X9,c("t1","t2"),sep="Parent=") %>% separate(t2, c("t3","t4"),sep=";") %>% separate(t4, c("t5","t6"),sep="GeneID:") %>% separate(t6,c("t7")) %>% select(t3,t7) %>% unique()

# #from chapter4 stuff (featureCounts full)
chick_ID_pull <- read_tsv("/Volumes/SeaBlue/rna-seq-2018/chicken_files/timChick_geneGrep.txt",col_names = F) %>% dplyr::select(X9) %>%
separate(X9,c("IDraw","GeneIDraw"),sep=";") %>% separate(IDraw,c("toss","ID"),sep="=") %>%
separate(GeneIDraw,c("tos","GeneIDraw2"),sep="GeneID:") %>% separate(GeneIDraw2,c("GeneID","ts")) %>% dplyr::select(-toss,-tos,-ts)


chick <- read_tsv(file="ck_featureCounts_4.3.19_Rmatrix.txt")
chick_FT <- full_join(chick,chick_ID_pull, by=c("Geneid"="ID")) 
emu <- read_tsv(file="emu_featureCounts_7.12.18.Rmatrix.txt")

temp1 <- full_join(rosetta,chick_FT,by=c("galGal"="GeneID")) %>% filter(!is.na(Geneid)) %>% dplyr::select(-Geneid)
full.features <- full_join(temp1,emu,by=c("droNov"="Geneid")) 

full.features <- full.features %>% filter(!is.na(geneName), !is.na(ck18FL1), !is.na(`707K`))

# rosetta.merge <- left_join(rosetta,gff, by=c("galGal" = "t7")) 
# rosetta.merge %>% filter(duplicated(galGal))

# start by identifying the geneNames for these duplicated IDs
dup.chick.geneName <- full.features %>% filter(duplicated(full.features$galGal)) %>% dplyr::select(geneName)
# and grab all rows containing these geneNames
dup.chick.full.df <- full.features %>% filter(geneName %in% dup.chick.geneName$geneName)
# make a simple vector of the geneNames
geneList <- as.vector(unique(dup.chick.geneName$geneName))

# and write a loop that returns the row for each gene that has the highest number of transcript counts across the chicken libraries.
b <- list() # initialize a list
for (gene in geneList){ # for every "split" gene
  a <- dup.chick.full.df %>% filter(geneName==gene) # grab the rows that contain that gene
  a <- a %>% mutate(sum.col = rowSums(a[,7:18])) # sum the counts for the chicken columns
  MAX = 0    # then we want to identify and return the row with the highest count, so start with max 0
  for (i in 1:nrow(a)) { # for every split copy of a given gene
    x <- as.integer(a[i,] %>% dplyr::select(sum.col)) #return the rowsum
    if (x > MAX) MAX <- x #if that rowsum is greater than 0, it becomes the new max
  }
  b[[gene]] <- a %>% filter(sum.col == MAX) #once all the rows have been checked, add the highest to the list
}
c <- bind_rows(b)
# now we want to delete all duplicated rows from the original dataframe
full.features <- full.features %>% filter(!geneName %in% dup.chick.full.df$geneName)
# and add back just the entry with the highest in the highest count
full.features.final <- bind_rows(full.features,c)

# look to see if there are any duplicates left
table(duplicated(full.features.final$galGal)) # nice - no more duplicated galGals
table(duplicated(full.features.final$geneName)) # there are a couple genes that have the same name however

write_tsv("emu_chick_count_April.txt",x=full.features.final %>% select(-sum.col))
HH18_Counts_4reps <- read.table("emu_chick_count_April.txt", header=T, row.names=3, sep = "\t")

#Load the counts table - Needed to add sep = "\t" for it to work
#HH18_Counts_4reps <- read.table("countTable_hh18_chick_emu_July2018.txt", header=T, row.names=3, sep = "\t")

#cut out flank libraries
HH18_Counts_4reps <- HH18_Counts_4reps %>% dplyr::select(ck18FL1,ck18FL2,ck18FL3,ck18FL4,ck18HL1,ck18HL2,ck18HL3,ck18HL4,X583FL,X595FL,X707FL,X710FL,X583HL,X595HL,X707HL,X710HL)

#Load the metadata - this has been edited to only include the Fl and Hl rows
HH18_Counts_4reps.meta <- read.table("HH18_meta_FL_HL_4reps_noFlank.txt", header = T, row.names = 1)
 
#read the counts table with searchable genes
#HH18_Counts_4reps.genes <- read.table("countTable_hh18_chick_emu_July2018.txt", header=T, sep = "\t")

#create the dds with species and species:tissue interactions
HH18_dds <- DESeqDataSetFromMatrix(countData = HH18_Counts_4reps, colData = HH18_Counts_4reps.meta, design =  ~ Species + Species:Tissue)

#run DESeq
HH18_dds <- DESeq(HH18_dds)

#Results for Species:Tissue interaction, but they are confusing since the log2FoldChange could mean different things (you'll see what I mean below)...
resMF.Species.Tissue <- results(HH18_dds,contrast=list("SpeciesChick.TissueHindlimb", "SpeciesEmu.TissueHindlimb"))

#So, if we also look at the direction of change within species between FL and HL, we can think about things the way we did a few months ago
resMF.ChickFLHL <- results(HH18_dds,name="SpeciesChick.TissueHindlimb")
resMF.EmuFLHL <- results(HH18_dds,name="SpeciesEmu.TissueHindlimb")

#Here I'm creating tbl dfs for each result set and renaming some columns, then merging together into a final tbl df 
df.resMF.Species.Tissue <- as.data.frame(resMF.Species.Tissue) %>% tbl_df() %>% mutate(geneName=rownames(resMF.Species.Tissue))
# the log2FoldChange and padj columns coming from the intra-species comparison will have a c. and an e. added to the start for ease of understanding
df.resMF.ChickFLHL <- as.data.frame(resMF.ChickFLHL ) %>% tbl_df() %>% mutate(geneName=rownames(resMF.ChickFLHL)) %>% dplyr::rename(c.padj=padj,c.log2FoldChange = log2FoldChange)
df.resMF.EmuFLHL <- as.data.frame(resMF.EmuFLHL ) %>% tbl_df() %>% mutate(geneName=rownames(resMF.EmuFLHL)) %>% dplyr::rename(e.log2FoldChange=log2FoldChange, e.padj=padj)

#Join chicken and emu FLHL dfs first, then the Species:Tissue
resMF.FLHL <- full_join(df.resMF.ChickFLHL,df.resMF.EmuFLHL,by="geneName") %>% dplyr::select(geneName, c.log2FoldChange, c.padj, e.log2FoldChange,e.padj)
resMF.all <- full_join(df.resMF.Species.Tissue,resMF.FLHL, by="geneName") %>% dplyr::rename(galGal=geneName)
resMF.all.fixed <- left_join(resMF.all,rosetta) %>% select(-droNov)

write.table(resMF.all, file = "resMF.all_April2019_4reps.txt", sep = "\t")

atac <- read_tsv("hh18_annotate_Nov25.txt") %>% tbl_df() %>% arrange(Chromosome,Start) %>% dplyr::select(-X45) %>%
  mutate(C.count.FL=rowSums(.[,c("C_1_FL","C_2_FL","C_3_FL")]>0)) %>%
  mutate(C.count.HL=rowSums(.[,c("C_1_HL","C_2_HL","C_3_HL")]>0)) %>%
  mutate(C.count.K=rowSums(.[,c("C_1_K","C_2_K","C_3_K")]>0)) %>%
  mutate(E.count.FL=rowSums(.[,c("E_1_FL","E_2_FL","E_3_FL")]>0)) %>%
  mutate(E.count.HL=rowSums(.[,c("E_1_HL","E_2_HL","E_3_HL")]>0)) %>%
  mutate(E.count.K=rowSums(.[,c("E_1_K","E_2_K","E_3_K")]>0)) %>%
  mutate(C.count.total=rowSums(.[,c("C_1_FL","C_2_FL","C_3_FL","C_1_HL","C_2_HL","C_3_HL","C_1_K","C_2_K","C_3_K")]>0)) %>%
  mutate(E.count.total=rowSums(.[,c("E_1_FL","E_2_FL","E_3_FL","E_1_HL","E_2_HL","E_3_HL","E_1_K","E_2_K","E_3_K")]>0)) %>%
  mutate(TotalActive=rowSums(.[,c("St16_whole_H3K27ac","St16_whole_H3K4me1","St21_limb_H3K27ac","St21_limb_H3K4me1","St21_whole_H3K27ac","St21_whole_H3K4me1","St32_limb_H3K27ac","St32_limb_H3K4me1","St32_whole_H3K27ac","St32_whole_H3K4me1")]> 0)) %>%
  mutate(LimbActive=rowSums(.[,c("St21_limb_H3K27ac","St21_limb_H3K4me1","St32_limb_H3K27ac","St32_limb_H3K4me1")]> 0)) %>%
  mutate(TotalRepressed=rowSums(.[,c("St16_whole_H3K27me3","St21_limb_H3K27me3","St21_whole_H3K27me3","St32_limb_H3K27me3","St32_whole_H3K27me3")]> 0)) %>% 
  rowwise() %>% 
  mutate(EverActive=max(St16_whole_H3K27ac,St16_whole_H3K4me1,St21_limb_H3K27ac,St21_limb_H3K4me1,St21_whole_H3K27ac,St21_whole_H3K4me1,St32_limb_H3K27ac,St32_limb_H3K4me1,St32_whole_H3K27ac,St32_whole_H3K4me1)) %>%
  mutate(EverRepressed=max(St16_whole_H3K27me3,St21_limb_H3K27me3,St21_whole_H3K27me3,St32_limb_H3K27me3,St32_whole_H3K27me3)) %>%
  ungroup()

closest <- read_tsv("final_closest_hh18.txt",col_names = F) %>% tbl_df() %>% dplyr::select(X1:X3,X7) %>% dplyr::rename(Chromosome=X1, Start=X2, Stop=X3, geneName=X7)
atac.gene <- full_join(atac,closest) %>% mutate(Size=Stop-Start) %>% dplyr::select(1:3,Size,geneName,everything())

# write.table(atac.gene, file = "Young_Grayson_2019_ATAC-seq_Matrix.txt", sep = "\t")

ATAC.rna <- left_join(atac.gene,resMF.all.fixed)

emuFLmissing.all <- ATAC.rna %>% 
  filter(C_strict_FL>0,C_strict_HL>0,C_strict_flank==0,E_strict_HL>0,E_strict_FL==0,E_strict_flank==0) %>% 
  filter(C.count.FL==3,C.count.HL==3,C.count.K==0,E.count.FL==0,E.count.HL==3,E.count.K==0) %>% mutate(c.log2FoldChange = round(c.log2FoldChange,5),c.padj=formatC(c.padj, format = "e", digits = 5),e.log2FoldChange = round(e.log2FoldChange,5),e.padj=formatC(e.padj, format = "e", digits = 5),CNEE=round(CNEE, digits = 5))

#write_tsv(emuFLmissing.all %>% dplyr::select(1:5,CNEE,EverActive,EverRepressed,66:69),"SupTable_emuFLmissing_all_hh18.txt",col_names=T)
#write_tsv(emuFLmissing.all %>% dplyr::select(1:5,CNEE,EverActive,EverRepressed),"SupTable_emuFLmissing_all_hh18_thin.txt",col_names=T)

emuFLmissing.no.interaction <- ATAC.rna %>% 
  filter(e.padj<0.1,c.padj>0.1) %>%
  filter(C_strict_FL>0,C_strict_HL>0,C_strict_flank==0,E_strict_HL>0,E_strict_FL==0,E_strict_flank==0) %>%
  filter(C.count.FL==3,C.count.HL==3,C.count.K==0,E.count.FL==0,E.count.HL==3,E.count.K==0) %>% mutate(c.log2FoldChange = round(c.log2FoldChange,5),c.padj=formatC(c.padj, format = "e", digits = 5),e.log2FoldChange = round(e.log2FoldChange,5),e.padj=formatC(e.padj, format = "e", digits = 5),CNEE=round(CNEE, digits = 5))

#write_tsv(emuFLmissing.no.interaction %>% dplyr::select(1:5,CNEE,EverActive,EverRepressed,66:69),"SupTable_emuFLmissing_all_hh18_DE.txt",col_names=T)
#write_tsv(emuFLmissing.no.interaction %>% dplyr::select(1:5,CNEE,EverActive,EverRepressed),"SupTable_emuFLmissing_all_hh18__DE_thin.txt",col_names=T)



emuFLmissing.interaction <- ATAC.rna %>% 
  filter(padj < 0.1, e.padj<0.1,c.padj>0.1) %>%
  filter(C_strict_FL>0,C_strict_HL>0,C_strict_flank==0,E_strict_HL>0,E_strict_FL==0,E_strict_flank==0) %>%
  filter(C.count.FL==3,C.count.HL==3,C.count.K==0,E.count.FL==0,E.count.HL==3,E.count.K==0)



# filtering the total peak dataset for numbers for the venn diagrams for figure

colnames(atac)

nrow(atac) # 128,411 peaks total 

# peaks that are strict and present across all three bio reps for at least one tissue - 119,638
nrow(atac %>%
        filter(E_strict_FL != 0 & E.count.FL==3 | E_strict_HL != 0 & E.count.HL==3 | E_strict_flank != 0 & E.count.K==3 | C_strict_FL != 0 & C.count.FL==3 | C_strict_HL != 0 & C.count.HL==3 | C_strict_flank != 0 & C.count.K==3))

# let's move forward with these as atac.clean  
atac.clean <- atac %>%
  filter(E_strict_FL != 0 & E.count.FL==3 | E_strict_HL != 0 & E.count.HL==3 | E_strict_flank != 0 & E.count.K==3 | C_strict_FL != 0 & C.count.FL==3 | C_strict_HL != 0 & C.count.HL==3 | C_strict_flank != 0 & C.count.K==3)

a <- atac.clean %>% 
  filter(E_strict_FL == 0, E_strict_HL == 0, E_strict_flank == 0) %>% 
  filter(E.count.FL==0,E.count.HL==0,E.count.K==0)
b <- atac.clean %>% 
  filter(C_strict_FL == 0, C_strict_HL == 0, C_strict_flank == 0) %>% 
  filter(C.count.FL==0,C.count.HL==0,C.count.K==0)
c <- atac.clean %>%
  filter(E_strict_FL != 0 | E_strict_HL != 0 | E_strict_flank != 0 | C_strict_FL != 0 | C_strict_HL != 0 | C_strict_flank != 0) %>%
  filter(C.count.total > 0, E.count.total > 0)

View(anti_join(atac.clean,a) %>% anti_join(.,b) %>% anti_join(.,c))
# one strange peak that has strict emu flank but not other emu tissues...i think we should drop that one
atac.clean <- atac.clean %>% filter(Start != "96939814") # that did it
# peaks that are strict and present across all three bio reps for at least one tissue - 119,638

nrow(atac.clean %>%
        filter(E_strict_FL != 0 & E.count.FL==3 | E_strict_HL != 0 & E.count.HL==3 | E_strict_flank != 0 & E.count.K==3 | C_strict_FL != 0 & C.count.FL==3 | C_strict_HL != 0 & C.count.HL==3 | C_strict_flank != 0 & C.count.K==3))




# peaks that are chicken specific - 37,118 
nrow(atac.clean %>% 
        filter(E_strict_FL == 0, E_strict_HL == 0, E_strict_flank == 0) %>% 
        filter(E.count.FL==0,E.count.HL==0,E.count.K==0))

# peaks that are emu specific - 23,814
nrow(atac.clean %>% 
        filter(C_strict_FL == 0, C_strict_HL == 0, C_strict_flank == 0) %>% 
        filter(C.count.FL==0,C.count.HL==0,C.count.K==0))

# peaks that are shared to at least some extent
atac.clean <- atac.clean %>% mutate(C.count.total = C.count.FL + C.count.HL + C.count.K) %>% 
  mutate(E.count.total = E.count.FL + E.count.HL + E.count.K)

nrow(atac.clean %>%
        filter(E_strict_FL != 0 | E_strict_HL != 0 | E_strict_flank != 0 | C_strict_FL != 0 | C_strict_HL != 0 | C_strict_flank != 0) %>%
        filter(C.count.total > 0, E.count.total >0))




# peaks that are shared in all tissues in both species - 5,563 
nrow(atac.clean %>% 
        filter(C_strict_FL != 0, C_strict_HL != 0, C_strict_flank != 0, E_strict_FL != 0, E_strict_HL != 0, E_strict_flank != 0) %>% 
        filter(C.count.FL==3,C.count.HL==3,C.count.K==3,E.count.FL==3,E.count.HL==3,E.count.K==3))

# peaks that are shared in all emu tissues (no chick) - 3,474
nrow(atac.clean %>% 
        filter(C_strict_FL == 0, C_strict_HL == 0, C_strict_flank == 0, E_strict_FL != 0, E_strict_HL != 0, E_strict_flank != 0) %>% 
        filter(C.count.FL==0,C.count.HL==0,C.count.K==0,E.count.FL==3,E.count.HL==3,E.count.K==3))

# peaks that are shared in all chick tissues (no emu) - 11,139
nrow(atac.clean %>% 
        filter(C_strict_FL != 0, C_strict_HL != 0, C_strict_flank != 0, E_strict_FL == 0, E_strict_HL == 0, E_strict_flank == 0) %>% 
        filter(C.count.FL==3,C.count.HL==3,C.count.K==3,E.count.FL==0,E.count.HL==0,E.count.K==0))

# peaks that are limb shared in chicken (not in flank or in emu) - 681
nrow(atac.clean %>% 
        filter(E_strict_FL == 0, E_strict_HL == 0, E_strict_flank == 0, C_strict_FL != 0, C_strict_HL != 0, C_strict_flank == 0) %>% 
        filter(E.count.FL==0,E.count.HL==0,E.count.K==0,C.count.FL==3,C.count.HL==3,C.count.K==0))

# peaks that are in forelimb and flank in chicken (no hl or any emu) - 282
nrow(atac.clean %>% 
        filter(E_strict_FL == 0, E_strict_HL == 0, E_strict_flank == 0, C_strict_FL != 0, C_strict_HL == 0, C_strict_flank != 0) %>% 
        filter(E.count.FL==0,E.count.HL==0,E.count.K==0,C.count.FL==3,C.count.HL==0,C.count.K==3))

# peaks that are in hindlimb and flank in chicken (no fl or any emu) - 155
nrow(atac.clean %>% 
        filter(E_strict_FL == 0, E_strict_HL == 0, E_strict_flank == 0, C_strict_FL == 0, C_strict_HL != 0, C_strict_flank != 0) %>% 
        filter(E.count.FL==0,E.count.HL==0,E.count.K==0,C.count.FL==0,C.count.HL==3,C.count.K==3))

# peaks that are limb shared in emu (not in flank or in chicken) - 132
nrow(atac.clean %>% 
        filter(C_strict_FL == 0, C_strict_HL == 0, C_strict_flank == 0, E_strict_FL != 0, E_strict_HL != 0, E_strict_flank == 0) %>% 
        filter(C.count.FL==0,C.count.HL==0,C.count.K==0,E.count.FL==3,E.count.HL==3,E.count.K==0))

# peaks that are in forelimb and flank in emu (no hl or any chick) - 418
nrow(atac.clean %>% 
        filter(C_strict_FL == 0, C_strict_HL == 0, C_strict_flank == 0, E_strict_FL != 0, E_strict_HL == 0, E_strict_flank != 0) %>% 
        filter(C.count.FL==0,C.count.HL==0,C.count.K==0,E.count.FL==3,E.count.HL==0,E.count.K==3))

# peaks that are in hindlimb and flank in emu (no fl or any chick) - 85
nrow(atac.clean %>% 
        filter(C_strict_FL == 0, C_strict_HL == 0, C_strict_flank == 0, E_strict_FL == 0, E_strict_HL != 0, E_strict_flank != 0) %>% 
        filter(C.count.FL==0,C.count.HL==0,C.count.K==0,E.count.FL==0,E.count.HL==3,E.count.K==3))

# peaks that are in non-proliferative tissues (both flanks, emu fl) - 4
nrow(atac.clean %>% 
        filter(C_strict_FL == 0, C_strict_HL == 0, C_strict_flank != 0, E_strict_FL != 0, E_strict_HL == 0, E_strict_flank != 0) %>% 
        filter(C.count.FL==0,C.count.HL==0,C.count.K==3,E.count.FL==3,E.count.HL==0,E.count.K==3))

# peaks that are in non-proliferative flank with emu hl control (both flanks, emu hl) - 0
nrow(atac.clean %>% 
        filter(C_strict_FL == 0, C_strict_HL == 0, C_strict_flank != 0, E_strict_FL == 0, E_strict_HL != 0, E_strict_flank != 0) %>% 
        filter(C.count.FL==0,C.count.HL==0,C.count.K==3,E.count.FL==0,E.count.HL==3,E.count.K==3))

# peaks that are in non-proliferative flank with chick fl control (both flanks, chick fl) - 0
nrow(atac.clean %>% 
        filter(C_strict_FL != 0, C_strict_HL == 0, C_strict_flank != 0, E_strict_FL == 0, E_strict_HL == 0, E_strict_flank != 0) %>% 
        filter(C.count.FL==3,C.count.HL==0,C.count.K==3,E.count.FL==0,E.count.HL==0,E.count.K==3))

# peaks that are in non-proliferative flank with chick hl control (both flanks, chick hl) - 0
nrow(atac.clean %>% 
        filter(C_strict_FL == 0, C_strict_HL != 0, C_strict_flank != 0, E_strict_FL == 0, E_strict_HL == 0, E_strict_flank != 0) %>% 
        filter(C.count.FL==0,C.count.HL==3,C.count.K==3,E.count.FL==0,E.count.HL==0,E.count.K==3))

# SOLO LIBRARIES

# peaks are only in chick FL - 632
nrow(atac.clean %>% 
        filter(C_strict_FL != 0, C_strict_HL == 0, C_strict_flank == 0, E_strict_FL == 0, E_strict_HL == 0, E_strict_flank == 0) %>% 
        filter(C.count.FL==3,C.count.HL==0,C.count.K==0,E.count.FL==0,E.count.HL==0,E.count.K==0))

# peaks are only in chick HL - 426
nrow(atac.clean %>% 
        filter(C_strict_FL == 0, C_strict_HL != 0, C_strict_flank == 0, E_strict_FL == 0, E_strict_HL == 0, E_strict_flank == 0) %>% 
        filter(C.count.FL==0,C.count.HL==3,C.count.K==0,E.count.FL==0,E.count.HL==0,E.count.K==0))

# peaks are only in chick flank - 995
nrow(atac.clean %>% 
        filter(C_strict_FL == 0, C_strict_HL == 0, C_strict_flank != 0, E_strict_FL == 0, E_strict_HL == 0, E_strict_flank == 0) %>% 
        filter(C.count.FL==0,C.count.HL==0,C.count.K==3,E.count.FL==0,E.count.HL==0,E.count.K==0))


# peaks are only in emu FL - 329
nrow(atac.clean %>% 
        filter(C_strict_FL == 0, C_strict_HL == 0, C_strict_flank == 0, E_strict_FL != 0, E_strict_HL == 0, E_strict_flank == 0) %>% 
        filter(C.count.FL==0,C.count.HL==0,C.count.K==0,E.count.FL==3,E.count.HL==0,E.count.K==0))

# peaks are only in emu HL - 2250
nrow(atac.clean %>% 
        filter(C_strict_FL == 0, C_strict_HL == 0, C_strict_flank == 0, E_strict_FL == 0, E_strict_HL != 0, E_strict_flank == 0) %>% 
        filter(C.count.FL==0,C.count.HL==0,C.count.K==0,E.count.FL==0,E.count.HL==3,E.count.K==0))

# peaks are only in emu flank- 174
nrow(atac.clean %>% 
        filter(C_strict_FL == 0, C_strict_HL == 0, C_strict_flank == 0, E_strict_FL == 0, E_strict_HL == 0, E_strict_flank != 0) %>% 
        filter(C.count.FL==0,C.count.HL==0,C.count.K==0,E.count.FL==0,E.count.HL==0,E.count.K==3))


# peaks that are flank specific in both species - 0 
nrow(atac.clean %>% 
        filter(C_strict_FL == 0, C_strict_HL == 0, C_strict_flank != 0, E_strict_FL == 0, E_strict_HL == 0, E_strict_flank != 0) %>% 
        filter(C.count.FL==0,C.count.HL==0,C.count.K==3,E.count.FL==0,E.count.HL==0,E.count.K==3))

# peaks that are in limbs (only) in both species - 3
nrow(atac.clean %>% 
        filter(C_strict_FL != 0, C_strict_HL != 0, C_strict_flank == 0, E_strict_FL != 0, E_strict_HL != 0, E_strict_flank == 0) %>% 
        filter(C.count.FL==3,C.count.HL==3,C.count.K==0,E.count.FL==3,E.count.HL==3,E.count.K==0))
View(atac %>% 
       filter(C_strict_FL != 0, C_strict_HL != 0, C_strict_flank == 0, E_strict_FL != 0, E_strict_HL != 0, E_strict_flank == 0) %>% 
       filter(C.count.FL==3,C.count.HL==3,C.count.K==0,E.count.FL==3,E.count.HL==3,E.count.K==0))

# peaks that are in proliferating limbs (only) in both species - 59
nrow(atac.clean %>% 
        filter(C_strict_FL != 0, C_strict_HL != 0, C_strict_flank == 0, E_strict_FL == 0, E_strict_HL != 0, E_strict_flank == 0) %>% 
        filter(C.count.FL==3,C.count.HL==3,C.count.K==0,E.count.FL==0,E.count.HL==3,E.count.K==0))


# peaks that are in reversed limbs (both chick, but only emu forelimb) in both species - 0
nrow(atac.clean %>% 
        filter(C_strict_FL != 0, C_strict_HL != 0, C_strict_flank == 0, E_strict_FL != 0, E_strict_HL == 0, E_strict_flank == 0) %>% 
        filter(C.count.FL==3,C.count.HL==3,C.count.K==0,E.count.FL==3,E.count.HL==0,E.count.K==0))





# 182 peaks that are shared in all limbs (no flank)
atac %>% filter(C_strict_FL != 0, C_strict_HL != 0, C_strict_flank == 0, E_strict_FL != 0, E_strict_HL != 0, E_strict_flank == 0) 
# 863 peaks that are shared in proliferating tissues
atac %>% filter(C_strict_FL != 0, C_strict_HL != 0, C_strict_flank == 0, E_strict_FL == 0, E_strict_HL != 0, E_strict_flank == 0) 

# what if we want non-messy clean (no bio reps with 1 or 2 counts)

# atac.clean.noMess <- atac.clean %>% filter(C.count.FL == 3 | C.count.FL == 0, 
#                                           C.count.HL == 3 | C.count.HL == 0,
#                                           C.count.K == 3 | C.count.K == 0,
#                                           E.count.FL == 3 | E.count.FL == 0, 
#                                           E.count.HL == 3 | E.count.HL == 0,
#                                           E.count.K == 3 | E.count.K == 0)

atac.clean.noMess <- atac.clean %>% filter(C_strict_FL > 0 & C.count.FL == 3 | C_strict_FL ==0 & C.count.FL == 0, 
                                           C_strict_HL > 0 & C.count.HL == 3 | C_strict_HL ==0 & C.count.HL == 0,
                                           C_strict_flank > 0 & C.count.K == 3 | C_strict_flank ==0 & C.count.K == 0,
                                           E_strict_FL > 0 & E.count.FL == 3 | E_strict_FL ==0 & E.count.FL == 0, 
                                           E_strict_HL > 0 & E.count.HL == 3 | E_strict_HL ==0 & E.count.HL == 0,
                                           E_strict_flank > 0 & E.count.K == 3 | E_strict_flank ==0 & E.count.K == 0)


#27,669 fall into this category

# peaks that are chicken specific - 14310 
nrow(atac.clean.noMess %>% 
        filter(E_strict_FL == 0, E_strict_HL == 0, E_strict_flank == 0) %>% 
        filter(E.count.FL==0,E.count.HL==0,E.count.K==0))


# peaks that are emu specific - 7162
nrow(atac.clean.noMess %>% 
        filter(C_strict_FL == 0, C_strict_HL == 0, C_strict_flank == 0) %>% 
        filter(C.count.FL==0,C.count.HL==0,C.count.K==0))

# peaks that are shared - 6197
nrow(atac.clean.noMess %>% filter(C_strict_FL != 0 & E_strict_FL !=0 |
                                     C_strict_FL != 0 & E_strict_HL !=0 |
                                     C_strict_FL != 0 & E_strict_flank !=0 |
                                     C_strict_HL != 0 & E_strict_FL !=0 |
                                     C_strict_HL != 0 & E_strict_HL !=0 |
                                     C_strict_HL != 0 & E_strict_flank !=0 |
                                     C_strict_flank != 0 & E_strict_FL !=0 |
                                     C_strict_flank != 0 & E_strict_HL !=0 |
                                     C_strict_flank != 0 & E_strict_flank !=0))

### need to carry out venn diagram calculations for both emu and chick separtately (for clean and no mess)

# atac.clean has 119,637 peaks  
# note that these are not STRICT - probably for supplements - only need some libraries
nrow(atac.clean %>% filter(C.count.total>0,E.count.total==0)) # chick only 37,118
nrow(atac.clean %>% filter(E.count.total>0,C.count.total==0)) # emu only 23,814
nrow(atac.clean %>% filter(E.count.total>0,C.count.total>0)) # shared 58,705

nrow(atac.clean %>% filter(C.count.FL > 0 | C.count.HL > 0 | C.count.K >0)) #95823 for all chick including shared
nrow(atac.clean %>% filter(E.count.FL > 0 | E.count.HL > 0 | E.count.K >0)) #82519 for all emu including shared
# 95823+82519-58705 (shared) = 119637

# chick "clean"
nrow(atac.clean %>% filter(C.count.FL > 0, C.count.HL == 0, C.count.K == 0)) # chick FL only 3791
nrow(atac.clean %>% filter(C.count.FL > 0, C.count.HL > 0, C.count.K == 0))# chick FL HL 6964
nrow(atac.clean %>% filter(C.count.FL > 0, C.count.HL == 0, C.count.K > 0)) # chick FL flank 5239
nrow(atac.clean %>% filter(C.count.FL == 0, C.count.HL > 0, C.count.K > 0))# chick HL flank 4026
nrow(atac.clean %>% filter(C.count.FL == 0, C.count.HL == 0, C.count.K > 0))# chick flank only 4424
nrow(atac.clean %>% filter(C.count.FL == 0, C.count.HL > 0, C.count.K == 0))# chick HL only 3420
nrow(atac.clean %>% filter(C.count.FL > 0, C.count.HL > 0, C.count.K > 0)) # chick all 67959

# emu "clean"
nrow(atac.clean %>% filter(E.count.FL > 0, E.count.HL == 0, E.count.K == 0)) # emu FL only 3716
nrow(atac.clean %>% filter(E.count.FL > 0, E.count.HL > 0, E.count.K == 0))# emu FL HL 8032
nrow(atac.clean %>% filter(E.count.FL > 0, E.count.HL == 0, E.count.K > 0)) # emu FL flank 4883
nrow(atac.clean %>% filter(E.count.FL == 0, E.count.HL > 0, E.count.K > 0))# emu HL flank 5368
nrow(atac.clean %>% filter(E.count.FL == 0, E.count.HL == 0, E.count.K > 0))# emu flank only 2482
nrow(atac.clean %>% filter(E.count.FL == 0, E.count.HL > 0, E.count.K == 0))# emu HL only 11578
nrow(atac.clean %>% filter(E.count.FL > 0, E.count.HL > 0, E.count.K > 0)) # emu all 46460


# chick "clean no mess"
nrow(atac.clean.noMess %>% filter(C.count.FL > 0, C.count.HL == 0, C.count.K == 0)) # chick FL only 657
nrow(atac.clean.noMess %>% filter(C.count.FL > 0, C.count.HL > 0, C.count.K == 0))# chick FL HL 765
nrow(atac.clean.noMess %>% filter(C.count.FL > 0, C.count.HL == 0, C.count.K > 0)) # chick FL flank 304
nrow(atac.clean.noMess %>% filter(C.count.FL == 0, C.count.HL > 0, C.count.K > 0))# chick HL flank 166
nrow(atac.clean.noMess %>% filter(C.count.FL == 0, C.count.HL == 0, C.count.K > 0))# chick flank only 1019
nrow(atac.clean.noMess %>% filter(C.count.FL == 0, C.count.HL > 0, C.count.K == 0))# chick HL only 461
nrow(atac.clean.noMess %>% filter(C.count.FL > 0, C.count.HL > 0, C.count.K > 0)) # chick all 17135

# emu "clean no mess"
nrow(atac.clean.noMess %>% filter(E.count.FL > 0, E.count.HL == 0, E.count.K == 0)) # emu FL only 351 - differs from 329 bc allows chick peaks
nrow(atac.clean.noMess %>% filter(E.count.FL > 0, E.count.HL > 0, E.count.K == 0))# emu FL HL 207
nrow(atac.clean.noMess %>% filter(E.count.FL > 0, E.count.HL == 0, E.count.K > 0)) # emu FL flank 460
nrow(atac.clean.noMess %>% filter(E.count.FL == 0, E.count.HL > 0, E.count.K > 0))# emu HL flank 109
nrow(atac.clean.noMess %>% filter(E.count.FL == 0, E.count.HL == 0, E.count.K > 0))# emu flank only 179
nrow(atac.clean.noMess %>% filter(E.count.FL == 0, E.count.HL > 0, E.count.K == 0))# emu HL only 2938
nrow(atac.clean.noMess %>% filter(E.count.FL > 0, E.count.HL > 0, E.count.K > 0)) # emu all 9115


#View(atac.clean %>% filter(C_strict_flank==0,C_strict_FL>0,C_strict_HL>0,E_strict_flank==0,E_strict_FL==0,E_strict_HL>0))
