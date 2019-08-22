library(tidyverse)

# how much LO do we have just from emu?
setwd("/Volumes/SeaBlue/rna-seq-2018/")

# load intersect files and parse them for emu, tinamou, rhea and NCBI ostrich
emu <- read.table(file="LO/droNov_intersect.bed",header=F, sep="\t") %>% tbl_df() %>% separate(V8, c("test1","gene"),sep="gene=") %>% separate(gene,c("gene","rest"),sep=";") %>% separate(V4,c("garbage","genetemp"),sep="ID=") %>% separate(genetemp,c("droNov","garbage2"),sep=";") %>% separate(test1, c("garbage","test1"),sep="GeneID:") %>% separate(test1, c("GeneID"))

LO.list <- list(emu)
names(LO.list)<-c("emu")  #super important to be able to name the output

for (i in 1:length(LO.list)){
  palaeo.unique <- LO.list[i][[1]] %>% select(droNov,gene,GeneID) %>% unique()
  palaeo.dup1 <- palaeo.unique %>% select(droNov) %>% filter(duplicated(.))
  palaeo.dup2 <- palaeo.unique %>% select(gene) %>% filter(duplicated(.))
  palaeo.dup3 <- palaeo.unique %>% select(GeneID) %>% filter(duplicated(.))
  palaeo.final <- palaeo.unique %>% filter(!droNov %in% palaeo.dup1$droNov, !gene %in% palaeo.dup2$gene, !GeneID %in% palaeo.dup3$GeneID)
  assign(paste0(names(LO.list[i]),  "_singleLO"),palaeo.final) 
}
# well that ain't bad - 7819.

# and HOG?
setwd("/Volumes/SeaBlue/rna-seq-2018/HOG")
emu_singleHOG_file <- read.table(file="emu_HOG_singles_newGFF.txt",header=T,sep="\t") %>% tbl_df() %>% dplyr::rename(droNov=gene,galGal=GeneID) %>% separate(droNov,c("droNov")) %>% select(geneName,droNov,galGal)

duplicate.hogs <- read.table("/Volumes/SeaBlue/rna-seq-2018/HOG/duplicated_hog_analysis/new.hog.candidates.midpoint.oldPruned.basal",header=T) %>% tbl_df() %>% select(gene,droNov,galGal) %>% dplyr::rename(geneName=gene)

all.hogs <- bind_rows(duplicate.hogs,emu_singleHOG_file)


rosetta_2sp <- full_join(all.hogs %>% mutate(from="HOG",galGal=as.character(galGal)),emu_singleLO %>% mutate(from="LO") %>% dplyr::rename(galGal=GeneID),by="galGal")
table(!is.na(rosetta_2sp$from.x),!is.na(rosetta_2sp$from.y)) # we have 7392 shared between LO and HOGs   
#     FALSE TRUE
#FALSE     0  605
#TRUE   4830 7392

rosetta_2sp %>% filter(!is.na(droNov.x),!is.na(droNov.y),droNov.y!=droNov.x) # 8 mismatches of the 7392 that overlap.  not bad.
rosetta_2sp %>% filter(!is.na(droNov.x),!is.na(droNov.y),geneName!=gene) #this is due to different geneNames in Tim's gff versus the NCBI gff

# going to do the same thing as with the 5 species rosetta, which is take the HOG gene and forget the LO gene.

emu_singleLO <- emu_singleLO %>% filter("PDCD7"!=gene,"LOC415976"!=gene,"TTC9"!=gene,"MMP11"!=gene,"SLC14A2"!=gene,"CCDC94"!=gene,"CCDC85A"!=gene,"DNOV00015978"!=droNov) #can't filter out C11H16ORF7 for some reason, so used its droNov code

# we have a problem with slightly different gene names in the two versions of the GFF file, so we will use the ones from Tim's file (replacing those from the normal liftover) 
newGFF.cds <- read.table("/Volumes/SeaBlue/rna-seq-2018/HOG/duplicated_hog_analysis/gcf_galGal4_CDSgrep_3columns.txt",header = T) %>% tbl_df()

emu_singleLO.GFFnames <- left_join(emu_singleLO %>% select(-gene),newGFF.cds %>% mutate(GeneID=as.character(GeneID)) %>% select(-XP.NP), by="GeneID") %>% unique()

rosetta_2sp <- full_join(all.hogs %>% mutate(from="HOG",galGal=as.character(galGal)),emu_singleLO.GFFnames %>% mutate(from="LO") %>% dplyr::rename(galGal=GeneID),by="galGal")

rosetta_2sp %>% filter(!is.na(droNov.x),!is.na(droNov.y),droNov.y!=droNov.x) # that fixed it.
rosetta_2sp %>% filter(!is.na(droNov.x),!is.na(droNov.y),geneName!=gene) # perfect


rosetta_2sp <- bind_rows(all.hogs %>% select(geneName,droNov,galGal) %>% mutate(galGal=as.character(galGal)),emu_singleLO.GFFnames %>% dplyr::rename(galGal=GeneID,geneName=gene)) %>% unique() #12514
rosetta_2sp.NA.droNov <- rosetta_2sp %>% filter(is.na(geneName)) # these 165 genes are either not present in the Tim gff or do not have an annotated CDS 
# there were only 43 in the 5 species comparisson, but that's to be expected.

# cut our losses on those:
rosetta_2sp <- rosetta_2sp %>% filter(!is.na(geneName))

rosetta_2sp.dup.droNov <- rosetta_2sp %>% filter(duplicated(droNov)) # now there are just 17 duplicated droNov entries 
#all the ones that appeared as duplicates in the 5 species analysis are present here too

rosetta_2sp.dup.galGal <- rosetta_2sp %>% filter(duplicated(galGal)) #0

rosetta_2sp <- rosetta_2sp %>% filter(!droNov %in% rosetta_2sp.dup.droNov$droNov) %>% filter(!is.na(geneName)) #12315 left after everything

# not a lot of mismatches here, but some.

rosetta_2sp %>% filter(duplicated(geneName)) # 10 expected due to same name but different GeneID in different places of the genome.
rosetta_2sp %>% filter(duplicated(droNov))
rosetta_2sp %>% filter(duplicated(galGal))

#write.table(file="emu_chick_species_GeneID_rosetta_12315.txt",rosetta_2sp,col.names = T,row.names = F,quote = F)
