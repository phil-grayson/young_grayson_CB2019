library(tidyverse)
library(ape)
library(phytools)

### ---------- LOAD AND PARSE NECESSARY FILES ---------- ###

setwd("/Volumes/SeaBlue/rna-seq-2018/HOG")
hog2.0 <- read.table(file = "new_hog_list.txt",sep="\t") %>% tbl_df() %>% filter(V4=="galGal" | V4=="droNov" | V4=="notPer" | V4=="strCam" | V4=="rheAme")

newGFF <- read.delim(file = "~/../../Volumes/SeaBlue/rna-seq-2018/chicken_files/grep_CDS.gff",sep="\t",header=F) %>% tbl_df()

newGFF.cds <- newGFF %>% filter(V3=="CDS") %>% select(V9) %>% separate(V9, c("Garbage","Rest"),sep="GeneID:") %>% select(-Garbage) %>% separate(Rest, c("tempID","Rest"),sep="Genbank:") %>% separate(tempID,c("GeneID","Garbage"))  %>% separate(Rest,c("XP.NP"),remove=F,sep="[,;.]") %>% separate(Rest,c("drop","tempName"),sep="gene=") %>% separate(tempName,c("gene","rest"),sep=";") %>% select(GeneID, XP.NP, gene) %>% distinct() %>% mutate(GeneID = as.factor(GeneID),XP.NP=as.factor(XP.NP),gene=as.factor(gene))

# write.table(file="gcf_galGal4_CDSgrep_3columns.txt", newGFF.cds, sep="\t",quote = FALSE,col.names = T,row.names = F) 

setwd("/Volumes/SeaBlue/rna-seq-2018/HOG/duplicated_hog_analysis/2500_trees/")

### ---------- BELOW ARE A NUMBER OF POSSIBLE TREE EDITORS TO GET THE MONOPHYLETIC GROUPS ---------- ###
# We have three origins that the trees can be called monophyletic from:
#   1 - original (which looks for groups of 5 using subtrees from the midpoint rooted or unrooted trees)
#   2 - pruned (which looks through the tree once the 5 species groups have been removed)
#   3 - basal (which looks at the basal-most split and checks if there is a group of 5 there)

### ---------- MIDPOINT WITH UPDATED PRUNED (SHOULD NOT NEED BASAL - BASAL ADDED AND PROVIDES NOTHING NEW) ---------- ###

files <- list.files(pattern = "*.txt", full.names=TRUE, recursive=FALSE)

multi_hog_parser <- function(x) {
  lst <- list()
  n=0
  for (t in x){
  tree <- read.tree(file=t)
  tree <- midpoint.root(tree)
  tiplist <- c()
  subs <-subtrees(tree)
  for (i in 1:length(subs)) {
    if (length(subs[[i]]$tip.label) == 5) {
      tips = c(subs[[i]]$tip.label)
      if (nrow(hog2.0 %>% filter(V2 %in% tips) %>% select(V4) %>% unique()) == 5){
        tiplist <- c(tiplist, subs[[i]]$tip.label)
        n <- n + 1
        tempdf <- subs[[i]]$tip.label %>% tbl_df() %>% rename(XP.NP = value) %>% 
          left_join(.,hog2.0 %>% rename(ID=V3),by=c("XP.NP" = "V2")) %>% select(ID,V4,V1) %>%
          rename(species = V4, HOG=V1) %>% spread(species,ID) %>% mutate(origin = "original")
        lst[[n]] <- tempdf
      }
    }
  }
  if (length(tiplist)>0 & length(tiplist) < length(tree$tip.label)){
    pruned <- drop.tip(tree, tiplist, trim.internal = TRUE, subtree = FALSE,
                       root.edge = 0, rooted = is.rooted(tree), collapse.singles = TRUE,interactive = FALSE)
    subs <-subtrees(pruned)
    for (i in 1:length(subs)) {
      if (length(subs[[i]]$tip.label) == 5) {
        tips = c(subs[[i]]$tip.label)
        if (nrow(hog2.0 %>% filter(V2 %in% tips) %>% select(V4) %>% unique()) == 5){
          tiplist <- c(tiplist, subs[[i]]$tip.label)
          n <- n + 1
          tempdf <- subs[[i]]$tip.label %>% tbl_df() %>% rename(XP.NP = value) %>% 
            left_join(.,hog2.0,by=c("XP.NP" = "V2")) %>% select(V3,V4,V1) %>%
            rename(species = V4, HOG=V1) %>% spread(species,V3) %>% mutate(origin = "pruned")
          lst[[n]] <-tempdf
        }
      }
    }
  }
  if(length(tiplist)>0 & length(tiplist) < length(tree$tip.label)){
    final.pruned <- drop.tip(tree, tiplist, trim.internal = TRUE, subtree = FALSE,
                             root.edge = 0, rooted = is.rooted(tree), collapse.singles = TRUE,interactive = FALSE)
    if (length(final.pruned$tip.label) > 5){
      drops.for.basal.isolation <- extract.clade(final.pruned, length(final.pruned$tip.label)+2 , root.edge = 0, collapse.singles = TRUE, interactive = FALSE)$tip.label
      basal.tree <- drop.tip(final.pruned, drops.for.basal.isolation, trim.internal = TRUE, subtree = FALSE, root.edge = 0, rooted = is.rooted(tree), collapse.singles = TRUE,interactive = FALSE)
      if (length(basal.tree$tip.label) == 5) {
        if (nrow(hog2.0 %>% filter(V2 %in% basal.tree$tip.label) %>% select(V4) %>% unique()) == 5){
          n <- n + 1
          tempdf <- subs[[i]]$tip.label %>% tbl_df() %>% rename(XP.NP = value) %>% 
            left_join(.,hog2.0,by=c("XP.NP" = "V2")) %>% select(V3,V4,V1) %>%
            rename(species = V4, HOG=V1) %>% spread(species,V3) %>% mutate(origin = "basal")
          lst[[n]] <-tempdf
            }
          }
  }
  }
  }
  return(bind_rows(lst) %>% left_join(.,newGFF.cds %>% select(-XP.NP) %>% unique(),by=c("galGal" = "GeneID")) %>% select(HOG,gene,droNov,galGal,notPer,rheAme,strCam,origin))
}

new.hog.candidates.midpoint.newpruned.basal <- multi_hog_parser(files)
system.time(new.hog.candidates.midpoint.newpruned.basal <- multi_hog_parser(files))
write.table(file="new.hog.candidates.midpoint.oldPruned.basal", new.hog.candidates.midpoint.newpruned.basal, sep="\t",quote = FALSE,col.names = T,row.names = F) 

