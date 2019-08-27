#Analysis of HH10-HH13 LPM samples of the chick and emu

#Chick
ck_LPM_FL <- read.table("ck_LPM_FL.txt", header = T, sep = "\t")
ck_LPM_FL_meta <- read.table("ck_LPM_FL_meta.txt")

all(colnames(ck_LPM_FL) %in% rownames(ck_LPM_FL_meta))
ck_LPM_dds <- DESeqDataSetFromMatrix(countData = ck_LPM_FL, colData = ck_LPM_FL_meta, design = ~ tissue)
ck_LPM_dds_esf <- estimateSizeFactors(ck_LPM_dds)
rlog_ck_LPM <- rlog(ck_LPM_dds_esf, blind = TRUE)
plotPCA(rlog_ck_LPM, intgroup=c("tissue"))

#ggplot PCA
Ck_FL_pcaData <- plotPCA(rlog_ck_LPM, intgroup=c("tissue"), returnData=TRUE)
percentVar <- round(100 * attr(Ck_FL_pcaData, "percentVar"))
ggplot(Ck_FL_pcaData, aes(PC1, PC2, color=tissue)) +
  scale_color_manual(values = c("#0033FF", "#33FF33", "#FF66CC")) +
  theme_classic() +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme(axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25)))



#Emu
emu_LPM_FL <- read.table("emu_LPM_FL.txt", header = T, sep = "\t")
emu_LPM_FL_meta <- read.table("LPM_Emu_FL_meta.txt",header = T, row.names = 1)

all(colnames(emu_LPM_FL) %in% rownames(emu_LPM_FL_meta))
emu_LPM_FL_dds <- DESeqDataSetFromMatrix(countData = emu_LPM_FL, colData = emu_LPM_FL_meta, design = ~tissue)
emu_LPM_FL_dds_esf <- estimateSizeFactors(emu_LPM_FL_dds)
rlog_emu_FL_LPM <- rlog(emu_LPM_FL_dds_esf, blind = TRUE)
plotPCA(rlog_emu_FL_LPM, intgroup=c("tissue"))

#ggplot_Emu
emu_FL_pcaData <- plotPCA(rlog_emu_FL_LPM, intgroup=c("tissue"), returnData=TRUE)
percentVar <- round(100 * attr(emu_FL_pcaData, "percentVar"))
ggplot(emu_FL_pcaData, aes(PC1, PC2, color=tissue)) +
  scale_color_manual(values = c("#0033FF", "#33FF33", "#FF66CC")) +
  theme_classic() +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme(axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25)))


#DESeq2 analysis of HH18 chick and emu samples

#read the tables into R
ck_HH18 <- read.table("ck_featureCounts_4.3.19_Rmatrix.txt", header = T, sep = "\t")
emu_HH18 <- read.table("emu_featureCounts_Rmatrix.txt", header = T, sep = "\t")

#Add an X to the emu library identifier in the emu counts table so that 583FL becomes X583FL.  This will allow for subsequent steps.


# Load the rosetta
rosetta <- read.table("Rosetta_EmuChick_Sept2018_12338total.txt", header=T,sep = "\t")

#Associate rosetta info to the chick counts based on geneID
ck_rosetta <- full_join(ck_HH18,rosetta,by=c("Geneid"="geneID"))


#Associate the emu counts into the ck_rosetta table (needed to change geneID.y to gene because the new rosetta looks different)
emu_ck_Rosetta <- full_join(ck_rosetta,emu_HH18,by=c("gene"="Geneid"))

HH18_Counts_4.2019 <- emu_ck_Rosetta %>% filter(!is.na(geneName)) %>% dplyr::select(geneName,ck18FL1,ck18FL2,ck18FL3,ck18FL4,ck18HL1,ck18HL2,ck18HL3,ck18HL4,ck18K1,ck18k2,ck18K3,ck18K4,X583FL,X583HL,X583K,X595FL,X595HL,X595K,X707FL,X707HL,X707K,X710FL,X710HL) %>% tibble::column_to_rownames(var="geneName")


# Load the counts table
HH18_Counts_4.2019 <- read.table("emu_chick_count_April.txt", header=T, row.names=3, sep = "\t")

#for now, cut out flank dplyR
colnames(HH18_Counts_4.2019)

HH18_Counts_4.2019_FLvHL <- HH18_Counts_4.2019 %>% select(-X583K,-X595K,-X707K,-ck18K1,-ck18k2,-ck18K3,-ck18K4)

# Load the metadata
HH18_Counts_4.2019_FLvHL_meta <- read.table("HH18_4.2019_FLvHL_meta.txt", header = T, row.names = 1)

#remove gene identifier columns
colnames(HH18_Counts_4.2019_FLvHL)
HH18_Counts_4.2019_FLvHL_DESeq2 <- HH18_Counts_4.2019_FLvHL %>% select(-droNov,-geneName)



library(DESeq2)


#create the dds with species and species:tissue interactions

HH18_FvH_dds_4.2019 <- DESeqDataSetFromMatrix(countData = HH18_Counts_4.2019_FLvHL_DESeq2, colData = HH18_Counts_4.2019_FLvHL_meta, design =  ~ species + species:tissue)

#PCA
pcaData <- plotPCA(rlog.dds.nonNormal, intgroup=c("Tissue","Species","Stage"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Species:Stage, shape=Tissue)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()


rlog_HH18_FvH_dds_4.2019 <- rlog(HH18_FvH_dds_4.2019, blind = TRUE)

plotPCA(rlog_HH18_FvH_dds_4.2019, intgroup=c("tissue","species","tissue_species"))

PCA_data <- plotPCA(rlog_HH18_FvH_dds_4.2019, intgroup=c("tissue","species","tissue_species") ,returnData=T)
percentVar <- round(100 * attr(PCA_data, "percentVar"))

ggplot(PCA, aes(PC1,PC2, color=tissue_species)) + geom_point(size=3) +
  geom_text(aes(label=name)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

#run DESeq
HH18_FLvHL_DESeq_4.2019 <- DESeq(HH18_FvH_dds_4.2019)

#Results for Species:Tissue interaction, but they are confusing since the log2FoldChange could mean different things 
resMF.Species.Tissue_4.2019 <- results(HH18_FLvHL_DESeq_4.2019,contrast=list("specieschick.tissueHL", "speciesemu.tissueHL"))
#Ck_ FLvHL
resMF.ChickFLHL_4.2019 <- results(HH18_FLvHL_DESeq_4.2019,name="specieschick.tissueHL")
resMF.EmuFLHL_4.2019 <- results(HH18_FLvHL_DESeq_4.2019,name="speciesemu.tissueHL")

#Here I'm creating tbl dfs for each result set and renaming some columns, then merging together into a final tbl df 
df.resMF.Species.Tissue_4.2019 <- as.data.frame(resMF.Species.Tissue_4.2019) %>% tbl_df() %>% mutate(geneName=rownames(resMF.Species.Tissue_4.2019))


# the log2FoldChange and padj columns coming from the intra-species comparison will have a c. and an e. added to the start for ease of understanding  Note: need to add "dplyr::" in front of rename because DESeq has a rename function also. stupid

df.resMF.ChickFLHL_4.2019 <- as.data.frame(resMF.ChickFLHL_4.2019 ) %>% tbl_df() %>% mutate(geneName=rownames(resMF.ChickFLHL_4.2019)) %>% dplyr::rename(c.log2FoldChange=log2FoldChange, c.padj=padj)


df.resMF.EmuFLHL_4.2019 <- as.data.frame(resMF.EmuFLHL_4.2019 ) %>% tbl_df() %>% mutate(geneName=rownames(resMF.EmuFLHL_4.2019)) %>% dplyr::rename(e.log2FoldChange=log2FoldChange, e.padj=padj)

#Join chicken and emu FLHL dfs first, then the Species:Tissue
resMF.FLHL_4.2019 <- full_join(df.resMF.ChickFLHL_4.2019,df.resMF.EmuFLHL_4.2019,by="geneName") %>% dplyr::select(geneName, c.log2FoldChange, c.padj, e.log2FoldChange,e.padj)
resMF.all4.2019 <- full_join(df.resMF.Species.Tissue_4.2019,resMF.FLHL_4.2019, by="geneName")


#Adding the genenames to the resMF
resMF.all.genenames_4.2019 <- HH18_Counts_4.2019_FLvHL %>% tibble:: rownames_to_column() %>% dplyr::rename(name=geneName,geneName=rowname) %>% select(geneName,droNov,name)

#joining the geneame table to the resMF
HH18_FvH_resMF.all_4.2019 <- full_join(resMF.all.genenames_4.2019,resMF.all4.2019)

write.table(HH18_FvH_resMF.all_4.2019, file = "/Users/johnyoung/Dropbox/John_Mac/Documents/Tabin_Lab/Emu/Emu_Sequences/ForEmu_RNAseq/re-analysis_2.1.18/emu_chick_HH18/Data/HH18_FvH_resMF.all_4.2019.txt", sep="\t", quote=F, col.names=NA)

#Adding the thresholds and cutoffs

padj.cutoff <- 0.05
lfc <- 2

threshold <-HH18_FvH_resMF.all_4.2019$e.padj < padj.cutoff & HH18_FvH_resMF.all_4.2019$c.padj > padj.cutoff & abs(HH18_FvH_resMF.all_4.2019$e.log2FoldChange) > lfc & abs(HH18_FvH_resMF.all_4.2019$c.log2FoldChange) < lfc

test <- HH18_FvH_resMF.all_4.2019


test$unique <-HH18_FvH_resMF.all_4.2019$e.padj < padj.cutoff & HH18_FvH_resMF.all_4.2019$c.padj > padj.cutoff & abs(HH18_FvH_resMF.all_4.2019$e.log2FoldChange) > lfc & abs(HH18_FvH_resMF.all_4.2019$c.log2FoldChange) < lfc

test$shared <-HH18_FvH_resMF.all_4.2019$e.padj < padj.cutoff & HH18_FvH_resMF.all_4.2019$c.padj < padj.cutoff & abs(HH18_FvH_resMF.all_4.2019$e.log2FoldChange) > lfc & abs(HH18_FvH_resMF.all_4.2019$c.log2FoldChange) > lfc

test$not_sig <-HH18_FvH_resMF.all_4.2019$e.padj > padj.cutoff | HH18_FvH_resMF.all_4.2019$c.padj > padj.cutoff | abs(HH18_FvH_resMF.all_4.2019$e.log2FoldChange) < lfc | abs(HH18_FvH_resMF.all_4.2019$c.log2FoldChange) < lfc



test <- test %>% unite(threshold,unique,shared,not_sig,sep="_")

test$Threshold <- threshold

test$targetGenes <- test$name %in% c("TBX4","TBX5","SALL1","FGF10","NKX2-5","PTCH2","FGF13","FGF8","FGF19","FGF4","NDNF","FGF20","FGF18","FGF3","ISL1","HOXD13","HOXD12","PITX1","NOX4","CALC4","HOXC5","HOXC12","HOXC13","HOXA2")

#add more packages
library(ggplot2)
library(ggbeeswarm)
library(RColorBrewer)
library(pheatmap)
library(ggrepel)
library(calibrate)

# Removing the "NA" values
install.packages("naniar")
library(naniar)
test <- test %>% replace_with_na(replace = list(threshold = c ("FALSE_NA_NA","NA_FALSE_TRUE","NA_NA_NA")))
table(test$threshold)


#Building volcano plots


#Figure 3 Volcano Plot
ggplot(test) +
  geom_point(aes(x=e.log2FoldChange, y=-log10(e.padj), colour=threshold)) + geom_text_repel(size = 3,check_overlap=F,aes(x=e.log2FoldChange, y=-log10(e.padj),label=ifelse(targetGenes==TRUE,as.character(name),'')),hjust=2, vjust=2)  + scale_color_manual(values = c("#666666", "#FF0000", "#0033FF")) + 
  theme_classic() +
  xlim(c(-6,6)) +
  ylim(c(0,50)) +
  ggtitle('Emu Forelimb vs Hindlimb') +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  geom_hline(yintercept = 1.43, colour="#999999", linetype="dashed") + geom_vline(xintercept = 2, colour="#999999", linetype="dashed") + geom_vline(xintercept = -2, colour="#999999", linetype="dashed") + # Add p-adj value cutoff line
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25)))


#Supplemental full plot
ggplot(test) +
  geom_point(aes(x=e.log2FoldChange, y=-log10(e.padj), colour=threshold)) + geom_text_repel(size = 3,check_overlap=F,aes(x=e.log2FoldChange, y=-log10(e.padj),label=ifelse(targetGenes==TRUE,as.character(name),'')),hjust=2, vjust=2)  + scale_color_manual(values = c("#666666", "#FF0000", "#0033FF")) + 
  theme_classic() +
  xlim(c(-10,10)) +
  ylim(c(0,200)) +
  ggtitle('Emu Forelimb vs Hindlimb') +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  geom_hline(yintercept = 1.43, colour="#999999", linetype="dashed") + geom_vline(xintercept = 2, colour="#999999", linetype="dashed") + geom_vline(xintercept = -2, colour="#999999", linetype="dashed") + # Add p-adj value cutoff line
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25)))

# generating tables of diffentially expressed genes

emu_unique <- HH18_FvH_resMF.all_4.2019 %>% filter(HH18_FvH_resMF.all_4.2019$e.padj < padj.cutoff & HH18_FvH_resMF.all_4.2019$c.padj > padj.cutoff & abs(HH18_FvH_resMF.all_4.2019$e.log2FoldChange) > lfc & abs(HH18_FvH_resMF.all_4.2019$c.log2FoldChange) < lfc)


write.table(emu_unique, file = "/Users/johnyoung/Dropbox/John_Mac/Documents/Tabin_Lab/Emu/Emu_Sequences/ForEmu_RNAseq/re-analysis_2.1.18/emu_chick_HH18/Data/emu_unique.txt", sep="\t", quote=F, col.names=NA)










