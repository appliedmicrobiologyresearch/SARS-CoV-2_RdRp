##script including all the draft figures for the resistance paper:
setwd("/Volumes/AMBR$/Mari/COVID1/ResistancePaper/GSAID/")
outdir <- "/Volumes/AMBR$/Mari/COVID1/ResistancePaper/GSAID/COVGAP7d/"
dir.create(outdir)
library("ggplot2")
#library(ggExtra)
library(ggspatial)

all <- read.csv("All.GSAID.COVGAP.reports.tab", sep = "\t", stringsAsFactors = F, row.names = NULL, header=T)
head(all)
#trim the header in the middle of the file 
length(which(all$X.CHROM=="X.CHROM"))
all <- all[-which(all$X.CHROM=="X.CHROM"),]
dim(all)
meta <- read.csv("metadata_2020-08-06_09-41.tsv", sep="\t", stringsAsFactors = F, row.names = NULL, header=T)
head(meta)
meta$strain <- gsub("/","-",meta$strain)
head(meta)
#now let's add some metadata features to the table 
all$DATE <- meta$date[match(all$SAMPLEID,meta$strain)]
all$SEX <- meta$sex[match(all$SAMPLEID,meta$strain)]
all$AGE <- meta$age[match(all$SAMPLEID,meta$strain)]
all$PANGOLIN <- meta$pangolin_lineage[match(all$SAMPLEID,meta$strain)]
all$REGION <- meta$region_exposure[match(all$SAMPLEID,meta$strain)]
all$COUNTRY <- meta$country_exposure[match(all$SAMPLEID,meta$strain)]
all$DIVISION <- meta$division_exposure[match(all$SAMPLEID,meta$strain)]
all$LOCATION <- meta$location[match(all$SAMPLEID,meta$strain)]
head(all)
all

#exploration-> to be removed
table(all$Remdesivir_binding_affecting)
table(all$Mpro_pockets)
table(all$Chloroquine_binding_affecting)
table(all$FiveFU_binding_affecting)
table(all$Mpro)
table(all$RdRp_NSP12)
all[which(grepl("Yes",all$Remdesivir_binding_affecting)==T),]
#need to assess in general how many sequences show one or more mutation per protein/pocket.
basin <- unique(all$SAMPLEID)
infoRemde <- c()
for (id in basin){
  #id <- "Australia-VIC2093-2020"
  single <- all[which(all$SAMPLEID==id),]
  r <- table(single$Remdesivir_binding_affecting)
  dim(r)
  d <- c(id,names(r), unname(r), "p", "p")
  if (dim(r)>1){
    infoRemde <- rbind(infoRemde,d)
    }else{}
}
infoRemde ## no samples carrying both mutations (pockets or codons)
dim(infoRemde)
#now only the samples showing more than 1 mutation in either pockets
infoRemde1 <- c()
for (id in basin){
  #id <- "Australia-VIC2093-2020"
  single <- all[which(all$SAMPLEID==id),]
  r <- table(single$Remdesivir_binding_affecting)
  dim(r)
  d <- c(id,names(r), unname(r), "p", "p")
  if (dim(r)>2){
    infoRemde1 <- rbind(infoRemde1,d)
  }else{}
}
infoRemde1
all[which(all$SAMPLEID=="England-20124001702-2020"),]
table(all$Remdesivir_binding_affecting)
#now for lopi
table(all$Remdesivir_binding_affecting)
infoMpro <- c()
for (id in basin){
  #id <- "Australia-VIC2093-2020"
  single <- all[which(all$SAMPLEID==id),]
  r <- table(single$Mpro_pockets)
  dim(r)
  d <- c(id,names(r),unname(r), "p", "p", "p")
  if (dim(r)>1){
    infoMpro <- rbind(infoMpro,d)
  }else{}
}
infoMpro
#now for Mpro as well the actual samples that show at least two or more mutation in Mpro
infoMpro1 <- c()
for (id in basin){
  #id <- "Australia-VIC2093-2020"
  single <- all[which(all$SAMPLEID==id),]
  r <- table(single$Mpro_pockets)
  dim(r)
  d <- c(id,names(r),unname(r), "p", "p", "p", "p")
  if (dim(r)>2){
    infoMpro1 <- rbind(infoMpro1,d)
  }else{}
}
infoMpro1
dim(seqcount)
write.table(infoRemde, paste(outdir,"Info_remdesivir_mutationsall.tab", sep=""), sep="\t", col.names = F, row.names = F, quote = F)
write.table(infoRemde1, paste(outdir,"Info_remdesivir_mutationsmultiple.tab", sep=""), sep="\t", col.names = F, row.names = F, quote = F)
write.table(infoMpro, paste(outdir,"Info_Mpro_mutationsall.tab", sep=""), sep="\t", col.names = F, row.names = F, quote = F)
write.table(infoMpro1, paste(outdir,"Info_Mpro_mutationsmultiple.tab", sep=""), sep="\t", col.names = F, row.names = F, quote = F)
##now some numbers on the progresses get so far, 
sink(paste(outdir,("General_stats.txt"),sep=""))
print("Total mutation events")
dim(all)
print("Classification of mutation events")
table(all$ANNOTATION)
print("which in percentage is..")
format(table(all$ANNOTATION)/dim(all)[1]*100, scientific = F)
print("Total amount of genomes")
length(unique(all$SAMPLEID))
print("Total amount of lineages")
length(unique(all$PANGOLIN))
print("####RdRP")
print("Total amount of mutations affecting RdRp")
dim(all[which(all$RdRp_NSP12=="Yes"),])
print("..distributed across how many genomes:")
length(unique(all$SAMPLEID[which(all$RdRp_NSP12=="Yes")]))
print("Amount of genomes showing at least ONE mutation in remdesivir relevant pockets:")
dim(infoRemde)[1]
print("To know more about which sample and other specs, please consult Info_remdesivir_mutationsall.tab ")
print("Amount of genomes showing at least TWO mutation in remdesivir relevant pockets:")
dim(infoRemde1)[1]
print("To know more about which sample and other specs, please consult Info_remdesivir_mutationsmultiple.tab ")
print("###Mpro")
print("Total amount of mutations affecting Mpro")
dim(all[which(all$Mpro=="Yes"),])
print("..distributed across how many genomes:")
length(unique(all$SAMPLEID[which(all$Mpro=="Yes")]))
print("Amount of genomes showing at least ONE mutation in Mpro relevant pockets:")
dim(infoMpro)
print("To know more about which sample and other specs, please consult Info_Mpro_mutationsall.tab ")
print("Amount of genomes showing at least TWO mutation in Mpro relevant pockets:")
dim(infoMpro1)
print("To know more about which sample and other specs, please consult Info_Mpro_mutationsmultiple.tab ")
sink()

##Todo list:
#-spatial map with the world and a zoom on europe on Remdesivir, (for all mutation, in case Madlen finds something for selection) - Alr. there
#-temporal map of mutation occurrence over time (again for all), superimpose the actual cumulative count -Alr. there
#-a event frequency map of mutations, across all the genome. total and divided into syn/nonsyn -to do
#-same as above, but focussed on the pockets. divided syn/nonsym - Alr. there
#-diversity positive/negative overt time - to modify

#Now it's important that in seqcount we have a column stating not just positive or negative, but which mutation they have, 
#this time including multiple, if that is the case. useful to put out stats of co-occurrences in the paper.
##now a dataframe per sample, not per mutation, mutation later
seqcount <- data.frame(ID=unique(all$SAMPLEID))
seqcount$Country <- all$COUNTRY[match(seqcount$ID, all$SAMPLEID)]
seqcount$Sex <- all$SEX[match(seqcount$ID, all$SAMPLEID)]
seqcount$Age <- all$AGE[match(seqcount$ID, all$SAMPLEID)]
seqcount$Pangolin <- all$PANGOLIN[match(seqcount$ID, all$SAMPLEID)]
head(seqcount)

#First the proteins as a whole:
#Mpro
YM <- all$SAMPLEID[which(all$Mpro=="Yes")]
seqcount$Mpro <- "No"
seqcount$Mpro[match(YM,seqcount$ID)] <- "Yes"
#RdRP
YR <- all$SAMPLEID[which(all$RdRp_NSP12=="Yes")]
seqcount$RdRP <- "No"
seqcount$RdRP[match(YR,seqcount$ID)] <- "Yes"
#now the single pockets
#Remdesivir
#F480
YP480 <- all$SAMPLEID[which(all$Remdesivir_binding_affecting=="Yes_pocket_F480")]
YC480 <- all$SAMPLEID[which(all$Remdesivir_binding_affecting=="Yes_pocket_codon_F480")]
seqcount$Pocket_R_F480 <- "No"
seqcount$Pocket_R_F480[match(YP480,seqcount$ID)] <- "Pocket"
seqcount$Pocket_R_F480[match(YC480,seqcount$ID)] <- "Codon"
#V557
YP557 <- all$SAMPLEID[which(all$Remdesivir_binding_affecting=="Yes_pocket_V557")]
YC557 <- all$SAMPLEID[which(all$Remdesivir_binding_affecting=="Yes_pocket_codon_V557")]
seqcount$Pocket_R_V557 <- "No"
seqcount$Pocket_R_V557[match(YP557,seqcount$ID)] <- "Pocket"
seqcount$Pocket_R_V557[match(YC557,seqcount$ID)] <- "Codon"
#Mpro
#Lopi1
YPL1 <- all$SAMPLEID[which(all$Mpro_pockets=="Yes_Pocket_Lopinavir1")]
YC1L1 <- all$SAMPLEID[which(all$Mpro_pockets=="Yes_Pocket_Lopinavir1_Thr24")]
YC2L1 <- all$SAMPLEID[which(all$Mpro_pockets=="Yes_Pocket_Lopinavir1_Thr26")]
seqcount$Pocket_L_Lopi1 <- "No"
seqcount$Pocket_L_Lopi1[match(YPL1,seqcount$ID)] <- "Pocket"
seqcount$Pocket_L_Lopi1[match(YC1L1,seqcount$ID)] <- "CodonThr24"
seqcount$Pocket_L_Lopi1[match(YC2L1,seqcount$ID)] <- "CodonThr26"
#Lopi2
YPL2 <- all$SAMPLEID[which(all$Mpro_pockets=="Yes_Pocket_Lopinavir2")]
YCL2 <- all$SAMPLEID[which(all$Mpro_pockets=="Yes_Pocket_Lopinavir2_Asn119")]
seqcount$Pocket_L_Lopi2 <- "No"
seqcount$Pocket_L_Lopi2[match(YPL2,seqcount$ID)] <- "Pocket"
seqcount$Pocket_L_Lopi2[match(YCL2,seqcount$ID)] <- "CodonAsn119"
#GC373
YPGC <- all$SAMPLEID[which(all$Mpro_pockets=="Yes_PocketGC373")]
YC1GC <- all$SAMPLEID[which(all$Mpro_pockets=="Yes_PocketGC373_Cys145")]
YC2GC <- all$SAMPLEID[which(all$Mpro_pockets=="Yes_PocketGC373_Glu166")]
YC3GC <- all$SAMPLEID[which(all$Mpro_pockets=="Yes_PocketGC373_Gly143")]
YC4GC <- all$SAMPLEID[which(all$Mpro_pockets=="Yes_PocketGC373_His163")]
YC5GC <- all$SAMPLEID[which(all$Mpro_pockets=="Yes_PocketGC373_Ser144")]

seqcount$Pocket_L_GC373 <- "No"
seqcount$Pocket_L_GC373[match(YPGC,seqcount$ID)] <- "Pocket"
seqcount$Pocket_L_GC373[match(YC1GC,seqcount$ID)] <- "Codon_Cys145"
seqcount$Pocket_L_GC373[match(YC2GC,seqcount$ID)] <- "Codon_Glu166"
seqcount$Pocket_L_GC373[match(YC1GC,seqcount$ID)] <- "Codon_Gly143"
seqcount$Pocket_L_GC373[match(YC2GC,seqcount$ID)] <- "Codon_His163"
seqcount$Pocket_L_GC373[match(YC1GC,seqcount$ID)] <- "Codon_Ser144"

#Virulence
YC1V <- all$SAMPLEID[which(all$Mpro_pockets=="Yes_Virulence_Ala285")]
seqcount$Pocket_LV_Ala285 <- "No"
seqcount$Pocket_LV_Ala285[match(YC1V,seqcount$ID)] <- "Codon_Ala285"

#5FU
YP615 <- all$SAMPLEID[which(all$FiveFU_binding_affecting=="Yes_pocket_M615")]
YC615 <- all$SAMPLEID[which(all$FiveFU_binding_affecting=="Yes_pocket_codon_M615")]
seqcount$Pocket_5FU_M615 <- "No"
seqcount$Pocket_5FU_M615[match(YP615,seqcount$ID)] <- "Pocket"
seqcount$Pocket_5FU_M615[match(YC615,seqcount$ID)] <- "Codon"

head(seqcount)
remde <- seqcount[which(seqcount$Pocket_R_F480 != "No" | seqcount$Pocket_R_V557 != "No"),]
lopilike <- seqcount[which(seqcount$Pocket_L_Lopi1 != "No" | seqcount$Pocket_L_Lopi2 != "No" | 
                         seqcount$Pocket_L_GC373 != "No" | seqcount$Pocket_LV_Ala285 != "No"),]
remdemut <- all[which(all$SAMPLEID %in% remde$ID ==T),]
length(unique(remdemut$SAMPLEID))
remdemutonly <- remdemut[which(remdemut$Remdesivir_binding_affecting != "No"),]
length(unique(remdemutonly$SAMPLEID))
lopimut <- all[which(all$SAMPLEID %in% lopilike$ID ==T),]
length(unique(lopimut$SAMPLEID))
lopimutonly <- lopimut[which(lopimut$Mpro_pockets != "No"),]
length(unique(lopimutonly$SAMPLEID))

#again some stats
sink(paste(outdir,"Particular_stats.tab"))
print("###Remdesivir")
print("Genomes presenting at least one mutation in the pocket of 480 or 557")
dim(remde)[1]
print("Of which the following only 480 pocket")
dim(remde[which(remde$Pocket_R_F480 != "No"),])[1]
remde[which(remde$Pocket_R_F480 != "No"),]
print("..of which the following impacts the exact 480 codon")
dim(remde[which(remde$Pocket_R_F480 == "Codon"),])[1]
remde[which(remde$Pocket_R_F480 == "Codon"),]
print("Here instead we have the mutations impacting the 557 pocket")
dim(remde[which(remde$Pocket_R_V557 != "No"),])[1]
remde[which(remde$Pocket_R_V557 != "No"),]
print("..of which the following impacts the exact 557 codon")
dim(remde[which(remde$Pocket_R_V557 == "Codon"),])[1]
remde[which(remde$Pocket_R_V557 == "Codon"),]
print("Pangolin lineage prevalence of genome presenting the mutation in either pockets")
table(remde$Pangolin)
print("Prevalence of mutations in samples carrying at least one mutation in either two pockets:")
table(remdemut$ANNOTATION)
print("..in percentage on the total")
format(table(remdemut$ANNOTATION)/dim(all)[1]*100, scientific = F)
print("..in percentage of the group..")
format(table(remdemut$ANNOTATION)/dim(remdemut)[1]*100, scientific = F)
print("Prevalence of mutations in either two pockets:")
table(remdemutonly$ANNOTATION)
print("..in percentage on the total")
format(table(remdemutonly$ANNOTATION)/dim(all)[1]*100, scientific = F)
print("..in percentage of the group..")
format(table(remdemutonly$ANNOTATION)/dim(remdemutonly)[1]*100, scientific = F)


print("###LopiLike")
print("Genomes presenting at least one mutation in the pocket lopi1, lopi2, gc373, virulence")
dim(lopilike)[1]
print("Of which the following on pocket1")
dim(lopilike[which(lopilike$Pocket_L_Lopi1 != "No"),])[1]
lopilike[which(lopilike$Pocket_L_Lopi1 != "No"),]
print("..of which the following impacts the exact codons")
dim(lopilike[which(grepl("Codon",lopilike$Pocket_L_Lopi1)==T),])[1]
lopilike[which(grepl("Codon",lopilike$Pocket_L_Lopi1)==T),]
print("the following on pocket2")
dim(lopilike[which(lopilike$Pocket_L_Lopi2 != "No"),])[1]
lopilike[which(lopilike$Pocket_L_Lopi2 != "No"),]
print("..of which the following impacts the exact codons")
dim(lopilike[which(grepl("Codon",lopilike$Pocket_L_Lopi2)==T),])[1]
lopilike[which(grepl("Codon",lopilike$Pocket_L_Lopi2)==T),]
print("the following on pocket gc373")
dim(lopilike[which(lopilike$Pocket_L_GC373 != "No"),])[1]
lopilike[which(lopilike$Pocket_L_GC373 != "No"),]
print("..of which the following impacts the exact codons")
dim(lopilike[which(grepl("Codon",lopilike$Pocket_L_GC373)==T),])[1]
lopilike[which(grepl("Codon",lopilike$Pocket_L_GC373)==T),]
print("the following on pocket virulence")
dim(lopilike[which(lopilike$Pocket_LV_Ala285 != "No"),])[1]
lopilike[which(lopilike$Pocket_LV_Ala285 != "No"),]
print("..of which the following impacts the exact codons")
dim(lopilike[which(grepl("Codon",lopilike$Pocket_LV_Ala285)==T),])[1]
lopilike[which(grepl("Codon",lopilike$Pocket_LV_Ala285)==T),]

print("Pangolin lineage prevalence of genome presenting the mutation in all lopilike pockets")
table(lopilike$Pangolin)
print("Prevalence of mutations in samples carrying at least a mutation in one pocket:")
table(lopimut$ANNOTATION)
print("..in percentage on the total")
format(table(lopimut$ANNOTATION)/dim(all)[1]*100, scientific = F)
print("..in percentage of the group..")
format(table(lopimut$ANNOTATION)/dim(remdemut)[1]*100, scientific = F)
print("Prevalence of mutations inside the pockets:")
table(lopimutonly$ANNOTATION)
print("..in percentage on the total")
format(table(lopimutonly$ANNOTATION)/dim(all)[1]*100, scientific = F)
print("..in percentage of the group..")
format(table(lopimutonly$ANNOTATION)/dim(lopimutonly)[1]*100, scientific = F)
sink()
## now we subset from *mut only the actual mutation of interest (excluding the other regions)
#summary of the datasets now available for the figures
####MAINS:
##ALL: the original report, each entry is a MUTATION 
##SEQCOUNT: here each entry is a GENOME, and custom columns indicate if there are interesting mutation 
####DERIVED:
##REMDE & LOPILIKE: subsets of seqcount, entry: GENOME, 
#defined as containing samples with AT LEAST ONE mutation in one pocket
##REMDEMUT & LOPIMUT: subsets of all with the genomes of remde and lopilike, each entry is a MUTATION: 
#it includes ALL the mutation appearing in remde and lopilike df
##REMDEMUTONLY & LOPIMUTONLY: subsets of remdemut and lopimut, each entry is a MUTATION, only those affecting the pockets.
head(remde)


