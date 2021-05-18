setwd("/Volumes/AMBR$/Mari/COVID1/ResistancePaper/GSAID/")
library("ggplot2")
#library(ggExtra)
library(ggspatial)

all <- read.csv("All.GSAID.reports.tab", sep = "\t", stringsAsFactors = F, row.names = NULL, header=T)
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
##now a dataframe per sample, not per mutation, mutation later
seqcount <- data.frame(ID=unique(all$SAMPLEID))
seqcount$Country <- all$COUNTRY[match(seqcount$ID, all$SAMPLEID)]
seqcount$Sex <- all$SEX[match(seqcount$ID, all$SAMPLEID)]
seqcount$Age <- all$AGE[match(seqcount$ID, all$SAMPLEID)]
seqcount$Pangolin <- all$PANGOLIN[match(seqcount$ID, all$SAMPLEID)]
head(seqcount)
length(seqcount$Country)
##now let's produce a dataframe containing a single entry per sequence, and label it whether this sequence contains any mutation for lopi, remde or 5FFU
#this gives an idea of how many patients and NOT mutations carry the resistance to LOP or REMDE or 5FFU. 
#attention: normally one patient carries max one mutation per drug, in case , f.e. we have one patient carrying 2 remdesivir mutation we will label it the same as the patient carrying 1 mutation for remde.
#we don't go for the loop, takes too long, here is a workaround
#three dataframe listing all the positive mutations
remde <- all[which(grepl("Yes",all$Remdesivir_binding_affecting)==TRUE & grepl("_M611", all$Remdesivir_binding_affecting)==FALSE),]
remdeNS <- remde[which(remde$ANNOTATION != "synonymous_variant"),]
lopi <- all[which(grepl("Yes",all$Lopinavir_binding_pocket)==TRUE),]
lopiNS <- lopi[which(lopi$ANNOTATION != "synonymous_variant"),]
fivefu <- all[which(grepl("_M611", all$Remdesivir_binding_affecting)==TRUE),]
fivefuNS <- fivefu[which(fivefu$ANNOTATION != "synonymous_variant"),]
dim(fivefu)  ### discrepancy in counts of positive between all and seqcount can be cause maybe a sigle sample harbours multiple mutations
table(fivefuNS$ANNOTATION)
##check annotation
remde[which(remde$PANGOLIN=="B.1.1.29"),]
remde[which(remde$POS==14877),]
?
########analyse now how the mutations are distributed in terms of missense synonymous and how many exactly affect the binding pockets.
#let's also take a chance to write some numbers down to be put in the manuscript
pock476 <- remde[which(remde$Remdesivir_binding_affecting=="Yes_F476"),]
pock553 <- remde[which(remde$Remdesivir_binding_affecting=="Yes_V553"),]
#solve this 3
un <- unique(remde$SAMPLEID)
length(which(duplicated(match(remde$SAMPLEID,un)))==TRUE)
remde$SAMPLEID[which(duplicated(match(remde$SAMPLEID,un)))]
#first some stats for all of them 
sink("Stats_remdesivir.txt")
print("total number of mutations occurring in the window area")
print(length(remde$SAMPLEID))
print("Total number of samples presenting at least one variant in one of the two binding pockets")
print(length(unique(remde$SAMPLEID)))
print("Number of samplesID that presents more than one mutation in the pockets, unique")
print(length(which(table(remde$SAMPLEID)>1)))
print("..being the following IDs..")
print(names(which(table(remde$SAMPLEID)>1)))
print("number of entries presenting more than one mutations, non unique")
print(length(remde$SAMPLEID[which(duplicated(match(remde$SAMPLEID,un)))]))
print("complete list of duplicated mutated: one entry means at least two mutation, 2 entry 3 mutatios and so on")
print(remde$SAMPLEID[which(duplicated(match(remde$SAMPLEID,un)))])
#now separate the pockets
print("numbers of different annotations across all remde subset")
print(table(remde$ANNOTATION))
print("numbers which AA is involved remde subset")
print(table(remde$AA_INVOLVED))
print("numbers of different annotations across pocket 480 subset")
print(table(pock476$ANNOTATION))
print("numbers which AA is involved 480 subset")
print(table(pock476$AA_INVOLVED))
print("numbers of different annotations across pocket 557 subset")
print(table(pock553$ANNOTATION))
print("numbers which AA is involved 557 subset")
table(pock553$AA_INVOLVED)
sink()
table(pock553$AA_INVOLVED)
library("ComplexHeatmap")
library(cowplot)
p <- ggplot(pock553, aes(x=POS,y=SAMPLEID, fill=ANNOTATION))+
  geom_tile()+
  scale_fill_manual(values = c("#68983e",
                               "#8255b9",
                               "#c36331",
                               "#bc4775"), name= "Variant annotation", labels=c("Non synonymous","Stop gained","Stop lost","Synonymous"))+
  theme(
    panel.background = element_blank(),
    panel.grid = element_line(size = 0.001, linetype = 'dashed',
                              colour = "ghostwhite"), 
    legend.key = element_rect (fill = "white"),
    legend.text = element_text (size = 20),
    plot.title = element_text(face = "bold", hjust = 1,vjust = 1, size = 30),
    legend.title = element_text (size = 20, face = "bold"),
    axis.text = element_text(size =20),
    axis.line = element_line(color= "black", size =0.6),
    axis.text.x = element_blank(),#text(size = 18, angle = 45,vjust= 1,hjust=1),
    axis.text.y = element_blank(),#text(size = 20),
    axis.title = element_blank(),#text(size =22, face = "bold"),
    legend.position="top"
  )
p
marginal1 <- qplot(pock553$POS)+
  scale_y_reverse()+
  xlab("Genome position")+
  ylab("Counts")+
  theme(
    panel.background = element_blank(),
    panel.grid = element_line(size = 0.001, linetype = 'dashed',
                              colour = "ghostwhite"), 
    legend.key = element_rect (fill = "white"),
    legend.text = element_text (size = 20),
    plot.title = element_text(face = "bold", hjust = 1,vjust = 1, size = 30),
    legend.title = element_text (size = 20, face = "bold"),
    axis.text = element_text(size =20),
    axis.line = element_line(color= "black", size =0.6),
    axis.text.x = element_text(size = 18, angle = 45,vjust= 1,hjust=1),
    axis.text.y = element_blank(),#text(size = 20),
    axis.title = element_blank(),#text(size =22, face = "bold"),
    legend.position="top"
  )
marginal1
head(pock476)
length(unique(pock476$COUNTRY))
marginal2 <- ggplot(pock553, aes( x=X.CHROM,y=SAMPLEID, fill=COUNTRY))+
  geom_tile()+
  scale_fill_manual(values=c("#ea317c",
                             "#5ebd31",
                             "#743dcd",
                             "#88a734",
                             "#d24adf",
                             "#509b3d",
                             "#9664f1",
                             "#bfa626",
                             "#567df1",
                             "#ed3c1b",
                             "#598fe7",
                             "#e38924",
                             "#654ebe",
                             "#40af6a",
                             "#cc46ba",
                             "#ab9036",
                             "#5a66c4",
                             "#d95726",
                             "#b67de3",
                             "#d88840",
                             "#a25bb1",
                             "#9e5314",
                             "#cd4997",
                             "#e17b56",
                             "#ca4972",
                             "#d43f34",
                             "#008f65",
                             "#b54637"), name="Country")+
  theme(
    panel.background = element_blank(),
    panel.grid = element_line(size = 0.001, linetype = 'dashed',
                              colour = "ghostwhite"), 
    legend.key = element_rect (fill = "white"),
    legend.text = element_text (size = 20),
    plot.title = element_text(face = "bold", hjust = 1,vjust = 1, size = 30),
    legend.title = element_text (size = 20, face = "bold"),
    axis.text = element_text(size =20),
    axis.line = element_line(color= "black", size =0.6),
    axis.text.x = element_blank(),#(size = 18, angle = 45,vjust= 1,hjust=1),
    axis.text.y = element_blank(),#text(size = 20),
    axis.title = element_blank(),#text(size =22, face = "bold"),
    legend.position="right"
  )
marginal2
library(ggExtra)
library(ggpubr)
pt <- ggarrange(p, marginal2, marginal1, ncol=2, nrow=2,  heights = c(2, 0.7), align = "v")
pt<- ggarrange(p, marginal2,marginal1, ncol=2, nrow=2, common.legend = FALSE)
pt
ggsave (plot=pt, filename="pock553heat.pdf", width=45, height=20, device="pdf", units = "in", limitsize=FALSE)

table(pock476$COUNTRY)
pock476[which(pock476$POS==14877),]
#from here we extract the sampleid, and we match it to samplecount ONLY NON SYNONYMOUS.
head(remde)
head(seqcount)
seqcount$RemdeRes <- "No"
seqcount$LopiRes <- "No"
seqcount$FiveFURes <- "No"
head(remde)
seqcount$RemdeRes <- remdeNS$Remdesivir_binding_affecting[match(seqcount$ID, remdeNS$SAMPLEID)]
seqcount$LopiRes <- lopiNS$Lopinavir_binding_pocket[match(seqcount$ID, lopiNS$SAMPLEID)]
seqcount$FiveFURes <- fivefuNS$Remdesivir_binding_affecting[match(seqcount$ID, fivefuNS$SAMPLEID)]
head(seqcount)
#uniform the names for sex
table(seqcount$Sex)
seqcount$Sex <- gsub("FEmale","Female",seqcount$Sex)
seqcount$Sex <- gsub("\\?","Unknown",seqcount$Sex)
seqcount$Sex <- gsub("Woman","Female",seqcount$Sex)
seqcount$Sex <- gsub("unknwon","Unknown",seqcount$Sex)
seqcount$Sex <- gsub("Unkown","Unknown",seqcount$Sex)
seqcount$Sex <- gsub("U","Unknown",seqcount$Sex)
seqcount$Sex <- gsub("Unknownnknown","Unknown",seqcount$Sex)
seqcount$Sex[which(is.na(seqcount$Sex))==TRUE] <- "Unknown"
#uniform the name for age
unique(seqcount$Age)
seqcount$Age <- gsub("\\?","Unknown", seqcount$Age)
seqcount$Age[which(is.na(seqcount$Age)==TRUE)] <- "Unknown"

#check
seqcount[which(seqcount$ID=="Australia-SAP319-2020"),]
#turn the NA into "No"
seqcount$RemdeRes[which(is.na(seqcount$RemdeRes)==TRUE)] <- "No"
seqcount$LopiRes[which(is.na(seqcount$LopiRes)==TRUE)] <- "No"
seqcount$FiveFURes[which(is.na(seqcount$FiveFURes)==TRUE)] <- "No"

#now let's see the counts/percentage per country of drug positive an S at the beginning indecates the dataframe comes from seqcounts
Sremde <- seqcount[which(seqcount$RemdeRes != "No"),]
Slopi <- seqcount[which(seqcount$LopiRes != "No"),]
Sfivefu <- seqcount[which(seqcount$FiveFURes != "No"),]
dim(Sfivefu)

#now dataframe of frequences per country
Scountrydf <- data.frame(Country=names(table(seqcount$Country)), 
                         Totoalseqs=unname(table(seqcount$Country)))[,-2]
#write.table(seqcount, "Seqcount.tab", sep = "\t", quote = F, row.names = F, col.names = T)
colnames(Scountrydf) <- c("Country","TotalSeqs")
stransr <- data.frame(Country=names(table(Sremde$Country)), RemdeposSeq= unname((table(Sremde$Country))))
stransl <- data.frame(Country=names(table(Slopi$Country)), LopiposSeq= unname((table(Slopi$Country))))
stransf <- data.frame(Country=names(table(Sfivefu$Country)), FiveposSeq= unname((table(Sfivefu$Country))))

Scountrydf$RemdeposC <- stransr$RemdeposSeq.Freq[match(Scountrydf$Country,stransr$Country)]
Scountrydf$RemdeposP <- Scountrydf$RemdeposC/Scountrydf$TotalSeqs*100
Scountrydf$LopiposC <- stransl$LopiposSeq.Freq[match(Scountrydf$Country,stransl$Country)]
Scountrydf$LopiposP <- Scountrydf$LopiposC/Scountrydf$TotalSeqs*100
Scountrydf$FiveFUposC <- stransf$FiveposSeq.Freq[match(Scountrydf$Country,stransf$Country)]
Scountrydf$FiveFUposP <- Scountrydf$FiveFUposC/Scountrydf$TotalSeqs*100
head(Scountrydf)
#the non match means they were not among the positives, so we put them to 0
Scountrydf[is.na(Scountrydf)] <- 0
###########################PANGOLIN
################same procedure for PANGOLIN
Spangodf <- data.frame(Pango=names(table(seqcount$Pangolin)), 
                         Totoalseqs=unname(table(seqcount$Pangolin)))[,-2]
colnames(Spangodf) <- c("Pango","TotalSeqs")

stransrP <- data.frame(Pango=names(table(Sremde$Pangolin)), RemdeposSeq= unname((table(Sremde$Pangolin))))
stranslP <- data.frame(Pango=names(table(Slopi$Pangolin)), LopiposSeq= unname((table(Slopi$Pangolin))))
stransfP <- data.frame(Pango=names(table(Sfivefu$Pangolin)), FiveposSeq= unname((table(Sfivefu$Pangolin))))


Spangodf$RemdeposC <- stransrP$RemdeposSeq.Freq[match(Spangodf$Pango,stransrP$Pango)]
Spangodf$RemdeposP <- Spangodf$RemdeposC/Spangodf$TotalSeqs*100
Spangodf$LopiposC <- stranslP$LopiposSeq.Freq[match(Spangodf$Pango,stranslP$Pango)]
Spangodf$LopiposP <- Spangodf$LopiposC/Spangodf$TotalSeqs*100
Spangodf$FiveFUposC <- stransfP$FiveposSeq.Freq[match(Spangodf$Pango,stransfP$Pango)]
Spangodf$FiveFUposP <- Spangodf$FiveFUposC/Spangodf$TotalSeqs*100
head(Spangodf)


#the non match means they were not among the positives, so we put them to 0
Spangodf[is.na(Spangodf)] <- 0
######SEX
################same procedure for sex
Ssexdf <- data.frame(Sex=names(table(seqcount$Sex)), 
                     Totoalseqs=unname(table(seqcount$Sex)))[,-2]
colnames(Ssexdf) <- c("Sex","TotalSeqs")

stransrS <- data.frame(Sex=names(table(Sremde$Sex)), RemdeposSeq= unname((table(Sremde$Sex))))
stranslS <- data.frame(Sex=names(table(Slopi$Sex)), LopiposSeq= unname((table(Slopi$Sex))))
stransfS <- data.frame(Sex=names(table(Sfivefu$Sex)), FiveposSeq= unname((table(Sfivefu$Sex))))


Ssexdf$RemdeposC <- stransrS$RemdeposSeq.Freq[match(Ssexdf$Sex,stransrS$Sex)]
Ssexdf$RemdeposP <- Ssexdf$RemdeposC/Ssexdf$TotalSeqs*100
Ssexdf$LopiposC <- stranslS$LopiposSeq.Freq[match(Ssexdf$Sex,stranslS$Sex)]
Ssexdf$LopiposP <- Ssexdf$LopiposC/Ssexdf$TotalSeqs*100
Ssexdf$FiveFUposC <- stransfS$FiveposSeq.Freq[match(Ssexdf$Sex,stransfS$Sex)]
Ssexdf$FiveFUposP <- Ssexdf$FiveFUposC/Ssexdf$TotalSeqs*100
head(Ssexdf)
sum(Ssexdf$TotalSeqs)
#the non match means they were not among the positives, so we put them to 0
Ssexdf[is.na(Ssexdf)] <- 0
######AGE
################same procedure for page
Sagedf <- data.frame(Age=names(table(seqcount$Age)), 
                     Totoalseqs=unname(table(seqcount$Age)))[,-2]
colnames(Sagedf) <- c("Age","TotalSeqs")

stransrA <- data.frame(Age=names(table(Sremde$Age)), RemdeposSeq= unname((table(Sremde$Age))))
stranslA <- data.frame(Age=names(table(Slopi$Age)), LopiposSeq= unname((table(Slopi$Age))))
stransfA <- data.frame(Age=names(table(Sfivefu$Age)), FiveposSeq= unname((table(Sfivefu$Age))))


Sagedf$RemdeposC <- stransrA$RemdeposSeq.Freq[match(Sagedf$Age,stransrA$Age)]
Sagedf$RemdeposP <- Sagedf$RemdeposC/Sagedf$TotalSeqs*100
Sagedf$LopiposC <- stranslA$LopiposSeq.Freq[match(Sagedf$Age,stranslA$Age)]
Sagedf$LopiposP <- Sagedf$LopiposC/Sagedf$TotalSeqs*100
Sagedf$FiveFUposC <- stransfA$FiveposSeq.Freq[match(Sagedf$Age,stransfA$Age)]
Sagedf$FiveFUposP <- Sagedf$FiveFUposC/Sagedf$TotalSeqs*100
sum(Sagedf$TotalSeqs)
#linear regression with the positive count
lmRemcount <- lm(RemdeposC~Age, data=Sagedf)
summary(lmRemcount)
#the non match means they were not among the positives, so we put them to 0
Sagedf[is.na(Sagedf)] <- 0
Sagedf$Age[104] <- NA
Sagedf <- Sagedf[order(as.numeric(as.character(Sagedf$Age))),]
Sagedf$Age <- as.numeric(as.character(Sagedf$Age))

########now all the dfs together in a list, so that we can plot it all together in one loop
Sdflist <- list(Country=Scountrydf,Pangolin=Spangodf,Sex=Ssexdf,Age=Sagedf)
###some quick check that everything adds up
##20 sequences discrepancy between age df and others because infos was not provided already from the metadata file 
#and instead we put age as Unknown to NA, we didn't do such a thing with other variables
setwd("NSMUT/")
for (w in 1:length(Sdflist)){
  G <- as.data.frame(Sdflist[[w]])
  listtag <- names(Sdflist)[w]
  print(paste("Total seqs sum for ", listtag, ":", sep = ""))
  print(sum(G$TotalSeqs))
  print(paste("Total seqs mean for ", listtag, ":", sep = ""))
  print(mean(G$TotalSeqs))
  
  print(paste("RemdeposC sum for ", listtag, ":", sep = ""))
  print(sum(G$RemdeposC))
  print(paste("RemdeposC mean for ", listtag, ":", sep = ""))
  print(mean(G$RemdeposC))
  
  print(paste("RemdeposP sum for ", listtag, ":", sep = ""))
  print(sum(G$RemdeposP))
  print(paste("RemdeposP mean for ", listtag, ":", sep = ""))
  print(mean(G$RemdeposP))
  
  print(paste("LopiposC sum for ", listtag, ":", sep = ""))
  print(sum(G$LopiposC))
  print(paste("LopiposC mean for ", listtag, ":", sep = ""))
  print(mean(G$LopiposC))
  
  print(paste("LopiposP sum for ", listtag, ":", sep = ""))
  print(sum(G$LopiposP))
  print(paste("LopiposP mean for ", listtag, ":", sep = ""))
  print(mean(G$LopiposP))
  
  print(paste("FiveFUposC sum for ", listtag, ":", sep = ""))
  print(sum(G$FiveFUposC))
  print(paste("FiveFUposC mean for ", listtag, ":", sep = ""))
  print(mean(G$FiveFUposC))
  
  print(paste("FiveFUposP sum for ", listtag, ":", sep = ""))
  print(sum(G$FiveFUposP))
  print(paste("FiveFUposP mean for ", listtag, ":", sep = ""))
  print(mean(G$FiveFUposP))
}
head(seqcount)
dim(seqcount)
library(reshape2)
library(ggplot2)
#install.packages("ggforce")
library(ggforce)
#ideas for plots
for (k in 1:length(Sdflist)){
  datf <- as.data.frame(Sdflist[[k]])
  #dim(datf)
  plottag <- names(Sdflist)[k]
  head(datf)
dfrR <- datf[,-c(2,5,6,7,8)]
dfrL <- datf[,-c(2,3,4,7,8)]
dfrF <- datf[,-c(2,3,4,5,6)]
#now absolute counts
dfrRA <- datf[,-c(4,5,6,7,8)]
dfrLA <- datf[,-c(3,4,6,7,8)]
dfrFA <- datf[,-c(3,4,5,6,8)]
#dfr$RemdeposP <- -1*dfr$RemdeposP
if (plottag=="Age"){
  mdfR <- melt(dfrR, id.vars="Age")
  mdfL <- melt(dfrL, id.vars="Age")
  mdfF <- melt(dfrF, id.vars="Age")
  
  mdfRA <- melt(dfrRA, id.vars="Age")
  mdfLA <- melt(dfrLA, id.vars="Age")
  mdfFA <- melt(dfrFA, id.vars="Age")
}else{
mdfR <- melt(dfrR)
mdfL <- melt(dfrL)
mdfF <- melt(dfrF)

mdfRA <- melt(dfrRA)
mdfLA <- melt(dfrLA)
mdfFA <- melt(dfrFA)

}
mdflist <- list(Remdesivir_relative=mdfR,Lopinavir_relative=mdfL, FiveF_Uracile_relative=mdfF, 
                Remdesivir_absolute=mdfRA,Lopinavir_absolute=mdfLA, FiveF_Uracile_absolute=mdfFA)
names(mdflist)
for (i in 1:length(mdflist)){
  D <- as.data.frame(mdflist[[i]])
  tag <- names(mdflist)[i]
print("Plottag is")
print(plottag)
print("tag is")
print(tag)
if (grepl("relative", tag)==TRUE){
  breaks <- c(5,10,15,20,25,30,100,300,400)
  leglabs <- c("Absolute positive consensus counts","Percentage positive consensus")
  if (plottag =="Age"){
    breaks <- c(5,10,15)
  }
}else{
  breaks <- c(100,200,300,400,500,1000,5000,10000,15000, 20000)
  leglabs <- c("Absolute consensus counts","Absolute positive consensus counts")
  if (plottag =="Age"){
    breaks <- c(100,200,300,400,500)
  }
}
tit <- paste(tag, " resistance prevalence", sep = "")
#plot trial
print("tit is")
print(tit)
q <- ggplot(D,aes(x=D[,1], y=value, fill=variable))+
  geom_bar(stat="identity",position="dodge")+
  xlab(plottag)+
  ylab("Resistance prevalence")+
  scale_fill_manual(name=tit,values = c("#FFA373","#50486D", "#464730"), labels=leglabs)+
  scale_y_sqrt(breaks=breaks)+
  ggtitle(tit)+
  theme(
    panel.background = element_blank(),
    panel.grid = element_line(size = 0.001, linetype = 'dashed',
                              colour = "ghostwhite"), 
    legend.key = element_rect (fill = "white"),
    legend.text = element_text (size = 20),
    plot.title = element_text(face = "bold", hjust = 1,vjust = 1, size = 30),
    legend.title = element_text (size = 20, face = "bold"),
    axis.text = element_text(size =20),
    axis.line = element_line(color= "black", size =0.6),
    axis.text.x = element_text(size = 19, angle = 60,vjust= 1,hjust=1),
    axis.text.y = element_text(size = 24),
    axis.title = element_text(size =22, face = "bold"),
    legend.position="top"
    
  )
out <- paste(plottag,"_",tag,".pdf", sep="")
print("Out is")
print(out)

#pdf(out,width=40, height=12)
#q
#dev.off()
if(plottag=="Pangolin"){
  ggsave (plot=q, filename=out, width=30, height=15, device="pdf", units = "in", limitsize=FALSE)
}else{
  ggsave (plot=q, filename=out, width=30, height=20, device="pdf", units = "in", limitsize=FALSE)
}
  }
}
#now anova on the country plot to see which one stands out
head(all)
allNS <- all[which(all$ANNOTATION!="synonymous_variant"),]
allNS$REMC <- 0
allNS$REMC[which(grepl("Yes",allNS$Remdesivir_binding_affecting)==TRUE & grepl("611",allNS$Remdesivir_binding_affecting)==FALSE)] <- 1
table(allNS$REMC)
aov <- aov(REMC~COUNTRY, data=allNS)
summary(aov)
t <- TukeyHSD(aov)
head(t[[1]][3])
head(t)
t[[1]][,3]
anova(aov)
plot(aov,2)
aov.res <- residuals(aov)
shapiro.test(x=aov.res)
kr <- kruskal.test(Country~RemdeposC, data= Sdflist$Country)
#no significance predictor among countries for remdepos
#some follow ups: lineage B.1.1.29 (100% Remde)-B.1.111 (100% Lopi)-B.1.5.8 -B.1.6 (both 100% 5FU) show all to have resistance mutations
#let's see more of them

seqcount[which(seqcount$Pangolin=="B.1.1.29"),] #all from Wales (Remde)
seqcount[which(seqcount$Pangolin=="B.1.5.5"),] #all from Wales (Remde)

allNS[which(allNS$PANGOLIN=="B.1.1.29"),]
seqcount[which(seqcount$Pangolin=="B.1.111"),] #1 Venezuela, 2 USA, 11 England, 11 Colombia, 2 Chile, 21 Canada (Lopi)
seqcount[which(seqcount$Pangolin=="B.1.5.8"),] #1 Belgium, 1 Denmark, 19 USA Texas (5FU)
seqcount[which(seqcount$Pangolin=="B.1.6"),] #1 Australia, 1 Austria, 21 Belgium, 4 DR. Congo, 2 England (5FU)

########now we add the date information and we plot cumulative date plot, following same strategy of the other df
head(seqcount)
seqcount$Date <- all$DATE[match(seqcount$ID, all$SAMPLEID)]
dim(seqcount)
library("reshape2")
brdate <- colsplit(string = seqcount$Date, pattern= "-", names=c("Year","Month","Day"))
table(brdate$Year, useNA = "always")
#uniform the entries
brdate$Day[which(brdate$Day=="")] <- NA
brdate$Day[which(brdate$Day=="XX")] <- NA
seqcount <- cbind(seqcount,brdate)
head(seqcount)
library(lubridate)
seqcount$Week <- lubridate::week(ymd(seqcount$Date))
seqcount$Week <- paste(seqcount$Year,seqcount$Week, sep ="_")
head(seqcount)
unique(seqcount$Week)
seqcount$Week <- factor(seqcount$Week, levels=c(
  "2019_52","2020_1","2020_2","2020_3","2020_4","2020_5","2020_6",
  "2020_7","2020_8","2020_9","2020_10","2020_11","2020_12","2020_13",
  "2020_14","2020_15","2020_16","2020_17","2020_18","2020_19","2020_20",
  "2020_21","2020_22","2020_23","2020_24","2020_25",
  "2020_26","2020_27","2020_28","2020_NA",
  "NA_NA"
))
head(seqcount)
#new dataframe just for the date
SremdeD <- seqcount[which(seqcount$RemdeRes != "No"),]
SlopiD <- seqcount[which(seqcount$LopiRes != "No"),]
SfivefuD <- seqcount[which(seqcount$FiveFURes != "No"),]
table(seqcount$Week)
Sdatedf <- data.frame(Date=names(table(seqcount$Week)), 
                     Totoalseqs=unname(table(seqcount$Week)))[,-2]
colnames(Sdatedf) <- c("Date","TotalSeqs")

stransrD <- data.frame(Date=names(table(SremdeD$Week)), RemdeposSeq= unname((table(SremdeD$Week))))
stranslD <- data.frame(Date=names(table(SlopiD$Week)), LopiposSeq= unname((table(SlopiD$Week))))
stransfD <- data.frame(Date=names(table(SfivefuD$Week)), FiveposSeq= unname((table(SfivefuD$Week))))


Sdatedf$RemdeposC <- stransrD$RemdeposSeq.Freq[match(Sdatedf$Date,stransrD$Date)]
Sdatedf$RemdeposP <- Sdatedf$RemdeposC/Sdatedf$TotalSeqs*100
Sdatedf$LopiposC <- stranslD$LopiposSeq.Freq[match(Sdatedf$Date,stranslD$Date)]
Sdatedf$LopiposP <- Sdatedf$LopiposC/Sdatedf$TotalSeqs*100
Sdatedf$FiveFUposC <- stransfD$FiveposSeq.Freq[match(Sdatedf$Date,stransfD$Date)]
Sdatedf$FiveFUposP <- Sdatedf$FiveFUposC/Sdatedf$TotalSeqs*100
head(Sdatedf)
#now we try cumulative
Sdatedf$TotalSeqsCum <- cumsum(Sdatedf$TotalSeqs)
Sdatedf$RemdeposCCum <- cumsum(Sdatedf$RemdeposC)
Sdatedf$LopiposCCum <- cumsum(Sdatedf$LopiposC)
Sdatedf$FiveFUposCCum <- cumsum(Sdatedf$FiveFUposC)
Sdatedf$RemdeposCCumP <- Sdatedf$RemdeposCCum/Sdatedf$TotalSeqsCum*100
Sdatedf$LopiposCCumP <- Sdatedf$LopiposCCum/Sdatedf$TotalSeqsCum*100
Sdatedf$FiveFUposCCumP <- Sdatedf$FiveFUposCCum/Sdatedf$TotalSeqsCum*100
head(Sdatedf)
#get a plot
Sdateplot <- Sdatedf[,-c(2:8)] 
head(Sdateplot)
Sdateplot$Date1 <- factor(Sdateplot$Date, levels=c(
  "2019_52","2020_1","2020_2","2020_3","2020_4","2020_5","2020_6",
  "2020_7","2020_8","2020_9","2020_10","2020_11","2020_12","2020_13",
  "2020_14","2020_15","2020_16","2020_17","2020_18","2020_19","2020_20",
  "2020_21","2020_22","2020_23","2020_24","2020_25",
  "2020_26","2020_27","2020_28","2020_NA",
  "NA_NA"
)) 
summary(Sdateplot)
head(Sdateplot)
#let's see the mean and the variance for remdesivir over time
#to be rethought
#Sdateplot$RemdeposCCumMean <- Sdateplot$RemdeposCCum
head(Sdateplot)
marginalb <- ggplot(Sdateplot,aes(x=Date1, y=TotalSeqsCum))+
  geom_col()+
  scale_y_reverse()+
  xlab("Week")+
  ylab("Counts")+
  theme(
    panel.background = element_blank(),
    panel.grid = element_line(size = 0.001, linetype = 'dashed',
                              colour = "ghostwhite"), 
    legend.key = element_rect (fill = "white"),
    legend.text = element_text (size = 20),
    plot.title = element_text(face = "bold", hjust = 1,vjust = 1, size = 30),
    legend.title = element_text (size = 20, face = "bold"),
    axis.text = element_text(size =20),
    axis.line = element_line(color= "black", size =0.6),
    axis.text.x = element_text(size = 18, angle = 45,vjust= 1,hjust=1),
    axis.text.y = element_blank(),#text(size = 20),
    axis.title = element_blank(),#text(size =22, face = "bold"),
    legend.position="top"
  )

marginalb
mdf <- melt(Sdateplot[,-c(1,6,7,8)], id.vars="Date1")
mdf1 <- melt(Sdateplot[,-c(1:5)], id.vars="Date1") 
qi <-ggplot(mdf1, aes(x=Date1, y= value, color=variable))+
  geom_line(aes(group=variable), size=1.9)+
  #scale_y_sqrt()+
  xlab("Weeks")+
  ylab("Percentage of positive genomes")+
  scale_color_manual(name="Percentage of Positive genomes:",values=c("#73c8a1","#dd827f","#beb563"),labels=c("Remdesivir", "Lopinavir",  "5-Fluoro Uracile"))+
  theme(
    panel.background = element_blank(),
    panel.grid = element_line(size = 0.001, linetype = 'dashed',
                              colour = "ghostwhite"), 
    legend.key = element_rect (fill = "white"),
    legend.text = element_text (size = 20),
    plot.title = element_text(face = "bold", hjust = 1,vjust = 1, size = 30),
    legend.title = element_text (size = 20, face = "bold"),
    axis.text = element_text(size =20),
    axis.line = element_line(color= "black", size =0.6),
    axis.text.x = element_blank(),#text(size = 18, angle = 45,vjust= 1,hjust=1),
    axis.text.y = element_blank(),#text(size = 20),
    axis.title.x = element_blank(),
    axis.title = element_blank(),#text(size =22, face = "bold"),
    legend.position="top"
  )+
  guides(color = guide_legend(override.aes = list(size = 2)))
qi
marginalb
qb<- ggarrange(qi, marginalb, ncol=1, nrow=2, common.legend = FALSE, align = NULL)
qb
ggsave (plot=qb, filename="Pecentage_resistance_overtime2.pdf", width=15, height=18, device="pdf", units = "in", limitsize=FALSE)

head(seqcount)
head(seqcount)
Sdateplot$Month <- seqcount$Month[match(Sdateplot$Date1, seqcount$Week)]
Sdateplot[,c(9,10)]
##NOw let's investigate the diversity or mutation per sequence, to see how much the positiveness to resistance 
#is standing out
library("SpadeR")
library("Matrix")
#starting df
df <- data.frame(Positions= allNS$POS, Samples= allNS$SAMPLEID)
df$Samplenumb <- as.numeric(df$Samples)
df$Counts <- 1
#j will be the position, i the numering of the sampleid
head(df)
df$Positions <- as.numeric(as.character(df$Positions))
length(df$Samples)
A <- sparseMatrix(df$Samplenumb, df$Positions, x = df$Counts, dimnames = list(unique(df$Samples), seq(1,29903, by=1)))
#A[A="."] <- 0
which(A[,241] >0)==TRUE
#B <- as.matrix(as.data.frame(as.matrix(A)))
A1 <- A[1:10000,1:3000]

abu_count_sample <- data.frame(Names= names(table(df$Samples)), Counts=unname(table(df$Samples)))[,-2]
abu_count_s <- data.frame(row.names = abu_count_sample$Names, Counts = abu_count_sample$Counts.Freq)
head(abu_count_s)
sum(abu_count_s$Counts)
#DiversityData$Abu_count
#div <- Diversity(abu_count_s,"abundance_freq_count", q=NULL)
#div$Simpson_diversity
head(A)
#A[1,]
#div <- Diversity(A[1,], "incidence_freq_count", q=NULL)
library("fossil")
#A11 <- as.matrix(as.data.frame(as.matrix(A)))
bip <- data.frame()
for (i in 1:nrow(A)){
  samname <- rownames(A)[i]
  ch2 <- chao2(A[i,], taxa.row = F)
  bop <- data.frame(Sample=samname, chao2=ch2)
  bip <- rbind(bip, bop)
}
head(bip)
head(abu_count_s)
head(seqcount)
seqcount$Chao2 <- bip$chao2[match(seqcount$ID, bip$Sample)]
head(seqcount)
class(seqcount$Chao2)
#we plot now the diversity of negative and positive samples
divplot <- seqcount[,-c(1:5,9:13)]
head(divplot)
qp <- ggplot(divplot, aes(x=RemdeRes, y=Chao2, color=RemdeRes))+
  geom_boxplot(size=1.5)+
  geom_jitter(width = 0.1)+
  ylim(0,200)+
  xlab("Resistant")+
  ylab("Mutation diversity (Chao2 incidence index)")+
  scale_color_manual(name="Remdesivir resistant",values=c("#FFA373","#50486D", "#464730"),labels=c("Non resistant", "Resistant"))+
  theme(
    panel.background = element_blank(),
    panel.grid = element_line(size = 0.001, linetype = 'dashed',
                              colour = "ghostwhite"), 
    legend.key = element_rect (fill = "white"),
    legend.text = element_text (size = 20),
    plot.title = element_text(face = "bold", hjust = 1,vjust = 1, size = 30),
    legend.title = element_text (size = 20, face = "bold"),
    axis.text = element_text(size =20),
    axis.line = element_line(color= "black", size =0.6),
    axis.text.x = element_text(size = 18, angle = 45,vjust= 1,hjust=1),
    axis.text.y = element_text(size = 20),
    axis.title = element_text(size =22, face = "bold")
    #legend.position="top"
  )+
  guides(color = guide_legend(override.aes = list(size = 2)))
ggsave (plot=qp, filename="Diversity_Remdesivir23.pdf", width=12, height=10, device="pdf", units = "in", limitsize=FALSE)
pdf("Hist_chao2diversity_disclaimer.pdf")
hist(divplot$Chao2,xlim = c(0,1000), ylim=c(0,25000),breaks = 3000)
dev.off()
seqcount[which(seqcount$Chao2>5000),]
divplot$RemdeRes[which(divplot$RemdeRes != "No")] <- "Yes"
divplot$LopiRes[which(divplot$LopiRes != "No")] <- "Yes"
divplot$FiveFURes[which(divplot$FiveFURes != "No")] <- "Yes"
t.test(Chao2~RemdeRes,data = divplot)
pairwise.t.test(divplot$Chao2,divplot$RemdeRes, data=divplot)
t.test(Chao2~LopiRes,data = divplot)

sink("Remdesivir_annexdiversity_T.tab")
t.test(Chao2~RemdeRes, data = divplot)
sink()
sink("Lopinavir_annexdiversity_T.tab")
t.test(Chao2~LopiRes, data = divplot)
sink()
sink("5Fluoro_Uracil_annexdiversity_T.tab")
t.test(Chao2~FiveFURes,data = divplot)
sink()

#two graphics
head(seqcount)
ggplot(seqcount, aes(x=Week, y=Chao2))+
  geom_point()+
  facet_grid(.~RemdeRes)+
  ylim(c(0,1000))

##
#Now the maps sequence based and also we will do the diversity per country
library("ggspatial")
library("sf")
library("rnaturalearthdata")
library("rnaturalearth")

#install.packages("rgeos")
library("rgeos")
#set the df and change the unclear name into standard names, so that we can fit them into the world map
Scountrydf$Official <- as.character(Scountrydf$Country)
Scountrydf$Official[which(Scountrydf$Official=="Gambia")] <- "Gambia, The"
Scountrydf$Official[which(Scountrydf$Official=="Russia")] <- "Russian Federation"
Scountrydf$Official[which(Scountrydf$Official=="Democratic Republic of the Congo")] <- "Congo, Dem. Rep."
Scountrydf$Official[which(Scountrydf$Official=="USA")] <- "United States of America"
Scountrydf$Official[which(Scountrydf$Official=="Kyrgyz Rep.")] <- "Kyrgyz Republic"
Scountrydf$Official[which(Scountrydf$Official=="Egypt")] <- "Egypt, Arab Rep."
Scountrydf$Official[which(Scountrydf$Official=="Iran")] <- "Iran, Islamic Rep."
Scountrydf$Official[which(Scountrydf$Official=="North Macedonia")] <- "Macedonia, FYR"
Scountrydf$Official[which(Scountrydf$Official=="Hong Kong")] <- "Hong Kong SAR, China"
Scountrydf$Official[which(Scountrydf$Official=="Venezuela")] <- "Venezuela, RB"
Scountrydf$Official[which(Scountrydf$Official=="South Korea")] <- "Korea, Rep."
Scountrydf$Official[which(Scountrydf$Official=="Slovakia")] <- "Slovak Republic"
#now we call the worldmap dataframe

world <- ne_countries(scale="medium", returnclass = "sf")



head(Scountrydf)
world$remde_perc <- Scountrydf$RemdeposP[match(world$name_sort, Scountrydf$Official)]
world$lopi_perc <- Scountrydf$LopiposP[match(world$name_sort, Scountrydf$Official)]
world$fiveFU_perc <- Scountrydf$FiveFUposP[match(world$name_sort, Scountrydf$Official)]
# and now plot
p <- ggplot(data = world) +
  geom_sf(color="slategrey", aes(fill=fiveFU_perc))+
  coord_sf(expand = FALSE)+
  xlab("Longitude") + ylab("Latitude") +
  scale_fill_viridis_c(option = "viridis", trans="sqrt", na.value = "white", 
                       name="5-Fluoro-Uracile resistance percentage")+
  theme(
    panel.background = element_rect(fill = "aliceblue")
  )
ggsave (p,file = "Seqcount_Map_FiveFUperc.pdf", width=18, height=12, units = "in", limitsize=FALSE)

#a zoom in europe
head(world)
p1 <- ggplot(data = world) +
  geom_sf(color="slategrey", aes(fill=fiveFU_perc))+
  #coord_sf(expand = FALSE)+
  xlab("Longitude") + ylab("Latitude") +
  scale_fill_viridis_c(option = "viridis", trans="sqrt", na.value = "white", 
                       name="5-Fluoro Uracile resistance percentage")+
  theme(
    panel.background = element_rect(fill = "aliceblue")
  )+
    coord_sf(xlim = c(-30,50), ylim = c(30,80), expand = FALSE)+
  annotation_scale(location = "tl", width_hint = 0.3) +
  annotation_north_arrow(location = "tr", which_north = "true", 
                         pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"),
                         style = north_arrow_fancy_orienteering) 
p1
ggsave (p1,file = "Seqcount_Map_FiveFUperc_ZoomEu.pdf", width=12, height=12, units = "in", limitsize=FALSE)
# now let's calculate the diversity per country, and then plot it the same way
df1 <- data.frame(Positions= all$POS, Countries= all$COUNTRY)
df1$Samplenumb <- as.numeric(df1$Countries)
head(df1)
df1$Counts <- 1
#j will be the position, i the numering of the sampleid
head(df1)
df1$Positions <- as.numeric(as.character(df1$Positions))
#df2 <- df1[-is.na(df1$Countries),]
df2 <- df1[-which(is.na(df1$Countries)==TRUE),]
df1[-is.na(df2$Countries),]
AC <- sparseMatrix(df2$Samplenumb, df2$Positions, x = df2$Counts, dimnames = list(unique(df2$Countries), seq(1,29903, by=1)))
head(AC)
bipC <- data.frame()
for (i in 1:nrow(AC)){
  samname <- rownames(AC)[i]
  ch1 <- chao1(AC[i,], taxa.row = F)
  bopC <- data.frame(Sample=samname, chao1=ch1)
  bipC <- rbind(bipC, bopC)
}
head(bipC)
#we get the name straigth, so that we can match it ti world df
bipC$Official <- Scountrydf$Official[match(bipC$Sample, Scountrydf$Country)]
world$chao1_div <- bipC$chao1[match(world$name_sort, bipC$Official)]
#now plot the world diversity map

p <- ggplot(data = world) +
  geom_sf(color="slategrey", aes(fill=chao1_div))+
  coord_sf(expand = FALSE)+
  xlab("Longitude") + ylab("Latitude") +
  scale_fill_viridis_c(option = "viridis", trans="sqrt", na.value = "white", 
                       name="Genome entropy (chao1 index)")+
  theme(
    panel.background = element_rect(fill = "aliceblue")
  )
p
ggsave (p,file = "Seqcount_Map_Diversity.pdf", width=18, height=10, units = "in", limitsize=FALSE)

#a zoom in europe
head(world)
p1 <- ggplot(data = world) +
  geom_sf(color="slategrey", aes(fill=chao1_div))+
  #coord_sf(expand = FALSE)+
  xlab("Longitude") + ylab("Latitude") +
  scale_fill_viridis_c(option = "viridis", trans="sqrt", na.value = "white", 
                       name="Genome entropy (chao1 index)")+
  theme(
    panel.background = element_rect(fill = "aliceblue")
  )+
  coord_sf(xlim = c(-30,50), ylim = c(30,80), expand = FALSE)+
  annotation_scale(location = "tl", width_hint = 0.3) +
  annotation_north_arrow(location = "tr", which_north = "true", 
                         pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"),
                         style = north_arrow_fancy_orienteering) 
p1
ggsave (p1,file = "Seqcount_Map_Diversity_ZoomEu.pdf", width=12, height=12, units = "in", limitsize=FALSE)



#  ylim(c(0,500))
#scale_y_continuous(breaks = pretty(mdf$value),labels = abs(pretty(mdf$value)))
  #theme_scientific()
?paste
ggplot(mdf, aes(x=Country, y=value))+
  geom_bar(stat="identity",position="identity")
  #scale_y_discrete(breaks=seq(-100,100,by=50),labels=abs(seq(-100,100,by=50)))+
  scale_x_log10()

######Now we investigate per mutation, that is, not counting the sample ids, but just the mutationn as a whole
Countryseq <- data.frame(Country=names(table(seqcount$Country)), Counts=unname(table(seqcount$Country)))
head(Countryseq)
length(unique(meta$strain))
table(all$Remdesivir_binding_affecting)

table(all$ANNOTATION)
remde <- all[which(grepl("Yes",all$Remdesivir_binding_affecting)==TRUE & grepl("_M611", all$Remdesivir_binding_affecting)==FALSE),]
lopi <- all[which(grepl("Yes",all$Lopinavir_binding_pocket)==TRUE),]
fivefu <- all[which(grepl("_M611", all$Remdesivir_binding_affecting)==TRUE),]

#we create a percentage of rmede resistance in percentage of the MUTATIONS AVAILABLE
Countrydf <- data.frame(Countries=names(table(all$COUNTRY)), Allmutations=unname(table(all$COUNTRY)))
Countrydf <- Countrydf[,-2]

transremde <- data.frame(Countries=names(table(remde$COUNTRY)), Resistant=unname(table(remde$COUNTRY)))
translopi <- data.frame(Countries=names(table(lopi$COUNTRY)), Resistant=unname(table(lopi$COUNTRY)))
trans5Ffu <- data.frame(Countries=names(table(fivefu$COUNTRY)), Resistant=unname(table(fivefu$COUNTRY)))
#percentages of remde resistant
Countrydf$RemResistantSeqs <- transremde$Resistant.Freq[match(Countrydf$Countries,transremde$Countries)]
Countrydf$RemResistantSeqs[which(is.na(Countrydf$RemResistantSeqs)==TRUE)] <- 0
Countrydf$RemdePerc <- Countrydf$RemResistantSeqs/Countrydf$Allmutations.Freq*100
#percentages of lopi resistant
Countrydf$LopiResistantSeqs <- translopi$Resistant.Freq[match(Countrydf$Countries,translopi$Countries)]
Countrydf$LopiResistantSeqs[which(is.na(Countrydf$LopiResistantSeqs)==TRUE)] <- 0
Countrydf$LopiPerc <- Countrydf$LopiResistantSeqs/Countrydf$Allmutations.Freq*100
#percentages of 5ffu resistant
Countrydf$FiveFUResistantSeqs <- trans5Ffu$Resistant.Freq[match(Countrydf$Countries,trans5Ffu$Countries)]
Countrydf$FiveFUResistantSeqs[which(is.na(Countrydf$FiveFUResistantSeqs)==TRUE)] <- 0
Countrydf$FiveFUPerc <- Countrydf$FiveFUResistantSeqs/Countrydf$Allmutations.Freq*100
Countrydf$SecCount <- Countryseq$Counts.Freq[match(Countrydf$Countries,Countryseq$Country)]
Countrydf
#plot it 
library(ggplot2)
library(reshape2)
install.packages("rgeos")
library("rnaturalearthdata")
world <- ne_countries(scale="medium", returnclass = "sf")
which(is.na(match(Countrydf$Countries,world$name_sort))==TRUE)
Countrydf$Official[which(is.na(match(Countrydf$Official,world$name_sort))==TRUE)]
Countrydf$Official <- as.character(Countrydf$Countries)
world[which(grepl("Korea", world$name_sort)==TRUE),]
Countrydf$Official[which(Countrydf$Official=="Gambia")] <- "Gambia, The"
Countrydf$Official[which(Countrydf$Official=="Russia")] <- "Russian Federation"
Countrydf$Official[which(Countrydf$Official=="Democratic Republic of the Congo")] <- "Congo, Dem. Rep."
Countrydf$Official[which(Countrydf$Official=="USA")] <- "United States of America"
Countrydf$Official[which(Countrydf$Official=="Kyrgyz Rep.")] <- "Kyrgyz Republic"
Countrydf$Official[which(Countrydf$Official=="Egypt")] <- "Egypt, Arab Rep."
Countrydf$Official[which(Countrydf$Official=="Iran")] <- "Iran, Islamic Rep."
Countrydf$Official[which(Countrydf$Official=="North Macedonia")] <- "Macedonia, FYR"
Countrydf$Official[which(Countrydf$Official=="Hong Kong")] <- "Hong Kong SAR, China"
Countrydf$Official[which(Countrydf$Official=="Venezuela")] <- "Venezuela, RB"
Countrydf$Official[which(Countrydf$Official=="South Korea")] <- "Korea, Rep."
Countrydf$Official[which(Countrydf$Official=="Slovakia")] <- "Slovak Republic"

world$remde_perc <- Countrydf$RemdePerc[match(world$name_sort, Countrydf$Official)]
world$lopi_perc <- Countrydf$LopiPerc[match(world$name_sort, Countrydf$Official)]
world$fiveFU_perc <- Countrydf$FiveFUPerc[match(world$name_sort, Countrydf$Official)]

library("ggspatial")
library("sf")
ggplot(data = world) +
  geom_sf(aes(fill=remde_lopi))+
  scale_fill_viridis_c(option = "plasma", trans="sqrt", na.value = "white")
p <- ggplot(data = world) +
  geom_sf(color="slategrey", aes(fill=remde_perc))+
  coord_sf(expand = FALSE)+
  xlab("Longitude") + ylab("Latitude") +
  scale_fill_viridis_c(option = "viridis", trans="sqrt", na.value = "white", name="Remdesivir resistance percentage")+
  theme(
    panel.background = element_rect(fill = "aliceblue")
  )

ggsave (p,file = "map_remde1.pdf", width=16, height=12, units = "in", limitsize=FALSE)

r <- ggplot(data = world) +
  geom_sf(color="slategrey", aes(fill=lopi_perc))+
  coord_sf(expand = FALSE)+
  xlab("Longitude") + ylab("Latitude") +
  scale_fill_viridis_c(option = "viridis", trans="sqrt", na.value = "white", name="Lopinavir resistance percentage")+
  theme(
    panel.background = element_rect(fill = "aliceblue")
  )

ggsave (r,file = "map_lopi1.pdf", width=16, height=12, units = "in", limitsize=FALSE)

s <- ggplot(data = world) +
  geom_sf(color="slategrey", aes(fill=fiveFU_perc))+
  coord_sf(expand = FALSE)+
  xlab("Longitude") + ylab("Latitude") +
  scale_fill_viridis_c(option = "viridis", trans="sqrt", na.value = "white", name="5-FluoroU resistance percentage")+
  theme(
    panel.background = element_rect(fill = "aliceblue")
  )

ggsave (s,file = "map_5fu1.pdf", width=16, height=12, units = "in", limitsize=FALSE)



#ggsave(p,"map_remde.pdf", width = 16, height = 12)
  #annotation_scale(location = "bl", width_hint = 0.5) +
  #annotation_north_arrow(location = "bl", which_north = "true", 
  #                       pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
  #                       style = north_arrow_fancy_orienteering)
  


head(Countrydf)
#melt
percdf <- Countrydf[,-c(2,3,5,7)]
head(percdf)
mCD <- melt(percdf)
head(mCD)
ggplot(mCD)+
  #geom_col()+
  geom_sf()
  #theme(
  #  axis.text.x = element_text(angle=45,vjust=1,hjust=1),
  #)+
  #facet_grid(variable~.)
##world map
install.packages("BiocManager")
BiocManager::install("rnaturalearth")
library(rnaturalearth)
library()
install.packages("googleway","ggspatial","libwgeom","rnaturalearthdata")
ggplot(data=world)+
  geom_sf()

##now let's see overtime taking the date
datedf <- data.frame(Dates=names(table(all$DATE)), Availablemut=unname(table(all$DATE)))
datedf <- datedf[,-2]
head(datedf)
dim(datedf)
#datedf1 <- cbind(colsplit(string = datedf$Dates,pattern = "-",names = c("Year","Month","Day")),datedf$Availablemut.Freq)
transdremde <- data.frame(Dates=names(table(remde$DATE)), Resistant=unname(table(remde$DATE)))
transdlopi <- data.frame(Dates=names(table(lopi$DATE)), Resistant=unname(table(lopi$DATE)))
transd5Ffu <- data.frame(Dates=names(table(fivefu$DATE)), Resistant=unname(table(fivefu$DATE)))
#percentages of remde resistant
datedf$RemResistantMut <- transdremde$Resistant.Freq[match(datedf$Dates,transdremde$Dates)]
datedf$RemResistantMut[which(is.na(datedf$RemResistantMut)==TRUE)] <- 0
datedf$RemdePerc <- datedf$RemResistantMut/datedf$Availablemut.Freq*100
#percentages of lopi resistant
datedf$LopiResistantMut <- transdlopi$Resistant.Freq[match(datedf$Dates,transdlopi$Dates)]
datedf$LopiResistantMut[which(is.na(datedf$LopiResistantMut)==TRUE)] <- 0
datedf$LopiPerc <- datedf$LopiResistantMut/datedf$Availablemut.Freq*100
#percentages of 5ffu resistant
datedf$FiveFUResistantMut <- transd5Ffu$Resistant.Freq[match(datedf$Dates,transd5Ffu$Dates)]
datedf$FiveFUResistantMut[which(is.na(datedf$FiveFUResistantMut)==TRUE)] <- 0
datedf$FiveFUPerc <- datedf$FiveFUResistantMut/datedf$Availablemut.Freq*100
datedf

#plot it 
library(ggplot2)
library(reshape2)
head(Countrydf)
#melt
percdf <- datedf[,-c(2,3,5,7)]
perddf <- cbind(colsplit(string = datedf$Dates,pattern = "-",names = c("Year","Month","Day")),percdf[,c(2,3,4)])
#filter out the untraceable date

perddf$Month <- gsub("12","December",perddf$Month)
perddf$Month <- gsub("1","January",perddf$Month)
perddf$Month <- gsub("2","February",perddf$Month)
perddf$Month <- gsub("3","March",perddf$Month)
perddf$Month <- gsub("4","April",perddf$Month)
perddf$Month <- gsub("5","May",perddf$Month)
perddf$Month <- gsub("6","June",perddf$Month)
perddf$Month <- gsub("7","July",perddf$Month)
head(perddf)
perddf1 <- perddf[-which(is.na(perddf$Month==TRUE)),-c(1,3)]
perddf2 <- perddf

perddf2 <- perddf2[order(perddf2$Year,perddf2$Month,perddf2$Day),]
perddf2$pasted <- paste(perddf$Year,perddf$Month,perddf$Day, sep ="-")
perddf2 <- perddf2[-which(is.na(perddf2$Month)==TRUE),-c(1,2,3)]
head(perddf2)
mD <- melt(perddf1)
head(mD)
mdd <- mD[which(mD$value>0),]
mdd
unique(mD$Month)
mD$Months <- factor(mD$Month, levels=c("December","January","February","March","April","May","June","July"))
q <- ggplot(mD, aes(x=Months, y=value))+
  geom_col()+
  
  theme(
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  )+
  ylab("Mutation percentage over all detected mutations")+
  facet_grid(variable~.)
ggsave (q,file = "Mutations_overTime.pdf", width=8, height=10, units = "in", limitsize=FALSE)
ggsave(q,"Mutations_overTime.pdf", width=8, device = "pdf", height = 10)
ggplot(mD, aes(x=Month, y=value))+
  geom_col()+
  
  theme(
    axis.text.x = element_text(angle=45,vjust=1,hjust=1),
    xlab("Percentage of Mutations")
  )+
  facet_grid(variable~.)

