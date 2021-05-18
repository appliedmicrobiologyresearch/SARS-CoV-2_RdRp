#surrounding mutation exploring: the idea is to check the number of events in a give n position  
#loading all the NS mutations

setwd("/Volumes/AMBR$/Mari/COVID1/ResistancePaper/GSAID/COVGAP7d/")
all <- read.csv("All.mutations.filtered.tab", quote = "", sep= "\t", row.names = NULL, header = T, stringsAsFactors = F)
length(unique(all$SAMPLEID))
stg_mut <- all[which(all$ANNOTATION=="stop_gained" & all$RdRp_NSP12=="Yes"),]
stgt_samp <- unique(stg_mut$SAMPLEID)
allnst <- all[which(all$SAMPLEID %in% stgt_samp==F),]
length(unique(allnst$SAMPLEID))
length(stgt_samp)
allns <- allnst[which(allnst$ANNOTATION=="missense_variant"),]
length(allns$SAMPLEID[which(allns$Remdesivir_binding_affecting!="No")])
dim(all)
dim(allns)
head(allns)
table(allns$Remdesivir_binding_affecting)

#all the genomes showing at least one non-synonymous mutation in F480 pocket 
pos480 <- allns$SAMPLEID[which(allns$Remdesivir_binding_affecting=="Yes_pocket_F480" | 
                        allns$Remdesivir_binding_affecting=="Yes_pocket_codon_F480")]

#all the genomes showing at least one non-synonymous mutation in V557 pocket
pos557 <- allns$SAMPLEID[which(allns$Remdesivir_binding_affecting=="Yes_pocket_V557" | 
                        allns$Remdesivir_binding_affecting=="Yes_pocket_codon_V557")]

#all the genomes showing at least one mutation in both pockets
both <- pos557[match(pos480,pos557)]
both <- both[which(complete.cases(both)==T)] #NAs out

#now a genome label (which will be our read-out phenotype) --> it labels the genome, not the mutation
allns$GenomeLabelRemdesivir <- "Negative"
allns$GenomeLabelRemdesivir[which(allns$SAMPLEID%in%pos480==T)] <- "F480_Positive"
allns$GenomeLabelRemdesivir[which(allns$SAMPLEID%in%pos557==T)] <- "V557_Positive"
allns$GenomeLabelRemdesivir[which(allns$SAMPLEID%in%both==T)] <- "F480_V557_Both_Positive"
allns[which(allns$GenomeLabelRemdesivir=="F480_V557_Both_Positive"),] ##double check
#doublecheck that the label is in the right position
unique(allns$SAMPLEID[which(allns$GenomeLabelRemdesivir=="Negative")])

allns[which(allns$SAMPLEID=="Australia-NSW04-2020"),]
##now lets consider the lineages of any positive genome
#subset all the samples and mutations showing at least one positive
allpos <- allns[which(allns$GenomeLabelRemdesivir !="Negative"),]
allneg <- allns[which(allns$GenomeLabelRemdesivir =="Negative"),]
#for the glm fitting, gotta take the positions that are present in both, positive, and negative genomes, 
#and belong to the same lineage, we have to have a position 
#the grepping must go from positive to negative, as we are not interested in position present only in the negative but not in the positive, 
#we are only interested in the other way around
pospos <- unique(allpos$POS)
posneg <- unique(allneg$POS)
match(pospos,posneg) #positions shared between positive and negative,NA mean the ones only in the positive, 
#we have to treat them differently inserting a zero in their negative frequency

#doublecheck
pospos[20] #that should be in the neg
posneg[960] #shall be the same result
#negative doublecheck
pospos[40]
posneg[which(posneg==40)]
#now instead i take all the negative samples that belong to the same lineages of the positive
sort(pospos)
allneg_pos <- allneg[which(allneg$POS%in%pospos==T),]
lineages <- unique(allpos$PANGOLIN)  #define the lineages that show some positive
allneg_poslineage <- allneg_pos[which(allneg_pos$PANGOLIN%in%lineages),]

#now let's see how many mutations 

#count dataframe
statdf <- rbind(allpos,allneg_poslineage)
sort(unique(statdf$POS))
allpos
head(statdf)
length(unique(statdf$SAMPLEID))
statdf$Lab <- "Negative" 
#now that we have at least one mutation per phenotype, 
statdf$Lab <- "Negative" 
statdf$Lab[which(statdf$GenomeLabelRemdesivir != "Negative")] <- "Positive" 

library(tidyr)
library(dplyr)
p <- statdf %>%
  select(POS,PANGOLIN,Lab)%>%
  group_by(POS,PANGOLIN, Lab)%>%
  count()
p
#let's spread it so to add missing values, in here in fact the 0s are meaningful
g <- p %>% 
  spread(value = n, key=Lab,  fill=0)
g
#now we re-gathering again, to format it for glm models
p1 <- g %>% 
  gather(data=g, Phenotype,Negative, Positive)
p1
p1$revfreq <- length(unique(statdf$SAMPLEID))-p1$Phenotype
colnames(p1) <- c("Position", "Lineage", "Phenotype", "Frequency", "Reverse_Frequency")
p1$Position <- as.factor(p1$Position)
p1$Lineage <- as.factor(p1$Lineage)
p1$Phenotype <- as.factor(p1$Phenotype)

#doublecheck that EACH position has at least two entries, positive and negative, even if one of them has frequency 0
head(p1)
table(p1$Position, useNA = "ifany")<2
p1[which(p1$Position==320),]
allns[which(allns$POS==320),]

#to loop through it, we need to make sure that each subset has same position, same lineage, both phenotypes


#we do a subset to see how it works in small scale 
length(p1$Position)
#y <- rbind(p1[which(p1$Position==282 | p1$Position==285 | p1$Position==320),])
#tail(y)
#y[which(y$Position==282),]
r <- unique(p1$Position)
p1[which(p1$Position==320),]
d <- data.frame()
for (j in 1:length(r)){
  subset <- p1[which(p1$Position==r[j]),]
  if (length(unique(subset$Lineage))<2){
    i <- glm(cbind(Frequency,Reverse_Frequency)~Phenotype, data=subset, family = "binomial")
    lin <- "single_lineage: "
  }else{
    i <- glm(cbind(Frequency,Reverse_Frequency)~Phenotype+Lineage, data=subset, family = "binomial")
    lin <- "multilineage: "
  }
  lineages <- paste(unique(subset$Lineage), collapse= ", " )
  linlab <- paste(lin,lineages, sep ="")
  o <- summary(i)
  pval <- data.frame(Position=r[j], Lineage_prevalence=linlab,Feature=names(o$coefficients[,4][2]),p_value=unname(o$coefficients[,4])[2])
  d <- rbind(d,pval)
}
o$coefficients
d
d$adj_pvalBY <- p.adjust(d$p_value, n=length(d$p_value), method="BY")
d$SigBY <- ""
d$SigBY[which(d$adj_pvalBY<0.05)] <- "*"
d$SigBY[which(d$adj_pvalBY<0.01)] <- "**"
d$SigBY[which(d$adj_pvalBY<0.001)] <- "***"
d$adj_pvalBH <- p.adjust(d$p_value, n=length(d$p_value), method="BH")
d$SigBH <- ""
d$SigBH[which(d$adj_pvalBH<0.05)] <- "*"
d$SigBH[which(d$adj_pvalBH<0.01)] <- "**"
d$SigBH[which(d$adj_pvalBH<0.001)] <- "***"


sink("Summary_glm_Compensatory_mut_noST.txt")
print(o)
sink()
head(p1)
head(statdf)
p1Pos <- p1[which(p1$Phenotype=="Positive"),]
p1Neg <- p1[which(p1$Phenotype=="Negative"),]
p1Cpos <- p1Pos %>%
  group_by(Position)%>%
  summarise(Frequency = sum(Frequency))
p1Cneg <- p1Neg %>%
  group_by(Position)%>%
  summarise(Frequency = sum(Frequency))


d$OriginalPoscounts <- "None" 
d$OriginalPoscounts <- p1Cpos$Frequency[match(d$Position,p1Cpos$Position)] 
d$OriginalNegcounts <- "None"
d$OriginalNegcounts <- p1Cneg$Frequency[match(d$Position,p1Cneg$Position)] 
d$Position <- as.numeric(as.character(d$Position))
d$PosNegratio <- d$OriginalPoscounts/d$OriginalNegcounts
#first of, we consider only the significantly correlated mutations,
#then we establish if they are associated to negative or positive phenotype
#loophole
#d <- read.csv("Glm_Compensatory_mut_adjusted_noSTG.txt", header=T, row.names = NULL, sep = "\t", stringsAsFactors = F)
#head(d)
allsig <- d[which(grepl("\\*",d$SigBY)==T),]
hist(allsig$PosNegratio)
#define the percentiles of posneg ratio for all the mutations 
#--> highest percentile=associated with positive phenotype
#--> lowest percentile=associated with negative phenotype
allsig_q <- quantile(allsig$PosNegratio, probs= seq(0,1, 0.05), type=8)
#now we narrow down the mutations associated with positive pheno as the the one > 95th percentile
n95 <- unname(allsig_q[20])
n05 <- unname(allsig_q[2])
posassoc <- allsig[which(allsig$PosNegratio>n95),]
negassoc <- allsig[which(allsig$PosNegratio<n05),]
posassoc$Position
negassoc$Position
dim(d)
dim(allsig)
#isolate the mutations happening in rdrp domain
#rdrpd <- allsig[which(allsig$Position>13442 & allsig$Position < 16236),]
#define the percentiles of posneg ratio for the rdrp mutations
#rdrp_q <- quantile(rdrpd$PosNegratio, probs= seq(0,1, 0.05), type=8)
#nr95 <- unname(rdrp_q[20])
#nr05 <- unname(rdrp_q[2])
#posassocR <- rdrpd[which(rdrpd$PosNegratio>nr95),]
#negassocR <- rdrpd[which(rdrpd$PosNegratio<nr05),]
#allns[which(allns$POS==15857),]
#allns[which(allns$POS==16210),]
#small resuslts parenthesis
#length(unique(allns$SAMPLEID))
#unique(allns$Remdesivir_binding_affecting)
sort(unique(allns[which(allns$Remdesivir_binding_affecting!="No"),which(colnames(allns)=="DATE")]))
sort(unique(allns[,which(colnames(allns)=="DATE")]))
allns[which(allns$DATE=="2020-01-20"),]
gen <- read.csv("/Volumes/AMBR$/Mari/COVID1/ResistancePaper/GSAID/COVGAP7d/All.genomes.seqcount.tab", 
                quote="", header = T, sep = "\t", stringsAsFactors = F)
head(gen)

#the ones enriched in the negative samples
rdrpd[which(rdrpd$OriginalNegcounts>0 & rdrpd$adj_pvalBY<0.05),]
#the ones enriched in the positive samples
rdrpd[which(rdrpd$OriginalPoscounts>0 & rdrpd$adj_pvalBY<0.05),]

d[which(d$OriginalPoscounts>0 & d$adj_pvalBY < 0.05 ),]
dim(p1)
write.table(d, "Glm_Compensatory_mut_adjusted_noSTG.txt", quote = F, col.names = T, row.names = F, sep = "\t")
write.table(posassoc, "Glm_Compensatory_POSASSOC_adjusted_noSTG.txt", quote = F, col.names = T, row.names = F, sep = "\t")
write.table(negassoc, "Glm_Compensatory_NEGASSOC_adjusted_noSTG.txt", quote = F, col.names = T, row.names = F, sep = "\t")
samplefinder <- function(posit) {
  IDs <- allns$SAMPLEID[which(allns$POS==posit)]
  DF <- allns[which(allns$SAMPLEID %in% IDs==T),]
  return(DF)
}
posassoc$Position
slimdf <- data.frame()
for (i in 1:length(posassoc$Position)){
  df <- samplefinder(posassoc$Position[i])
  df$POSAssociatedMutation_GLM <- posassoc$Position[i]
  slimdf <- rbind(slimdf,df)
}

slimdf1 <- data.frame()
for (i in 1:length(negassoc$Position)){
  df1 <- samplefinder(negassoc$Position[i])
  df1$NEGAssociatedMutation_GLM <- negassoc$Position[i]
  slimdf1 <- rbind(slimdf1,df1)
}


slimdf1

write.table(slimdf, "Glm_Compensatory_POSASSOC_adjusted_noSTG__SAMPLES.txt", quote = F, col.names = T, row.names = F, sep = "\t")
write.table(slimdf1, "Glm_Compensatory_NEGASSOC_adjusted_noSTG__SAMPLES.txt", quote = F, col.names = T, row.names = F, sep = "\t")



#now make the pic
head(d)
library(ggplot2)
pl <- ggplot(allsig, aes(x=Position,y=PosNegratio ))+
  geom_col(position = "dodge", colour="black",size=0.6)+
  scale_x_continuous(breaks = c(0,2000,4000,6000, 8000, 10000, 12000,14000, 16000, 18000,20000, 22000, 24000,26000, 28000, 30000))+
  scale_y_continuous(breaks=c(0,0.025,0.05,0.075,0.1,0.125,0.15, 0.175))+
  geom_hline(yintercept=n95, linetype="dashed", size=0.2)+
  ylab("Positive/Negative ratio")+
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = n95, ymax = Inf, fill = "cornflowerblue", alpha = .5, color = NA)+
  geom_hline(yintercept=n05, linetype="dashed", size=0.2)+
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = n05, fill = "darkgreen", alpha = .5, color = NA)+
  theme(
    panel.background = element_blank(),
    panel.grid = element_line(size = 0.001, linetype = 'dashed',
                              colour = "ghostwhite"), 
    legend.key = element_rect (fill = "white"),
    legend.text = element_text (size = 20),
    plot.title = element_text(face = "bold", hjust = 1,vjust = 1, size = 30),
    legend.title = element_text (size = 20, face = "bold"),
    axis.text = element_blank(),#(size =20),
    axis.line = element_line(color= "black", size =0.6),
    axis.text.x = element_blank(),#(size = 18, angle = 45,vjust= 1,hjust=1),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size =22, face = "bold"),
    axis.title.x = element_blank(),
    legend.position="right"
  )
  #xlim(c(0,30000))
pl
head(allsig)
plo <- ggplot(allsig, aes(x=Position, y=adj_pvalBY))+
  geom_col(position = "dodge", colour="black",size=0.6)+
  scale_y_reverse()+
  scale_x_continuous(breaks = c(0,2000,4000,6000, 8000, 10000, 12000,14000, 16000, 18000,20000, 22000, 24000,26000, 28000, 30000))+
  xlab("Mutated positions in the Genome")+
  ylab("Adjusted pvalue -Benjamini Yekutieli-")+
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
    axis.text.x = element_text(size = 20, angle = 45,vjust= 1,hjust=1),
    axis.text.y = element_text(size = 20),
    axis.title = element_text(size =22, face = "bold"),
    legend.position="right"
  )
plo  

library(ggpubr)
library(ggExtra)

#now combine them
pt<- ggarrange(pl, plo, ncol=1, nrow=2, common.legend = T, legend = "top", align = "v")
ggsave (plot=pt, filename="CompensatoryMutationsNEW.pdf", width=20, height=12, device="pdf", units = "in", limitsize=FALSE)

