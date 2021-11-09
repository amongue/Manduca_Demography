###redoing Manduca demography for G3 reviewers using JHU as our baseline.

#Part 1. Structure.
#complication: Evanno method now thinks K=4 is most probable..let's see what this works out to be:
setwd("/Users/Andrew/Documents/Manduca_Demography/")
probs<-read.table("JHU_logs.txt",header=F)
colnames(probs)<-c("K","Prob")
boxplot((-1*probs$Prob)~probs$K,outline=F)
plot(x=probs$K,y=(-1*probs$Prob),pch=16,ylab="L",xlab="K")
#likelihoods don't stall out cleanly as in the Evanno paper, but 2->3 is still the biggest likelihood jump




probs<-read.table("JHU_no12.txt",header=F)
colnames(probs)<-c("K","Prob")
boxplot((-1*probs$Prob)~probs$K,outline=F)
plot(x=probs$K,y=(-1*probs$Prob),pch=16,ylab="L",xlab="K")



#Figure 1 
k39<- t(as.matrix(read.table("/Users/Andrew/Documents/Manduca_Demography/clumpak_2/9strap_K3.qopt")))
k32<- t(as.matrix(read.table("/Users/Andrew/Documents/Manduca_Demography/struct/2strap_K3_J.qopt")))
#labeling by sample
k39l<-k39
colnames(k39l)<-c("S32","S33","S34","S35","S36","S37","S38","S39","S40","S42","S44","S45","A36","A70",
"A71","A76","A78","A82","A84","A85","WK2","GK3","LK5","KC1")

#par(mfrow=c(1,1))
barplot(k39l,space=0,main="Tobacco hornworm population structure (K=3)",ylab="Admixture proportion",col=c("indianred3","royalblue3","khaki3"),
 las=1, cex.lab=1.1,cex.axis=1.5,cex.main=2,las=2)
abline(v=12,col="black",lwd=3.5)
abline(v=20,col="black",lwd=3.5)
mtext(side = 1, line = 3.5, at = 6, "NC", cex = 2)
mtext(side = 1, line = 3.5, at = 16, "AZ", cex = 2)
mtext(side = 1, line = 3.5, at = 22, "KS", cex = 2)
mtext(side = 3, line = 0.75, at = -0.6, "a", cex = 2)




#PCA plots
#Let's combine with the PCA data
library("RcppCNPy")
setwd("/Users/Andrew/Documents/Manduca_Demography/")
PCA_test_cov<-npyLoad("PCA_sel_JHU.cov.npy")
#barplot(PCA_test_cov)

tpc<-prcomp(PCA_test_cov)
colvec<-c(rep("indianred3",12),rep("khaki3",8),rep("royalblue3",4))
par(mfrow=c(1,2))
plot(tpc$rotation[,1],tpc$rotation[,2],xlab="PCA1",ylab="PCA2",main="Principle component clustering of samples",col=colvec,pch=16,cex=1.4,las=1, cex.lab=1.3,cex.axis=1.2,cex.main=1.4)
text(-0.04243051,-0.332,"KS",cex=1.8)
text(-0.118,0.1,"NC",cex=1.8)
text(0.2,-0.13,"AZ",cex=1.8)
text(0.3,0.18,"AZ",cex=1.8)
text(0.07,-0.32,"AZ",cex=1.8)
mtext(side = 3, line = 0.5, at = -0.2, "c", cex = 1.7)
#basically no change with the new reference, that's good!



PCA_test_cov<-npyLoad("PCA_sel_JHU_no12.cov.npy")
#barplot(PCA_test_cov)

tpc<-prcomp(PCA_test_cov)
colvec<-c(rep("indianred3",12),rep("khaki3",8),rep("royalblue3",4))
plot(tpc$rotation[,1],tpc$rotation[,2],xlab="PCA1",ylab="PCA2",main="Principle component clustering of samples
 excluding Chromosome 12",col=colvec,pch=16,cex=1.2,las=1, cex.lab=1.3,cex.axis=1.2,cex.main=1.4)
text(-0.158,0.1,"NC",cex=1.8)
text(0.2,0.13,"AZ",cex=1.8)
text(0.07,-0.32,"KS",cex=1.8)
mtext(side = 3, line = 0.5, at = -0.25, "c", cex = 1.7)
#removing the inversion only makes KS and AZ more distinct, this is good!


#
#read in our chromosome assignments
library("data.table")
setwd("/Users/Andrew/Documents/Manduca_Demography/")
mstrans<-read.csv("msexta_JHU_assembly_chr_anchor.csv",header=T,stringsAsFactors=F)
mstabr<-as.data.table(cbind(mstrans$Scaffold,mstrans$Synteny_liftover_Chr))
colnames(mstabr)<-c("Scaffold","Chr")

#short aside, but we know we want to subset out the Z and chr12 from our structure analyses
#msts<-mstabr[which(mstabr$Chr!="Z"),]
#msts<-msts[which(msts$Chr!="12"),]
#write.table(msts$Scaffold,"JHU_noZ_no12.txt",quote=F,row.names=F,col.names=F)

#read Fst data in
anf<-fread("AZ.NC.10kbfst.txt",stringsAsFactors=F, header=T)
kaf<-fread("KS.AZ.10kbfst.txt",stringsAsFactors=F,header=T)
knf<-fread("KS.NC.10kbfst.txt",stringsAsFactors=F, header=T)

#add chr info
anf<-merge(anf,mstabr,by="Scaffold",all.x=T)
anf[which(anf$Scaffold=="HiC_scaffold_31"),6]<-"2"
anf[which(anf$Scaffold=="HiC_scaffold_25"),6]<-"26"
anf$Chr<-factor(anf$Chr,levels=c("Z","2","3","4","5","6","7","8","9","10","11","12","13","14","15",
"16","17","18","19","20","21","22","23","24","25","26","27","28","0","NA"))

anfo<-anf[order(anf$Chr),]
kaf<-merge(kaf,mstabr,by="Scaffold",all.x=T)
kaf[which(kaf$Scaffold=="HiC_scaffold_31"),6]<-"2"
kaf[which(kaf$Scaffold=="HiC_scaffold_25"),6]<-"26"
kaf$Chr<-factor(kaf$Chr,levels=c("Z","2","3","4","5","6","7","8","9","10","11","12","13","14","15",
"16","17","18","19","20","21","22","23","24","25","26","27","28","0","NA"))

kafo<-kaf[order(kaf$Chr),]
knf<-merge(knf,mstabr,by="Scaffold",all.x=T)
knf[which(knf$Scaffold=="HiC_scaffold_31"),6]<-"2"
knf[which(knf$Scaffold=="HiC_scaffold_25"),6]<-"26"
knf$Chr<-factor(knf$Chr,levels=c("Z","2","3","4","5","6","7","8","9","10","11","12","13","14","15",
"16","17","18","19","20","21","22","23","24","25","26","27","28","0","NA"))
knfo<-knf[order(knf$Chr),]


#Genome-wide Fst

plot(1:nrow(anfo),anfo$Fst,col=rgb(255,0,0,max=255,alpha=15),ylim=c(0,1),xlab="",ylab=expression('F'[ST])
,pch=20, xaxt='n',las=1, main = "Population differentiation across the genome",cex.main=2,xlim=c(0,nrow(anfo)+14000))
points(1:nrow(kafo),kafo$Fst,col=rgb(0,0,255,max=255,alpha=15),pch=18)
points(1:nrow(knfo),knfo$Fst,col=rgb(0,0,0,max=255,alpha=15),pch=15)
segments(x0=0,x1=217379,y0=0.5,y1=0.5, lty=3)
#delimit chromosomes
scafbreaks<-match(unique(knfo$Chr),knfo$Chr)
#we have one too many breaks because both 0 and NA are unplaced
scafbreaks[30]<-417379
abline(v=scafbreaks)
scaflab<-scafbreaks
for(i in 1:length(unique(knfo$Chr))-1)
{scaflab[i]<-((scafbreaks[i+1]-scafbreaks[i])/2) + scafbreaks[i]}
text(c("Z","2","3","4","5","6","7","8","9","10","11","12","13","14","15",
"16","17","18","19","20","21","22","23","24","25","26","27","28","",""),x=scaflab,y=0.95)
text("unplaced",x=223991,y=0.95)
legend(218991,0.6,col=c("black","blue","red"),c("NC-KS","KS-AZ","NC-AZ"),pch=c(15,18,20))
mtext("Genomic location by chromosome",side=1,line=1,cex=1.6)

####Z chr Fst

knfoz<-knfo[which(knfo$Chr=="Z"),]
anfoz<-anfo[which(anfo$Chr=="Z"),]
kafoz<-kafo[which(kafo$Chr=="Z"),]

plot(1:nrow(anfoz),anfoz$Fst,col="red",ylim=c(0,1),xlab="Genomic location by chromosome",ylab=expression('F'[ST])
,pch=20, xaxt='n',las=1, main = "Evidence of a single inversion on Z chromosome in Arizona",cex.main=2)
points(1:nrow(kafoz),kafoz$Fst,col="blue",pch=18)
points(1:nrow(knfoz),knfoz$Fst,col="black",pch=15)
scafbreaks<-match(unique(knfoz$Scaffold),knfoz$Scaffold)
abline(v=scafbreaks)

####chr 12 Fst

knfoz<-knfo[which(knfo$Chr=="12"),]
anfoz<-anfo[which(anfo$Chr=="12"),]
kafoz<-kafo[which(kafo$Chr=="12"),]

plot(1:nrow(anfoz),anfoz$Fst,col="red",ylim=c(0,1),xlab="Genomic location by chromosome",ylab=expression('F'[ST])
,pch=20, xaxt='n',las=1, main = "Evidence of a single inversion on chromosome 12 in Arizona",cex.main=2)
points(1:nrow(kafoz),kafoz$Fst,col="blue",pch=18)
points(1:nrow(knfoz),knfoz$Fst,col="black",pch=15)
scafbreaks<-match(unique(knfoz$Scaffold),knfoz$Scaffold)
#abline(v=scafbreaks)
abline(v=4600)


#new analysis: is the chr12 inversion enriched for any class of genes
#we're doing this with the old assembly, because the annotations are probably better and the RNAseq is all done it as well
inv12<-c("scaffold00042","scaffold00095","scaffold00119","scaffold00158","scaffold00203","scaffold00240","scaffold00264","scaffold00300","scaffold00335",
"scaffold00460","scaffold00466","scaffold00571","scaffold00616","scaffold00633","scaffold00790","scaffold00888","scaffold00921","scaffold00966",
"scaffold01138","scaffold01198","scaffold01282","scaffold01317","scaffold01366","scaffold01482","scaffold01552")
zkmanmaster<-read.csv("ZENKOKU_m_sexta_genome_mat.csv",header=T, stringsAsFactors=F)
expressmaster<-read.csv("Msexta_pan_RNA_specificity_update.csv",header=T, stringsAsFactors=F)

omni<-merge(zkmanmaster,expressmaster,by="Gene",all.x=T,all.y=T)
#this has all our stage and tissue SPMs
#but we need to bin things if we're going to do independence tests
omni[which(omni$SPM.Larva>0.70),22]<-"Larva"
omni[which(omni$SPM.Pupa>0.70),22]<-"Pupa"
omni[which(omni$SPM.Adult>0.70),22]<-"Adult"
omni[which(omni$SPM.Adult>0.20 & omni$SPM.Pupa>0.20 & omni$SPM.Larva>0.20),22]<-"Lifelong"
omni[which(omni$SPM.Head>0.70),26]<-"Head"
omni[which(omni$SPM.Fat>0.70),26]<-"Fat"
omni[which(omni$SPM.Midgut>0.70),26]<-"Midgut"
omni[which(omni$SPM.Malphig>0.70),26]<-"Malphig"
omni[which(omni$SPM.Muscle>0.70),26]<-"Muscle"
omni[which(omni$SPM.Testes>0.70),26]<-"Testes"
omni[which(omni$SPM.Ovary>0.70),26]<-"Ovary"
omni[which(omni$SPM.Anten>0.70),26]<-"Antennae"
omni[which(omni$SPM.Male>0.70),19]<-"Male"
omni[which(omni$SPM.Female>0.70),19]<-"Female"
omni[which(omni$SPM.Male<0.70&omni$SPM.Male>0.30),19]<-"Unbiased"


omni12<-omni[which(omni$Scaffold %in% inv12),]
omninot<-omni[which(!(omni$Scaffold %in% inv12)),]

table(omninot$Stage.Specificity)
table(omni12$Stage.Specificity)
invtest<-rbind(table(omninot$Stage.Specificity),table(omni12$Stage.Specificity))
chisq.test(invtest)
#X-squared = 2.0597, df = 3, p-value = 0.5601

#the same but for tissue?

table(omninot$Tissue.Specificity)
table(omni12$Tissue.Specificity)
invtest<-rbind(table(omninot$Tissue.Specificity),table(omni12$Tissue.Specificity))
chisq.test(invtest)
#X-squared = 7.4478, df = 7, p-value = 0.3838
#nothing doing...
#one last time for sex-bias
table(omninot$Sex.Specificity)
table(omni12$Sex.Specificity)
invtest<-rbind(table(omninot$Sex.Specificity),table(omni12$Sex.Specificity))
chisq.test(invtest)

#The Z

 
 #contiguous Z
azz<-as.data.frame(fread("/Users/Andrew/Documents/Manduca_Demography/AZ.JHU.Z.freqs.abr.txt"))
azzf<-azz[which(azz$n_chr=="16"),]
azzf$F2<-((azzf$F2*16)-1)/14
azzf<-azzf[which(azzf$F2>.2),]
ncz<-as.data.frame(fread("/Users/Andrew/Documents/Manduca_Demography/NC.JHU.Z.freqs.abr.txt"))
nczf<-ncz[which(ncz$n_chr=="24"),]
nczf<-nczf[which(nczf$F2>0.2),]

 anfoz<-anfo[which(anfo$Scaffold=="HiC_scaffold_19"),]
knfoz<-knfo[which(knfo$Scaffold=="HiC_scaffold_19"),]
kafoz<-kafo[which(kafo$Scaffold=="HiC_scaffold_19"),]

ncld<-as.data.frame(fread("/Users/Andrew/Documents/Manduca_Demography/Msex_50bpJ_NC.Z.LD.txt"))
colnames(ncld)<-c("Scaffold","Start","Stop","n","r2","Dist")

azld<-as.data.frame(fread("/Users/Andrew/Documents/Manduca_Demography/Msex_50bpJ_AZ.Z.LD.txt"))
colnames(azld)<-c("Scaffold","Start","Stop","n","r2","Dist")

par(mfrow=c(3,1))
#B L T R
par(mai=c(0.05,1,1.2,0.05))
plot(anfoz$midPos,anfoz$Fst,col="red",ylim=c(0,1),xlab="",ylab=expression('F'[ST])
,pch=20, xaxt='n',las=1, main = "A putative Arizona-specific inversion on the Z chromosome",cex.main=2,type='l')
points(kafoz$midPos,kafoz$Fst,col="blue",pch=18,type='l')
points(knfoz$midPos,knfoz$Fst,col="black",pch=15,type='l')
#abline(h=0.5, lty=3)
abline(v=c( 13706400, 14710588),lwd=1.5,lty=3)
legend(20000000, 0.95,col=c("black","blue","red"),c("NC-KS","KS-AZ","NC-AZ"),pch=c(15,15,15))
par(mai=c(1.2,1,0.05,0.05))
plot(azzf$POS,azzf$F2,pch=16,col=rgb(255, 0, 0, max = 255, alpha = 1),las=1,ylab="Allele frequency",xaxt='n',
xlab="Z chromosome (HiC_scaffold_19) position (Mb)",ylim=c(0.2,1))
points(nczf$POS,nczf$F2,pch=16,col=rgb(0, 0, 0, max = 255, alpha = 1))
abline(v=c( 13706400, 14710588),lwd=1.5,lty=3)
axis(1, at=c(1000000,2000000,3000000,4000000,5000000,6000000,7000000,8000000,9000000,10000000,11000000,12000000,13000000,14000000
,15000000,16000000,17000000,18000000,19000000,20000000,21000000), labels=seq(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21)))
legend(20000000,0.9,c("NC","AZ"),pch=16,col=c("black","red"))

plot(ncld$Start,ncld$r2,col=rgb(0, 0, 0, max = 255, alpha = 1))
points(azld$Start,azld$r2,col=rgb(255, 0, 0, max = 255, alpha = 1))

#per the fast Z paper
X-squared = 1.3172, df = 2, p-value = 0.5176

#scaff460:
#          Female     Male Unbiased 
##whole Z      55      295      177
#scaff460       1        3        3
#scaff22inv     2       12       21
#z-minus test  52      280      153
#z-lump+        3       15       24


zp<-c(52,280,153)
zm<-c(3,15,24)
ztest<-rbind(zp,zm)

chisq.test(ztest)
#data:  ztest
#X-squared = 11.362, df = 2, p-value = 0.00341

zp<-c(53,283,156)
zm<-c(2,12,21)
ztest<-rbind(zp,zm)

chisq.test(ztest)
#data:  ztest
#X-squared = 11.74, df = 2, p-value = 0.002823


omni2[which(omni2$Scaffold=="scaffold00022"),c(1,35)]
#msex2.01413 should be the last we include, so:
#          Gene Sex.Bias
#Msex2.01379   Female
#Msex2.01380 Unbiased
#Msex2.01381 Unbiased
#Msex2.01382 Unbiased
#Msex2.01383 Unbiased
#Msex2.01384 Unbiased
#Msex2.01385 Unbiased
#Msex2.01386   Female
#Msex2.01387 Unbiased
#Msex2.01388 Unbiased
#Msex2.01389 Unbiased
#Msex2.01390 Unbiased
#Msex2.01391 Unbiased
#Msex2.01392 Unbiased
#Msex2.01393 Unbiased
#Msex2.01394 Unbiased
#Msex2.01395     Male
#Msex2.01396     Male
#Msex2.01397 Unbiased
#Msex2.01398 Unbiased
#Msex2.01399 Unbiased
#Msex2.01400     Male
#Msex2.01401     Male
#Msex2.01402     Male
#Msex2.01403 Unbiased
#Msex2.01404     Male
#Msex2.01405     Male
#Msex2.01406 Unbiased
#Msex2.01407     Male
#Msex2.01408 Unbiased
#Msex2.01409 Unbiased
#Msex2.01410     Male
#Msex2.01411     Male
#Msex2.01412     Male
#Msex2.01413     Male



#not
#Msex2.01414     Male
#Msex2.01415 Unbiased
#Msex2.01416 Unbiased
#Msex2.01417     Male
#Msex2.01418 Unbiased
#Msex2.01419 Unbiased
#Msex2.01420     Male
#Msex2.01421     Male
#Msex2.01422 Unbiased
#Msex2.01423 Unbiased
#Msex2.01424 Unbiased
#Msex2.01425     Male

#plotting the Z inversion:

fst10kb<-fread("/Users/Andrew/Documents/Manduca_Demography/3_pop_10kb_slidingwindow.txt",header=T, sep="\t",stringsAsFactors=F)
dumx<-seq(1,33040,1)

#let's try to reorder things to be vaguely in chromosomal order
colnames(fst10kb)[2]<-"Scaffold"
zkmanmaster<-read.csv("Manduca_genome_matrix_for_pub.csv",header=T, stringsAsFactors=F)
zksub<-as.data.frame(cbind(zkmanmaster$Scaffold,zkmanmaster$Chr,zkmanmaster$ZorA))
zksub<-unique(zksub[,1:3])
zksub[,2]<-as.numeric(as.character(zksub[,2]))
colnames(zksub)<-c("Scaffold","Chr","ZorA")
zksub[which(zksub$Chr=="1"),2]<-"01"
zksub[which(zksub$ZorA=="Z"),2]<-"01"
zksub[which(zksub$Chr=="2"),2]<-"02"
zksub[which(zksub$Chr=="3"),2]<-"03"
zksub[which(zksub$Chr=="4"),2]<-"04"
zksub[which(zksub$Chr=="5"),2]<-"05"
zksub[which(zksub$Chr=="6"),2]<-"06"
zksub[which(zksub$Chr=="7"),2]<-"07"
zksub[which(zksub$Chr=="8"),2]<-"08"
zksub[which(zksub$Chr=="9"),2]<-"09"
#using the hi_C assembly to resolve some of the unplaced/misplaced
#scaffold00460
#scaffold00419 
#scaffold01602
#scaffold01636
zksub[which(zksub$Scaffold=="scaffold00460"),2]<-"12"
zksub[which(zksub$Scaffold=="scaffold00419"),2]<-"19"
zksub[which(zksub$Scaffold=="scaffold01602"),2]<-"11"
zksub[which(zksub$Scaffold=="scaffold01636"),2]<-"13"
#and more
zksub[which(zksub$Scaffold=="scaffold01482"),2]<-"12"
zksub[which(zksub$Scaffold=="scaffold01563"),2]<-"12"
zksub[which(zksub$Scaffold=="scaffold01354"),2]<-"22"
zksub[which(zksub$Scaffold=="scaffold01731"),2]<-"24"
nonation<-fst10kb[which(fst10kb$Chr=="0"),]
head(nonation[order(-nonation$Fst.KS.AZ),],20)
#scaffold01482
#scaffold01563
#scaffold01354
#scaffold00546
#scaffold01731

fst10kb<-merge(fst10kb,zksub, by = "Scaffold",all.x=T)

reord<- fst10kb[order(fst10kb$Chr),]
#reord<-fst10kb[with(fst10kb, order(Chr)),]

fst10kb<-reord


scaf22<-fst10kb[which(fst10kb$Scaffold=="scaffold00022"),]

zfreqs<-fread("/Users/Andrew/Documents/Manduca_Demography/az.22genes.txt",header=T)
zfreqs$AF<-as.numeric(zfreqs$AF)
zfreqs$POS<-as.numeric(zfreqs$POS)
#we need to correct for hemizygosity
#AF assumes 16 alleles...so AF*(16/14)
zfreqs$AF<-((zfreqs$AF*16)-1)/14
#of course this messes up fixed alts
zfreqs[which(zfreqs$AF>1),5]<-1

par(mfrow=c(2,1))
plot(scaf22$midPos,scaf22$Fst.NC.AZ ,col="red",type="l",
ylim=c(0,1),ylab="",lty=3,xlab="Scaffold00022 position (bp)", las=1, cex.lab=1.2,
main="")
mtext(side=2, line = 2.3, at = 0.5, expression("F"["ST"]), cex=1.35)
lines(scaf22$midPos,scaf22$Fst.KS.AZ ,col="blue",lty=2)
lines(scaf22$midPos,scaf22$Fst.NC.KS ,col="black",lty=1)
mtext(side=3, line =1, at =-296000, "Differentiated portions of Z chromosome", cex=2.3)
legend(1385150, 0.89, legend=c("NC", "KS","AZ"),
       col=c("black", "blue", "red"), lty=1:3)

plot(zfreqs$POS,zfreqs$AF,col=rgb(255, 0, 0, max = 255, alpha = 1),ylim=c(0,1))
abline(h=13/14)


#chr12 
ld12<-as.data.frame(fread("/Users/Andrew/Documents/Manduca_Demography/LD_1kb_chr12.geno.ld",header=T,stringsAsFactors=F))




#Chr 12
#okay the scaffold with the stop-loss (00381) aligns to:
#                JHU
#1:  HiC_scaffold_12
#2:  HiC_scaffold_30
#3: HiC_scaffold_600

#but only 600 and 30 are syntenic to 13...

 
ncld<-as.data.frame(fread("/Users/Andrew/Documents/Manduca_Demography/Msex_50bpJ_NC.chr12LD.txt"))
colnames(ncld)<-c("Scaffold","Start","Stop","n","r2","Dist")

azld<-as.data.frame(fread("/Users/Andrew/Documents/Manduca_Demography/Msex_50bpJ_AZ.chr12LD.txt"))
colnames(azld)<-c("Scaffold","Start","Stop","n","r2","Dist")

anfo13<-anfo[which(anfo$Scaffold=="HiC_scaffold_12"),]
knfo13<-knfo[which(knfo$Scaffold=="HiC_scaffold_12"),]
kafo13<-kafo[which(kafo$Scaffold=="HiC_scaffold_12"),]
par(mfrow=c(2,1))
plot(1:nrow(anfo13),anfo13$Fst,col="red",ylim=c(0,1),xlab="Genomic location by chromosome",ylab=expression('F'[ST])
,pch=20, xaxt='n',las=1, main = "Population differentiation across the genome",cex.main=2,type='l')
points(1:nrow(kafo13),kafo13$Fst,col="blue",pch=18,type='l')
points(1:nrow(knfo13),knfo13$Fst,col="black",pch=15,type='l')
abline(h=0.5, lty=3)
abline(v=c(5414046/2000, 13233866/2000),lwd=2)
 
plot(ncld$Start,ncld$r2,col=rgb(0, 0, 0, max = 255, alpha = 1))
points(azld$Start,azld$r2,col=rgb(255, 0, 0, max = 255, alpha = 1))

#alternate view:

#left window
l12<-as.data.frame(azld[which(azld$Stop<4000000),])
i12<-azld[which(azld$Stop>7000000 & azld$Stop<11000000),]
r12<-azld[which(azld$Stop>14000000),]

i12n<-ncld[which(ncld$Stop>7000000 & ncld$Stop<11000000),]

par(mfrow=c(2,2))
boxplot(r2~Dist,data=l12,main="Non-inverted left",outline=F)
boxplot(r2~Dist,data=i12,main="Inversion",outline=F)
boxplot(r2~Dist,data=r12,main="Non-inverted right",outline=F)
boxplot(r2~Dist,data=i12n,main="Inversion but NC",outline=F)

#let's see if it works with shared snps

sld<-merge(ncld,azld,by=c("Scaffold","Start","Stop"))
r2diff<-sld$r2.y - sld$r2.x

plot(sld$Stop,r2diff,col=rgb(0, 0, 0, max = 255, alpha = 5), xlab="Chromosome position",ylab="AZ - NC r2")
abline(v=c(5414046, 13233866),lty=3)

plot(i12$Stop,i12$r2,col=rgb(255, 0, 0, max = 255, alpha = 5))
points(i12n$Stop,i12n$r2,col=rgb(0, 0, 0, max = 255, alpha = 5))


####the Z LD

ncld<-as.data.frame(fread("/Users/Andrew/Documents/Manduca_Demography/Msex_50bpJ_NC.Z.LD.txt"))
colnames(ncld)<-c("Scaffold","Start","Stop","n","r2","Dist")

azld<-as.data.frame(fread("/Users/Andrew/Documents/Manduca_Demography/Msex_50bpJ_AZ.Z.LD.txt"))

colnames(azld)<-c("Scaffold","Start","Stop","n","r2","Dist")
sld<-merge(ncld,azld,by=c("Scaffold","Start","Stop"))
r2diff<-sld$r2.y - sld$r2.x

plot(sld$Stop,r2diff,col=rgb(0, 0, 0, max = 255, alpha = 5), xlab="Chromosome position",ylab="AZ - NC r2")
abline(v=c(5414046, 13233866),lty=3)



#plotting Pi for AA, AB, BB

aa<-as.data.frame(fread("/Users/Andrew/Documents/Manduca_Demography/AZ_AA_pi.sites.txt"))
ab<-as.data.frame(fread("/Users/Andrew/Documents/Manduca_Demography/AZ_AB_pi.sites.txt"))
bb<-as.data.frame(fread("/Users/Andrew/Documents/Manduca_Demography/NC_BB_pi.sites.txt"))

library("zoo")
aam<-zoo(aa$PI)
abm<-zoo(ab$PI)
bbm<-zoo(bb$PI)

aamm<-rollapply(aam, width = 600, by = 150, FUN = mean, align = "left")
abmm<-rollapply(abm, width = 600, by = 150, FUN = mean, align = "left")
bbmm<-rollapply(bbm, width = 600, by = 150, FUN = mean, align = "left")




par(mfrow=c(1,1))
plot(1:length(aamm),aamm,col="blue",pch=16,cex=0.8, ylim=c(0,0.6),ylab=expression(pi),xlab="Position along chromosome", xaxt='n')
points(1:length(abmm),abmm,col="orange",pch=16,cex=0.8)
points(1:length(bbmm),bbmm,col="gray",pch=16,cex=0.8)


par(mfrow=c(3,1))
par(mai=c(0,0.5,0.5,0.05))
plot(aa$POS,aa$PI,pch=16,col=rgb(0, 0, 0, max = 255, alpha = 1),main="AA",ylab=expression(pi),las=1,xlim=c(4000000,15000000))
abline(v=c(5414046, 13233866),lty=3)
rect(xleft=9257027, ybottom=0, xright= 9340442,ytop=1, density = 100, col = rgb(255, 0, 0, max = 255, alpha = 50))
plot(ab$POS,ab$PI,pch=16,col=rgb(0, 0, 0, max = 255, alpha = 1),main="AB",ylab=expression(pi),las=1,xlim=c(4000000,15000000))
abline(v=c(5414046, 13233866),lty=3)
par(mai=c(0.4,0.5,0.5,0.05))
rect(xleft=9257027, ybottom=0, xright= 9340442,ytop=1, density = 100, col = rgb(255, 0, 0, max = 255, alpha = 50))
plot(bb$POS,bb$PI,pch=16,col=rgb(0, 0, 0, max = 255, alpha = 1),main="BB",ylab=expression(pi),las=1,xlim=c(4000000,15000000))
abline(v=c(5414046, 13233866),lty=3)
rect(xleft=9257027, ybottom=0, xright= 9340442,ytop=1, density = 100, col = rgb(255, 0, 0, max = 255, alpha = 50))

#pseudogene zoom in

par(mfrow=c(3,1))
par(mai=c(0,0.5,0.5,0.05))
plot(aa$POS,aa$PI,pch=16,col=rgb(0, 0, 0, max = 255, alpha = 30),main="AA",ylab=expression(pi),las=1,xlim=c(8000000,10000000))
abline(v=c(5414046, 13233866),lty=3)
rect(xleft=9257027, ybottom=0, xright= 9340442,ytop=1, density = 100, col = rgb(255, 0, 0, max = 255, alpha = 50))
plot(ab$POS,ab$PI,pch=16,col=rgb(0, 0, 0, max = 255, alpha = 30),main="AB",ylab=expression(pi),las=1,xlim=c(8000000,10000000))
abline(v=c(5414046, 13233866),lty=3)
par(mai=c(0.4,0.5,0.5,0.05))
rect(xleft=9257027, ybottom=0, xright= 9340442,ytop=1, density = 100, col = rgb(255, 0, 0, max = 255, alpha = 50))
plot(bb$POS,bb$PI,pch=16,col=rgb(0, 0, 0, max = 255, alpha = 30),main="BB",ylab=expression(pi),las=1,xlim=c(8000000,10000000))
abline(v=c(5414046, 13233866),lty=3)
rect(xleft=9257027, ybottom=0, xright= 9340442,ytop=1, density = 100, col = rgb(255, 0, 0, max = 255, alpha = 50))



###same logic for the Z
aa<-as.data.frame(fread("/Users/Andrew/Documents/Manduca_Demography/AZ_AA_Z.pi.txt"))
bb<-as.data.frame(fread("/Users/Andrew/Documents/Manduca_Demography/NC_BB_Z.pi.txt"))
par(mfrow=c(2,1))
par(mai=c(0,0.5,0.5,0.05))
plot(aa$POS,aa$PI,pch=16,col=rgb(0, 0, 0, max = 255, alpha = 1),main="AA",ylab=expression(pi),las=1,xaxt='n', xlim=c(12000000,16000000))
abline(v=c( 13706400, 14710588),lwd=1.5,lty=3)
par(mai=c(0.4,0.5,0.5,0.05))
plot(bb$POS,bb$PI,pch=16,col=rgb(0, 0, 0, max = 255, alpha = 1),main="BB",ylab=expression(pi),las=1,xlim=c(12000000,16000000))
abline(v=c( 13706400, 14710588),lwd=1.5,lty=3)


#Chr19
af19<-as.data.frame(fread("AZ.19.frq.txt"))
af19n<-as.data.frame(fread("NC.J.19.abr.txt"))
cov19<-as.data.frame(fread("A70J19cov.txt"))
colnames(cov19)<-c("Scaff",'S1"',"S2","cov")


anfo13<-anfo[which(anfo$Scaffold=="HiC_scaffold_21"),]
knfo13<-knfo[which(knfo$Scaffold=="HiC_scaffold_21"),]
kafo13<-kafo[which(kafo$Scaffold=="HiC_scaffold_21"),]

#Figure S6
par(mfrow=c(1,1))
plot(anfo13$midPos,anfo13$Fst,col="red",ylim=c(0,1),xlab="Position along chromosome (bp)",ylab=expression('F'[ST])
,pch=20,las=1, main = "Population differentiation across chromosome 19",cex.main=2,type='l')
#,xlim=c(500000,2000000)
points(kafo13$midPos,kafo13$Fst,col="blue",pch=18,type='l')
points(knfo13$midPos,knfo13$Fst,col="black",pch=15,type='l')
abline(h=0.5, lty=3)
#abline(v=c(695215,1045369))

#looks like two small inversions near each other starting around 650kb to 900kb or 1.07Mb and another from 1.3Mb to 1.75Mb maybe? but no interesting genes within


#Chr 12 pseudogene

#
scaf381<-fst10kb[which(fst10kb$chr=="scaffold00381"),]
layout(matrix(c(1,1,1,2,3,4),2,3,byrow=T))
plot(scaf381$midPos,scaf381$Fst.NC.AZ ,col="red",type="l",
ylim=c(0,1),ylab="",lty=3,xlab="Scaffold00381 position (bp)", las=1, cex.lab=1.2, cex.main=2.3,
main="A Misassembly reveals a gene of interest on Chromosome 12")
mtext(side=2, line = 2.3, at = 0.5, expression("F"["ST"]), cex=1.35)
lines(scaf381$midPos,scaf381$Fst.KS.AZ ,col="blue",lty=2)
lines(scaf381$midPos,scaf381$Fst.NC.KS ,col="black",lty=1)
rect(xleft=245015,ybottom=-10.1,xright=248158,ytop=91.1,col = rgb(0.5,0.5,0.5,1/3))
#add synteny information here:
rect(xleft=0,ybottom=-10.1,xright=223392,ytop=0.4,col = rgb(0.5,0.5,0.5,1/4))
rect(xleft=223392,ybottom=-10.1,xright=348158,ytop=0.4,col = rgb(1,0,0,1/4))
text(100000, 0.25, "Syntenic to Chr13")
text(264000, 0.25, "Syntenic to Chr12")
legend(38150, .88, legend=c("NC-KS", "KS-AZ","NC-AZ"),
       col=c("black", "blue", "red"), lty=1:3)
text(227000,0.8, "Msex2.10493->", cex=1.5)
mtext(side =2, line =1,"a", at=1.15, cex=1.8, las=1)

plot(-1,-1, xlim=c(0,12),ylim=c(0,12),xaxt='n',yaxt='n',ylab="",xlab="")
abline(v=6)	
abline(h=4)
abline(h=8)
mtext(side=3, line = 2.5, at = 3, "Non-synonymous", cex=0.9)
mtext(side=3, line = 1, at = 3, "Polymorphisms",cex=0.9)
mtext(side=3, line = 2.5, at = 9, "Synonymous", cex=0.9)
mtext(side=3, line = 1, at = 9, "Polymorphisms",cex=0.9)
mtext(side=2, line = 1, at = 10, "NC", las=1)	
mtext(side=2, line = 1, at = 6, "KS", las=1)	
mtext(side=2, line = 1, at = 2, "AZ", las=1)	
text(3, 10, "2", cex =3)		
text(9, 10, "1", cex =3)		
text(3, 6, "1", cex =3)		
text(9, 6, "0", cex =3)	
text(3, 2, "27", cex =3)		
text(9, 2, "34", cex =3)	
mtext(side =2, line =1,"b", at=13.8, cex=1.8, las=1)

library("data.table")


azef<-fread("/Users/Andrew/Documents/Manduca_Demography/msextafilt_az.genes.txt",header=T,stringsAsFactors=F)
ksef<-fread("/Users/Andrew/Documents/Manduca_Demography/msextafilt_ks.genes.txt",header=T,stringsAsFactors=F)
az_gene<-azef[which(azef$GeneName=="Msex2.10493"),]
ks_gene<-ksef[which(ksef$GeneName=="Msex2.10493"),]
az_persnp<-fread("/Users/Andrew/Documents/Manduca_Demography/gene_10493_eff_az.txt",header=T,stringsAsFactors=F)

plot(az_persnp$POS,(az_persnp$AF),ylab="Alt allele frequency",xlab="Scaffold position",
xlim=c(245000,248180),las=1,cex.lab=1.2, pch=16, main="Msex2.10493 SNPs in Arizona",cex.main=1.3)
abline(v=247732, col="red",lty=3)
abline(v=245015)
abline(v=248158)
rect(xleft=245015,ybottom=-0.1,xright=245795 ,ytop=1.1,col = rgb(0.5,0.5,0.5,1/4))
rect(xleft=245955,ybottom=-0.1,xright=246175 ,ytop=1.1,col = rgb(0.5,0.5,0.5,1/4))
rect(xleft=246704,ybottom=-0.1,xright=246818 ,ytop=1.1,col = rgb(0.5,0.5,0.5,1/4))
rect(xleft=247555,ybottom=-0.1,xright=247874 ,ytop=1.1,col = rgb(0.5,0.5,0.5,1/4))
mtext(side =2, line =1,"c", at=1.15, cex=1.8, las=1)

cov381<-fread("/Users/Andrew/Documents/Manduca_Demography/scaf381_raws.txt",header=T, sep="\t",stringsAsFactors=F)
colnames(cov381)<-c("Scaffold","Start","End","S32", "S33", "S34", "S35", "S36", "S37", "S38", "S39", "S40", "S42", "S44", "S45", "Q6", "A36", "A70", "A71", "A76", "A78", "A82", "A84", "A85", "WK2", "GK3", "LK5", "KC1")

attach(cov381)
nccov<-(S32+S33+S33+S34+S35+S36+S37+S38+S39+S40+S42+S44+S45)/12
azcov<-(A36 + A70 + A71 + A76 +A78 + A82 + A84 +A85)/8
kscov<-(WK2 + GK3 + LK5 + KC1)/4

plot((cov381$Start+25),nccov,col="black",xlim=c(245000,248180),ylim=c(0,72),ylab="Depth of coverage",xlab="Scaffold position",pch=16,main="Per population coverage of Msex2.10493", type="l",lwd=2,las=1,cex.lab=1.2, cex.main=1.3)
lines((cov381$Start+25),kscov,col="blue", lwd=2, lty=2)
lines((cov381$Start+25),azcov,col="red", lwd=2, lty=3)
#rect(xleft=245015,ybottom=-10.1,xright=248158,ytop=91.1,col = rgb(0.5,0.5,0.5,1/4))
abline(v=247732, col="red",lty=3)
abline(v=245015)
abline(v=248158)
points(az_persnp$POS,y=rep(0,length(az_persnp$POS)),pch=16)
legend(245150, 73, legend=c("NC", "KS","AZ"),
       col=c("black", "blue", "red"), lty=1:3)
mtext(side =2, line =1,"d", at=82, cex=1.8, las=1)

#end pseudogene