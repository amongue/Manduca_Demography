#chr12LD...can we plot it?
library("data.table")
setwd("/Users/Andrew/Documents/Manduca_Demography/")
a12ld<-as.data.table(fread("LD_1kb_chr12abr.txt",header=F))
colnames(a12ld)<-c("Start","Stop","LD")
#what if we subset to the longest LD comps?
a12lab<-a12ld[!duplicated(a12ld$Stop, fromLast =T),]



n12ld<-as.data.table(fread("LD_1kb_chr12abrnc.txt",header=F))
colnames(n12ld)<-c("Start","Stop","LD")
n12lab<-n12ld[!duplicated(n12ld$Stop, fromLast =T),]

k12ld<-as.data.table(fread("LD_1kb_chr12abrks.txt",header=F))
colnames(k12ld)<-c("Start","Stop","LD")
k12lab<-k12ld[!duplicated(k12ld$Stop, fromLast =T),]

library('zoo')

aldm<-rollmedian(na.omit(a12lab$LD),9999)
nldm<-rollmedian(na.omit(n12lab$LD),9999)
kldm<-rollmedian(na.omit(k12lab$LD),9999)

nx<-seq(1:length(nldm))*(length(aldm)/length(nldm))
kx<-seq(1:length(kldm))*(length(aldm)/length(kldm))


plot(seq(1:length(aldm)),aldm,ylim=c(0,1),type='l',lwd=4, col="red", las=1,ylab="Linkage disequilibrium",main="Linkage across chromosome 12",xaxt='n',xlab="Position along chromosome")
points(nx,nldm,col='black',type='l',lwd=4)
points(kx,kldm,col='blue',type='l',lwd=4)  
legend(2, 0.98, c("NC","KS","AZ"),col=c("black","blue","red"),pch=15)

abline(v=c(5414046/17279714*length(aldm), 13233866/17279714*length(aldm)),lwd=2,lty=3)

length(aldm)
length(nldm)
length(kldm)






####and pop gen:

#########now manduca
msld<-as.data.frame(fread("/Users/Andrew/Documents/Manduca_Demography/FastZ/Msex50bp_full_abr.txt",header=F,stringsAsFactors=F))
colnames(msld)<-c("Scaffold","r2","dist")

#lazily copying to make sure we're workign with the same scaffold assignments

#sanity check: what does coverage claim to be z-linked?
setwd("/Users/Andrew/Documents/Manduca_Demography")
library("data.table")

msamp<-fread("S35_median_cov.txt",header=F,stringsAsFactors=F)
colnames(msamp)<-c("Scaffold","Cov")
m_norm<-msamp$Cov/mean(msamp$Cov)
fsamp<-fread("LK5_median_cov.txt",header=F,stringsAsFactors=F)
colnames(fsamp)<-c("Scaffold","Cov")
f_norm<-fsamp$Cov/mean(fsamp$Cov)


cov_trans<-log((m_norm/f_norm),base=2)

cov_check<-as.data.frame(cbind(msamp$Scaffold,cov_trans),stringsAsFactors=F)
colnames(cov_check)<-c("Scaffold","Log2Cov")
zkmanmaster<-read.csv("ZENKOKU_m_sexta_genome_mat.csv",header=T, stringsAsFactors=F)
zk_hack<-as.data.frame(cbind(zkmanmaster$Scaffold,zkmanmaster$ZorA,zkmanmaster$Chromosome),stringsAsFactors=F)
colnames(zk_hack)<-c("Scaffold","ZorA","Chromosome")
scaf_pred <- zk_hack[match(unique(zk_hack$Scaffold), zk_hack$Scaffold),]


zcheck<-as.data.frame(merge(cov_check,scaf_pred,by="Scaffold"))
zcheck$Log2Cov<-as.numeric(as.character(zcheck$Log2Cov))

scaf_length<-fread("Manduca_OGS2.bed",header=F,stringsAsFactors=F)
scaf_length<-cbind(scaf_length[,1],scaf_length[,3])
colnames(scaf_length)<-c("Scaffold","Length")
zcheck<-merge(zcheck,scaf_length,by="Scaffold")

head(zcheck,1)
#per the Neo-Z, we'll only consider scaffolds over the N90 length
zcheck2<-zcheck[which(zcheck$Length>46400),]
put_z<-zcheck2[which(zcheck2$Log2Cov > 0.7),]
put_z <- put_z[is.finite(put_z$Log2Cov),]
put_no<-zcheck2[which(zcheck2$Log2Cov < 0.7),]
put_no <- put_no[is.finite(put_no$Log2Cov),]
put_no[which(put_no$Log2Cov==min(put_no$Log2Cov)),]
#incidentally, there's a few strongly female biased scaffolds
#        Scaffold          Log2Cov ZorA Chromosome Length
#419 scaffold00419 -3.13415    A          0                227848
#502 scaffold00505    -Inf         A          0                 164272
#815 scaffold00850    -Inf         A          0                   62167
#939 scaffold01000    -Inf         A          0                   46784

table(put_z$ZorA,put_z$Chromosome)
#     0  1 13 25
#  A  7  0  1  1
#  Z  0 27  0  0
#so we're recovering 9 additional Z linked scaffolds, 7 of which were previously unplaced...feels good. chr13 and 25 might merit investigation


table(put_no$ZorA,put_no$Chromosome)
#      0  10  11  12  13  14  15  16  17  18  19   2  20  21  22  23  24  25  26  27  28   3   4   5   6   7   8   9
#  A 157  26  61  34  30  17  35  23  23  26  31   6  24  19  42  22  47  27  24  16  22  20  30  28  17  23  19  24
#furthermore, we've got no erroneously Z-linked things


hist(zcheck2$Log2Cov,breaks=30)
abline(v=0.7)

#okay what if we update Z-linkage based on coverage
new_z_scafs<-unique(put_z$Scaffold)

zk_ishin<-zkmanmaster

zk_ishin[which(zk_ishin$Scaffold%in%new_z_scafs),14]<-"Z"
table(zk_ishin$ZorA)
table(zkmanmaster$ZorA)
#from this we pick up 43 genes
#the sum of newly recovered Z-linked sequence is 2,128,494

zkmanmaster<-zk_ishin
zkm<-as.data.table(cbind(zkmanmaster$Scaffold,zkmanmaster$ZorA))
zkm<-unique(zkm[,1:2])
colnames(zkm)<-c("Scaffold","Loc")


##working line

#KS params
median(na.omit(zkmanmaster$pN.KS/zkmanmaster$Non.sites))
azn<-zkmanmaster$pN.KS/zkmanmaster$Non.sites
sd(azn[is.finite(azn)])



median(na.omit(zkmanmaster$pS.KS/zkmanmaster$Syn.sites))
azn<-zkmanmaster$pS.KS/zkmanmaster$Syn.sites
sd(azn[is.finite(azn)])

median(na.omit((zkmanmaster$pN.KS/zkmanmaster$Non.sites)/(zkmanmaster$pN.KS/zkmanmaster$Syn.sites)))
azn<-(zkmanmaster$pN.KS/zkmanmaster$Non.sites)/(zkmanmaster$pS.KS/zkmanmaster$Syn.sites)
sd(azn[is.finite(azn)])



#0x
a0x<-fread("/Users/Andrew/Documents/Manduca_Demography/KS_0x_auto.thetawin.txt.pestPG",header=T,stringsAsFactors=F)
a0p<-a0x$Tajima
mean(a0p)
#-0.03669585
median(a0p)
#0
sd(a0p)
#0.2988597

a0x<-fread("/Users/Andrew/Documents/Manduca_Demography/KS_0x_auto_unfold_pi.txt",header=T,stringsAsFactors=F)
a0p<-exp(a0x$Pairwise)
mean(a0p)
#0.002876751
median(a0p)
#4.265988e-08
sd(a0p)
#0.03419486




#4x
a0x<-fread("/Users/Andrew/Documents/Manduca_Demography/KS_4x_auto.thetawin.txt.pestPG",header=T,stringsAsFactors=F)
a0p<-a0x$Tajima
mean(a0p)
#-0.0148472
median(a0p)
#0
sd(a0p)
#0.367562


a0x<-fread("/Users/Andrew/Documents/Manduca_Demography/KS_4x_auto_unfold_pi.txt",header=T,stringsAsFactors=F)
a0p<-exp(a0x$Pairwise)
mean(a0p)
#0.02035638
median(a0p)
#3.374637e-07
sd(a0p)
#0.09044618


#AZ params

median(na.omit(zkmanmaster$pN.AZ/zkmanmaster$Non.sites))
sd(is.finite(zkmanmaster$pN.AZ/zkmanmaster$Non.sites))
azn<-zkmanmaster$pN.AZ/zkmanmaster$Non.sites
sd(azn[is.finite(azn)])


median(na.omit(zkmanmaster$pS.AZ/zkmanmaster$Syn.sites))
azn<-zkmanmaster$pS.AZ/zkmanmaster$Syn.sites
sd(azn[is.finite(azn)])

median(na.omit((zkmanmaster$pN.AZ/zkmanmaster$Non.sites)/(zkmanmaster$pN.AZ/zkmanmaster$Syn.sites)))
azn<-(zkmanmaster$pN.AZ/zkmanmaster$Non.sites)/(zkmanmaster$pS.AZ/zkmanmaster$Syn.sites)
sd(azn[is.finite(azn)])

#0x
a0x<-fread("/Users/Andrew/Documents/Manduca_Demography/AZ_0x_auto.thetawin.txt.pestPG",header=T,stringsAsFactors=F)
a0p<-a0x$Tajima
mean(a0p)
#-0.06731429
median(a0p)
#0
sd(a0p)
#0.3537572

a0x<-fread("/Users/Andrew/Documents/Manduca_Demography/AZ_0x_auto_unfold_pi.txt",header=T,stringsAsFactors=F)
a0p<-exp(a0x$Pairwise)
mean(a0p)
#0.00297649
median(a0p)
#7.947924e-08
sd(a0p)
0.03171317


#4x
a0x<-fread("/Users/Andrew/Documents/Manduca_Demography/AZ_4x_auto.thetawin.txt.pestPG",header=T,stringsAsFactors=F)
a0p<-a0x$Tajima
mean(a0p)
#-0.02589311
median(a0p)
#0
sd(a0p)
#0.4056982


a0x<-fread("/Users/Andrew/Documents/Manduca_Demography/AZ_4x_auto_unfold_pi.txt",header=T,stringsAsFactors=F)
a0p<-exp(a0x$Pairwise)
mean(a0p)
#0.02099016
median(a0p)
#4.322602e-07
sd(a0p)
0.08562383

#nc check
median(na.omit((zkmanmaster$pN/zkmanmaster$Non.sites)/(zkmanmaster$pN/zkmanmaster$Syn.sites)))
azn<-(zkmanmaster$pN/zkmanmaster$Non.sites)/(zkmanmaster$pS/zkmanmaster$Syn.sites)
sd(azn[is.finite(azn)])




