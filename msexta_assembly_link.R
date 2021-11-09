#Against my better judgment, I'm curious how the new M. sexta assembly fairs with chromosomes
library("data.table")
setwd("/Users/Andrew/Documents/Manduca_Demography/")



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


hist(zcheck2$Log2Cov,breaks=60, xlim = c(-0.7,1.55),xlab="Log2(M:F coverage)",main="Coveraged-based Z-linkage in the Carolina sphinx moth",ylab="Number of scaffolds")
abline(v=0.75)
text(-0.3,100, "Autosomal")
text(1.3,100, "Z-linked")


#okay what if we update Z-linkage based on coverage
new_z_scafs<-unique(put_z$Scaffold)

zk_ishin<-zkmanmaster

zk_ishin[which(zk_ishin$Scaffold%in%new_z_scafs),14]<-"Z"

#write.csv(zk_ishin,"ishin_m_sexta_genome_mat.csv",row.names=F,quote=F)

zki_abridged<-unique(zk_ishin[,4:5])


jhuman<-fread("satsuma_summary.chained.out",header=F)
colnames(jhuman)<-c("Scaffold", "S1", "S2", "JHU", "S3", "S4", "X", "Strand")

#stoploss search
#the OG assembly put the peak at 245000 to 250000 on scaffold 381
j381<-jhuman[which(jhuman$Scaffold=="scaffold00381"),]
plot(-1,-1,xlim=c(240000,260000),ylim=c(0,1))
segments(x0=j381$S1,y0=0.5,x1=j381$S2,y1=0.5)
j381[which(j381$S1>240000),]
#we should be looking for Scaffold 12, 9,241,951 - 9,261,972
#supplemental figure###
plot(-1,-1,xlim=c(000,310000),ylim=c(0.4,0.6),yaxt='n',xlab="Scaffold00381 position",main="A missassembly misplaced the stop-loss variant in the old assembly",ylab="")
j3811<-j381[which(j381$JHU=="HiC_scaffold_12"),]
j3812<-j381[which(j381$JHU=="HiC_scaffold_30"),]
j3813<-j381[which(j381$JHU=="HiC_scaffold_600"),]
segments(x0=j3811$S1,y0=0.506,x1=j3811$S2,y1=0.506,col="red",lwd=7)
segments(x0=j3812$S1,y0=0.5,x1=j3812$S2,y1=0.5,col="black",lwd=7)
segments(x0=j3813$S1,y0=0.5,x1=j3813$S2,y1=0.5,col="black",lwd=7)
text(100000,0.53,"Syntenic to Hi_C_scaffold_30 
(chr 13)")
text(255000,0.53,"Syntenic to Hi_C_scaffold_12
 (chr 12)",col="red")
 rect(xleft=245000,ybottom=0.44,xright=248000,ytop=0.51,col=rgb(1,1,1,alpha=0.5))
text(188000,0.46,"Location of Msex2.10493 ->")

chromassign<-merge(jhuman,zki_abridged,by="Scaffold")
msexta_truth<-table(chromassign$JHU,chromassign$Chromosome)
#write.csv(msexta_truth,"msexta_hic_chroms.csv",row.names=T,quote=F)

#460 is actually part of 12 invert by this count as well


#now let's make a simplfying assumption:
ms_simp<-colnames(msexta_truth)[apply(msexta_truth,1,which.max)]
mss<-cbind(msexta_truth,ms_simp)
mss <- cbind(rownames(mss), data.frame(mss, row.names=NULL))
colnames(mss)[1]<-"Scaffold"
colnames(mss)[31]<-"Synteny_liftover_Chr"


jbed<-fread("moth_canu_nanopolish_racon1.FINAL.bed",header=F)
colnames(jbed)<-c("Scaffold", "Start", "Length")

Mst<-(as.data.frame(merge(mss,jbed,by="Scaffold"),stringsAsFactors=F))
#Mst$Synteny_liftover_Chr<-as.numeric(as.character(Mst$Synteny_liftover_Chr))
Mst$Synteny_liftover_Chr<-numeric(as.character(Mst$Synteny_liftover_Chr))
Mst[which(Mst$Synteny_liftover_Chr=="1"),]
Mst[which(Mst$Synteny_liftover_Chr=="1"),31]<-"Z"
#write.csv(Mst,"msexta_JHU_assembly_chr_anchor.csv",row.names=T,quote=F)