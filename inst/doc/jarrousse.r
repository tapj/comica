library(BiotypeR)
library(ggplot2)
library(ade4)
library(vegan)
library(gridExtra)
library(reshape2)

user=Sys.info()[["user"]]
setwd(paste("/Users/",user,"/Dropbox/GitHub/comica/inst/extdata", sep=""))

# preprocessing

comica=comica_metadata=read.csv2("comica_metadata.csv", row.names=1)
comica_metadata=comica_metadata[-which(comica$Date=="200903"),]


seq_mice_metadata=read.csv2("seq_metadata.csv",row.names=1)
seq_mice_metadata=seq_mice_metadata[-which(seq_mice_metadata$comment=="donneur KJ28, 32 ans"),]

otu=read.csv2("OTU_Table.csv",row.names=1, check.names=F)
tax=read.csv2("Taxonomy.csv")



# microbiota analysis

donor=levels(as.factor(as.character(seq_mice_metadata$donor)))

donor.dd=data.frame(row.names=donor, time_corrected=rep(0,length(donor)), donor=donor, donor_status=c("cancer","cancer","normal","normal"))

seq_metadata=rbind(seq_mice_metadata[,7:9],donor.dd)

genus=apply(otu[,row.names(seq_metadata)], 2, tapply, tax$genus, sum)
#genus=apply(otu, 2, tapply, tax$genus, sum)
genus=as.data.frame(prop.table(as.matrix(genus),2))
genus.jsd=dist.JSD(genus)
genus.pco=dudi.pco(as.dist(as.matrix(genus.jsd)[-c(40:43),-c(40:43)]), scannf=F, nf=3)
s.label(genus.pco$li)

dd=data.frame(pco=genus.pco$li, seq_metadata[-c(40:43),])

ggplot(dd, aes(x=pco.A1, y=pco.A2)) +
geom_point(size=5, aes(colour=donor_status, shape=donor)) +
xlab(paste("PCo1",round(genus.pco$eig[1]/sum(genus.pco$eig)*100,2 ),"%")) +
ylab(paste("PCo2",round(genus.pco$eig[2]/sum(genus.pco$eig)*100,2 ),"%")) + 
scale_shape_manual("Donors", labels=c("HK-JAN10","HL-JUN10","HN-JAN09","HN-JUN10"), values=15:18) +
 scale_color_discrete(name="Health status") -> jsd_pcoa_plot


donor.jsd=NULL

for(i in donor) {

tmp=as.matrix(genus.jsd)[seq_metadata$donor==i,seq_metadata$donor==i]

n=colnames(tmp)
tmp=tmp[,i]

tps=seq_metadata[n,"time_corrected"]

dd=data.frame(n=rep(i,length(tmp)), time=tps, jsd=tmp)

donor.jsd=rbind(donor.jsd,dd)


}

dd=data.frame(

	median=tapply(donor.jsd$jsd, donor.jsd$time, median),
	q25=tapply(donor.jsd$jsd, donor.jsd$time, quantile, 0.25),
	q75=tapply(donor.jsd$jsd, donor.jsd$time, quantile, 0.75),
	time=levels(as.factor(donor.jsd$time))

)


ggplot(donor.jsd, aes(x=time, y=jsd)) + 
#geom_point(aes(colour=n)) +
geom_line(data=dd, aes(x=as.numeric(as.character(time)), y=median)) +
geom_line(data=dd, aes(x=as.numeric(as.character(time)), y=q25) , linetype=2 ) +
geom_line(data=dd, aes(x=as.numeric(as.character(time)), y=q75) , linetype=2 ) + 
xlab("Sampling days") + ylab("Mice microbiota JSD distance\nfrom donors") -> jsd_baseline_plot


genus.denoized=noise.removal(genus[,-c(40:43)], percent=1)

genus.corr=data.frame(
Corr1=cor(t(genus.denoized), genus.pco$li[,1]),
Corr2=cor(t(genus.denoized), genus.pco$li[,2])
)[-1,]

genus.corr=genus.corr[which(abs(genus.corr[,1]) > 0.5 |  abs(genus.corr[,2]) > 0.5), ]

ggplot(genus.corr, aes(x=Corr1, y=Corr2)) +
geom_text(label=row.names(genus.corr), fontface=3, size=3) +
xlim(-1,1) + ylim(-1,1) + 
geom_segment(xend=0, yend=0,linetype=2) -> genus_plot

figure1_plot=grid.arrange(jsd_pcoa_plot, arrangeGrob(jsd_baseline_plot, genus_plot, ncol=1), ncol=2, widths=c(2,1))
                    


					
					
# clinical analysis
group=as.factor(comica_metadata$Donnor_status:comica_metadata$conditions)
comica.tmp=comica_metadata
#group=as.factor(as.character(group[-which(comica_metadata$Date=="200903")]))

#pdf(); for(i in c(6:10,15:29)) {boxplot(comica.tmp[,i]~group, main=names(comica.tmp[i]))}; dev.off()

comica.tmp=comica_metadata[,c(6:10,15:29)]
for(i in 1:20) { idx=attr(na.omit(comica.tmp[,i]),"na.action"); comica.tmp[idx,i]<-mean(na.omit(comica.tmp[,i]))}


library(ade4)
pca=dudi.pca(comica.tmp[,c(3:4,6:18)], scannf=F, nf=3)
scatter(pca)

comica.bca=bca(pca, fac=group, scannf=F, nf=3)
plot(comica.bca, yax=2)
dev.new()
plot(comica.bca, yax=3)


# partial co-inertia between clinical and microbiota


#comica.bca.partial=comica.bca
#comica.bca.partial$ls=comica.bca$ls[which(comica_metadata$Microbiota_sequencing=="yes"),]

pca.partial=pca

pca.partial$tab=pca$tab[which(comica_metadata$Microbiota_sequencing=="yes"),]
pca.partial$li=pca$li[which(comica_metadata$Microbiota_sequencing=="yes"),]
pca.partial$l1=pca$l1[which(comica_metadata$Microbiota_sequencing=="yes"),]
pca.partial$lw=pca$lw[which(comica_metadata$Microbiota_sequencing=="yes")]

seq_names=row.names(seq_mice_metadata[(as.character(seq_mice_metadata$Mice_ID) %in% row.names(pca.partial$li)) & seq_mice_metadata$time_corrected==42, ])


genus.pco.partial=genus.pco

genus.pco.partial$tab=genus.pco$tab[seq_names,]
genus.pco.partial$li=genus.pco$li[seq_names,]
genus.pco.partial$l1=genus.pco$l1[seq_names,]
genus.pco.partial$lw=pca$lw[which(comica_metadata$Microbiota_sequencing=="yes")]

partial.coi=coinertia(genus.pco.partial,pca.partial, scannf=F, nf=2)

partial.coi.microbiota.lX=as.matrix(genus.pco$tab) %*% as.matrix(partial.coi$c1)
partial.coi.clinical.lY=as.matrix(pca$tab) %*% as.matrix(partial.coi$l1)


normalise.w <- function(X, w) {
        f2 <- function(v) sqrt(sum(v * v * w))
        norm <- apply(X, 2, f2)
        X <- sweep(X, 2, norm, "/")
        return(X)
    }

	
partial.coi.microbiota.mX=normalise.w(partial.coi.microbiota.lX, genus.pco$lw)
partial.coi.clinical.mY=normalise.w(partial.coi.clinical.lY, pca$lw)

	
	
dd.clin=data.frame(partial.coi.clinical.mY,
		comica_metadata[row.names(partial.coi.clinical.mY),4:5],
		type=rep("clinical",dim(partial.coi.clinical.mY)[1]),
		Mice_ID=row.names(partial.coi.clinical.mY),
		time=rep("42", dim(partial.coi.clinical.mY)[1]))

dd.microbiota=data.frame(partial.coi.microbiota.mX,
				Donnor_status=seq_mice_metadata[row.names(partial.coi.microbiota.mX),9],
				conditions=rep("NaCl",dim(partial.coi.microbiota.mX)[1]),
				type=rep("microbiota",dim(partial.coi.microbiota.mX)[1]),
				Mice_ID=seq_mice_metadata[row.names(partial.coi.microbiota.mX),"Mice_ID"],
				time=as.character(seq_mice_metadata[row.names(partial.coi.microbiota.mX),"time_corrected"]))

names(dd.microbiota) = names(dd.clin)

dd=rbind(dd.clin, dd.microbiota)

dd.connected=cbind(dd[dd$Mice_ID %in% row.names(pca.partial$li) & dd$time=="42" & dd$type=="clinical"  ,1:2],dd[dd$Mice_ID %in% row.names(pca.partial$li) & dd$time=="42" & dd$type=="microbiota"  ,1:2])
names(dd.connected)=paste(c("c","c","m","m"), names(dd.connected), sep="_")

ggplot(dd, aes(x=RS1, y=RS2)) +  geom_segment(data=dd.connected, aes(x=c_RS1, y=c_RS2, xend=m_RS1, yend=m_RS2) ) +
 geom_point(data=dd, aes(shape=type, colour=Donnor_status, alpha=conditions), size=5) +
 scale_alpha_discrete(range = c(0.5, 0.9)) +
 xlab("PC1") + ylab("PC2") + xlim(-2,2) -> partial_coinertia_plot
 
 #+ xlim(-2,2)  
 
 
 

dd.microbiota=cbind(dd.microbiota,t(genus[,-c(40:43)]))
 
dd.clin=cbind(dd.clin, comica.tmp[,c(3:4,6:18)])


ggplot(dd.microbiota, aes(x=RS1, y=Bacteroides, colour=c("black")))+
 geom_point(alpha=0.8) +
 geom_point(aes(y=Coprococcus, colour="red"), alpha=0.8 )  +
 scale_y_log10("Microbial\nrel. abundance\n ") +
 scale_colour_manual("",values=c("black","red"), labels=c("Bacteroides","Coprococcus")) +
 xlab("PC1") + xlim(-2,2) -> coinertia_bacteria_fit_plot

ggplot(comica, aes(x=as.factor(conditions:Donnor_status), y=FCA_cm_colon)) +  
geom_boxplot(aes(fill=as.factor(conditions:Donnor_status))) + geom_jitter(position = position_jitter(w = 0.1)) +
#ggplot(dd.clin, aes(x=as.factor(conditions:Donnor_status), y=FCA_cm_colon)) +
 xlab("") + ylab("ACFs per cm") -> ACF_boxplot
 
 
### ACF rarefaction per cm ###

comica.acf=na.omit(comica[c(4,5,6,7,8)])
group=as.factor(comica.acf$Donnor_status: comica.acf$conditions)
comica.acf$group=group

rarefaction_comica= function(group, per_cm, len) {

	group=as.factor(group)
	nb.group=length(levels(group))
	len.max=round(max(tapply(len, group, sum)))
	
	results=matrix(nr=len.max, nc=nb.group)
	colnames(results)=levels(group)
	#rownames(results)=1:len.max
	
	for(i in levels(group)) {
	
	idx=which(group ==i)
	
	cm.group=rep(per_cm[idx], round(len[idx]))
	
		for (s in seq(1, round(sum(len[idx])), 1)) {
						
						
			m=sum(sample(cm.group,s, replace=TRUE))
			results[s,i]=m
			
	
		}
	
	}

	return(data.frame(len=1:len.max,results))


} 

r=rarefaction_comica(comica.acf$group, comica.acf$FCA_cm_colon, comica.acf$lg_colon)

r=NULL

for( i in 1:100) {
	r=rbind(r,rarefaction_comica(comica.acf$group, comica.acf$FCA_cm_colon, comica.acf$lg_colon))
}

ggplot(data=melt(r, id="len"), aes(x=len, y=value, group=variable)) + geom_point(aes(colour=variable))

ggplot(data=melt(r, id="len"), aes(x=len, y=value, group=variable))  +
 geom_quantile(aes(colour=variable), quantiles=c(0.25,0.75), linetype=2) +
 geom_quantile(aes(colour=variable), quantiles=c(0.50), size=1.5) +
 xlab("accumulated\ncolon length observed (cm)") + ylab("nb of ACFs") -> ACF_curve_plot

 ###############################################################


 
 dd.clin.melt=melt(cbind(dd.clin[,c(1,3,4)],apply(dd.clin[,11:18],2,scale)), id=c("RS1","Donnor_status","conditions"))
 
 ggplot(dd.clin.melt, aes(x=RS1, y=value)) + geom_point(aes(colour=variable))
 
 dd.clin.melt=cbind( dd.clin.melt, t(sapply(strsplit(as.character(dd.clin.melt$variable), split="_"), function(x) x)))
 names(dd.clin.melt)=c(names(dd.clin.melt)[1:5],"organ","gene","raw")
 
 
ggplot(dd.clin.melt, aes(x=RS1, y=value)) +
 geom_point(aes(colour=Donnor_status, shape=organ), alpha=0.8) +
 geom_smooth() +
 xlab("PC1") + ylab("scale normalized\ntranscriptomic assays\n(ELF3, HES1, KLF4\nand MATH1 genes)") + xlim(-2,2) -> transcriptomic_plot
 
 # 

pdf("Figure1_microbiota.pdf", h=9, w=14)
figure1_plot <- grid.arrange(jsd_pcoa_plot+theme_bw(), arrangeGrob(jsd_baseline_plot+theme_bw(), genus_plot+theme_bw(), ncol=1), ncol=2, widths=c(2,1))
dev.off()

pdf("Figure2_ACFs.pdf", h=5, w=10)
figure2_plot <- grid.arrange(ACF_boxplot+theme_bw()+theme(legend.position="none"), ACF_curve_plot+theme_bw()+theme(legend.title=element_blank(),legend.position=c(0.2, 0.8)), ncol=2)
dev.off()

pdf("Figure3_clinical.pdf", h=9, w=14)
figure3_plot <- grid.arrange(partial_coinertia_plot+theme_bw()+theme(legend.position = c(0.1, 0.2)), arrangeGrob(coinertia_bacteria_fit_plot+theme_bw(),  transcriptomic_plot+theme_bw() ,  ncol=1), ncol=2, widths=c(4,3))
dev.off()          

figureS1_transcripto_boxplot <-ggplot(dd.clin.melt, aes(y=value, x=as.factor(Donnor_status:conditions))) + geom_boxplot(aes(fill=as.factor(Donnor_status:conditions))) + facet_grid(gene~organ) + theme_bw()
 
 
























