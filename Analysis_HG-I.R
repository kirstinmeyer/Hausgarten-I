##Data analysis for HG I time series

############################################################################
#Load packages needed for analysis
library(vegan)
library(ggplot2)
library(devtools)
library(pairwiseAdonis)
library(car)
library(dbstats)
library(ggrepel)
library(ggpubr)
library(tidyverse)
library(MASS)
library(ggpattern)

#Read in files for all data
apath<-file.choose() #EGIV_all_species
Allspec<-read.csv(apath,row.names=1) 

bpath<-file.choose() #EGIV_year
Year<-read.csv(bpath,row.names=0) 

cpath<-file.choose() #EGIV_dropstones
Dropston<-read.csv(cpath,row.names=1)
HGI.env[,11]<-HGI.env$ChlA/HGI.env$CPE*100
HGI.env.5<-matrix(,nrow=450,ncol=11)
HGI.env.5<-data.frame(HGI.env.5)
HGI.env.5[1:90,]<-HGI.env[1,]
HGI.env.5[91:180,]<-HGI.env[6,]
HGI.env.5[181:270,]<-HGI.env[11,]
HGI.env.5[271:360,]<-HGI.env[16,]
HGI.env.5[361:450,]<-HGI.env[21,]
colnames(HGI.env.5)<-c("ChlA","CPE","Protein","NAO","AO","Bot_temp","Bot_sal","Surf_temp","Surf_sal","C.org","PropChl")
HGI.env.5<-HGI.env.5[,2:11]

HGI.env.1<-matrix(,nrow=450,ncol=11)
HGI.env.1<-data.frame(HGI.env.1)
HGI.env.1[1:90,]<-HGI.env[16,]
HGI.env.1[91:180,]<-HGI.env[17,]
HGI.env.1[181:270,]<-HGI.env[18,]
HGI.env.1[271:360,]<-HGI.env[20,]
HGI.env.1[361:450,]<-HGI.env[21,]
colnames(HGI.env.1)<-c("ChlA","CPE","Protein","NAO","AO","Bot_temp","Bot_sal","Surf_temp","Surf_sal","C.org","PropChl")
HGI.env.1<-HGI.env.1[,2:11]

############################################################################
#Univariate analyses for 5-year dataset
#Total abundance
HGI.5year$abund<-rowSums(HGI.5year)
abund.aov<-aov(HGI.5year$abund~as.character(HGI.years.5), data=HGI.5year)
summary(abund.aov)
abund.aov.factor <- aov(HGI.5year$abund~as.factor(HGI.years.5), data=HGI.5year)
abund.aov.tukey<-TukeyHSD(abund.aov.factor)
abund.aov.tukey

#Species richness
for (j in 1:450) {HGI.5year$rich[j]<-sum(HGI.5year[j,1:12]>0, na.rm=TRUE)} 
rich.aov<-aov(HGI.5year$rich~as.character(HGI.years.5), data=HGI.5year)
summary(rich.aov)
rich.aov.factor <- aov(HGI.5year$rich~as.factor(HGI.years.5), data=HGI.5year)
rich.aov.tukey<-TukeyHSD(rich.aov.factor)
rich.aov.tukey

#Bar graphs
#Creating standard error function
stand_err <- function(data){
  sd(data)/sqrt(length(data))
}

#Total abundance
abund_sum <- aggregate(HGI.5year$abund~HGI.years.5, FUN = mean)
colnames(abund_sum) = c("Year","mean_density") 
abund_sum_se <- aggregate(HGI.5year$abund~HGI.years.5, FUN = stand_err)
abund_sum$se <- abund_sum_se$`HGI.5year$abund`
abund.bar.5<- ggplot(data=abund_sum, aes(x=factor(Year),y=mean_density)) +
  geom_bar(stat="identity", fill="gray", alpha=0.7)  +
  geom_errorbar(aes(ymin=mean_density-se, ymax=mean_density+se, width=.2)) + 
  theme_classic()+ xlab("Year") + ylab(expression(paste(plain("Total individuals")," m"^plain(-2))))+
  theme(legend.text=element_text(size=12),legend.title=element_text(size=14))+
  theme(axis.text.x=element_text(size=12),axis.title.x=element_text(size=14))+
  theme(axis.text.y=element_text(size=12),axis.title.y=element_text(size=14))
abund.bar.5

#Species richness
rich_sum <- aggregate(HGI.5year$rich~HGI.years.5, FUN = mean)
colnames(rich_sum) = c("Year","mean_density") 
rich_sum_se <- aggregate(HGI.5year$rich~HGI.years.5, FUN = stand_err)
rich_sum$se <- rich_sum_se$`HGI.5year$rich`
rich.bar.5<- ggplot(data=rich_sum, aes(x=factor(Year),y=mean_density)) +
  geom_bar(stat="identity", fill="gray", alpha=0.7)  +
  geom_errorbar(aes(ymin=mean_density-se, ymax=mean_density+se, width=.2)) + 
  theme_classic()+ xlab("Year") + ylab(expression(paste(plain("Species richness")," m"^plain(-2))))+
  theme(legend.text=element_text(size=12),legend.title=element_text(size=14))+
  theme(axis.text.x=element_text(size=12),axis.title.x=element_text(size=14))+
  theme(axis.text.y=element_text(size=12),axis.title.y=element_text(size=14))
rich.bar.5

############################################################################
##Univariate regression for 5-year dataset

#Find best-fit model for total abundance
fit<-lm(HGI.5year$abund~CPE+Protein+NAO+AO+Bot_temp+Bot_sal+Surf_temp+Surf_sal+C.org+PropChl,data=HGI.env.5)
step<-stepAIC(fit,direction="both")
anova(step)
RsquareAdj(step)

#Find best-fit model for species richness
fit<-lm(HGI.5year$rich~CPE+Protein+NAO+AO+Bot_temp+Bot_sal+Surf_temp+Surf_sal+C.org+PropChl,data=HGI.env.5)
step<-stepAIC(fit,direction="both")
anova(step)
RsquareAdj(step)

############################################################################
##Multivariate analyses for 5-year dataset
#Bray-Curtis similarity matrix
HGI.5year.dist <-vegdist(HGI.5year, method="bray")

#PERMANOVA main test
HGI.years.5<-as.data.frame(HGI.years.5)
colnames(HGI.years.5)<-c("Year")
HGI.5year.perm <-adonis2(HGI.5year~HGI.years.5,data=HGI.years.5, permutations = 999, method="bray")
HGI.5year.perm

#PERMANOVA pairwise post hoc test
pairwise.adonis(HGI.5year,HGI.years.5)

#Simper analysis
HGI.simper<-simper(HGI.5year, HGI.years.5, permutations=999, parellel=1)
HGI.simper

#Find best-fit model for multivariate data
m0 <- rda(HGI.5year~AO, HGI.env.5)
m1 <- rda(HGI.5year~AO + ., HGI.env.5)
m <- ordistep(m0, scope = list(lower=m0, upper=m1))
m
anova(m)
RsquareAdj(m)
dbRDA.5<-capscale(HGI.5year~AO+Surf_temp+C.org+Protein,HGI.env.5)
summary(dbRDA.5)

#Make dbRDA plot
env.5 = as.data.frame(dbRDA.5$CCA$biplot)
env.5["env"] = c("AO","Surf. temp.","C-org","PP")
bio.5=as.data.frame(dbRDA.5$Ybar)
bio.5$Year = HGI.years.5$Year
Fig5A<-ggplot() +
  geom_point(data=bio.5,aes(x=Dim1,y=Dim2,color=as.factor(Year),shape=as.factor(Year)),size=2,stroke=1)+
  theme_classic()+xlab("CAP1 (79.9% of variation)")+ylab("CAP2 (0.2% of variation)")+
  geom_segment(data=env.5,aes(x=0,y=0,xend=CAP1,yend=CAP2),arrow=arrow(length=unit(0.2,"cm"))) +
  geom_text_repel(data=env.5,aes(x=CAP1,y=CAP2,label=c("AO","Surf. temp.","C-org","PP")),nudge_y=-0.05,color="black",size=4)+
  scale_color_manual(values=c("2002"="#332288","2007"="#117733",
                              "2012"="#AA4499","2017"="#88CCEE",
                              "2022"="#888888"))+
  scale_shape_manual(values=c("2002"=0,"2007"=1,"2012"=2,"2017"=5,"2022"=6))+
  guides(color=guide_legend("Year"),shape=guide_legend("Year"))+
  theme(legend.text=element_text(size=12),legend.title=element_text(size=14))+
  theme(axis.text.x=element_text(size=12),axis.title.x=element_text(size=14))+
  theme(axis.text.y=element_text(size=12),axis.title.y=element_text(size=14))
Fig5A

############################################################################
#Univariate analyses for 1-year dataset
#Total abundance
HGI.1year$abund<-rowSums(HGI.1year)
abund.aov<-aov(HGI.1year$abund~as.character(HGI.years.1), data=HGI.1year)
summary(abund.aov)
abund.aov.factor <- aov(HGI.1year$abund~as.factor(HGI.years.1), data=HGI.1year)
abund.aov.tukey<-TukeyHSD(abund.aov.factor)
abund.aov.tukey

#Species richness
for (j in 1:450) {HGI.1year$rich[j]<-sum(HGI.1year[j,1:19]>0, na.rm=TRUE)} 
rich.aov<-aov(HGI.1year$rich~as.character(HGI.years.1), data=HGI.1year)
summary(rich.aov)
rich.aov.factor <- aov(HGI.1year$rich~as.factor(HGI.years.1), data=HGI.1year)
rich.aov.tukey<-TukeyHSD(rich.aov.factor)
rich.aov.tukey

#Bar graphs
#Total abundance
abund_sum <- aggregate(HGI.1year$abund~HGI.years.1, FUN = mean)
colnames(abund_sum) = c("Year","mean_density") 
abund_sum_se <- aggregate(HGI.1year$abund~HGI.years.1, FUN = stand_err)
abund_sum$se <- abund_sum_se$`HGI.1year$abund`
abund.bar.1<- ggplot(data=abund_sum) +
  geom_bar_pattern(aes(x=factor(Year),y=mean_density),stat="identity",fill="gray",alpha=0.7,pattern="stripe",pattern_fill="white",pattern_color="white")+
  geom_errorbar(aes(x=factor(Year),ymin=mean_density-se, ymax=mean_density+se, width=.2)) + 
  theme_classic()+ xlab("Year") + ylab(expression(paste(plain("Total individuals")," m"^plain(-2))))+
  theme(legend.text=element_text(size=12),legend.title=element_text(size=14))+
  theme(axis.text.x=element_text(size=12),axis.title.x=element_text(size=14))+
  theme(axis.text.y=element_text(size=12),axis.title.y=element_text(size=14))
abund.bar.1

#Species richness
rich_sum <- aggregate(HGI.1year$rich~HGI.years.1, FUN = mean)
colnames(rich_sum) = c("Year","mean_density") 
rich_sum_se <- aggregate(HGI.1year$rich~HGI.years.1, FUN = stand_err)
rich_sum$se <- rich_sum_se$`HGI.1year$rich`
rich.bar.1<- ggplot(data=rich_sum) +
  geom_bar_pattern(aes(x=factor(Year),y=mean_density),stat="identity", fill="gray", alpha=0.7,pattern="stripe",pattern_fill="white",pattern_color="white")  +
  geom_errorbar(aes(x=factor(Year),ymin=mean_density-se, ymax=mean_density+se, width=.2)) + 
  theme_classic()+ xlab("Year") + ylab(expression(paste(plain("Species richness")," m"^plain(-2))))+
  theme(legend.text=element_text(size=12),legend.title=element_text(size=14))+
  theme(axis.text.x=element_text(size=12),axis.title.x=element_text(size=14))+
  theme(axis.text.y=element_text(size=12),axis.title.y=element_text(size=14))
rich.bar.1


############################################################################
##Univariate regression for 1-year dataset

#Find best-fit model for total abundance
fit<-lm(HGI.1year$abund~CPE+Protein+NAO+AO+Bot_temp+Bot_sal+Surf_temp+Surf_sal+C.org+PropChl,data=HGI.env.1)
step<-stepAIC(fit,direction="both")
anova(step)
RsquareAdj(step)

#Find best-fit model for species richness abundance
fit<-lm(HGI.1year$rich~CPE+Protein+NAO+AO+Bot_temp+Bot_sal+Surf_temp+Surf_sal+C.org+PropChl,data=HGI.env.1)
step<-stepAIC(fit,direction="both")
anova(step)
RsquareAdj(step)

############################################################################
##Multivariate analyses for 1-year dataset
#Bray-Curtis similarity matrix
HGI.1year.dist <-vegdist(HGI.1year, method="bray")

#PERMANOVA main test
HGI.years.1<-as.data.frame(HGI.years.1)
colnames(HGI.years.1)<-c("Year")
HGI.1year.perm <-adonis2(HGI.1year~HGI.years.1,data=HGI.years.1, permutations = 999, method="bray")
HGI.1year.perm

#PERMANOVA pairwise post hoc test
pairwise.adonis(HGI.1year[,1:19],HGI.years.1)

#Simper analysis
HGI.simper<-simper(HGI.1year,group=HGI.years.1, permutations=999)
HGI.simper

#Find best-fit model for data
m0 <- rda(HGI.1year~AO,HGI.env.1)
m1 <- rda(HGI.1year~AO + ., HGI.env.1)
m <- ordistep(m0, scope = list(lower=m0, upper=m1))
m
anova(m)
RsquareAdj(m)
dbRDA.1<-capscale(HGI.1year~AO+NAO+Protein,HGI.env.1)
summary(dbRDA.1)

#Make dbRDA plot
env.1 = as.data.frame(dbRDA.1$CCA$biplot)
env.1["env"] = c("AO","NAO","PP")
bio.1=as.data.frame(dbRDA.1$Ybar)
bio.1$Year = HGI.years.1$Year
Fig5B<-ggplot() +
  geom_point(data=bio.1,aes(x=Dim1,y=Dim2,color=as.factor(Year),shape=as.factor(Year)),size=2,stroke=1)+
  theme_classic()+xlab("CAP1 (20.3% of variation)")+ylab("CAP2 (5.1% of variation)")+
  geom_segment(data=env.1,aes(x=0,y=0,xend=CAP1,yend=CAP2),arrow=arrow(length=unit(0.2,"cm"))) +
  geom_text_repel(data=env.1,aes(x=CAP1,y=CAP2,label=env.1$env),nudge_y=-0.05,color="black",size=4)+
  scale_color_manual(values=c("2017"="#88CCEE","2018" = "#CC6677",
                              "2019"="#DDCC77", "2021"="#44AA99",
                              "2022"="#888888"))+
  scale_shape_manual(values=c("2017"=5,"2018"=0,"2019"=1,"2021"=2,"2022"=6))+
  guides(color=guide_legend("Year"),shape=guide_legend("Year"))+
  theme(legend.text=element_text(size=12),legend.title=element_text(size=14))+
  theme(axis.text.x=element_text(size=12),axis.title.x=element_text(size=14))+
  theme(axis.text.y=element_text(size=12),axis.title.y=element_text(size=14))
Fig5B

############################################################################
##Species accumulation curves
yr.2002<-HGI.5year[1:90,1:12]
yr.2007<-HGI.5year[91:180,1:12]
yr.2012<-HGI.5year[181:270,1:12]
yr.2017<-HGI.5year[271:360,1:12]
yr.2018<-HGI.1year[91:180,1:12]
yr.2019<-HGI.1year[181:270,1:12]
yr.2021<-HGI.1year[271:360,1:12]
yr.2022<-HGI.1year[361:450,1:12]
yr.2002.spa<-specaccum(yr.2002,permutations=100)
yr.2007.spa<-specaccum(yr.2007,permutations=100)
yr.2012.spa<-specaccum(yr.2012,permutations=100)
yr.2017.spa<-specaccum(yr.2017,permutations=100)
yr.2018.spa<-specaccum(yr.2018,permutations=100)
yr.2019.spa<-specaccum(yr.2019,permutations=100)
yr.2021.spa<-specaccum(yr.2021,permutations=100)
yr.2022.spa<-specaccum(yr.2022,permutations=100)
SPAC<-data.frame(yr.2002.spa$sites,yr.2002.spa$richness,yr.2007.spa$richness,
                 yr.2012.spa$richness,yr.2017.spa$richness,yr.2018.spa$richness,
                 yr.2019.spa$richness,yr.2021.spa$richness,yr.2022.spa$richness,
                 yr.2002.spa$sd*0.196,yr.2007.spa$sd*0.196,
                 yr.2012.spa$sd*0.196,yr.2017.spa$sd*0.196,yr.2018.spa$sd*0.196,
                 yr.2019.spa$sd*0.196,yr.2021.spa$sd*0.196,yr.2022.spa$sd*0.196)
colnames(SPAC)<-c("Photos","Rich02","Rich07","Rich12","Rich17","Rich18","Rich19","Rich21","Rich22","CI02","CI07","CI12","CI17","CI18","CI19","CI21","CI22")

#Species accumulation figure
Specaccumfig<-ggplot(data=SPAC,aes(x=Photos))+
  geom_line(aes(y=Rich02,color="2002"),linewidth=1)+
  geom_ribbon(aes(x=Photos,ymin=(Rich02-CI02),ymax=(Rich02+CI02)),fill="#332288",alpha=0.5)+
  geom_line(aes(y=Rich07,color="2007"),linewidth=1)+
  geom_ribbon(aes(x=Photos,ymin=(Rich07-CI07),ymax=(Rich07+CI07)),fill="#117733",alpha=0.5)+
  geom_line(aes(y=Rich12,color="2012"),linewidth=1)+
  geom_ribbon(aes(x=Photos,ymin=(Rich12-CI12),ymax=(Rich12+CI12)),fill="#AA4499",alpha=0.5)+
  geom_line(aes(y=Rich17,color="2017"),linewidth=1)+
  geom_ribbon(aes(x=Photos,ymin=(Rich17-CI17),ymax=(Rich17+CI17)),fill="#88CCEE",alpha=0.5)+
  geom_line(aes(y=Rich18,color="2018"),linewidth=1)+
  geom_ribbon(aes(x=Photos,ymin=(Rich18-CI18),ymax=(Rich18+CI18)),fill="#CC6677",alpha=0.5)+
  geom_line(aes(y=Rich19,color="2019"),linewidth=1)+
  geom_ribbon(aes(x=Photos,ymin=(Rich19-CI19),ymax=(Rich19+CI19)),fill="#DDCC77",alpha=0.5)+
  geom_line(aes(y=Rich21,color="2021"),linewidth=1)+
  geom_ribbon(aes(x=Photos,ymin=(Rich21-CI21),ymax=(Rich21+CI21)),fill="#44AA99",alpha=0.5)+
  geom_line(aes(y=Rich22,color="2022"),linewidth=1)+
  geom_ribbon(aes(x=Photos,ymin=(Rich22-CI22),ymax=(Rich22+CI22)),fill="#888888",alpha=0.5)+
  theme_classic()+xlab("No. photos") + ylab("Species richness")+
  scale_color_manual(values=c("2002"="#332288","2007"="#117733",
                              "2012"="#AA4499","2017"="#88CCEE",
                              "2018"="#CC6677","2019"="#DDCC77",
                              "2021"="#44AA99","2022"="#888888"))+
  guides(color=guide_legend("Year"))+
  theme(legend.text=element_text(size=12),legend.title=element_text(size=14))+
  theme(axis.text.x=element_text(size=12),axis.title.x=element_text(size=14))+
  theme(axis.text.y=element_text(size=12),axis.title.y=element_text(size=14))
Specaccumfig

#Saving Image
ggsave("Fig3_Specaccum.png", width = 7, height = 4, units=c("in"), dpi=300)

############################################################################
#Graph environmental data
hpath<-file.choose()#Environmental_data
HGI.env.fig<-read.csv(hpath)
HGI.env.fig[,13]<-HGI.env.fig$ChlA/HGI.env.fig$CPE
Envfig<-ggplot(data=HGI.env.fig,aes(x=Year))+
  geom_point(aes(y=Protein,color="Protein",shape="Protein"),size=3)+
  geom_point(aes(y=C.org,color="Organic carbon",shape="Organic carbon"),size=4)+
  scale_color_manual(values=c("Protein"="firebrick","Organic carbon"="#DDCC77"))+
  scale_shape_manual(values=c("Protein"=17,"Organic carbon"=18))+
  theme_classic()+
  guides(color=guide_legend("Parameter"),shape=guide_legend("Parameter"))+
  scale_y_continuous(name="Protein (mg/cm3), C-org (%)")+
  theme(legend.text=element_text(size=12),legend.title=element_text(size=14))+
  theme(axis.text.x=element_text(size=12),axis.title.x=element_text(size=14))+
  theme(axis.text.y=element_text(size=12),axis.title.y=element_text(size=13))+
  theme(legend.position="right")
Envfig
#ggsave("Fig5_Environmental.png", width = 6, height = 4, units=c("in"), dpi=300)

HGI.env.fig[,12]<-HGI.env.fig$ChlA/HGI.env.fig$CPE*100
ChlPhafig<-ggplot(data=HGI.env.fig,aes(x=Year))+
  geom_point(aes(y=V12,color="% Chlorophyll",shape="% Chlorophyll"),size=3)+
  geom_point(aes(y=CPE,color="Total pigments",shape="Total pigments"),size=3)+
  theme_classic()+
  scale_color_manual(values=c("% Chlorophyll"="gray67","Total pigments"="darkblue"))+
  scale_shape_manual(values=c("Total pigments"=16,"% Chlorophyll"=17))+
  guides(color=guide_legend("Parameter"),shape=guide_legend("Parameter"))+
  theme(legend.text=element_text(size=12),legend.title=element_text(size=14))+
  theme(axis.text.x=element_text(size=12),axis.title.x=element_text(size=14))+
  theme(axis.text.y=element_text(size=12),axis.title.y=element_text(size=13))+
  theme(legend.position="right")+
  ylim(0,60)+
  ylab("CPE (\U00B5g/cm3), % Chl A")
ChlPhafig
#ggsave("Fig6_ChlPha.png", width = 6, height = 5, units=c("in"), dpi=300)

CPEfig<-ggplot(data=HGI.env.fig,aes(x=factor(Year),y=V12))+
  geom_bar(stat="identity",fill="black",alpha=0.7)+
  theme_classic()+ xlab("Year")+ylab("% Chl A")
CPEfig

Pigm<-matrix(,nrow=42,ncol=3)
Pigm<-data.frame(Pigm)
Pigm[,1]<-HGI.env.fig$Year
Pigm[1:21,2]<-HGI.env.fig$ChlA
Pigm[22:42,2]<-HGI.env.fig$CPE-HGI.env.fig$ChlA
Pigm[1:21,3]<-c("ChlA")
Pigm[22:42,3]<-c("Rest")
colnames(Pigm)<-c("Year","value","parameter")

Chlfig<-ggplot(data=Pigm,aes(fill=parameter,x=factor(Year),y=value))+
  geom_bar(position=position_stack(reverse=T),stat="identity")+
  scale_fill_manual(values=c("ChlA"="darkgoldenrod1","Rest"="cornflowerblue"))+
  theme_classic()+ xlab("Year") + ylab(expression(paste(plain("Pigment (\U00B5g/cm3)"))))+
  theme(legend.position="right")+
  theme(legend.text=element_text(size=12),legend.title=element_text(size=14))+
  theme(axis.text.x=element_text(size=12),axis.title.x=element_text(size=14))+
  theme(axis.text.y=element_text(size=12),axis.title.y=element_text(size=14))
Chlfig
ggsave("Fig6_ChlPha.png", width = 6, height = 4, units=c("in"), dpi=300)

Envfig2<-ggplot(data=HGI.env.fig,aes(x=Year))+
  geom_line(aes(y=NAO,color="N. Atl. Oscill."),linewidth=1)+
  geom_point(aes(y=NAO,color="N. Atl. Oscill.",shape="N. Atl. Oscill."),size=3)+
  geom_line(aes(y=AO,color="Arctic Oscill."),linewidth=1)+
  geom_point(aes(y=AO,color="Arctic Oscill.",shape="Arctic Oscill."),size=3)+
  geom_smooth(aes(y=AO),method=lm,se=F,linetype="dashed",color="#88CCEE")+
  theme_classic()+ylab("Index value")+
  scale_color_manual(values=c("N. Atl. Oscill."="#CC6677","Arctic Oscill."="#88CCEE"))+
  scale_shape_manual(values=c("N. Atl. Oscill."=16,"Arctic Oscill."=15))+
  guides(shape=guide_legend("Index"),color=guide_legend("Index"))+
  theme(legend.text=element_text(size=12),legend.title=element_text(size=14))+
  theme(axis.text.x=element_text(size=12),axis.title.x=element_text(size=14))+
  theme(axis.text.y=element_text(size=12),axis.title.y=element_text(size=13))+
  theme(legend.position="right")
Envfig2
#ggsave("Fig7_NAOAO.png", width = 6, height = 4, units=c("in"), dpi=300)

############################################################################
##Combine panels into figures
SNfig<-ggarrange(abund.bar.5,abund.bar.1,rich.bar.5,rich.bar.1,ncol=2,nrow=2)
SNfig+annotate("text",x=0.1,y=0.98,label="A",size=5)+annotate("text",x=0.615,y=0.98,label="B",size=5)+
  annotate("text",x=0.1,y=0.48,label="C",size=5)+annotate("text",x=0.6,y=0.48,label="D",size=5)+
  annotate("text",x=0.13,y=0.71,label="a",size=4)+annotate("text",x=0.21,y=0.69,label="a",size=4)+
  annotate("text",x=0.29,y=0.94,label="b",size=4)+annotate("text",x=0.365,y=0.93,label="b",size=4)+
  annotate("text",x=0.442,y=0.835,label="c",size=4)+
  annotate("text",x=0.643,y=0.944,label="a",size=4)+annotate("text",x=0.718,y=0.88,label="b",size=4)+
  annotate("text",x=0.796,y=0.88,label="b",size=4)+annotate("text",x=0.868,y=0.87,label="b",size=4)+
  annotate("text",x=0.946,y=0.841,label="c",size=4)+
  annotate("text",x=0.125,y=0.278,label="a",size=4)+annotate("text",x=0.203,y=0.289,label="a",size=4)+
  annotate("text",x=0.285,y=0.43,label="b",size=4)+annotate("text",x=0.361,y=0.436,label="b",size=4)+
  annotate("text",x=0.44,y=0.323,label="c",size=4)+
  annotate("text",x=0.624,y=0.445,label="a",size=4)+annotate("text",x=0.706,y=0.38,label="b",size=4)+
  annotate("text",x=0.784,y=0.43,label="ac",size=4)+annotate("text",x=0.862,y=0.413,label="c",size=4)+
  annotate("text",x=0.944,y=0.35,label="d",size=4)
ggsave("Fig3_TotalS_TotalN.png",width=7,height=7,units=c("in"),dpi=300)

Trophicfig<-ggarrange(troph.bar.5,troph.bar.1,common.legend=T,legend="right",ncol=2,nrow=1)
Trophicfig+annotate("text",x=0.1,y=0.96,label="A",size=5)+annotate("text",x=0.49,y=0.96,label="B",size=5)
ggsave("Fig4_Trophic.png",width=10,height=4,units=c("in"),bg="white",dpi=300)

dbRDAfig<-ggarrange(Fig5A,Fig5B,ncol=2,nrow=1,labels=c("A","B"))
dbRDAfig
ggsave("Fig5_dbRDA.png",width=10,height=4,units=c("in"),bg="white",dpi=300)

Allfactors<-ggarrange(ChlPhafig,Envfig,Envfig2,ncol=1,nrow=3,common.legend=F,labels=c("A","B","C"))
Allfactors
ggsave("Fig6_Allenv.png",width=6,height=9,units=c("in"),bg="white",dpi=300)

############################################################################
############################################################################
##SUPPLEMENTARY MATERIAL

dpath<-file.choose() #Long-term_trophic
HGI.trophic.5<-read.csv(dpath,row.names=1)

#Suspension feeders
sus.aov<-aov(HGI.trophic.5$SF~as.character(HGI.years.5), data=HGI.trophic.5)
summary(sus.aov)
sus.aov.factor <- aov(HGI.trophic.5$SF~as.factor(HGI.years.5), data=HGI.trophic.5)
sus.aov.tukey<-TukeyHSD(sus.aov.factor)
sus.aov.tukey
fit<-lm(HGI.trophic.5$SF~CPE+Protein+NAO+AO+Bot_temp+Bot_sal+Surf_temp+Surf_sal+C.org+PropChl,data=HGI.env.5)
step<-stepAIC(fit,direction="both")
anova(step)
RsquareAdj(step)

#Deposit feeders
dep.aov<-aov(HGI.trophic.5$DF~as.character(HGI.years.5), data=HGI.trophic.5)
summary(dep.aov)
dep.aov.factor <- aov(HGI.trophic.5$DF~as.factor(HGI.years.5), data=HGI.trophic.5)
dep.aov.tukey<-TukeyHSD(dep.aov.factor)
dep.aov.tukey
fit<-lm(HGI.trophic.5$DF~CPE+Protein+NAO+AO+Bot_temp+Bot_sal+Surf_temp+Surf_sal+C.org+PropChl,data=HGI.env.5)
step<-stepAIC(fit,direction="both")
anova(step)
RsquareAdj(step)

#Predator/Scavengers
pre.aov<-aov(HGI.trophic.5$PS~as.character(HGI.years.5), data=HGI.trophic.5)
summary(pre.aov)
pre.aov.factor <- aov(HGI.trophic.5$PS~as.factor(HGI.years.5), data=HGI.trophic.5)
pre.aov.tukey<-TukeyHSD(pre.aov.factor)
pre.aov.tukey
fit<-lm(HGI.trophic.5$PS~CPE+Protein+NAO+AO+Bot_temp+Bot_sal+Surf_temp+Surf_sal+C.org+PropChl,data=HGI.env.5)
step<-stepAIC(fit,direction="both")
anova(step)
RsquareAdj(step)

#Aggregate data for stacked bar
sus_sum <- aggregate(HGI.trophic.5$SF~HGI.years.5, FUN = mean)
colnames(sus_sum) = c("Year","mean_density") 
sus_sum_se <- aggregate(HGI.trophic.5$SF~HGI.years.5, FUN = stand_err)
sus_sum$se <- sus_sum_se$`HGI.trophic.5$SF`
sus_sum$Group<-c("Suspension feeders","Suspension feeders","Suspension feeders","Suspension feeders","Suspension feeders")

dep_sum <- aggregate(HGI.trophic.5$DF~HGI.years.5, FUN = mean)
colnames(dep_sum) = c("Year","mean_density") 
dep_sum_se <- aggregate(HGI.trophic.5$DF~HGI.years.5, FUN = stand_err)
dep_sum$se <- dep_sum_se$`HGI.trophic.5$DF`
dep_sum$Group<-c("Deposit feeders","Deposit feeders","Deposit feeders","Deposit feeders","Deposit feeders")

pre_sum <- aggregate(HGI.trophic.5$PS~HGI.years.5, FUN = mean)
colnames(pre_sum) = c("Year","mean_density") 
pre_sum_se <- aggregate(HGI.trophic.5$PS~HGI.years.5, FUN = stand_err)
pre_sum$se <- pre_sum_se$`HGI.trophic.5$PS`
pre_sum$Group<-c("Predators/Scavengers","Predators/Scavengers","Predators/Scavengers","Predators/Scavengers","Predators/Scavengers")

HGI.trophic.5$UK<-numeric(450)
uk_sum <- aggregate(HGI.trophic.5$UK~HGI.years.5, FUN = mean)
colnames(uk_sum) = c("Year","mean_density") 
uk_sum_se <- aggregate(HGI.trophic.5$UK~HGI.years.5, FUN = stand_err)
uk_sum$se <- uk_sum_se$`HGI.trophic.5$UK`
uk_sum$Group<-c("Unknown","Unknown","Unknown","Unknown","Unknown")

#Stacked bar graph for trophic groups
troph_sum <- rbind(sus_sum,dep_sum,pre_sum,uk_sum)
troph.bar.5<- ggplot(data=troph_sum, aes(fill=Group,x=factor(Year),y=mean_density))+
  geom_bar(position="fill",stat="identity")+
  scale_fill_manual(values=c("Suspension feeders"="darkgoldenrod1","Deposit feeders"="cornflowerblue","Predators/Scavengers"="chocolate4","Unknown"="gray"))+
  theme_classic()+ xlab("Year") + ylab(expression(paste(plain("Percent of community"))))+
  scale_y_continuous(labels=scales::percent)+
  theme(legend.position="right")+
  theme(legend.text=element_text(size=12),legend.title=element_text(size=14))+
  theme(axis.text.x=element_text(size=12),axis.title.x=element_text(size=14))+
  theme(axis.text.y=element_text(size=12),axis.title.y=element_text(size=14))
troph.bar.5

#Cerianthid
cer.aov<-aov(HGI.5year$Cerianthid~as.character(HGI.years.5$Year), data=HGI.5year)
summary(cer.aov)
cer.aov.factor <- aov(HGI.5year$Cerianthid~as.factor(HGI.years.5$Year), data=HGI.5year)
cer.aov.tukey<-TukeyHSD(cer.aov.factor)
cer.aov.tukey

#Jasminiera schaudinni
jas.aov<-aov(HGI.5year$Jasminiera~as.character(HGI.years.5$Year), data=HGI.5year)
summary(jas.aov)
jas.aov.factor <- aov(HGI.5year$Jasminiera~as.factor(HGI.years.5$Year), data=HGI.5year)
jas.aov.tukey<-TukeyHSD(jas.aov.factor)
jas.aov.tukey

#Bylgides groenlandicus
byl.aov<-aov(HGI.5year$Byglides~as.character(HGI.years.5$Year), data=HGI.5year)
summary(byl.aov)
byl.aov.factor <- aov(HGI.5year$Byglides~as.factor(HGI.years.5$Year), data=HGI.5year)
byl.aov.tukey<-TukeyHSD(byl.aov.factor)
byl.aov.tukey

#Bythocaris leucopis
byt.aov<-aov(HGI.5year$Bythocaris~as.character(HGI.years.5$Year), data=HGI.5year)
summary(byt.aov)
byt.aov.factor <- aov(HGI.5year$Bythocaris~as.factor(HGI.years.5$Year), data=HGI.5year)
byt.aov.tukey<-TukeyHSD(byt.aov.factor)
byt.aov.tukey

#Colossendeis proboscidea
col.aov<-aov(HGI.5year$Colossendeis~as.character(HGI.years.5$Year), data=HGI.5year)
summary(col.aov)
col.aov.factor <- aov(HGI.5year$Colossendeis~as.factor(HGI.years.5$Year), data=HGI.5year)
col.aov.tukey<-TukeyHSD(col.aov.factor)
col.aov.tukey

#Nymphon macronyx
nym.aov<-aov(HGI.5year$Nymphon~as.character(HGI.years.5$Year), data=HGI.5year)
summary(nym.aov)
nym.aov.factor <- aov(HGI.5year$Nymphon~as.factor(HGI.years.5$Year), data=HGI.5year)
nym.aov.tukey<-TukeyHSD(nym.aov.factor)
nym.aov.tukey

#Bathybiaster vexillifer
bat.aov<-aov(HGI.5year$Bathybiaster~as.character(HGI.years.5$Year), data=HGI.5year)
summary(bat.aov)
bat.aov.factor <- aov(HGI.5year$Bathybiaster~as.factor(HGI.years.5$Year), data=HGI.5year)
bat.aov.tukey<-TukeyHSD(bat.aov.factor)
bat.aov.tukey

#Ophiocten gracilis
oph.aov<-aov(HGI.5year$Ophiocten~as.character(HGI.years.5$Year), data=HGI.5year)
summary(oph.aov)
oph.aov.factor <- aov(HGI.5year$Ophiocten~as.factor(HGI.years.5$Year), data=HGI.5year)
oph.aov.tukey<-TukeyHSD(oph.aov.factor)
oph.aov.tukey

#Elpidia heckeri
elp.aov<-aov(HGI.5year$Elpidia~as.character(HGI.years.5$Year), data=HGI.5year)
summary(elp.aov)
elp.aov.factor <- aov(HGI.5year$Elpidia~as.factor(HGI.years.5$Year), data=HGI.5year)
elp.aov.tukey<-TukeyHSD(elp.aov.factor)
elp.aov.tukey

#Lycodes squamiventer
lyc.aov<-aov(HGI.5year$Lycodes~as.character(HGI.years.5$Year), data=HGI.5year)
summary(lyc.aov)
lyc.aov.factor <- aov(HGI.5year$Lycodes~as.factor(HGI.years.5$Year), data=HGI.5year)
lyc.aov.tukey<-TukeyHSD(lyc.aov.factor)
lyc.aov.tukey

#Lycodonus flagellicauda
lyco.aov<-aov(HGI.5year$Lycodonus~as.character(HGI.years.5$Year), data=HGI.5year)
summary(lyco.aov)
lyco.aov.factor <- aov(HGI.5year$Lycodonus~as.factor(HGI.years.5$Year), data=HGI.5year)
lyco.aov.tukey<-TukeyHSD(lyco.aov.factor)
lyco.aov.tukey

#Mohnia mohni
moh.aov<-aov(HGI.5year$Mohnia~as.character(HGI.years.5$Year), data=HGI.5year)
summary(moh.aov)
moh.aov.factor <- aov(HGI.5year$Mohnia~as.factor(HGI.years.5$Year), data=HGI.5year)
moh.aov.tukey<-TukeyHSD(moh.aov.factor)
moh.aov.tukey

##Bar graphs
#Suspension feeders
sus_sum <- aggregate(HGI.trophic.5$SF~HGI.years.5$Year, FUN = mean)
colnames(sus_sum) = c("Year","mean_density") 
sus_sum_se <- aggregate(HGI.trophic.5$SF~HGI.years.5$Year, FUN = stand_err)
sus_sum$se <- sus_sum_se$`HGI.trophic.5$SF`
sus_sum$Group<-c("Suspension feeders","Suspension feeders","Suspension feeders","Suspension feeders","Suspension feeders")
sus.bar<- ggplot(data=sus_sum, aes(x=factor(Year),y=mean_density)) +
  geom_bar(stat="identity", fill="purple3", alpha=0.7)  +
  geom_errorbar(aes(ymin=mean_density-se, ymax=mean_density+se, width=.2)) + 
  theme_classic()+ xlab("Year") + ylab(expression(paste(plain("Suspension feeders")," m"^plain(-2))))
sus.bar
ggsave("Bar5_SF.png", width = 4, height = 4, units=c("in"), dpi=300)

#Deposit feeders
dep_sum <- aggregate(HGI.trophic.5$DF~HGI.years.5$Year, FUN = mean)
colnames(dep_sum) = c("Year","mean_density") 
dep_sum_se <- aggregate(HGI.trophic.5$DF~HGI.years.5$Year, FUN = stand_err)
dep_sum$se <- dep_sum_se$`HGI.trophic.5$DF`
dep_sum$Group<-c("Deposit feeders","Deposit feeders","Deposit feeders","Deposit feeders","Deposit feeders")
dep.bar<- ggplot(data=dep_sum, aes(x=factor(Year),y=mean_density)) +
  geom_bar(stat="identity", fill="orangered1", alpha=0.7)  +
  geom_errorbar(aes(ymin=mean_density-se, ymax=mean_density+se, width=.2)) + 
  theme_classic()+ xlab("Year") + ylab(expression(paste(plain("Deposit feeders")," m"^plain(-2))))
dep.bar
ggsave("Bar5_DF.png", width = 4, height = 4, units=c("in"), dpi=300)

#Predator/Scavengers
pre_sum <- aggregate(HGI.trophic.5$PS~HGI.years.5$Year, FUN = mean)
colnames(pre_sum) = c("Year","mean_density") 
pre_sum_se <- aggregate(HGI.trophic.5$PS~HGI.years.5$Year, FUN = stand_err)
pre_sum$se <- pre_sum_se$`HGI.trophic.5$PS`
pre_sum$Group<-c("Predators/Scavengers","Predators/Scavengers","Predators/Scavengers","Predators/Scavengers","Predators/Scavengers")
pre.bar<- ggplot(data=pre_sum, aes(x=factor(Year),y=mean_density)) +
  geom_bar(stat="identity", fill="turquoise4", alpha=0.7)  +
  geom_errorbar(aes(ymin=mean_density-se, ymax=mean_density+se, width=.2)) + 
  theme_classic()+ xlab("Year") + ylab(expression(paste(plain("Predator/Scavengers")," m"^plain(-2))))
pre.bar
ggsave("Bar5_PS.png", width = 4, height = 4, units=c("in"), dpi=300)

#Cerianthid
cer_sum <- aggregate(HGI.5year$Cerianthid~HGI.years.5$Year, FUN = mean)
colnames(cer_sum) = c("Year","mean_density") 
cer_sum_se <- aggregate(HGI.5year$Cerianthid~HGI.years.5$Year, FUN = stand_err)
cer_sum$se <- cer_sum_se$`HGI.5year$Cerianthid`
cer.bar<- ggplot(data=cer_sum, aes(x=factor(Year),y=mean_density)) +
  geom_bar(stat="identity", fill="purple3", alpha=0.7)  +
  geom_errorbar(aes(ymin=mean_density-se, ymax=mean_density+se, width=.2)) + 
  theme_classic()+ xlab("Year") + ylab(expression(paste(plain("Cerianthids")," m"^plain(-2))))
cer.bar
ggsave("Bar5_Cerianthid.png", width = 4, height = 4, units=c("in"), dpi=300)

#Jasminiera schaudinni
jas_sum <- aggregate(HGI.5year$Jasminiera~HGI.years.5$Year, FUN = mean)
colnames(jas_sum) = c("Year","mean_density") 
jas_sum_se <- aggregate(HGI.5year$Jasminiera~HGI.years.5$Year, FUN = stand_err)
jas_sum$se <- jas_sum_se$`HGI.5year$Jasminiera`
jas.bar<- ggplot(data=jas_sum, aes(x=factor(Year),y=mean_density)) +
  geom_bar(stat="identity", fill="purple3", alpha=0.7)  +
  geom_errorbar(aes(ymin=mean_density-se, ymax=mean_density+se, width=.2)) + 
  theme_classic()+ xlab("Year") + ylab(expression(paste(italic("Jasminiera schaudinni")," m"^plain(-2))))
jas.bar
ggsave("Bar5_Jasminiera.png", width = 4, height = 4, units=c("in"), dpi=300)

#Bylgides groenlandicus
byl_sum <- aggregate(HGI.5year$Byglides~HGI.years.5$Year, FUN = mean)
colnames(byl_sum) = c("Year","mean_density") 
byl_sum_se <- aggregate(HGI.5year$Byglides~HGI.years.5$Year, FUN = stand_err)
byl_sum$se <- byl_sum_se$`HGI.5year$Byglides`
byl.bar<- ggplot(data=byl_sum, aes(x=factor(Year),y=mean_density)) +
  geom_bar(stat="identity", fill="turquoise4", alpha=0.7)  +
  geom_errorbar(aes(ymin=mean_density-se, ymax=mean_density+se, width=.2)) + 
  theme_classic()+ xlab("Year") + ylab(expression(paste(italic("Bylgides groenlandicus")," m"^plain(-2))))
byl.bar
ggsave("Bar5_Bylgides.png", width = 4, height = 4, units=c("in"), dpi=300)

#Bythocaris leucopis
byt_sum <- aggregate(HGI.5year$Bythocaris~HGI.years.5$Year, FUN = mean)
colnames(byt_sum) = c("Year","mean_density") 
byt_sum_se <- aggregate(HGI.5year$Bythocaris~HGI.years.5$Year, FUN = stand_err)
byt_sum$se <- byt_sum_se$`HGI.5year$Bythocaris`
byt.bar<- ggplot(data=byt_sum, aes(x=factor(Year),y=mean_density)) +
  geom_bar(stat="identity", fill="turquoise4", alpha=0.7)  +
  geom_errorbar(aes(ymin=mean_density-se, ymax=mean_density+se, width=.2)) + 
  theme_classic()+ xlab("Year") + ylab(expression(paste(italic("Bythocaris leucopis")," m"^plain(-2))))
byt.bar
ggsave("Bar5_Bythocaris.png", width = 4, height = 4, units=c("in"), dpi=300)

#Colossendeis proboscidea
col_sum <- aggregate(HGI.5year$Colossendeis~HGI.years.5$Year, FUN = mean)
colnames(col_sum) = c("Year","mean_density") 
col_sum_se <- aggregate(HGI.5year$Colossendeis~HGI.years.5$Year, FUN = stand_err)
col_sum$se <- col_sum_se$`HGI.5year$Colossendeis`
col.bar<- ggplot(data=col_sum, aes(x=factor(Year),y=mean_density)) +
  geom_bar(stat="identity", fill="orangered1", alpha=0.7)  +
  geom_errorbar(aes(ymin=mean_density-se, ymax=mean_density+se, width=.2)) + 
  theme_classic()+ xlab("Year") + ylab(expression(paste(italic("Colossendeis proboscidea")," m"^plain(-2))))
col.bar
ggsave("Bar5_Colossendeis.png", width = 4, height = 4, units=c("in"), dpi=300)

#Nymphon macronyx
nym_sum <- aggregate(HGI.5year$Nymphon~HGI.years.5$Year, FUN = mean)
colnames(nym_sum) = c("Year","mean_density") 
nym_sum_se <- aggregate(HGI.5year$Nymphon~HGI.years.5$Year, FUN = stand_err)
nym_sum$se <- nym_sum_se$`HGI.5year$Nymphon`
nym.bar<- ggplot(data=nym_sum, aes(x=factor(Year),y=mean_density)) +
  geom_bar(stat="identity", fill="turquoise4", alpha=0.7)  +
  geom_errorbar(aes(ymin=mean_density-se, ymax=mean_density+se, width=.2)) + 
  theme_classic()+ xlab("Year") + ylab(expression(paste(italic("Nymphon macronyx")," m"^plain(-2))))
nym.bar
ggsave("Bar5_Nymphon.png", width = 4, height = 4, units=c("in"), dpi=300)

#Bathybiaster vexillifer  
bat_sum <- aggregate(HGI.5year$Bathybiaster~HGI.years.5$Year, FUN = mean)
colnames(bat_sum) = c("Year","mean_density") 
bat_sum_se <- aggregate(HGI.5year$Bathybiaster~HGI.years.5$Year, FUN = stand_err)
bat_sum$se <- bat_sum_se$`HGI.5year$Bathybiaster`
bat.bar<- ggplot(data=bat_sum, aes(x=factor(Year),y=mean_density)) +
  geom_bar(stat="identity", fill="turquoise4", alpha=0.7)  +
  geom_errorbar(aes(ymin=mean_density-se, ymax=mean_density+se, width=.2)) + 
  theme_classic()+ xlab("Year") + ylab(expression(paste(italic("Bathybiaster vexillifer")," m"^plain(-2))))
bat.bar
ggsave("Bar5_Bathybiaster.png", width = 4, height = 4, units=c("in"), dpi=300)

#Ophiocten gracilis
oph_sum <- aggregate(HGI.5year$Ophiocten~HGI.years.5$Year, FUN = mean)
colnames(oph_sum) = c("Year","mean_density") 
oph_sum_se <- aggregate(HGI.5year$Ophiocten~HGI.years.5$Year, FUN = stand_err)
oph_sum$se <- oph_sum_se$`HGI.5year$Ophiocten`
oph.bar<- ggplot(data=oph_sum, aes(x=factor(Year),y=mean_density)) +
  geom_bar(stat="identity", fill="orangered1", alpha=0.7)  +
  geom_errorbar(aes(ymin=mean_density-se, ymax=mean_density+se, width=.2)) + 
  theme_classic()+ xlab("Year") + ylab(expression(paste(italic("Ophiocten gracilis")," m"^plain(-2))))
oph.bar
ggsave("Bar5_Ophiocten.png", width = 4, height = 4, units=c("in"), dpi=300)

#Elpidia heckeri
elp_sum <- aggregate(HGI.5year$Elpidia~HGI.years.5$Year, FUN = mean)
colnames(elp_sum) = c("Year","mean_density") 
elp_sum_se <- aggregate(HGI.5year$Elpidia~HGI.years.5$Year, FUN = stand_err)
elp_sum$se <- elp_sum_se$`HGI.5year$Elpidia`
elp.bar<- ggplot(data=elp_sum, aes(x=factor(Year),y=mean_density)) +
  geom_bar(stat="identity", fill="orangered1", alpha=0.7)  +
  geom_errorbar(aes(ymin=mean_density-se, ymax=mean_density+se, width=.2)) + 
  theme_classic()+ xlab("Year") + ylab(expression(paste(italic("Elpidia heckeri")," m"^plain(-2))))
elp.bar
ggsave("Bar5_Elpidia.png", width = 4, height = 4, units=c("in"), dpi=300)

#Lycodes squamiventer
lyc_sum <- aggregate(HGI.5year$Lycodes~HGI.years.5$Year, FUN = mean)
colnames(lyc_sum) = c("Year","mean_density") 
lyc_sum_se <- aggregate(HGI.5year$Lycodes~HGI.years.5$Year, FUN = stand_err)
lyc_sum$se <- lyc_sum_se$`HGI.5year$Lycodes`
lyc.bar<- ggplot(data=lyc_sum, aes(x=factor(Year),y=mean_density)) +
  geom_bar(stat="identity", fill="turquoise4", alpha=0.7)  +
  geom_errorbar(aes(ymin=mean_density-se, ymax=mean_density+se, width=.2)) + 
  theme_classic()+ xlab("Year") + ylab(expression(paste(italic("Lycodes squamiventer")," m"^plain(-2))))
lyc.bar
ggsave("Bar5_Lycodes.png", width = 4, height = 4, units=c("in"), dpi=300)

#Lycodonus flagellicauda
lyco_sum <- aggregate(HGI.5year$Lycodonus~HGI.years.5$Year, FUN = mean)
colnames(lyco_sum) = c("Year","mean_density") 
lyco_sum_se <- aggregate(HGI.5year$Lycodonus~HGI.years.5$Year, FUN = stand_err)
lyco_sum$se <- lyco_sum_se$`HGI.5year$Lycodonus`
lyco.bar<- ggplot(data=lyco_sum, aes(x=factor(Year),y=mean_density)) +
  geom_bar(stat="identity", fill="turquoise4", alpha=0.7)  +
  geom_errorbar(aes(ymin=mean_density-se, ymax=mean_density+se, width=.2)) + 
  theme_classic()+ xlab("Year") + ylab(expression(paste(italic("Lycodonus flagellicauda")," m"^plain(-2))))
lyco.bar
ggsave("Bar5_Lycodonus.png", width = 4, height = 4, units=c("in"), dpi=300)

#Mohnia mohni
moh_sum <- aggregate(HGI.5year$Mohnia~HGI.years.5$Year, FUN = mean)
colnames(moh_sum) = c("Year","mean_density") 
moh_sum_se <- aggregate(HGI.5year$Mohnia~HGI.years.5$Year, FUN = stand_err)
moh_sum$se <- moh_sum_se$`HGI.5year$Mohnia`
moh.bar<- ggplot(data=moh_sum, aes(x=factor(Year),y=mean_density)) +
  geom_bar(stat="identity", fill="turquoise4", alpha=0.7)  +
  geom_errorbar(aes(ymin=mean_density-se, ymax=mean_density+se, width=.2)) + 
  theme_classic()+ xlab("Year") + ylab(expression(paste(italic("Mohnia mohni")," m"^plain(-2))))
moh.bar
ggsave("Bar5_Mohnia.png", width = 4, height = 4, units=c("in"), dpi=300)

############################################################################

epath<-file.choose() #Short-term_trophic
HGI.trophic.1<-read.csv(epath,row.names=1)

#Suspension feeders
sus.aov<-aov(HGI.trophic.1$SF~as.character(HGI.years.1), data=HGI.trophic.1)
summary(sus.aov)
sus.aov.factor <- aov(HGI.trophic.1$SF~as.factor(HGI.years.1), data=HGI.trophic.1)
sus.aov.tukey<-TukeyHSD(sus.aov.factor)
sus.aov.tukey
fit<-lm(HGI.trophic.1$SF~CPE+Protein+NAO+AO+Bot_temp+Bot_sal+Surf_temp+Surf_sal+C.org+PropChl,data=HGI.env.1)
step<-stepAIC(fit,direction="both")
anova(step)
RsquareAdj(step)

#Deposit feeders
dep.aov<-aov(HGI.trophic.1$DF~as.character(HGI.years.1), data=HGI.trophic.1)
summary(dep.aov)
dep.aov.factor <- aov(HGI.trophic.1$DF~as.factor(HGI.years.1), data=HGI.trophic.1)
dep.aov.tukey<-TukeyHSD(dep.aov.factor)
dep.aov.tukey
fit<-lm(HGI.trophic.1$DF~CPE+Protein+NAO+AO+Bot_temp+Bot_sal+Surf_temp+Surf_sal+C.org+PropChl,data=HGI.env.1)
step<-stepAIC(fit,direction="both")
anova(step)
RsquareAdj(step)

#Predator/Scavengers
pre.aov<-aov(HGI.trophic.1$PS~as.character(HGI.years.1), data=HGI.trophic.1)
summary(pre.aov)
pre.aov.factor <- aov(HGI.trophic.1$PS~as.factor(HGI.years.1), data=HGI.trophic.1)
pre.aov.tukey<-TukeyHSD(pre.aov.factor)
pre.aov.tukey
fit<-lm(HGI.trophic.1$PS~CPE+Protein+NAO+AO+Bot_temp+Bot_sal+Surf_temp+Surf_sal+C.org+PropChl,data=HGI.env.1)
step<-stepAIC(fit,direction="both")
anova(step)
RsquareAdj(step)

#Aggregate data for stacked bar
sus_sum <- aggregate(HGI.trophic.1$SF~HGI.years.1, FUN = mean)
colnames(sus_sum) = c("Year","mean_density") 
sus_sum_se <- aggregate(HGI.trophic.1$SF~HGI.years.1, FUN = stand_err)
sus_sum$se <- sus_sum_se$`HGI.trophic.1$SF`
sus_sum$Group<-c("Suspension feeders","Suspension feeders","Suspension feeders","Suspension feeders","Suspension feeders")

dep_sum <- aggregate(HGI.trophic.1$DF~HGI.years.1, FUN = mean)
colnames(dep_sum) = c("Year","mean_density") 
dep_sum_se <- aggregate(HGI.trophic.1$DF~HGI.years.1, FUN = stand_err)
dep_sum$se <- dep_sum_se$`HGI.trophic.1$DF`
dep_sum$Group<-c("Deposit feeders","Deposit feeders","Deposit feeders","Deposit feeders","Deposit feeders")

pre_sum <- aggregate(HGI.trophic.1$PS~HGI.years.1, FUN = mean)
colnames(pre_sum) = c("Year","mean_density") 
pre_sum_se <- aggregate(HGI.trophic.1$PS~HGI.years.1, FUN = stand_err)
pre_sum$se <- pre_sum_se$`HGI.trophic.1$PS`
pre_sum$Group<-c("Predators/Scavengers","Predators/Scavengers","Predators/Scavengers","Predators/Scavengers","Predators/Scavengers")

uk_sum <- aggregate(HGI.trophic.1$UK~HGI.years.1, FUN = mean)
colnames(uk_sum) = c("Year","mean_density") 
uk_sum_se <- aggregate(HGI.trophic.1$UK~HGI.years.1, FUN = stand_err)
uk_sum$se <- uk_sum_se$`HGI.trophic.1$UK`
uk_sum$Group<-c("Unknown","Unknown","Unknown","Unknown","Unknown")

troph_sum <- rbind(sus_sum,dep_sum,pre_sum,uk_sum)
troph.bar.1<- ggplot(data=troph_sum, aes(fill=Group,x=factor(Year),y=mean_density))+
  geom_bar(position="fill",stat="identity")+
  scale_fill_manual(values=c("Suspension feeders"="darkgoldenrod1","Deposit feeders"="cornflowerblue","Predators/Scavengers"="chocolate4","Unknown"="grey"))+
  theme_classic()+ xlab("Year") + ylab(expression(paste(plain("Percent of community"))))+
  scale_y_continuous(labels=scales::percent)+
  guides(color=guide_legend("Group"))+
  theme(legend.text=element_text(size=12),legend.title=element_text(size=14))+
  theme(axis.text.x=element_text(size=12),axis.title.x=element_text(size=14))+
  theme(axis.text.y=element_text(size=12),axis.title.y=element_text(size=14))
troph.bar.1

#Cerianthid
cer.aov<-aov(HGI.1year$Cerianthid~as.character(HGI.years.1$Year), data=HGI.1year)
summary(cer.aov)
cer.aov.factor <- aov(HGI.1year$Cerianthid~as.factor(HGI.years.1$Year), data=HGI.1year)
cer.aov.tukey<-TukeyHSD(cer.aov.factor)
cer.aov.tukey

#Jasminiera schaudinni
jas.aov<-aov(HGI.1year$Jasminiera~as.character(HGI.years.1$Year), data=HGI.1year)
summary(jas.aov)
jas.aov.factor <- aov(HGI.1year$Jasminiera~as.factor(HGI.years.1$Year), data=HGI.1year)
jas.aov.tukey<-TukeyHSD(jas.aov.factor)
jas.aov.tukey

#Bylgides groenlandicus
byl.aov<-aov(HGI.1year$Byglides~as.character(HGI.years.1$Year), data=HGI.1year)
summary(byl.aov)
byl.aov.factor <- aov(HGI.1year$Byglides~as.factor(HGI.years.1$Year), data=HGI.1year)
byl.aov.tukey<-TukeyHSD(byl.aov.factor)
byl.aov.tukey

#Bythocaris leucopis
byt.aov<-aov(HGI.1year$Bythocaris~as.character(HGI.years.1$Year), data=HGI.1year)
summary(byt.aov)
byt.aov.factor <- aov(HGI.1year$Bythocaris~as.factor(HGI.years.1$Year), data=HGI.1year)
byt.aov.tukey<-TukeyHSD(byt.aov.factor)
byt.aov.tukey

#Colossendeis proboscidea
col.aov<-aov(HGI.1year$Colossendeis~as.character(HGI.years.1$Year), data=HGI.1year)
summary(col.aov)
col.aov.factor <- aov(HGI.1year$Colossendeis~as.factor(HGI.years.1$Year), data=HGI.1year)
col.aov.tukey<-TukeyHSD(col.aov.factor)
col.aov.tukey

#Nymphon macronyx
nym.aov<-aov(HGI.1year$Nymphon~as.character(HGI.years.1$Year), data=HGI.1year)
summary(nym.aov)
nym.aov.factor <- aov(HGI.1year$Nymphon~as.factor(HGI.years.1$Year), data=HGI.1year)
nym.aov.tukey<-TukeyHSD(nym.aov.factor)
nym.aov.tukey

#Bathybiaster vexillifer
bat.aov<-aov(HGI.1year$Bathybiaster~as.character(HGI.years.1$Year), data=HGI.1year)
summary(bat.aov)
bat.aov.factor <- aov(HGI.1year$Bathybiaster~as.factor(HGI.years.1$Year), data=HGI.1year)
bat.aov.tukey<-TukeyHSD(bat.aov.factor)
bat.aov.tukey

#Ophiocten gracilis
oph.aov<-aov(HGI.1year$Ophiocten~as.character(HGI.years.1$Year), data=HGI.1year)
summary(oph.aov)
oph.aov.factor <- aov(HGI.1year$Ophiocten~as.factor(HGI.years.1$Year), data=HGI.1year)
oph.aov.tukey<-TukeyHSD(oph.aov.factor)
oph.aov.tukey

#Elpidia heckeri
elp.aov<-aov(HGI.1year$Elpidia~as.character(HGI.years.1$Year), data=HGI.1year)
summary(elp.aov)
elp.aov.factor <- aov(HGI.1year$Elpidia~as.factor(HGI.years.1$Year), data=HGI.1year)
elp.aov.tukey<-TukeyHSD(elp.aov.factor)
elp.aov.tukey

#Lycodes squamiventer
lyc.aov<-aov(HGI.1year$Lycodes~as.character(HGI.years.1$Year), data=HGI.1year)
summary(lyc.aov)
lyc.aov.factor <- aov(HGI.1year$Lycodes~as.factor(HGI.years.1$Year), data=HGI.1year)
lyc.aov.tukey<-TukeyHSD(lyc.aov.factor)
lyc.aov.tukey

#Mohnia mohni
moh.aov<-aov(HGI.1year$Mohnia~as.character(HGI.years.1$Year), data=HGI.1year)
summary(moh.aov)
moh.aov.factor <- aov(HGI.1year$Mohnia~as.factor(HGI.years.1$Year), data=HGI.1year)
moh.aov.tukey<-TukeyHSD(moh.aov.factor)
moh.aov.tukey

#Pontaster tennuispinus
pon.aov<-aov(HGI.1year$Pontaster~as.character(HGI.years.1$Year), data=HGI.1year)
summary(pon.aov)
pon.aov.factor <- aov(HGI.1year$Pontaster~as.factor(HGI.years.1$Year), data=HGI.1year)
pon.aov.tukey<-TukeyHSD(pon.aov.factor)
pon.aov.tukey

#Amphipoda/Isopoda
amp.aov<-aov(HGI.1year$Amphipoda_Isopoda~as.character(HGI.years.1$Year), data=HGI.1year)
summary(amp.aov)
amp.aov.factor <- aov(HGI.1year$Amphipoda_Isopoda~as.factor(HGI.years.1$Year), data=HGI.1year)
amp.aov.tukey<-TukeyHSD(amp.aov.factor)
amp.aov.tukey

#Sabellid polychaete
sab.aov<-aov(HGI.1year$Sabellid~as.character(HGI.years.1$Year), data=HGI.1year)
summary(sab.aov)
sab.aov.factor <- aov(HGI.1year$Sabellid~as.factor(HGI.years.1$Year), data=HGI.1year)
sab.aov.tukey<-TukeyHSD(sab.aov.factor)
sab.aov.tukey

#Bathyarca freilei
fre.aov<-aov(HGI.1year$Bathyarca~as.character(HGI.years.1$Year), data=HGI.1year)
summary(fre.aov)
fre.aov.factor <- aov(HGI.1year$Bathyarca~as.factor(HGI.years.1$Year), data=HGI.1year)
fre.aov.tukey<-TukeyHSD(fre.aov.factor)
fre.aov.tukey

#White long-tentacled anemone
ane.aov<-aov(HGI.1year$Anemone~as.character(HGI.years.1$Year), data=HGI.1year)
summary(ane.aov)
ane.aov.factor <- aov(HGI.1year$Anemone~as.factor(HGI.years.1$Year), data=HGI.1year)
ane.aov.tukey<-TukeyHSD(ane.aov.factor)
ane.aov.tukey

#Skinny sabellid polychaete
ski.aov<-aov(HGI.1year$Sabellid_skinny~as.character(HGI.years.1$Year), data=HGI.1year)
summary(ski.aov)
ski.aov.factor <- aov(HGI.1year$Sabellid_skinny~as.factor(HGI.years.1$Year), data=HGI.1year)
ski.aov.tukey<-TukeyHSD(ski.aov.factor)
ski.aov.tukey

#Gastropod
gas.aov<-aov(HGI.1year$Gastropod~as.character(HGI.years.1$Year), data=HGI.1year)
summary(gas.aov)
gas.aov.factor <- aov(HGI.1year$Gastropod~as.factor(HGI.years.1$Year), data=HGI.1year)
gas.aov.tukey<-TukeyHSD(gas.aov.factor)
gas.aov.tukey

#Liparid
lip.aov<-aov(HGI.1year$Liparid~as.character(HGI.years.1$Year), data=HGI.1year)
summary(lip.aov)
lip.aov.factor <- aov(HGI.1year$Liparid~as.factor(HGI.years.1$Year), data=HGI.1year)
lip.aov.tukey<-TukeyHSD(lip.aov.factor)
lip.aov.tukey

#Suspension feeders
sus_sum <- aggregate(HGI.trophic.1$SF~HGI.years.1, FUN = mean)
colnames(sus_sum) = c("Year","mean_density") 
sus_sum_se <- aggregate(HGI.trophic.1$SF~HGI.years.1, FUN = stand_err)
sus_sum$se <- sus_sum_se$`HGI.trophic.1$SF`
sus_sum$Group<-c("Suspension feeders","Suspension feeders","Suspension feeders","Suspension feeders","Suspension feeders")
sus.bar<- ggplot(data=sus_sum, aes(x=factor(Year),y=mean_density)) +
  geom_bar_pattern(stat="identity", fill="purple3", alpha=0.7,,pattern="stripe",pattern_fill="white",pattern_color="white")  +
  geom_errorbar(aes(x=factor(Year),ymin=mean_density-se, ymax=mean_density+se, width=.2)) + 
  theme_classic()+ xlab("Year") + ylab(expression(paste(plain("Suspension feeders")," m"^plain(-2))))
sus.bar
ggsave("Bar1_SF.png", width = 4, height = 4, units=c("in"), dpi=300)

#Deposit feeders
dep_sum <- aggregate(HGI.trophic.1$DF~HGI.years.1, FUN = mean)
colnames(dep_sum) = c("Year","mean_density") 
dep_sum_se <- aggregate(HGI.trophic.1$DF~HGI.years.1, FUN = stand_err)
dep_sum$se <- dep_sum_se$`HGI.trophic.1$DF`
dep_sum$Group<-c("Deposit feeders","Deposit feeders","Deposit feeders","Deposit feeders","Deposit feeders")
dep.bar<- ggplot(data=dep_sum, aes(x=factor(Year),y=mean_density)) +
  geom_bar_pattern(stat="identity", fill="orangered1", alpha=0.7,,pattern="stripe",pattern_fill="white",pattern_color="white")  +
  geom_errorbar(aes(x=factor(Year),ymin=mean_density-se, ymax=mean_density+se, width=.2)) + 
  theme_classic()+ xlab("Year") + ylab(expression(paste(plain("Deposit feeders")," m"^plain(-2))))
dep.bar
ggsave("Bar1_DF.png", width = 4, height = 4, units=c("in"), dpi=300)

#Predator/Scavengers
pre_sum <- aggregate(HGI.trophic.1$PS~HGI.years.1, FUN = mean)
colnames(pre_sum) = c("Year","mean_density") 
pre_sum_se <- aggregate(HGI.trophic.1$PS~HGI.years.1, FUN = stand_err)
pre_sum$se <- pre_sum_se$`HGI.trophic.1$PS`
pre_sum$Group<-c("Predators/Scavengers","Predators/Scavengers","Predators/Scavengers","Predators/Scavengers","Predators/Scavengers")
pre.bar<- ggplot(data=pre_sum, aes(x=factor(Year),y=mean_density)) +
  geom_bar_pattern(stat="identity", fill="turquoise4", alpha=0.7,,pattern="stripe",pattern_fill="white",pattern_color="white")  +
  geom_errorbar(aes(x=factor(Year),ymin=mean_density-se, ymax=mean_density+se, width=.2)) + 
  theme_classic()+ xlab("Year") + ylab(expression(paste(plain("Predator/Scavengers")," m"^plain(-2))))
pre.bar
ggsave("Bar1_PS.png", width = 4, height = 4, units=c("in"), dpi=300)

#Cerianthid
cer_sum <- aggregate(HGI.1year$Cerianthid~HGI.years.1, FUN = mean)
colnames(cer_sum) = c("Year","mean_density") 
cer_sum_se <- aggregate(HGI.1year$Cerianthid~HGI.years.1, FUN = stand_err)
cer_sum$se <- cer_sum_se$`HGI.1year$Cerianthid`
cer.bar<- ggplot(data=cer_sum, aes(x=factor(Year),y=mean_density)) +
  geom_bar_pattern(stat="identity", fill="purple3", alpha=0.7,,pattern="stripe",pattern_fill="white",pattern_color="white")+
  geom_errorbar(aes(x=factor(Year),ymin=mean_density-se, ymax=mean_density+se, width=.2)) + 
  theme_classic()+ xlab("Year") + ylab(expression(paste(plain("Cerianthids")," m"^plain(-2))))
cer.bar
ggsave("Bar1_Cerianthid.png", width = 4, height = 4, units=c("in"), dpi=300)

#Jasminiera schaudinni
jas_sum <- aggregate(HGI.1year$Jasminiera~HGI.years.1, FUN = mean)
colnames(jas_sum) = c("Year","mean_density") 
jas_sum_se <- aggregate(HGI.1year$Jasminiera~HGI.years.1, FUN = stand_err)
jas_sum$se <- jas_sum_se$`HGI.1year$Jasminiera`
jas.bar<- ggplot(data=jas_sum, aes(x=factor(Year),y=mean_density)) +
  geom_bar_pattern(stat="identity", fill="purple3", alpha=0.7,,pattern="stripe",pattern_fill="white",pattern_color="white")  +
  geom_errorbar(aes(x=factor(Year),ymin=mean_density-se, ymax=mean_density+se, width=.2)) + 
  theme_classic()+ xlab("Year") + ylab(expression(paste(italic("Jasminiera schaudinni")," m"^plain(-2))))
jas.bar
ggsave("Bar1_Jasminiera.png", width = 4, height = 4, units=c("in"), dpi=300)

#Bylgides groenlandicus
byl_sum <- aggregate(HGI.1year$Bylgides~HGI.years.1, FUN = mean)
colnames(byl_sum) = c("Year","mean_density") 
byl_sum_se <- aggregate(HGI.1year$Bylgides~HGI.years.1, FUN = stand_err)
byl_sum$se <- byl_sum_se$`HGI.1year$Bylgides`
byl.bar<- ggplot(data=byl_sum, aes(x=factor(Year),y=mean_density)) +
  geom_bar_pattern(stat="identity", fill="turquoise4", alpha=0.7,,pattern="stripe",pattern_fill="white",pattern_color="white")  +
  geom_errorbar(aes(x=factor(Year),ymin=mean_density-se, ymax=mean_density+se, width=.2)) + 
  theme_classic()+ xlab("Year") + ylab(expression(paste(italic("Bylgides groenlandicus")," m"^plain(-2))))
byl.bar
ggsave("Bar1_Bylgides.png", width = 4, height = 4, units=c("in"), dpi=300)

#Bythocaris leucopis
byt_sum <- aggregate(HGI.1year$Bythocaris~HGI.years.1, FUN = mean)
colnames(byt_sum) = c("Year","mean_density") 
byt_sum_se <- aggregate(HGI.1year$Bythocaris~HGI.years.1, FUN = stand_err)
byt_sum$se <- byt_sum_se$`HGI.1year$Bythocaris`
byt.bar<- ggplot(data=byt_sum, aes(x=factor(Year),y=mean_density)) +
  geom_bar_pattern(stat="identity", fill="turquoise4", alpha=0.7,,pattern="stripe",pattern_fill="white",pattern_color="white")  +
  geom_errorbar(aes(x=factor(Year),ymin=mean_density-se, ymax=mean_density+se, width=.2)) + 
  theme_classic()+ xlab("Year") + ylab(expression(paste(italic("Bythocaris leucopis")," m"^plain(-2))))
byt.bar
ggsave("Bar1_Bythocaris.png", width = 4, height = 4, units=c("in"), dpi=300)

#Colossendeis proboscidea
col_sum <- aggregate(HGI.1year$Colossendeis~HGI.years.1, FUN = mean)
colnames(col_sum) = c("Year","mean_density") 
col_sum_se <- aggregate(HGI.1year$Colossendeis~HGI.years.1, FUN = stand_err)
col_sum$se <- col_sum_se$`HGI.1year$Colossendeis`
col.bar<- ggplot(data=col_sum, aes(x=factor(Year),y=mean_density)) +
  geom_bar_pattern(stat="identity", fill="orangered1", alpha=0.7,,pattern="stripe",pattern_fill="white",pattern_color="white")  +
  geom_errorbar(aes(x=factor(Year),ymin=mean_density-se, ymax=mean_density+se, width=.2)) + 
  theme_classic()+ xlab("Year") + ylab(expression(paste(italic("Colossendeis proboscidea")," m"^plain(-2))))
col.bar
ggsave("Bar1_Colossendeis.png", width = 4, height = 4, units=c("in"), dpi=300)

#Nymphon macronyx
nym_sum <- aggregate(HGI.1year$Nymphon~HGI.years.1, FUN = mean)
colnames(nym_sum) = c("Year","mean_density") 
nym_sum_se <- aggregate(HGI.1year$Nymphon~HGI.years.1, FUN = stand_err)
nym_sum$se <- nym_sum_se$`HGI.1year$Nymphon`
nym.bar<- ggplot(data=nym_sum, aes(x=factor(Year),y=mean_density)) +
  geom_bar_pattern(stat="identity", fill="turquoise4", alpha=0.7,pattern="stripe",pattern_fill="white",pattern_color="white")  +
  geom_errorbar(aes(x=factor(Year),ymin=mean_density-se, ymax=mean_density+se, width=.2)) + 
  theme_classic()+ xlab("Year") + ylab(expression(paste(italic("Nymphon macronyx")," m"^plain(-2))))
nym.bar
ggsave("Bar1_Nymphon.png", width = 4, height = 4, units=c("in"), dpi=300)

#Bathybiaster vexillifer  
bat_sum <- aggregate(HGI.1year$Bathybiaster~HGI.years.1, FUN = mean)
colnames(bat_sum) = c("Year","mean_density") 
bat_sum_se <- aggregate(HGI.1year$Bathybiaster~HGI.years.1, FUN = stand_err)
bat_sum$se <- bat_sum_se$`HGI.1year$Bathybiaster`
bat.bar<- ggplot(data=bat_sum, aes(x=factor(Year),y=mean_density)) +
  geom_bar_pattern(stat="identity", fill="turquoise4", alpha=0.7,,pattern="stripe",pattern_fill="white",pattern_color="white")  +
  geom_errorbar(aes(x=factor(Year),ymin=mean_density-se, ymax=mean_density+se, width=.2)) + 
  theme_classic()+ xlab("Year") + ylab(expression(paste(italic("Bathybiaster vexillifer")," m"^plain(-2))))
bat.bar
ggsave("Bar1_Bathybiaster.png", width = 4, height = 4, units=c("in"), dpi=300)

#Ophiocten gracilis
oph_sum <- aggregate(HGI.1year$Ophiocten~HGI.years.1, FUN = mean)
colnames(oph_sum) = c("Year","mean_density") 
oph_sum_se <- aggregate(HGI.1year$Ophiocten~HGI.years.1, FUN = stand_err)
oph_sum$se <- oph_sum_se$`HGI.1year$Ophiocten`
oph.bar<- ggplot(data=oph_sum, aes(x=factor(Year),y=mean_density)) +
  geom_bar_pattern(stat="identity", fill="orangered1", alpha=0.7,,pattern="stripe",pattern_fill="white",pattern_color="white")  +
  geom_errorbar(aes(x=factor(Year),ymin=mean_density-se, ymax=mean_density+se, width=.2)) + 
  theme_classic()+ xlab("Year") + ylab(expression(paste(italic("Ophiocten gracilis")," m"^plain(-2))))
oph.bar
ggsave("Bar1_Ophiocten.png", width = 4, height = 4, units=c("in"), dpi=300)

#Elpidia heckeri
elp_sum <- aggregate(HGI.1year$Elpidia~HGI.years.1, FUN = mean)
colnames(elp_sum) = c("Year","mean_density") 
elp_sum_se <- aggregate(HGI.1year$Elpidia~HGI.years.1, FUN = stand_err)
elp_sum$se <- elp_sum_se$`HGI.1year$Elpidia`
elp.bar<- ggplot(data=elp_sum, aes(x=factor(Year),y=mean_density)) +
  geom_bar_pattern(stat="identity", fill="orangered1", alpha=0.7,,pattern="stripe",pattern_fill="white",pattern_color="white")  +
  geom_errorbar(aes(x=factor(Year),ymin=mean_density-se, ymax=mean_density+se, width=.2)) + 
  theme_classic()+ xlab("Year") + ylab(expression(paste(italic("Elpidia heckeri")," m"^plain(-2))))
elp.bar
ggsave("Bar1_Elpidia.png", width = 4, height = 4, units=c("in"), dpi=300)

#Lycodes squamiventer
lyc_sum <- aggregate(HGI.1year$Lycodes~HGI.years.1, FUN = mean)
colnames(lyc_sum) = c("Year","mean_density") 
lyc_sum_se <- aggregate(HGI.1year$Lycodes~HGI.years.1, FUN = stand_err)
lyc_sum$se <- lyc_sum_se$`HGI.1year$Lycodes`
lyc.bar<- ggplot(data=lyc_sum, aes(x=factor(Year),y=mean_density)) +
  geom_bar_pattern(stat="identity", fill="turquoise4", alpha=0.7,,pattern="stripe",pattern_fill="white",pattern_color="white")  +
  geom_errorbar(aes(x=factor(Year),ymin=mean_density-se, ymax=mean_density+se, width=.2)) + 
  theme_classic()+ xlab("Year") + ylab(expression(paste(italic("Lycodes squamiventer")," m"^plain(-2))))
lyc.bar
ggsave("Bar1_Lycodes.png", width = 4, height = 4, units=c("in"), dpi=300)

#Mohnia mohni
moh_sum <- aggregate(HGI.1year$Mohnia~HGI.years.1, FUN = mean)
colnames(moh_sum) = c("Year","mean_density") 
moh_sum_se <- aggregate(HGI.1year$Mohnia~HGI.years.1, FUN = stand_err)
moh_sum$se <- moh_sum_se$`HGI.1year$Mohnia`
moh.bar<- ggplot(data=moh_sum, aes(x=factor(Year),y=mean_density)) +
  geom_bar_pattern(stat="identity", fill="turquoise4", alpha=0.7,,pattern="stripe",pattern_fill="white",pattern_color="white")  +
  geom_errorbar(aes(x=factor(Year),ymin=mean_density-se, ymax=mean_density+se, width=.2)) + 
  theme_classic()+ xlab("Year") + ylab(expression(paste(italic("Mohnia mohni")," m"^plain(-2))))
moh.bar
ggsave("Bar1_Mohnia.png", width = 4, height = 4, units=c("in"), dpi=300)

#Pontaster tennuispinus
pon_sum <- aggregate(HGI.1year$Pontaster~HGI.years.1, FUN = mean)
colnames(pon_sum) = c("Year","mean_density") 
pon_sum_se <- aggregate(HGI.1year$Pontaster~HGI.years.1, FUN = stand_err)
pon_sum$se <- pon_sum_se$`HGI.1year$Pontaster`
pon.bar<- ggplot(data=pon_sum, aes(x=factor(Year),y=mean_density)) +
  geom_bar_pattern(stat="identity", fill="orangered1", alpha=0.7,,pattern="stripe",pattern_fill="white",pattern_color="white")  +
  geom_errorbar(aes(x=factor(Year),ymin=mean_density-se, ymax=mean_density+se, width=.2)) + 
  theme_classic()+ xlab("Year") + ylab(expression(paste(italic("Pontaster tennuispinus")," m"^plain(-2))))
pon.bar
ggsave("Bar1_Pontaster.png", width = 4, height = 4, units=c("in"), dpi=300)

#Amphipoda/Isopoda
amp_sum <- aggregate(HGI.1year$Amphipoda_Isopoda~HGI.years.1, FUN = mean)
colnames(amp_sum) = c("Year","mean_density") 
amp_sum_se <- aggregate(HGI.1year$Amphipoda_Isopoda~HGI.years.1, FUN = stand_err)
amp_sum$se <- amp_sum_se$`HGI.1year$Amphipoda_Isopoda`
amp.bar<- ggplot(data=amp_sum, aes(x=factor(Year),y=mean_density)) +
  geom_bar_pattern(stat="identity", fill="grey", alpha=0.7,,pattern="stripe",pattern_fill="white",pattern_color="white")  +
  geom_errorbar(aes(x=factor(Year),ymin=mean_density-se, ymax=mean_density+se, width=.2)) + 
  theme_classic()+ xlab("Year") + ylab(expression(paste(plain("Amphipoda/Isopoda")," m"^plain(-2))))
amp.bar
ggsave("Bar1_AmphiIsopoda.png", width = 4, height = 4, units=c("in"), dpi=300)

#Sabellid polychaete
sab_sum <- aggregate(HGI.1year$Sabellid~HGI.years.1, FUN = mean)
colnames(sab_sum) = c("Year","mean_density") 
sab_sum_se <- aggregate(HGI.1year$Sabellid~HGI.years.1, FUN = stand_err)
sab_sum$se <- sab_sum_se$`HGI.1year$Sabellid`
sab.bar<- ggplot(data=sab_sum, aes(x=factor(Year),y=mean_density)) +
  geom_bar_pattern(stat="identity", fill="purple3", alpha=0.7,,pattern="stripe",pattern_fill="white",pattern_color="white")  +
  geom_errorbar(aes(x=factor(Year),ymin=mean_density-se, ymax=mean_density+se, width=.2)) + 
  theme_classic()+ xlab("Year") + ylab(expression(paste(plain("Sabellid polychaetes")," m"^plain(-2))))
sab.bar
ggsave("Bar1_Sabellid.png", width = 4, height = 4, units=c("in"), dpi=300)

#Bathyarca frielei
fre_sum <- aggregate(HGI.1year$Bathyarca~HGI.years.1, FUN = mean)
colnames(fre_sum) = c("Year","mean_density") 
fre_sum_se <- aggregate(HGI.1year$Bathyarca~HGI.years.1, FUN = stand_err)
fre_sum$se <- fre_sum_se$`HGI.1year$Bathyarca`
fre.bar<- ggplot(data=fre_sum, aes(x=factor(Year),y=mean_density)) +
  geom_bar_pattern(stat="identity", fill="purple3", alpha=0.7,,pattern="stripe",pattern_fill="white",pattern_color="white")  +
  geom_errorbar(aes(x=factor(Year),ymin=mean_density-se, ymax=mean_density+se, width=.2)) + 
  theme_classic()+ xlab("Year") + ylab(expression(paste(italic("Bathyarca frielei")," m"^plain(-2))))
fre.bar
ggsave("Bar1_Bathyarca.png", width = 4, height = 4, units=c("in"), dpi=300)

#White long-tentacled anemone
ane_sum <- aggregate(HGI.1year$Anemone~HGI.years.1, FUN = mean)
colnames(ane_sum) = c("Year","mean_density") 
ane_sum_se <- aggregate(HGI.1year$Anemone~HGI.years.1, FUN = stand_err)
ane_sum$se <- ane_sum_se$`HGI.1year$Anemone`
ane.bar<- ggplot(data=ane_sum, aes(x=factor(Year),y=mean_density)) +
  geom_bar_pattern(stat="identity", fill="purple3", alpha=0.7,,pattern="stripe",pattern_fill="white",pattern_color="white")  +
  geom_errorbar(aes(x=factor(Year),ymin=mean_density-se, ymax=mean_density+se, width=.2)) + 
  theme_classic()+ xlab("Year") + ylab(expression(paste(plain("White long-tentacled anemone")," m"^plain(-2))))
ane.bar
ggsave("Bar1_Anemone.png", width = 4, height = 4, units=c("in"), dpi=300)

#Skinny worm
ski_sum <- aggregate(HGI.1year$Sabellid_skinny~HGI.years.1, FUN = mean)
colnames(ski_sum) = c("Year","mean_density") 
ski_sum_se <- aggregate(HGI.1year$Sabellid_skinny~HGI.years.1, FUN = stand_err)
ski_sum$se <- ski_sum_se$`HGI.1year$Sabellid_skinny`
ski.bar<- ggplot(data=ski_sum, aes(x=factor(Year),y=mean_density)) +
  geom_bar_pattern(stat="identity", fill="purple3", alpha=0.7,,pattern="stripe",pattern_fill="white",pattern_color="white")  +
  geom_errorbar(aes(x=factor(Year),ymin=mean_density-se, ymax=mean_density+se, width=.2)) + 
  theme_classic()+ xlab("Year") + ylab(expression(paste(plain("Skinny worms")," m"^plain(-2))))
ski.bar
ggsave("Bar1_Skinny.png", width = 4, height = 4, units=c("in"), dpi=300)

#Gastropod
gas_sum <- aggregate(HGI.1year$Gastropod~HGI.years.1, FUN = mean)
colnames(gas_sum) = c("Year","mean_density") 
gas_sum_se <- aggregate(HGI.1year$Gastropod~HGI.years.1, FUN = stand_err)
gas_sum$se <- gas_sum_se$`HGI.1year$Gastropod`
gas.bar<- ggplot(data=gas_sum, aes(x=factor(Year),y=mean_density)) +
  geom_bar_pattern(stat="identity", fill="grey", alpha=0.7,,pattern="stripe",pattern_fill="white",pattern_color="white")  +
  geom_errorbar(aes(x=factor(Year),ymin=mean_density-se, ymax=mean_density+se, width=.2)) + 
  theme_classic()+ xlab("Year") + ylab(expression(paste(plain("Gastropods")," m"^plain(-2))))
gas.bar
ggsave("Bar1_Gastropod.png", width = 4, height = 4, units=c("in"), dpi=300)

#Liparid fish
lip_sum <- aggregate(HGI.1year$Liparid~HGI.years.1, FUN = mean)
colnames(lip_sum) = c("Year","mean_density") 
lip_sum_se <- aggregate(HGI.1year$Liparid~HGI.years.1, FUN = stand_err)
lip_sum$se <- lip_sum_se$`HGI.1year$Liparid`
lip.bar<- ggplot(data=lip_sum, aes(x=factor(Year),y=mean_density)) +
  geom_bar_pattern(stat="identity", fill="turquoise4", alpha=0.7,,pattern="stripe",pattern_fill="white",pattern_color="white")  +
  geom_errorbar(aes(x=factor(Year),ymin=mean_density-se, ymax=mean_density+se, width=.2)) + 
  theme_classic()+ xlab("Year") + ylab(expression(paste(plain("Liparid fish")," m"^plain(-2))))
lip.bar
ggsave("Bar1_Liparid.png", width = 4, height = 4, units=c("in"), dpi=300)


#######################################################################
#######################################################################
############Trajectory analysis suggested by Reviewer 2################
library(vegclust)
library(ecotraj)

#Set up data for 5 year analysis
HGI.5year.mean<-aggregate(HGI.5year,by=list(HGI.years.5),FUN=sum)
HGI.5year.mean<-HGI.5year.mean[,2:13]
HGI.5year.mean.dist<-vegdist(HGI.5year.mean,method="bray")
HGI.sites<-rep(1,5)
HGI.years.5.short<-c(2002,2007,2012,2017,2022)
Long.traj<-trajectoryPCoA(HGI.5year.mean.dist,HGI.sites,HGI.years.5.short,lwd=2)
Long.traj
PCA1<-Long.traj$points[,1]
PCA2<-Long.traj$points[,2]
Long.traj.data<-data.frame(HGI.years.5.short,PCA1,PCA2)
colnames(Long.traj.data)=c("Year","PCAX","PCAY")

#Plot with ggplot
Traj.5year<-ggplot(data=Long.traj.data)+
  geom_point(aes(x=PCAX,y=PCAY,color=as.factor(Year),shape=as.factor(Year)),size=3,stroke=1.5)+
  geom_segment(aes(x=PCAX[1],y=PCAY[1],xend=PCAX[2],yend=PCAY[2]),arrow=arrow(length=unit(0.2,"cm")))+
  geom_segment(aes(x=PCAX[2],y=PCAY[2],xend=PCAX[3],yend=PCAY[3]),arrow=arrow(length=unit(0.2,"cm")))+
  geom_segment(aes(x=PCAX[3],y=PCAY[3],xend=PCAX[4],yend=PCAY[4]),arrow=arrow(length=unit(0.2,"cm")))+
  geom_segment(aes(x=PCAX[4],y=PCAY[4],xend=PCAX[5],yend=PCAY[5]),arrow=arrow(length=unit(0.2,"cm")))+
  geom_text_repel(aes(x=PCAX,y=PCAY,label=Year),nudge_y=0.01,color="black",size=4,segment.color="transparent")+
  scale_color_manual(values=c("2002"="#332288","2007"="#117733","2012"="#AA4499","2017"="#88CCEE","2022"="#888888"))+
  scale_shape_manual(values=c("2002"=0,"2007"=1,"2012"=2,"2017"=5,"2022"=6))+
  guides(color=guide_legend("Year"),shape=guide_legend("Year"))+
  theme_classic()+ xlab("PCA1 (96%)")+ylab("PCA2 (4%)")+
  theme(legend.position="none")+
  theme(axis.text.x=element_text(size=12),axis.title.x=element_text(size=14))+
  theme(axis.text.y=element_text(size=12),axis.title.y=element_text(size=14))
Traj.5year

#Factor analysis 5 year
HGI.fa<-factanal(HGI.5year,factors=1)
HGI.fa
apply(HGI.fa$loadings^2, 1, sum)

#Set up data for 1 year analysis
HGI.1year.mean<-aggregate(HGI.1year,by=list(HGI.years.1),FUN=sum)
HGI.1year.mean<-HGI.1year.mean[,2:22]
HGI.1year.mean.dist<-vegdist(HGI.1year.mean,method="bray")
HGI.sites<-rep(1,5)
HGI.years.1.short<-c(2017,2018,2019,2021,2022)
Short.traj<-trajectoryPCoA(HGI.1year.mean.dist,HGI.sites,HGI.years.1.short,lwd=2)
Short.traj
PCA1<-Short.traj$points[,1]
PCA2<-Short.traj$points[,2]
Short.traj.data<-data.frame(HGI.years.1.short,PCA1,PCA2)
colnames(Short.traj.data)=c("Year","PCAX","PCAY")

#Plot with ggplot
Traj.1year<-ggplot(data=Short.traj.data)+
  geom_point(aes(x=PCAX,y=PCAY,color=as.factor(Year),shape=as.factor(Year)),size=3,stroke=1.5)+
  geom_segment(aes(x=PCAX[1],y=PCAY[1],xend=PCAX[2],yend=PCAY[2]),arrow=arrow(length=unit(0.2,"cm")))+
  geom_segment(aes(x=PCAX[2],y=PCAY[2],xend=PCAX[3],yend=PCAY[3]),arrow=arrow(length=unit(0.2,"cm")))+
  geom_segment(aes(x=PCAX[3],y=PCAY[3],xend=PCAX[4],yend=PCAY[4]),arrow=arrow(length=unit(0.2,"cm")))+
  geom_segment(aes(x=PCAX[4],y=PCAY[4],xend=PCAX[5],yend=PCAY[5]),arrow=arrow(length=unit(0.2,"cm")))+
  geom_text_repel(aes(x=PCAX,y=PCAY,label=Year),nudge_y=0.009,color="black",size=4,segment.color="transparent")+
  scale_color_manual(values=c("2017"="#88CCEE","2018" = "#CC6677","2019"="#DDCC77", "2021"="#44AA99","2022"="#888888"))+
  scale_shape_manual(values=c("2017"=5,"2018"=0,"2019"=1,"2021"=2,"2022"=6))+
  guides(color=guide_legend("Year"),shape=guide_legend("Year"))+
  theme_classic()+ xlab("PCA1 (72%)")+ylab("PCA2 (23%)")+
  theme(legend.position="none")+
  theme(axis.text.x=element_text(size=12),axis.title.x=element_text(size=14))+
  theme(axis.text.y=element_text(size=12),axis.title.y=element_text(size=14))
Traj.1year

#Factor analysis 1 year
HGI.fa<-factanal(HGI.1year,factors=1)
HGI.fa
apply(HGI.fa$loadings^2, 1, sum)

#Combine trajectory plots
Trajfig<-ggarrange(Traj.5year,Traj.1year,ncol=2,nrow=1,labels=c("A","B"))
Trajfig
ggsave("Fig6_Trajectories.png",width=10,height=3,units=c("in"),bg="white",dpi=300)
