library(tidyverse)
library(quantmod)
#Clearing the Global enviroment:
remove(list = ls())
#Setting the directorate:
setwd("C:/Users/cvazque/Desktop/demos/")
tab<-read.csv("QiangP.csv", sep=",")%>%select(sample, chr)

tab<-separate(tab, chr, into=c("chro", "begin"), sep="[\\:]")%>%
  separate(begin, into=c("end", "begin"), sep = "[\\-]")%>%
  separate(begin, into=c("begin", "X"), sep = "[\\s]")%>%
  select(sample, begin, end)

#Switch coordinates##hg38
Sm<-c(105856987,105860436)
Sg3<-c(105772941,105774642)
Sg1<-c(105744600,105746938)
Sg2<-c(105645612,105647422)
Sg4<-c(105627601,105628640)
Sa1<-c(105709267,105712271)
Sa2<-c(105589201,105591275)
Se<- c(105602375,105603285)
#Switch coordinates entre puntos##hg38
SSm<-c(105856987:105860436)
SSg3<-c(105772941:105774642)
SSg1<-c(105744600:105746938)
SSg2<-c(105645612:105647422)
SSg4<-c(105627601:105628640)
SSa1<-c(105709267:105712271)
SSa2<-c(105589201:105591275)
SSe<- c(105602375:105603285)
#Create combined table:
tab6<-filter(tab, sample=="ATM106")%>%select(begin, end)%>%gather()%>%mutate("sample"=c("ATM106"))%>%group_by(sample, key, value)%>%summarise(n=n())
tab7<-filter(tab, sample=="ATM107")%>%select(begin, end)%>%gather()%>%mutate("sample"=c("ATM107"))%>%group_by(sample, key, value)%>%summarise(n=n())
tabc1<-filter(tab, sample=="Control1")%>%select(begin, end)%>%gather()%>%mutate("sample"=c("Control1"))%>%group_by(sample, key, value)%>%summarise(n=n())
tabL<-filter(tab, sample=="Lig4")%>%select(begin, end)%>%gather()%>%mutate("sample"=c("Lig4"))%>%group_by(sample, key, value)%>%summarise(n=n())
tabc2<-filter(tab, sample=="Control2")%>%select(begin, end)%>%gather()%>%mutate("sample"=c("Control2"))%>%group_by(sample, key, value)%>%summarise(n=n())
X<-rbind( tab6, tab7, tabc1, tabc2)
X.<-mutate(X, "type"= ifelse(grepl("Control", sample), "Control", "Case"))
#Breakpoints in the IgH locus:
ggplot(X.)+
  geom_rect(aes(xmin=Sm[1], xmax=Sm[2], ymin=0, ymax=74827), fill="grey95",alpha=0.3)+
  geom_rect(aes(xmin=Sg3[1], xmax=Sg3[2], ymin=0, ymax=74827),fill="grey95", alpha=0.3)+
  geom_rect(aes(xmin=Sg1[1], xmax=Sg1[2], ymin=0, ymax=74827), fill="grey95",alpha=0.3)+
  geom_rect(aes(xmin=Sg2[1], xmax=Sg2[2], ymin=0, ymax=74827),fill="grey95", alpha=0.3)+
  geom_rect(aes(xmin=Sg4[1], xmax=Sg4[2], ymin=0, ymax=74827), fill="grey95",alpha=0.3)+
  geom_rect(aes(xmin=Sa1[1], xmax=Sa1[2], ymin=0, ymax=74827),fill="grey95", alpha=0.3)+
  geom_rect(aes(xmin=Sa2[1], xmax=Sa2[2], ymin=0, ymax=74827), fill="grey95",alpha=0.3)+
  geom_histogram(aes(x=as.numeric(value), y=n, color=sample), stat="identity", alpha=0.5)+
  scale_color_manual(values=c("BRCA1"="coral", 
                              "AID"="darkgreen", 
                              "ATM106"="black",
                              "ATM107"="blue3",
                              "NIPBL"="orange",
                              "Control1"="purple",
                              "Control2"="red",
                              "Control3"="darkred"))+
  ggtitle("Coordinates of Breakpoints")+
  scale_y_log10()+
  geom_hline(yintercept = 0, color="black")+
  facet_grid(sample~.)+
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
        axis.title.y=element_blank(), axis.ticks.x = element_blank(),
        panel.background = element_rect(fill = "white"), plot.title = element_text(hjust=0.5, size=15),
        axis.text.y=element_text(size=5), panel.border = element_rect(fill=NA, color="grey75"), 
        legend.position = "none", strip.text.y = element_text(size=8, vjust=0.5))
#Plot only the breaks outside of the switch regions:
  #Remove the breakpoints occurring in the switch regions and create a new dataframe:
    #One vector with coordinates from the switch regions:
vector<- c(SSm, SSg3, SSg1, SSa1, SSe, SSa2, SSg4, SSg2)

df<-subset(X., !(value %in% vector))

    #Check in the plot if it is how it should be:
ggplot(df)+
  geom_rect(aes(xmin=Sm[1], xmax=Sm[2], ymin=0, ymax=400), fill="grey95",alpha=0.3)+
  geom_rect(aes(xmin=Sg3[1], xmax=Sg3[2], ymin=0, ymax=400),fill="grey95", alpha=0.3)+
  geom_rect(aes(xmin=Sg1[1], xmax=Sg1[2], ymin=0, ymax=400), fill="grey95",alpha=0.3)+
  geom_rect(aes(xmin=Sg2[1], xmax=Sg2[2], ymin=0, ymax=400),fill="grey95", alpha=0.3)+
  geom_rect(aes(xmin=Sg4[1], xmax=Sg4[2], ymin=0, ymax=400), fill="grey95",alpha=0.3)+
  geom_rect(aes(xmin=Sa1[1], xmax=Sa1[2], ymin=0, ymax=400),fill="grey95", alpha=0.3)+
  geom_rect(aes(xmin=Sa2[1], xmax=Sa2[2], ymin=0, ymax=400), fill="grey95",alpha=0.3)+
  geom_histogram(aes(x=as.numeric(value), y=n, color=sample), stat="identity", alpha=0.5)+
  scale_color_manual(values=c("BRCA1"="coral", 
                              "AID"="darkgreen", 
                              "ATM106"="black",
                              "ATM107"="blue3",
                              "NIPBL"="orange",
                              "Control1"="purple",
                              "Control2"="red",
                              "Control3"="darkred"))+
  ggtitle("Coordinates of Breakpoints")+
  scale_y_log10()+
  geom_hline(yintercept = 0, color="black")+
  facet_grid(sample~.)+
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
        axis.title.y=element_blank(), axis.ticks.x = element_blank(),
        panel.background = element_rect(fill = "white"), plot.title = element_text(hjust=0.5, size=15),
        axis.text.y=element_text(size=5), panel.border = element_rect(fill=NA, color="grey75"), 
        legend.position = "none", strip.text.y = element_text(size=8, vjust=0.5))
  #Plot the IgH locus in four parts:
    #92345 nt
    #Sm-Sg3
filter(df, value %in% c((Sm[1]+10000):Sg3[2]))%>%
ggplot()+
  geom_rect(aes(xmin=Sm[1], xmax=Sm[2], ymin=0, ymax=400), fill="grey95",alpha=0.3)+
  geom_rect(aes(xmin=Sg3[1], xmax=Sg3[2], ymin=0, ymax=400),fill="grey95", alpha=0.3)+
  #primer binding IgG3
  geom_rect(aes(xmin=105772166, xmax=105772194, ymin=0, ymax=3855),fill="lightpink", alpha=0.3)+
  #I3 promoter IgG3
  geom_rect(aes(xmin=105774999, xmax=105775994, ymin=0, ymax=3855),fill="cornsilk", alpha=0.3)+
  geom_histogram(aes(x=as.numeric(value), y=n, color=sample), stat="identity", alpha=0.5)+
  scale_color_manual(values=c("BRCA1"="coral", 
                              "AID"="darkgreen", 
                              "ATM106"="black",
                              "ATM107"="blue3",
                              "NIPBL"="orange",
                              "Control1"="purple",
                              "Control2"="red",
                              "Control3"="darkred"))+
  ggtitle("Coordinates of Breakpoints Sm to Sg3 - 92345 nt")+
  scale_y_log10()+
  geom_hline(yintercept = 0, color="black")+
  facet_grid(sample~.)+
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
        axis.title.y=element_blank(), axis.ticks.x = element_blank(),
        panel.background = element_rect(fill = "white"), plot.title = element_text(hjust=0.5, size=15),
        axis.text.y=element_text(size=5), panel.border = element_rect(fill=NA, color="grey75"), 
        legend.position = "none", strip.text.y = element_text(size=8, vjust=0.5))
    #Sg3-Sa1#
filter(df, value %in% c((Sg3[1]+5000):Sa1[2]))%>%
  ggplot()+
  geom_rect(aes(xmin=Sa1[1], xmax=Sa1[2], ymin=0, ymax=400), fill="grey95",alpha=0.3)+
  geom_rect(aes(xmin=Sg1[1], xmax=Sg1[2], ymin=0, ymax=400), fill="grey95",alpha=0.3)+
  geom_rect(aes(xmin=Sg3[1], xmax=Sg3[2], ymin=0, ymax=400),fill="grey95", alpha=0.3)+
  #primer binding IgG3
  geom_rect(aes(xmin=105772166, xmax=105772194, ymin=0, ymax=400),fill="lightpink", alpha=0.3)+#primer binding
  geom_rect(aes(xmin=105743825, xmax=105743853, ymin=0, ymax=400),fill="lightpink", alpha=0.3)+
  #I3 promoter IgG1:
  geom_rect(aes(xmin=105747472, xmax=105747984, ymin=0, ymax=400),fill="cornsilk", alpha=0.3)+
  #I3 promoter IgG3
  geom_rect(aes(xmin=105774999, xmax=105775994, ymin=0, ymax=400),fill="cornsilk", alpha=0.3)+
  #I3 promoter IgA1
  geom_rect(aes(xmin=105712828, xmax=105713605, ymin=0, ymax=400),fill="cornsilk", alpha=0.3)+
  geom_histogram(aes(x=as.numeric(value), y=n, color=sample), stat="identity", alpha=0.5)+
  scale_color_manual(values=c("BRCA1"="coral", 
                              "AID"="darkgreen", 
                              "ATM106"="black",
                              "ATM107"="blue3",
                              "NIPBL"="orange",
                              "Control1"="purple",
                              "Control2"="red",
                              "Control3"="darkred"))+
  ggtitle("Coordinates of Breakpoints Sg3 to Sa1 - 65170 nt")+
  scale_y_log10()+
  geom_hline(yintercept = 0, color="black")+
  facet_grid(sample~.)+
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
        axis.title.y=element_blank(), axis.ticks.x = element_blank(),
        panel.background = element_rect(fill = "white"), plot.title = element_text(hjust=0.5, size=15),
        axis.text.y=element_text(size=5), panel.border = element_rect(fill=NA, color="grey75"), 
        legend.position = "none", strip.text.y = element_text(size=8, vjust=0.5))
    #Sa1-Sg2#
filter(df, value %in% c((Sa1[1]+6500):Sg2[2]))%>%
  ggplot()+
  geom_rect(aes(xmin=Sa1[1], xmax=Sa1[2], ymin=0, ymax=400), fill="grey95",alpha=0.3)+
  geom_rect(aes(xmin=Sg2[1], xmax=Sg2[2], ymin=0, ymax=400), fill="grey95",alpha=0.3)+
  #I3 promoter IgA1
  geom_rect(aes(xmin=105712828, xmax=105713605, ymin=0, ymax=400),fill="cornsilk", alpha=0.3)+
  geom_histogram(aes(x=as.numeric(value), y=n, color=sample), stat="identity", alpha=0.5)+
  scale_color_manual(values=c("BRCA1"="coral", 
                              "AID"="darkgreen", 
                              "ATM106"="black",
                              "ATM107"="blue3",
                              "NIPBL"="orange",
                              "Control1"="purple",
                              "Control2"="red",
                              "Control3"="darkred"))+
  ggtitle("Coordinates of Breakpoints Sa1 to Sg2 - 68345 nt")+
  scale_y_log10()+
  geom_hline(yintercept = 0, color="black")+
  facet_grid(sample~.)+
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
        axis.title.y=element_blank(), axis.ticks.x = element_blank(),
        panel.background = element_rect(fill = "white"), plot.title = element_text(hjust=0.5, size=15),
        axis.text.y=element_text(size=5), panel.border = element_rect(fill=NA, color="grey75"), 
        legend.position = "none", strip.text.y = element_text(size=8, vjust=0.5))
    #Sg2-Sa2#
filter(df, value %in% c((Sg2[1]+6500):Sa2[2]))%>%
  ggplot()+
  geom_rect(aes(xmin=Sa2[1], xmax=Sa2[2], ymin=0, ymax=400), fill="grey95",alpha=0.3)+
  geom_rect(aes(xmin=Sg2[1], xmax=Sg2[2], ymin=0, ymax=400), fill="grey95",alpha=0.3)+
  geom_rect(aes(xmin=Sg4[1], xmax=Sg4[2], ymin=0, ymax=400), fill="grey95",alpha=0.3)+
  #I3 promoter
  geom_rect(aes(xmin=105591857, xmax=105592634, ymin=0, ymax=19144),fill="cornsilk", alpha=0.3)+
  #exon I
  geom_rect(aes(xmin=105628972, xmax=105629967, ymin=0, ymax=290),fill="cornsilk", alpha=0.3)+
  geom_histogram(aes(x=as.numeric(value), y=n, color=sample), stat="identity", alpha=0.5)+
  scale_color_manual(values=c("BRCA1"="coral", 
                              "AID"="darkgreen", 
                              "ATM106"="black",
                              "ATM107"="blue3",
                              "NIPBL"="orange",
                              "Control1"="purple",
                              "Control2"="red",
                              "Control3"="darkred"))+
  ggtitle("Coordinates of Breakpoints Sg2 to Sa2 - 60837 nt")+
  scale_y_log10()+
  geom_hline(yintercept = 0, color="black")+
  facet_grid(sample~.)+
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
        axis.title.y=element_blank(), axis.ticks.x = element_blank(),
        panel.background = element_rect(fill = "white"), plot.title = element_text(hjust=0.5, size=15),
        axis.text.y=element_text(size=5), panel.border = element_rect(fill=NA, color="grey75"), 
        legend.position = "none", strip.text.y = element_text(size=8, vjust=0.5))

##Plot the deletions##
setwd("C:/Users/cvazque/Desktop/demos/")
tab<-read.csv("deletions.csv", sep=",")

dfAT1<-filter(tab, grepl("06", donor))
dfAT2<-filter(tab, grepl("07", donor))
dfC1<-filter(tab, grepl("l 1", donor))
dfC2<-filter(tab, grepl("l 2", donor))

  #IF NECESSARY# Filter the read for a specific region in the genome#
  #Take into account that many columns have coordinates#
  #Vector from last coord of Sm to beginning of Sg3#
RPIG3<-c(105791037:105775994)
RPIG1<-c(105747984:105762941)
RPIA1<-c(105730148:105713605)
RPIA2<-c(105599634:105592634)
RPIG4<-c(105640972:105629967)

filter.coord<-function(df, vector){
  gg<-data.frame()
  for(i in (3:ncol(df))){
    ty<-data.frame()
    ty<-filter(df, str_detect(df[,i], paste(vector, collapse="|")))
    
    if(dim(ty)[1]>0){
      gg<-rbind(gg, ty)
    }
  
  }
  return(gg)
}

dfA1<-filter.coord(dfAT1, RPIA2)%>%unique()
dfA2<-filter.coord(dfAT2, RPIA2)%>%unique()
dfC1<-filter.coord(dfC1, RPIA2)%>%unique()
dfC2<-filter.coord(dfC2, RPIA2)%>%unique()

dfAT1<-filter(dfAT1, str_detect(dfAT1[,4], paste(RPIG1, collapse="|")))
dfAT2<-filter(dfAT2, str_detect(dfAT2[,4], paste(RPIG1, collapse="|")))


  #Randomly taking reads from the samples:
dfAT1<-dfA1[(sample(1:nrow(dfA1), 10)), ]
dfAT1<-mutate(dfAT1,"numb"=c(1:(nrow(dfAT1))), "ymin"=5*(numb-1)+2)
dfAT2<-dfA2[(sample(1:nrow(dfA2), 10)), ]
dfAT2<-mutate(dfAT2,"numb"=c(1:(nrow(dfAT2))), "ymin"=5*(numb-1)+2)
dfC1<-dfC1[(sample(1:nrow(dfC1), 10)), ]
dfC1<-mutate(dfC1,"numb"=c(1:(nrow(dfC1))), "ymin"=5*(numb-1)+2)
dfC2<-dfC2[(sample(1:nrow(dfC2), 10)), ]
dfC2<-mutate(dfC2,"numb"=c(1:(nrow(dfC2))), "ymin"=5*(numb-1)+2)

df500<- rbind(dfAT1, dfAT2, dfC1, dfC2)

demo5<-gather(df500,"type", "coor", -donor, -read.id, -numb, -ymin)%>%
  filter(coor != "")%>%select(-type)
demo5<-separate(demo5, coor, into=c("trash", "coor"), sep = "[\\:]")%>%
  separate(coor, into=c("coor1", "coor2"), sep="[\\-]")%>%
  separate(coor2, into=c("coor2", "trash"), sep = "[\\(]")%>%select(-trash)
demo5<-arrange(demo5, donor, numb)
bin5<-mutate(demo5, "deletion"=ifelse(demo5$numb == lead(demo5$numb, 1), lead(demo5$coor1),"NO"), "ymax"=ymin+2, "size"=as.numeric(deletion)-as.numeric(coor2))%>%
  mutate("case"=ifelse(grepl("Control", donor), "Control", "Case"))

    ##Plotting:
#filter(bin5, deletion != "NO", size > 0) %>%
  #filter(coor1 %in% c((Sm[1]+10000):(Sg3[2]-1000)))%>%
  #filter(deletion %in% c((Sm[1]+10000):(Sg3[2]-1000)))%>%
  ggplot(bin5)+
  #geom_rect(aes(xmin=as.numeric(coor2), xmax=as.numeric(deletion), ymin=ymin, ymax=ymax, color=donor),fill="black", alpha=0.5)+
  geom_rect(aes(xmin=Sg3[1], xmax=Sg3[2], ymin=0, ymax=50),fill="grey95", alpha=0.3)+
  geom_rect(aes(xmin=Sg1[1], xmax=Sg1[2], ymin=0, ymax=50), fill="grey95",alpha=0.3)+
  geom_rect(aes(xmin=Sg2[1], xmax=Sg2[2], ymin=0, ymax=50),fill="grey95", alpha=0.3)+
  geom_rect(aes(xmin=Sg4[1], xmax=Sg4[2], ymin=0, ymax=50), fill="grey95",alpha=0.3)+
  geom_rect(aes(xmin=Sa1[1], xmax=Sa1[2], ymin=0, ymax=50),fill="grey95", alpha=0.3)+
  geom_rect(aes(xmin=Sa2[1], xmax=Sa2[2], ymin=0, ymax=50),fill="grey95",alpha=0.3)+
  geom_rect(aes(xmin=Sm[1], xmax=Sm[2], ymin=0, ymax=50), fill="grey95",alpha=0.3)+
  #geom_rect(aes(xmin=106097612, xmax=106303092, ymin=0, ymax=50), fill="red",alpha=0.3)+
  #I3 promoter IgA2
  geom_rect(aes(xmin=105591857, xmax=105592634, ymin=0, ymax=50),fill="cornsilk", alpha=0.3)+
  #I3 promoter IgA1
  geom_rect(aes(xmin=105712828, xmax=105713605, ymin=0, ymax=50),fill="cornsilk", alpha=0.3)+
  #I3 promoter IgG1:
  geom_rect(aes(xmin=105747472, xmax=105747984, ymin=0, ymax=50),fill="cornsilk", alpha=0.3)+
  #I3 promoter IgG3
  geom_rect(aes(xmin=105774999, xmax=105775994, ymin=0, ymax=50),fill="cornsilk", alpha=0.3)+
    #exon I Ig4
  geom_rect(aes(xmin=105628972, xmax=105629967, ymin=0, ymax=50),fill="cornsilk", alpha=0.3)+
  ##
  ##geom_rect(aes(xmin=105640972, xmax=105629967, ymin=0, ymax=30),fill="darkgreen", alpha=0.3)+
  ##
  #geom_rect(aes(xmin=as.numeric(coor2), xmax=as.numeric(deletion), ymin=ymin, ymax=ymax, color=donor), fill="orange", alpha=0.5)+
  geom_rect(aes(xmin=as.numeric(coor1), xmax=as.numeric(coor2), ymin=ymin, ymax=ymax), color="grey65")+
  geom_histogram(aes(x=as.numeric(coor2)), bins=900, fill="red", alpha=0.5)+
  geom_histogram(aes(x=as.numeric(coor1)), bins=900, fill="red", alpha=0.5)+
  scale_color_manual(values=c("BRCA1"="coral", 
                              "Lig4"="black",
                              "AID"="orange", 
                              "ATM106"="turquoise4",
                              "ATM107"="blue3",
                              "NIPBL"="darkgreen",
                              "Control 1"="purple",
                              "Control 2"="red",
                              "Control 3"="darkred"))+
  ggtitle("Reads from 10 samples IgH locus")+
  geom_hline(yintercept = 0, color="black")+
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
        axis.title.y=element_blank(), axis.ticks.x = element_blank(),
        panel.background = element_rect(fill = "white"), plot.title = element_text(hjust=0.5, size=15),
        axis.text.y=element_text(size=5), panel.border = element_rect(fill=NA, color="grey75"), 
        strip.text.y = element_text(size=7.5,face="bold", vjust=0.5), legend.position="none")+
  facet_grid(donor~.)

#Estimate#
  dfAT1<-filter(tab, grepl("06", donor))
  dfAT2<-filter(tab, grepl("07", donor))
  dfC1<-filter(tab, grepl("l 1", donor))
  dfC2<-filter(tab, grepl("l 2", donor))
  #Get the nrows = nreads of these dataframes:
  dfA1<-filter.coord(dfAT1, RPIG1)%>%unique()
  dfA2<-filter.coord(dfAT2, RPIG1)%>%unique()
  dfC1<-filter.coord(dfC1, RPIG1)%>%unique()
  dfC2<-filter.coord(dfC2, RPIG1)%>%unique()
  
  nrow(dfA2)
  hh<-c()
  hh<-c(nrow(dfA1))
  hh<-c(hh,nrow(dfC2))
  
  #Get reads from here that have Sm:
  dffA2<- filter(dfC2, str_detect(dfC2[,4], paste(SSm, collapse="|")))
  ll<-c(0)
  ll<-c(ll, nrow(dffA2))
  
  #Make an excel and bring it up:
  setwd("C:/Users/cvazque/Desktop/demos/")
  df<-read.csv("reads.csv")
  ggplot(df)+
    geom_histogram(aes(x=donor, y=value, fill=class), alpha=0.5, stat="identity", position="dodge", color="black")+
    ylim(c(0,1000))+
    scale_y_log10()
    
  
  
####Getting the reads from FASTA####
  ##Create a vector with the names of the reads:
vector<-group_by(bin5, read.id, donor,numb)%>%
    summarise(n=n())%>%
  arrange(donor)

  ##Working with FASTA:
library(Biostrings)
library(seqinr)
#BiocManager::install("snpStats")
library(ShortRead)  
library(snpStats)
library(GenomicTools.fileHandler)
  #Set the directory where all the FASTA are:
setwd("S:/Clara VG/1 - Experiments/1 - MinION/1 - MinION runs/20190506_Qiang Pan/Results/fastq_pass/")
  #list.files() have the names of all files in the directory
a<-list.files()

    #Seeing if the loop will work##
#out<-c()
#fastafile<- readFastq(a[14])
#b<-attr(fastafile[1:5], "id")%>%as.data.frame()
#c<-attr(fastafile[1:5], "sread")%>%as.data.frame()
#d<-cbind(b,c)
#colnames(d)[1]<-"id"
#  gg<-filter(d, grepl(paste(vector$read.id, collapse="|"), id))%>%mutate("s"=paste(id,x))%>%select(s)
#  if(dim(gg)[1]> 0){
#    out<-c(out, gg)
#  }
#gg<-filter(d, grepl(vector$read.id[x],id))%>%mutate("s"=paste(id,x))%>%select(s)
#out<-c(out, gg)

    ####Loop#####
    ##Create an empty vector where the sequences matching with the reads will end up##
out<-c()
    ##You will read the fastq files from the folder of preference and will create
    ##dataframes that includes the read and the read-ID
    ##The you filter the reads each loop with the vector "vector" that contains the info about your reads
    ##If in a cycle the loop finds at least one match it will be included into the vector "out"
for(i in (1:length(a))){
  fastafile<- readFastq(a[i])
  b<-attr(fastafile, "id")%>%as.data.frame()
  c<-attr(fastafile, "sread")%>%as.data.frame()
  d<-cbind(b,c)
  colnames(d)[1]<-"id"
  gg<-filter(d, grepl(paste(vector$read.id, collapse="|"),id))%>%mutate("s"=paste(id,x))%>%select(s)
  
  if(dim(gg)[1]>0){
  out<-c(out, gg)
  }
}



#Get the sample to the right sequence#
setwd("S:/Clara VG/1 - Experiments/0 - Experiments 2020/CVG029_200526_QiangPanInserts/ATM breakpoints between switch regions/plots only 10/run 18_Sa2 control w data/")
lapply(out, write, "test4.txt", append=TRUE)
df<-read.csv("reads10x4.csv", header=FALSE)%>%separate(V1, into=c("ID", "sequence"), sep="QPrun1") 

select(vector, -n)%>%write.table( "filesamples.txt")
write.table(df, "file.txt")
