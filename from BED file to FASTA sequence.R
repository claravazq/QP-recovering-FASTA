library(tidyverse)
library(quantmod)
#Clearing the Global enviroment:
remove(list = ls())
##Setting the directorate and getting the file:
setwd("C:/Users/cvazque/Desktop/demos/")
tab<-read.csv("deletions.csv", sep=",")

##Filtering the df based by sample
dfAT1<-filter(tab, grepl("06", donor))
dfAT2<-filter(tab, grepl("07", donor))
dfC1<-filter(tab, grepl("l 1", donor))
dfC2<-filter(tab, grepl("l 2", donor))

##Function that allows to filter the reads in the df that contain certain coordinates given by a vector
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

##IF NECESSARY# Filter the read for a specific region in the genome#
  #Take into account that many columns have coordinates#
    #Vector containing desired coordinates#
RPIG3<-c(105791037:105775994)

##Usage of the function in the independent df
dfA1<-filter.coord(dfAT1, RPIA2)%>%unique()
dfA2<-filter.coord(dfAT2, RPIA2)%>%unique()
dfC1<-filter.coord(dfC1, RPIA2)%>%unique()
dfC2<-filter.coord(dfC2, RPIA2)%>%unique()

##Randomly take reads from the df:
dfAT1<-dfA1[(sample(1:nrow(dfA1), 10)), ]
dfAT1<-mutate(dfAT1,"numb"=c(1:(nrow(dfAT1))), "ymin"=5*(numb-1)+2)
dfAT2<-dfA2[(sample(1:nrow(dfA2), 10)), ]
dfAT2<-mutate(dfAT2,"numb"=c(1:(nrow(dfAT2))), "ymin"=5*(numb-1)+2)
dfC1<-dfC1[(sample(1:nrow(dfC1), 10)), ]
dfC1<-mutate(dfC1,"numb"=c(1:(nrow(dfC1))), "ymin"=5*(numb-1)+2)
dfC2<-dfC2[(sample(1:nrow(dfC2), 10)), ]
dfC2<-mutate(dfC2,"numb"=c(1:(nrow(dfC2))), "ymin"=5*(numb-1)+2)
  #Bind the independent df
df500<- rbind(dfAT1, dfAT2, dfC1, dfC2)

##Organize the df to work with it - 
  #separate the coordinates into two of the, and addition of a number to follow back the read:
demo5<-gather(df500,"type", "coor", -donor, -read.id, -numb, -ymin)%>%
  filter(coor != "")%>%select(-type)
demo5<-separate(demo5, coor, into=c("trash", "coor"), sep = "[\\:]")%>%
  separate(coor, into=c("coor1", "coor2"), sep="[\\-]")%>%
  separate(coor2, into=c("coor2", "trash"), sep = "[\\(]")%>%select(-trash)
demo5<-arrange(demo5, donor, numb)
bin5<-mutate(demo5, "deletion"=ifelse(demo5$numb == lead(demo5$numb, 1), lead(demo5$coor1),"NO"), "ymax"=ymin+2, "size"=as.numeric(deletion)-as.numeric(coor2))%>%
  mutate("case"=ifelse(grepl("Control", donor), "Control", "Case"))

##Plotting the 10 different samples form each donor:
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

##Estimate the number of valid reads in the total number of reads
#valid reads = contain Sm in the first coordinate#
#hh is a vector that contain all the information of all the independent df
  nrow(dfA2)
  hh<-c()
  hh<-c(nrow(dfA1))
  hh<-c(hh,nrow(dfC2))
  
  #Get reads from here that have Sm and store them in ll:
  dffA2<- filter(dfC2, str_detect(dfC2[,4], paste(SSm, collapse="|")))
  ll<-c(0)
  ll<-c(ll, nrow(dffA2))
  
  #Make an excel outside R and import it back:
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

#Save the information in the directory:
select(vector, -n)%>%write.table( "filesamples.txt")
write.table(df, "file.txt")
