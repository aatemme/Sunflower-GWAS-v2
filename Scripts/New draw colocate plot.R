


colocate<-rbind(sig.blocks,sug.blocks[sug.blocks$hapID%in%sig.blocks$hapID,])

colocate$sighap<-sighap_to_genomehap$sig.hap[match(colocate$hapID,sighap_to_genomehap$genome.hap)]


colocate$trait_env<-paste(colocate$trait,colocate$env,sep="_")

traits.per.block<-colocate %>% group_by(sighap) %>% summarise(trait_num=length(trait_env))

# single.trait.blocks<-colocate[colocate$hapID%in%traits.per.block$hapID[traits.per.block$trait_num==1],]
# 
# colocate<-colocate[!colocate$hapID%in%traits.per.block$hapID[traits.per.block$trait_num==1],]

colocate<-colocate %>% separate(sighap, sep= "_", c("chromosome","blocknum"),remove=F) %>%
                        arrange(chromosome, blocknum)

colocate<- colocate %>% group_by(chromosome) %>%
                        mutate(blockseq=match(blocknum,unique(blocknum)),beta.sign=sign(beta)) %>%
                        mutate(region=paste(formatC(as.numeric(chromosome),width=2, flag="0"), 
                                                               formatC(blockseq,width=2, flag="0"),sep="-"))

colocate$region<-factor(colocate$region)

chrom.borders<-colocate %>% group_by(chromosome)%>% summarise(bin=length(unique(region))) %>% arrange(as.integer(chromosome))
chrom.borders<-cumsum(chrom.borders$bin)
chrom.borders<-chrom.borders+0.5
chrom.borders<-chrom.borders[1:length(chrom.borders)-1]

baseplot<-ggplot(colocate,aes(x=region,y=trait,fill=as.factor(beta.sign)))

baseplot+geom_vline(xintercept=c(1:length(plot.data$region)),colour="darkgrey",linetype=3)+
  geom_vline(xintercept=chrom.borders,colour="black")+
  geom_tile(fill="white")+
  geom_tile(aes(alpha=pvalue),colour="black")+
  theme_minimal()+
  theme(axis.text.y = element_text(hjust = 0))+
  scale_fill_manual(values=c(colors[1],colors[2]))+
  scale_alpha_manual(values=c(1,0.1))+
  scale_x_discrete(drop=F)+
  theme_classic()+
  theme(axis.title.y=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5,hjust=1))+
  ggtitle("control")+theme(legend.position = "none")+theme(axis.title.x=element_blank())+
  facet_wrap(~env,nrow=3)






plot.data<-colocate[colocate$env=="water",]

# chrom.borders<-plot.data %>% group_by(chromosome)%>% summarise(bin=length(unique(region))) %>% arrange(as.integer(chromosome))
# chrom.borders<-cumsum(chrom.borders$bin)
# chrom.borders<-chrom.borders+0.5
# chrom.borders<-chrom.borders[1:length(chrom.borders)-1]

baseplot<-ggplot(plot.data,aes(x=region,y=trait,fill=as.factor(beta.sign)))

baseplot+geom_vline(xintercept=c(1:length(plot.data$region)),colour="darkgrey",linetype=3)+
  geom_vline(xintercept=chrom.borders,colour="black")+
  geom_tile(fill="white")+
  geom_tile(aes(alpha=pvalue),colour="black")+
  theme_minimal()+
  theme(axis.text.y = element_text(hjust = 0))+
  scale_fill_manual(values=c(colors[1],colors[2]))+
  scale_alpha_manual(values=c(1,0.1))+
  scale_x_discrete(drop=F)+
  theme_classic()+
  theme(axis.title.y=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5,hjust=1))+
  ggtitle("control")+theme(legend.position = "none")+theme(axis.title.x=element_blank())



