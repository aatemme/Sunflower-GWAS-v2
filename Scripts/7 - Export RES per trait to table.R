### script to make table of RES of traits
library(tidyverse)
library(data.table)

colocate<-read.table("Tables/Blocks/colocate_table.txt")

region.RES<-colocate %>% filter(pvalue=="significant") %>%
                            group_by(trait_env,region) %>% summarise(RES=max(RES))

trait.RES<-region.RES %>% group_by(trait_env) %>% summarise(RES=sum(RES), count=length(region)) %>% separate(trait_env,into=c("trait","env"),sep="_")

RES<-trait.RES %>% gather("metric","value", -c(trait,env)) %>% unite("spreader",env,metric,sep="_") %>% spread(spreader,value)

write.table(RES,"Tables/Blocks/RES.txt")
