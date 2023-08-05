
#The Effects of Cover Crops on Insect Biodiversity- A Meta-Analysis
---
title: "Project Cover Crop"
author: "Eugene Ohba"
date: "2023-08-05"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```


#Loads Package
```{r}
library(tidyverse)
library(sp)
library(sf)
library(metafor)
library(bitops)
library(RCurl)
library(Formula)
```

#Loads Data
```{r}
data=read.csv(file="Cover crop arthropod database - Data.csv") %>%  
  mutate(cover.mono=factor(cover.mono, levels=c("y", "n")),
         
         metric.simplified= factor(metric.simplified, levels=c("abundance", "diversity")), 
         
         arthropod.function= factor(arthropod.function, levels=c("predators", "pests")))

```

#Map of Studies
```{r}
map_data= data %>% 
  group_by(doi, location) %>% 
  summarize(lat=mean(lat), lon=mean(lon))

map_coordinates=map_data("world")
map_plot= ggplot() +
  geom_map(data=map_coordinates, map=map_coordinates, 
           aes(long, lat, map_id=region), fill="lightgreen")+ 
  geom_point(data=map_data, aes(x=lon, y=lat))+
  ylim(-55,100)+
  theme_classic()

map_plot
#Our dataset may be biased
```
![image](https://github.com/Eugenefut19/Meta-Analysis-on-Cover-Crops/assets/134546229/759f9121-ca48-4216-8709-3f8ac3a2afe6)


#Effect size
```{r}

data.metafor= escalc(data=data, 
                     n1i= n.treatment, 
                     n2i= n.control, 
                     m1i= treatment.mean, 
                     m2i = control.mean, 
                     sd1i= treatment.sd, 
                     sd2i= control.sd, 
                     measure= "ROM", 
                     append=T)
```

#Counts numebr of studies with mono
```{r}
c_data=data %>% 
  filter(cover.mono %in% c("y", "n")) %>%
  group_by(cover.mono) %>% 
  ggplot(aes(x=cover.mono)) +
  labs(y= "Count", x = "Studies with Monoculture")+
  ggtitle("Number of studies with Monoculture") +
  geom_bar()

c_data

```
![image](https://github.com/Eugenefut19/Meta-Analysis-on-Cover-Crops/assets/134546229/8b0c4b85-ec9b-4871-8c6a-c9a13f627897)

##Gives the count 
```{r}
data %>% 
  filter(data$cover.mono=="y") %>% 
  summarize(n=n())
```

#Diversity Summary

##Predator Diversity
```{r}
data.predator= data.metafor %>% 
  filter(arthropod.function=="predators")

predator.diversity= data.predator %>% 
  filter(metric.simplified=="diversity", vi <10)
  
predator.model.d= rma.mv(yi, vi, data=predator.diversity, random=~1|doi, method= "REML")
  
summary(predator.model.d)

forest(predator.model.d)
```
![image](https://github.com/Eugenefut19/Meta-Analysis-on-Cover-Crops/assets/134546229/b4616651-7353-45e8-a18e-11443e10d30a)

##Pest Diversity
```{r}
data.pest= data.metafor %>% 
  filter(arthropod.function=="pests")

pest.diversity= data.pest %>% 
  filter(metric.simplified=="diversity", vi <10)
  
predator.model.d= rma.mv(yi, vi, data=predator.diversity, random=~1|doi, method= "REML")
  
summary(predator.model.d)

forest(predator.model.d)
```
![image](https://github.com/Eugenefut19/Meta-Analysis-on-Cover-Crops/assets/134546229/e7d02dcd-704f-46f8-9a0a-9e203da50328)



##Diversity of All Insects
```{r}
#Provides summary of Insect Diversity model
#Creates a forest plot of Insect Diversity
data.diversity= data.metafor %>% 
  filter(metric.simplified=="diversity")

diversity.model.a= rma.mv(yi, vi, data=data.diversity, random=~1|doi, method= "REML")
  
summary(diversity.model.a)

forest(diversity.model.a, order="obs", cex=0.5)


mono.diversity.model= rma.mv(yi, vi, data=data.diversity, random=~1|doi, method= "REML")


```
![image](https://github.com/Eugenefut19/Meta-Analysis-on-Cover-Crops/assets/134546229/8ef4ff15-d606-4a8e-a554-df416f5f27cd)


#Abundance

##Predator Abundance
```{r}
predator.abundance= data.predator %>% 
  filter(metric.simplified=="abundance", vi <10)
  
predator.model.a= rma.mv(yi, vi, data=predator.abundance, random=~1|doi, method= "REML")
  
summary(predator.model.a)

forest(predator.model.a, cex=0.5, order="obs")
```
![image](https://github.com/Eugenefut19/Meta-Analysis-on-Cover-Crops/assets/134546229/7b163595-35cf-46c1-9eac-29df3b82edd6)


##Pest Abundance
```{r}
pest.abundance= data.pest %>% 
  filter(metric.simplified=="abundance", vi <10, vi>0)
  
pest.model.a= rma.mv(yi, vi, data=pest.abundance, random=~1|doi, method= "REML")
  
summary(pest.model.a)

forest(pest.model.a, cex=0.5, order="obs")

```
![image](https://github.com/Eugenefut19/Meta-Analysis-on-Cover-Crops/assets/134546229/9ccce26e-7e8a-4a79-b43c-fe6c366fac10)


#Mono vs Polyculture

##Mono vs Polycultrue in Pests Diversity
```{r}
mono.pests.metafor= rma.mv(yi, vi, data=predator.diversity, random=~1|doi, method= "REML", mod=~cover.mono)
summary(mono.pests.metafor)


forest(mono.pests.metafor$b, ci.lb=mono.pests.metafor$ci.lb, ci.ub=mono.pests.metafor$ci.ub, slab=c("Monoculture", "Polyculture"))
```
![image](https://github.com/Eugenefut19/Meta-Analysis-on-Cover-Crops/assets/134546229/5e16ab58-18b6-465a-b43a-6e4c0fb823c5)


##Mono vs Polyculture in Pests Abundance
```{r}
mono.pests.metafor= rma.mv(yi, vi, data=pest.abundance, random=~1|doi, method= "REML", mod=~cover.mono)
summary(mono.pests.metafor)

forest(mono.pests.metafor)

forest(mono.pests.metafor$b, ci.lb=mono.pests.metafor$ci.lb, ci.ub=mono.pests.metafor$ci.ub, slab=c("Monoculture", "Polyculture"))

```
![image](https://github.com/Eugenefut19/Meta-Analysis-on-Cover-Crops/assets/134546229/f8188bd1-4c45-4b58-a9ed-272e168c53df)


![image](https://github.com/Eugenefut19/Meta-Analysis-on-Cover-Crops/assets/134546229/89ba0fef-dad6-43a5-9d69-84be05945f28)


##Mono vs Polycultrue in Predator Diversity
```{r}
mono.predator.metafor= rma.mv(yi, vi, data=predator.diversity, random=~1|doi, method= "REML", mod=~cover.mono)
summary(mono.predator.metafor)


forest(mono.predator.metafor$b, ci.lb=mono.predator.metafor$ci.lb, ci.ub=mono.predator.metafor$ci.ub, slab=c("Monoculture", "Polyculture"))
```
![image](https://github.com/Eugenefut19/Meta-Analysis-on-Cover-Crops/assets/134546229/7c911cfe-f72d-43c7-8123-0c73f92dedd9)


##Mono vs Polyculture in Predator Abundance
```{r}
mono.predator.metafor= rma.mv(yi, vi, data=predator.abundance, random=~1|doi, method= "REML", mod=~cover.mono)
summary(mono.predator.metafor)

forest(mono.predator.metafor$b, ci.lb=mono.predator.metafor$ci.lb, ci.ub=mono.predator.metafor$ci.ub, slab=c("Monoculture", "Polyculture"))
```
![image](https://github.com/Eugenefut19/Meta-Analysis-on-Cover-Crops/assets/134546229/7064161b-ee18-41f8-93b8-76823b981056)

#Overall Insects Diversity Poly vs Mono
```{r}
diversity_overall= data.metafor %>% 
   filter(metric.simplified=="diversity", vi> 0)


diversity.metafor= rma.mv(yi, vi, data=diversity_overall, random=~1|doi, method= "REML", mod=~cover.mono)
summary(diversity.metafor)

forest(diversity.metafor$b, ci.lb=diversity.metafor$ci.lb, ci.ub=diversity.metafor$ci.ub, slab=c("Monoculture", "Polyculture"))
```
![image](https://github.com/Eugenefut19/Meta-Analysis-on-Cover-Crops/assets/134546229/bac715dd-db9c-44f3-824b-ef2c9c90f995)

#Synchronocity

##Synchronocity of Insect Predators
```{r}
sync.predator= rma.mv(yi, vi, data=data.predator, random= ~1|doi, method="REML")
yes.sync= rma.mv(yi, vi, data=data.predator, random= ~1|doi, method="REML", mod=~synchronous)
summary(yes.sync)

forest(yes.sync$b, ci.lb= yes.sync$ci.lb, ci.ub=yes.sync$ci.ub, slab= c("Synchrnous", "Not Synchronous"))
```
![image](https://github.com/Eugenefut19/Meta-Analysis-on-Cover-Crops/assets/134546229/c6c8b7e1-2d8b-4d82-a4c7-030467aef905)

#Other Visualizations

##Graphs Count of Insecticides
```{r}
insecticide.data= data %>% 
  filter(insecticide %in% c("y", "n")) %>% 
  group_by() %>% 
  ggplot(aes(x=insecticide))+
  geom_bar()
  
insecticide.data
```
![image](https://github.com/Eugenefut19/Meta-Analysis-on-Cover-Crops/assets/134546229/5f6fa22e-231b-4e92-b90a-f3993f21f3d2)

##Pie Chart of Different Crop Groups
```{r}
crop.group.plot= data %>%
  group_by(doi, crop.group) %>% 
  summarize(count=n()) %>% 
  group_by(crop.group) %>% 
  summarize(count=n()) %>%
  arrange(desc(crop.group)) %>% 
  mutate(percent=round(100* count/sum(count)), position=cumsum(percent)-0.5*(percent)) %>% 
  ggplot(aes(x="", y=percent, fill=crop.group))+   geom_bar(stat="identity", width=1)+
  coord_polar("y", start=0) +
  theme_void()+
  geom_text(aes(y=position, label=percent))

crop.group.plot
```
![image](https://github.com/Eugenefut19/Meta-Analysis-on-Cover-Crops/assets/134546229/6776efe3-9509-4936-95d4-fd3ab1386ede)

##Different types of Crops Grown
```{r}
crop.species.plot= data %>%
  group_by(doi, crop.species) %>% 
  summarise(count=n()) %>% 
  group_by(crop.species) %>% 
  summarize(count=n()) %>% 
  arrange(desc(crop.species)) %>% 
  mutate(percent=round(100* count/sum(count)), position=cumsum(percent)-0.5*(percent)) %>% 
  ggplot(aes(x="", y=percent, fill=crop.species))+   geom_bar(stat="identity", width=100)+
  coord_polar("y", start=0) +
  theme_void()+
  geom_text(aes(y=position, label=percent))

crop.species.plot
```
![image](https://github.com/Eugenefut19/Meta-Analysis-on-Cover-Crops/assets/134546229/2be46afe-27f2-4c68-81d8-e830683fc5c3)
**
