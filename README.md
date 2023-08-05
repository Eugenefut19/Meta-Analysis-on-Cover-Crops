
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
![image](https://github.com/Eugenefut19/Arthropod/assets/134546229/6cac845a-4120-48dc-b407-a80445a23474)

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
![image](https://github.com/Eugenefut19/Arthropod/assets/134546229/f11614d1-c398-48a5-b11c-cafac1624791)

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
![image](https://github.com/Eugenefut19/Arthropod/assets/134546229/74e71735-d76b-44fd-9269-8640651f74ff)

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
![image](https://github.com/Eugenefut19/Arthropod/assets/134546229/c0616519-4d0b-4a7c-abd6-aede35921735)


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
![image](https://github.com/Eugenefut19/Arthropod/assets/134546229/fc6f40c8-deac-4863-8476-c22a2b697eb4)


#Abundance

##Predator Abundance
```{r}
predator.abundance= data.predator %>% 
  filter(metric.simplified=="abundance", vi <10)
  
predator.model.a= rma.mv(yi, vi, data=predator.abundance, random=~1|doi, method= "REML")
  
summary(predator.model.a)

forest(predator.model.a, cex=0.5, order="obs")
```
![image](https://github.com/Eugenefut19/Arthropod/assets/134546229/1ba72280-ce30-4e4f-ab0a-ad0bddbd55e2)

##Pest Abundance
```{r}
pest.abundance= data.pest %>% 
  filter(metric.simplified=="abundance", vi <10, vi>0)
  
pest.model.a= rma.mv(yi, vi, data=pest.abundance, random=~1|doi, method= "REML")
  
summary(pest.model.a)

forest(pest.model.a, cex=0.5, order="obs")

```
![image](https://github.com/Eugenefut19/Arthropod/assets/134546229/4be4a995-a50d-4730-9003-3bb41ea0f558)

#Mono vs Polyculture

##Mono vs Polycultrue in Pests Diversity
```{r}
mono.pests.metafor= rma.mv(yi, vi, data=predator.diversity, random=~1|doi, method= "REML", mod=~cover.mono)
summary(mono.pests.metafor)


forest(mono.pests.metafor$b, ci.lb=mono.pests.metafor$ci.lb, ci.ub=mono.pests.metafor$ci.ub, slab=c("Monoculture", "Polyculture"))
```
![image](https://github.com/Eugenefut19/Arthropod/assets/134546229/27d9602e-7666-48b6-9c50-31542b916d3b)

##Mono vs Polyculture in Pests Abundance
```{r}
mono.pests.metafor= rma.mv(yi, vi, data=pest.abundance, random=~1|doi, method= "REML", mod=~cover.mono)
summary(mono.pests.metafor)

forest(mono.pests.metafor)

forest(mono.pests.metafor$b, ci.lb=mono.pests.metafor$ci.lb, ci.ub=mono.pests.metafor$ci.ub, slab=c("Monoculture", "Polyculture"))

```
![image](https://github.com/Eugenefut19/Arthropod/assets/134546229/d2d7b111-829d-4ecb-b9d3-16b52f424cc5)

![image](https://github.com/Eugenefut19/Arthropod/assets/134546229/62ddc189-48a7-4cfb-adb1-77d6cb23f045)

##Mono vs Polycultrue in Predator Diversity
```{r}
mono.predator.metafor= rma.mv(yi, vi, data=predator.diversity, random=~1|doi, method= "REML", mod=~cover.mono)
summary(mono.predator.metafor)


forest(mono.predator.metafor$b, ci.lb=mono.predator.metafor$ci.lb, ci.ub=mono.predator.metafor$ci.ub, slab=c("Monoculture", "Polyculture"))
```
![image](https://github.com/Eugenefut19/Arthropod/assets/134546229/54382982-90d8-4091-84b7-3fc0de267c28)

##Mono vs Polyculture in Predator Abundance
```{r}
mono.predator.metafor= rma.mv(yi, vi, data=predator.abundance, random=~1|doi, method= "REML", mod=~cover.mono)
summary(mono.predator.metafor)

forest(mono.predator.metafor$b, ci.lb=mono.predator.metafor$ci.lb, ci.ub=mono.predator.metafor$ci.ub, slab=c("Monoculture", "Polyculture"))
```
![image](https://github.com/Eugenefut19/Arthropod/assets/134546229/40d78141-93b8-41a1-85fa-913f614e4924)

#Overall Insects Diversity Poly vs Mono
```{r}
diversity_overall= data.metafor %>% 
   filter(metric.simplified=="diversity", vi> 0)


diversity.metafor= rma.mv(yi, vi, data=diversity_overall, random=~1|doi, method= "REML", mod=~cover.mono)
summary(diversity.metafor)

forest(diversity.metafor$b, ci.lb=diversity.metafor$ci.lb, ci.ub=diversity.metafor$ci.ub, slab=c("Monoculture", "Polyculture"))
```
![image](https://github.com/Eugenefut19/Arthropod/assets/134546229/32d6b838-bde3-4116-8063-88ac4f497c3b)

#Synchronocity

##Synchronocity of Insect Predators
```{r}
sync.predator= rma.mv(yi, vi, data=data.predator, random= ~1|doi, method="REML")
yes.sync= rma.mv(yi, vi, data=data.predator, random= ~1|doi, method="REML", mod=~synchronous)
summary(yes.sync)

forest(yes.sync$b, ci.lb= yes.sync$ci.lb, ci.ub=yes.sync$ci.ub, slab= c("Synchrnous", "Not Synchronous"))
```
![image](https://github.com/Eugenefut19/Arthropod/assets/134546229/5de216f7-df55-41fc-9587-b7a7294b6b3b)

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
![image](https://github.com/Eugenefut19/Arthropod/assets/134546229/cd725fef-c3ff-4988-8c38-65795615931f)

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
![image](https://github.com/Eugenefut19/Arthropod/assets/134546229/604843ea-95e3-41ea-98ce-3b7f46db821d)

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
![image](https://github.com/Eugenefut19/Arthropod/assets/134546229/4003dadb-01e7-4095-8c4a-81169cd17994)
**
