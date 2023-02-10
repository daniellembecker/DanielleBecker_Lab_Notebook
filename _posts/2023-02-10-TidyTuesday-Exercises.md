---
layout: post
title: TidyTuesday Exercises
date: '2023-02-10'
categories: Analysis
tags: R TidyTuesday
---
This post is a running log of miscellaneous coding and TidyTuesday exercises! 

# 13 February 2023 

We are doing a TidyTuesday exercise as a group activity for the Roberts Lab meeting today. We are analyzing the FeederWatch dataset with the goal that each person makes 1-2 fun plots. 

I was interested in 1) how bird abundance in Washington State varies by vegetation type and 2) how bird abundance responds to pet (dog and cat) abndance.

## Code  

Loading data for Tidy Tuesday lab meeting practice. 

#### *Setup*  

Load packages. 

```
library(tidyverse)
library(RColorBrewer)
```

Load data frame for FeederWatch

```
# Read in the data manually

feederwatch <- readr::read_csv('https://raw.githubusercontent.com/rfordatascience/tidytuesday/master/data/2023/2023-01-10/PFW_2021_public.csv')
site_data <- readr::read_csv('https://raw.githubusercontent.com/rfordatascience/tidytuesday/master/data/2023/2023-01-10/PFW_count_site_data_public_2021.csv')

```

### *Examine dataset*  

Look at feederwatch dataset. 

```
str(feederwatch)
```

Look at site dataset. 

```
str(site_data)
```

How does bird abundance relate to habitat type in Washington State? 

Plot: Distribution plot for bird abundance per unit effort for each habitat type 

How does bird abundance relate to pet abundance?  

Plot: Correlation plot/regression plot of total birds sighted across abundance of cats and dogs 

### *Preparing data* 

Next, prepare the data so that I can answer the above questions. 

We will need one dataframe that has site descriptions and bird counts for each site at each time point. 

### Manipulate feeder bird count data  

Calculate total birds sighted at each location over all time points normalized to the effort of observation (total birds normalized = birds / effort).  

```{r}
feeder_df<-feederwatch %>%
  select(loc_id, Month, Day, Year, subnational1_code, how_many, effort_hrs_atleast)%>%
  arrange(loc_id)%>%
  mutate(norm_count=how_many/effort_hrs_atleast)%>%
  group_by(loc_id, subnational1_code)%>%
  summarise(total_birds_norm=sum(how_many))

head(feeder_df)
```

Now manipulate the site data to be able to pull out cat and dog presence and types of trees. 

```{r}
site_df<-site_data%>%
  select(loc_id, evgr_trees_atleast, dcid_trees_atleast, fru_trees_atleast, cacti_atleast, cats, dogs)%>%
  gather(key="vegetation_type", value="tree_abundance", -loc_id, -cats, -dogs)%>% #gather 
  group_by(loc_id)%>%
  mutate(mean_cats=mean(cats, na.rm=TRUE))%>% #calculate mean dogs and cat abundance 
  group_by(loc_id)%>%
  mutate(mean_dogs=mean(dogs, na.rm=TRUE))%>%
  select(!cats)%>%
  select(!dogs)%>%
  unique()
```

Note that some sites have multiple habitat types. We will keep these for now. 

Now merge the two data frames together and subset by WA.  

```{r}
df<-left_join(feeder_df, site_df)

df<-df%>%
  rename("cats"="mean_cats", "dogs"="mean_dogs")%>%
  gather(key="animal", value="animal_abundance", cats:dogs)%>%
  mutate(vegetation_type=if_else(vegetation_type=="dcid_trees_atleast", "deciduous", 
                                 if_else(vegetation_type=="evgr_trees_atleast", "evergreen",
                                         if_else(vegetation_type=="fru_trees_atleast", "fruit", "NA"))))


df_WA<-df%>%
  filter(subnational1_code=="US-WA")
```

### Plots 

Plot effect of tree type and abundance on bird sightings.   

```{r}
plot1<-df_WA%>%
  filter(!vegetation_type=="NA")%>%
  filter(!vegetation_type=="cacti_atleast")%>%
  
  ggplot(aes(x=tree_abundance, y=total_birds_norm, group=vegetation_type, colour=vegetation_type))+
  #facet_wrap(~vegetation_type, ncol=3)+
  geom_point(aes(group=vegetation_type), alpha=0.2, position=position_jitterdodge(0.2))+
  stat_smooth(method="lm", show.legend=FALSE)+
  scale_color_brewer(palette="Dark2")+
  scale_fill_brewer(palette="Dark2")+
  theme_classic()+
  ylab("Birds (normalized to effort)")+
  xlab("Tree Abundance")+
  ggtitle("More birds with more trees, regardless of tree type!")+
  theme(
    legend.position="right",
    legend.title=element_blank()
  );plot1

ggsave("trees_birds.png", plot1, width=6, height=6)
```

Plot effects of pet abundance on bird sightings.  

```{r}
plot2<-df%>%
  ggplot(aes(x=animal_abundance, y=total_birds_norm, group=animal, colour=animal))+
  geom_smooth(aes(fill=animal))+
  scale_color_brewer(palette="Dark2")+
  scale_fill_brewer(palette="Dark2")+
  theme_classic()+
  ylab("Birds (normalized to effort)")+
  xlab("Animal Presence")+
  ggtitle("Birds are scarce when pets are sometimes present!")+
  theme(
    legend.position="right",
    legend.title=element_blank()
  );plot2 

ggsave("pets_birds.png", plot2, width=6, height=6)
```

[The full script can be found on GitHub here](https://github.com/AHuffmyer/TidyTuesday/blob/main/13Feb2023/script.Rmd)  

## Plot Results 

Bird abundance increases with increased vegetation, but there is no difference by type of trees.  

![](https://github.com/AHuffmyer/TidyTuesday/raw/main/13Feb2023/trees_birds.png)

Bird abundance is lowest when pets are *sometimes* present. 

![](https://github.com/AHuffmyer/TidyTuesday/raw/main/13Feb2023/pets_birds.png)




