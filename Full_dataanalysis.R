#### Full data analysis Sand lizards #####

rm(list=ls()) #Empty the global environment

# Load packages & set functions
library(MASS)
library(DescTools)
library(reshape2)
library(tidyverse)
library(RColorBrewer)
library(lubridate)
library(DHARMa)
library(emmeans)
library(lme4)
library(ggpubr)
library(EnvStats)
library(grid)
library(ggpmisc)

select <- dplyr::select #always use select() in dplyr

# Function to produce summary statistics
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymin <- ifelse(ymin < 0, 0, ymin)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

###############################################################################
# L. AGILIS TICK BURDEN

# Import data from github
lizards.complete <- read_csv(url("https://raw.githubusercontent.com/clarakoehler/sandlizards_TBD/main/sandlizards.csv"))

#remove incomplete case
l <- lizards.complete %>% filter(sex == "m" | sex=="f")

l <- l %>% #add variables day of the year and life stage/age (group)
  mutate(day = as.Date(date, format="%d/%m")) %>% 
  mutate(daynr = yday(day),
         month = as.factor(month(day)))

l %>% summarise(larval.burden = mean(larvae),
                nymphal.burden = mean(nymph),
                total.burden = mean(total),
                total.larvae = sum(larvae),
                total.nymphs = sum(nymph))

# Nagative-binomial GLM for individual tick burden
summary(m1 <- glm.nb(total ~ habitat + sex + weight + daynr, l))

# test model assumptions with DHARMa package
testDispersion(m1) # test dispersion of model
sim.res <- simulateResiduals(fittedModel = m1, plot = T) #simulate residuals

par(mfrow = c(1,1))
testDispersion(sim.res)
testZeroInflation(sim.res)
testQuantiles(sim.res)

# Estimated marginal means with emmeans package
emmeans(m1, "habitat", type="response")
emmeans(m1, "sex", type="response")

exp(0.089361) #back-transform weight
exp(0.002206) #back-transform daynr

drop1(m1, test="LRT") #drop variables for likelihood-ratio tests

## Plot lizard tick burden per habitat and sex
ggplot(data = l, aes(x = fct_reorder(habitat, total), y =total, fill=sex), color="black")+
  theme_linedraw(base_size = 14)+
  theme(panel.grid.major = element_blank())+
  stat_boxplot(geom = 'errorbar', width = 0.4, position = position_dodge(0.7))+
  geom_boxplot(outlier.shape = NA, width = 0.6, position = position_dodge(0.7))+
  ylab("Individual tick burden (larvae and nymphs)") + ylim(0,55)+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle=45, hjust=1, size=14), 
        axis.title.y = element_text(size=14, face="bold"),
        axis.text.y = element_text(size=14),
        legend.position="top",
        legend.text = element_text(size=14))+
  scale_fill_manual(name = element_blank(), labels = c("Female", "Male"),
                    values = c("#FF9999", "#9999FF"))+
  geom_point(position = position_jitterdodge(0.3), size= 0.9, alpha=0.8)

rm(lizards.complete, m1, sim.res, l)

###############################################################################
# SPECIES SPECIFIC TICK BURDEN & INFECTION PREVALENCE

#load the database
database <- read_csv(url("https://raw.githubusercontent.com/clarakoehler/sandlizards_TBD/main/database_complete.csv"))

species <- database %>% group_by(Species.scientific) %>%
  summarize(hosts.studied = sum(na.omit(Hosts.studied)),
            all.larvae = sum(na.omit(Larvae.on.hosts)),
            BBsl.Hosts.tested = sum(na.omit(BBsl.Hosts.tested)),
            BBsl.Hosts.positive = sum(na.omit(BBsl.Hosts.positive)),
            BBsl.Larvae.tested = sum(na.omit(BBsl.Larvae.tested)),
            BBsl.Larvae.positive = sum(na.omit(BBsl.Larvae.positive)),
            AP.Hosts.tested = sum(na.omit(AP.Hosts.tested)),
            AP.Hosts.positive = sum(na.omit(AP.Hosts.positive)),
            AP.Larvae.tested = sum(na.omit(AP.Larvae.tested)),
            AP.Larvae.positive = sum(na.omit(AP.Larvae.positive))) %>% 
  mutate(mean.larval.burden = round(all.larvae / hosts.studied, digits=3),
         host.prev = round(BBsl.Hosts.positive / BBsl.Hosts.tested, digits=3),
         prev.larvae = round(BBsl.Larvae.positive / BBsl.Larvae.tested, digits=3),
         host.prev.a = round(AP.Hosts.positive / AP.Hosts.tested, digits = 3),
         prev.larvae.a = round(AP.Larvae.positive / AP.Larvae.tested, digits=3))

#Add own data on L. agilis
dbl <- species %>%
  add_row(Species.scientific = "Lacerta agilis", hosts.studied = 100, 
          all.larvae = 1048, 
          mean.larval.burden = round(all.larvae/hosts.studied, digits=3),
          BBsl.Larvae.tested = 1031, BBsl.Larvae.positive = 2, 
          prev.larvae = round(BBsl.Larvae.positive/BBsl.Larvae.tested, digits=3),
          AP.Larvae.tested = 1031, AP.Larvae.positive = 5,
          prev.larvae.a = round(AP.Larvae.positive/AP.Larvae.tested, digits=3))

#add reservoir competence per species with cutoff deciding whether the prevalence 
#in hosts or the prevalence in feeding larvae is used
dbl <- dbl %>% 
  mutate(rc = ifelse(BBsl.Larvae.tested >= 20, prev.larvae,
                     ifelse(BBsl.Larvae.tested < 20 & BBsl.Hosts.tested >= 20, 
                            host.prev, NA)),
         rc.a = ifelse(AP.Larvae.tested >= 20, prev.larvae.a,
                       ifelse(AP.Larvae.tested < 20 & AP.Hosts.tested >= 20, 
                              host.prev.a, NA)))

# Add A. phagocytophilum results of Fabri et.al. 2021
mutate_cond <- function(.data, condition, ..., envir = parent.frame()) {
  condition <- eval(substitute(condition), .data, envir)
  .data[condition, ] <- .data[condition, ] %>% mutate(...)
  .data
}

dbl <- dbl %>% mutate_cond(Species.scientific == 'Dama dama', rc.a = 0.79)

dbl.reduced <- dbl %>% #keep only relevant columns
  select(Species.scientific, mean.larval.burden, rc, rc.a)

###############################################################################
# CALCULATE COMMUNITY COMPETENCE

#create dataframe containing species densities per habitat type
Species.scientific <- c("Lacerta agilis", "Apodemus sylvaticus", 
                        "Erinaceus europaeus", "Myodes glareolus", 
                        "Sciurus vulgaris", "Dama dama", "Fringilla coelebs",
                        "Parus major", "Troglodytes troglodytes", 
                        "Turdus merula", "Turdus philomelos")
Forest <- c(50, 700, 10, 100, 25, 60, 400, 300, 100, 100, 100)
Open_forest <- c(2000, 500, 10, 500, 0, 60, 100, 100, 100, 50, 50)
Dense_shrub <- c(500, 1500, 25, 1000, 0, 60, 100, 50, 100, 50, 50)
Open_shrub <- c(6000, 500, 10, 500, 0, 60, 50, 25, 50, 10, 10)
Grassland <- c(3000, 500, 5, 500, 0, 60, 50, 25, 50, 10, 10)
Bracken <- c(2000, 1500, 25, 1000, 0, 60, 50, 25, 50, 50, 50)

dens <- data.frame(Species.scientific, Forest, Open_forest, Dense_shrub, 
                 Open_shrub, Grassland, Bracken)

comp.sp <- dbl.reduced %>% #add densities to data.frame
  left_join(dens, by= "Species.scientific")

long <- comp.sp %>% #elongate data.frame
  pivot_longer(cols = (Forest : Bracken), names_to = "habitat", values_to = "density")

rm(dens, dbl.reduced, mutate_cond, comp.sp, Forest, Open_forest, Dense_shrub, Open_shrub,
   Grassland, Bracken, Species.scientific, database)

long <- long %>% #larvae fed per species per habitat type
  mutate(larvae.fed = mean.larval.burden * density) 

tmp <- long %>% group_by(habitat) %>% #sum of larvae fed per habitat type
  summarize(all.larvae.fed = sum(larvae.fed))

long <- long %>% left_join(tmp, by = "habitat") #add to df

long <- long %>% # larvae infected by each species per habitat
  mutate(inf.bsl = larvae.fed * rc, inf.ap = larvae.fed * rc.a)

tmp <- long %>% group_by(habitat) %>% 
  summarize(all.inf.bsl = sum(na.omit(inf.bsl)),
            all.inf.ap = sum(na.omit(inf.ap)))

long <- long %>% left_join(tmp, by = "habitat")

rm(tmp)

long <- long %>% #community competence
  mutate(cc.bsl = all.inf.bsl/all.larvae.fed,
         cc.ap = all.inf.ap/all.larvae.fed)

per.habitat <- distinct(long, habitat, .keep_all = TRUE)

# Realized community competence
per.habitat <- per.habitat %>% #scale down to per 100m
  mutate(absolute.cc.bsl = all.inf.bsl/10000,
         absolute.cc.ap = all.inf.ap/10000)

per.habitat <- per.habitat %>% 
  select(habitat, absolute.cc.bsl, absolute.cc.ap, cc.bsl, cc.ap)

long$habitat <- gsub('_', ' ', long$habitat) #remove underscores

# Stacked bar plot for species contribution to feeding/infecting larvae
tmp <- long %>% 
  mutate(ri.species = larvae.fed / all.larvae.fed, #relative importance of each species
         ri.infecting = inf.bsl / all.larvae.fed, 
         ri.infecting.a = inf.ap /all.larvae.fed) %>% 
  mutate(Species = case_when(Species.scientific == "Lacerta agilis" ~ Species.scientific,
                             Species.scientific == "Apodemus sylvaticus" ~ "Rodents",
                             Species.scientific == "Myodes glareolus" ~ "Rodents",
                             Species.scientific == "Erinaceus europaeus" ~ Species.scientific,
                             TRUE ~ "Others")) %>% #make groups
  mutate(Species = factor(Species, levels = c("Lacerta agilis", "Rodents", 
                                              "Erinaceus europaeus", "Others")),
         habitat = ordered(habitat, levels=c("Open shrub","Grassland", "Open forest", 
                                             "Bracken", "Dense shrub", "Forest")))

plot.frame <- tmp %>% group_by(habitat, Species) %>% #summarize data frame
  summarize(density = sum(density),
            ri.species = sum(ri.species),
            ri.infecting=sum(ri.infecting),
            ri.infecting.a = sum(ri.infecting.a, na.rm=TRUE))

# Stacked bar plot relative contribution to feeding larvae
feeding <- ggplot(plot.frame, aes(x=habitat), color="black")+ 
  geom_col(aes(y=ri.species, fill=Species), position=position_stack(), colour="black")+
  theme_classic(base_size=12)+
  scale_fill_manual(values=c("#79AF97CC","#DF8F44CC","#6A6599CC","#525252"), 
                    labels=c(expression(paste("Sand lizard", italic(" (L. agilis)"))), 
                             expression(paste("Rodents", italic(" (A. sylvaticus & M. glareolus)"))),
                             expression(paste("Hedgehog", italic(" (E. europaeus)"))), "Others")) +
  ylab("Proportion of larvae fed") +
  theme(axis.text.x=element_text(angle=40, hjust=1, size=14, color="black"), axis.title.x=element_blank(), 
        axis.title.y=element_text(face="bold", size=16), axis.text.y=element_text(size=14, color="black"))+
  theme(legend.text.align=0, legend.position="right", legend.spacing.y=unit(8, "mm"), 
        legend.background=element_blank(), legend.title=element_blank(),
        legend.text=element_text(size=14, color="black"))+
  guides(fill=guide_legend(byrow=TRUE))

print(feeding)

# Stacked bar plot relative contribution to infecting larvae with B. burgdorferi s.l.
bor <- ggplot(plot.frame, aes(x = habitat, y = ri.infecting, fill = Species), color="black") + 
  geom_col(position = position_stack(), colour="black")+
  theme_classic(base_size = 12)+
  scale_fill_manual(values=c("#79AF97CC","#DF8F44CC","#6A6599CC","#525252"), 
                    labels=c(expression(paste("Sand lizard ", italic("(L. agilis)"))), 
                             expression(paste("Rodents ", italic("(A. sylvaticus & M. glareolus)"))),
                             expression(paste("Hedgehog ", italic("(E. europaeus)"))),"Others")) +
  ylab("Proportion of larvae infected") + ylim(0,0.25) +
  theme(axis.text.x=element_text(angle=40, hjust=1, size=14, color="black"), 
        axis.title.x=element_blank(), 
        axis.title.y=element_text(face="bold", size=16), 
        axis.text.y=element_text(size=14, color="black"))+
  theme(legend.text.align=0, legend.position="right", legend.spacing.y=unit(3, "mm"), 
        legend.background=element_blank(), legend.title=element_blank(),
        legend.text=element_text(size=12))+
  guides(fill=guide_legend(byrow=TRUE))

print(bor)

# Stacked bar plot relative contribution to infecting larvae with A. phagocytophilum
ana <- ggplot(plot.frame, aes(x = habitat, y = ri.infecting.a, fill = Species), color="black") + 
  geom_col(position = position_stack(), colour="black")+
  theme_classic(base_size = 12)+
  scale_fill_manual(values=c("#79AF97CC","#DF8F44CC","#6A6599CC","#525252"), 
                    labels=c(expression(paste("Sand lizard ", italic("(L. agilis)"))), 
                             expression(paste("Rodents ", italic("(A. sylvaticus & M. glareolus)"))),
                             expression(paste("Hedgehog ", italic("(E. europaeus)"))),"Others")) +
  ylab("Proportion of larvae infected") + ylim(0,0.25) +
  theme(axis.text.x=element_text(angle=40, hjust=1, size=14, color="black"), 
        axis.title.x=element_blank(), 
        axis.title.y=element_text(face="bold", size=16), 
        axis.text.y=element_text(size=14, color="black"))+
  theme(legend.text.align=0, legend.position="right", legend.spacing.y=unit(3, "mm"), 
        legend.background=element_blank(), legend.title=element_blank(),
        legend.text=element_text(size=12))+
  guides(fill=guide_legend(byrow=TRUE))

ggarrange(feeding + rremove("legend") + rremove("x.text"), #combine plots
          feeding+ rremove("legend") + rremove("y.title") + rremove("x.text"), 
          bor + rremove("legend"), ana + rremove("legend") + rremove("y.title"),
          ncol=2, nrow=2, align = "hv")

rm(dbl, tmp, long, ana, bor, feeding, plot.frame)

###############################################################################
# ADD PCR RESULTS OF QUESTING TICKS

# Load data on questing ticks
plots <- read_csv(url("https://raw.githubusercontent.com/clarakoehler/sandlizards_TBD/main/questing_ticks.csv"))

plots <- plots %>% filter(habitat != "Short grass" & habitat !="High grass") %>% 
  select(date:t.transect, latitude, longitude)

pcr.all <- read_csv(url("https://raw.githubusercontent.com/clarakoehler/sandlizards_TBD/main/questingticks_PCR_results.csv"))

pcr.nymphs <- pcr.all[(pcr.all$Stadia == "N"), ] #remove adults

#select samples where the procedure worked
pcr.success <- pcr.nymphs %>% filter(is.na(Remarks) | Remarks != "DNA coloured" &
                             Remarks != "Gezogen/DNAcoloured")

# make new column with number of samples initially analyzed per plot
tmp <- pcr.nymphs %>% 
  select(PlotNr) %>% 
  count(as.character(PlotNr)) %>% 
  as_tibble() %>% 
  select(volgorde = `as.character(PlotNr)`, 
         analyzed = n)

plots <- plots %>% 
  mutate(volgorde = volgorde %>% as.character()) %>% 
  left_join(tmp, by = "volgorde")

rm(tmp)

tmp <- pcr.success %>% #add column with number of samples successfully analyzed
  select(PlotNr) %>% 
  count(as.character(PlotNr)) %>% 
  as_tibble() %>% 
  select(volgorde = `as.character(PlotNr)`, 
         successf.analyzed = n)

plots <- plots %>% 
  mutate(volgorde = volgorde %>% as.character()) %>% 
  left_join(tmp, by = "volgorde")

rm(tmp)

# B. burgdorferi s.l. & A. phagocytophilum 
tmp <- pcr.success %>% #count positive samples per plot
  group_by(PlotNr) %>%
  summarize(positives.b = sum(trans_borrelia),
            positives.a = sum(trans_anaplasma)) %>%
  as_tibble() %>% 
  select(volgorde=PlotNr, positives.b=positives.b, positives.a=positives.a)

plots <- plots %>% #add as column to main data
  mutate(volgorde = volgorde %>% as.character()) %>% 
  left_join(tmp, by = "volgorde")

rm(tmp)

plots <- plots %>% #add prevalence to plot data.frame
  mutate(prev.b = as.numeric(positives.b) / as.numeric(successf.analyzed),
         prev.a = as.numeric(positives.a) / as.numeric(successf.analyzed))

plots %>% #summary stats for report
  summarize(all.successf.analyzed = sum(na.omit(successf.analyzed)),
            positives.b.sum = sum(na.omit(positives.b)),
            positives.a.sum = sum(na.omit(positives.a)))

# Genotypes & Ecotypes
tmp <- pcr.success %>% #count positive samples per plot
  group_by(PlotNr) %>%
  summarize(afzelii = sum(genospecies=="afzelii"),
            garinii = sum(genospecies=="garinii"),
            valaisiana = sum(genospecies=="valaisiana"),
            spielmanii = sum(genospecies=="spielmanii"),
            bavariensis = sum(genospecies=="bavariensis"),
            sensustricto = sum(genospecies=="sensustricto"),
            eco1 = sum(eco1 == 1), eco2 = sum(eco2 == 1)) %>%
  as_tibble() %>% 
  select(volgorde=PlotNr, afzelii=afzelii, garinii=garinii, 
         valaisiana=valaisiana, spielmanii=spielmanii, bavariensis=bavariensis,
         sensustricto=sensustricto, eco1=eco1, eco2=eco2)

plots_genotypes <- plots %>% #add as column to main data
  mutate(volgorde = volgorde %>% as.character()) %>% 
  left_join(tmp, by = "volgorde")

rm(tmp)

plots_genotypes <- plots_genotypes %>% #number of successfully genotyped nymphs & prevalence per plot
  mutate(success.genotyped = afzelii + garinii + valaisiana + spielmanii + 
           bavariensis + sensustricto) %>% 
  mutate(fraction.afz = afzelii / success.genotyped,
         fraction.gar = garinii / success.genotyped,
         fraction.spi = valaisiana / success.genotyped,
         fraction.val = spielmanii / success.genotyped,
         fraction.ss = bavariensis / success.genotyped,
         fraction.bav = sensustricto / success.genotyped) %>% 
  mutate(prev.afz = (fraction.afz*positives.b)/successf.analyzed,
         prev.gar = (fraction.gar*positives.b)/successf.analyzed, 
         prev.spi = (fraction.spi*positives.b)/successf.analyzed,
         prev.val = (fraction.val*positives.b)/successf.analyzed,
         prev.ss = (fraction.ss*positives.b)/successf.analyzed,
         prev.bav = (fraction.bav*positives.b)/successf.analyzed)


plots_genotypes %>% #summary stats for report
  summarize(all.successf.genotyped = sum(na.omit(success.genotyped)),
            positives.afz = sum(na.omit(afzelii)),
            positives.gar = sum(na.omit(garinii)),
            positives.spi = sum(na.omit(spielmanii)),
            positives.ss = sum(na.omit(sensustricto)),
            eco1 = sum(na.omit(eco1)),
            eco2 = sum(na.omit(eco2)))

rm(pcr.all, pcr.nymphs, pcr.success, plots_genotypes)

###############################################################################
### Add CC of the surrounding habitat
###############################################################################
edges <- read_csv(url("https://raw.githubusercontent.com/clarakoehler/sandlizards_TBD/main/surrounding_habitat.csv"))

edges <- edges %>% filter(Struco > 20) %>% 
  mutate(volgorde=as.factor(volgorde),
         Struco=as.character(Struco))

tmp <- edges %>% group_by(volgorde) %>% 
  summarize(tot.area=sum(area_2))

edges <- edges %>%
  mutate(volgorde=as.character(volgorde)) %>% 
  left_join(tmp, by = "volgorde") %>% 
  mutate(volgorde=as.factor(volgorde))

edges <- edges %>% 
  mutate(struco = case_when(Struco == 10 ~ "Sand/Moss",
                            Struco == 20 ~ "Low herbaceous",
                            Struco == 40 ~ "Grassland",
                            Struco == 50 ~ "Open shrub", Struco == 61 ~ "Dense shrub",
                            Struco == 100 ~ "Open forest", Struco == 110 ~ "Forest",
                            Struco == 130 ~ "Bracken", TRUE ~ "Other"))  
edges <- edges %>%
  group_by(volgorde, struco) %>% 
  summarize(area_2=sum(area_2),
            tot.area=first(tot.area)) %>% 
  mutate(relative.area = area_2/tot.area)

edges <- edges %>% 
  mutate(cc.bsl = case_when(struco == "Forest" ~ 0.24,
                            struco == "Dense shrub" ~ 0.21,
                            struco == "Bracken" ~ 0.13,
                            struco == "Open forest" ~ 0.07,
                            struco == "Grassland" ~ 0.05,
                            struco == "Open shrub" ~ 0.03,
                            struco == "Other" ~ 0, TRUE ~ 0),
         cc.ap = case_when(struco == "Forest" ~ 0.04,
                           struco == "Dense shrub" ~ 0.04,
                           struco == "Bracken" ~ 0.02,
                           struco == "Open forest" ~ 0.02,
                           struco == "Grassland" ~ 0.01,
                           struco == "Open shrub" ~ 0.01,
                           struco == "Other" ~ 0, TRUE ~ 0))

tmp<- edges %>% group_by(volgorde) %>% 
  summarize(sur.cc.bsl = mean(cc.bsl),
            sur.cc.ap = mean(cc.ap))

plots <- plots %>% 
  left_join(tmp, by="volgorde")

rm(edges,tmp)


###############################################################################
# Density of Nymphs (DON)
###############################################################################

plots <- plots %>% # Create date object and add month & day variables
  mutate(date = date %>% as.Date(.,"%d/%m/%y")) %>%
  mutate(month = month(date) %>% as.factor(),
         daynr = yday(date),
         habitat = ordered(habitat, 
                           levels=c("Forest", "Dense shrub","Bracken", 
                                    "Open forest", "Open shrub", "Grassland"))) %>% 
  mutate(habitat = habitat %>% factor(ordered = FALSE))

# GLM for DON ####
summary(m1 <- glm.nb(n.transect ~ habitat + sd + daynr, data = plots))

VIF(m1) #variance inflation factor

# test model assumptions
testDispersion(m1)
sim.res <- simulateResiduals(fittedModel = m1, plot = T)

plotResiduals(sim.res, plots$habitat)
testDispersion(sim.res)
testZeroInflation(sim.res)
testQuantiles(sim.res)

# Test if negative binomial is the right distribution
1 - pchisq(summary(m1)$deviance, summary(m1)$df.residual)
#The model fits a negative binomial distribution if P > 0.05

# Estimated marginal means with emmeans package
emmeans(m1, "habitat", type="response")

exp(-0.024477)-1 #backtransform SD
exp(-0.018573)-1 #backtransform daynr

drop1(m1, test="LRT") #drop variables for likelihood-ratio tests

# Make estimates at mean temperature
est.DON <- data.frame(daynr = mean(plots$daynr),
                      sd = mean(plots$sd),
                      habitat = factor(1:6, levels = 1:6, 
                                       labels = levels(plots$habitat)))

est.DON$DON <- predict(m1, est.DON, type = "response")
est.DON

plots <- plots %>% 
  left_join(est.DON[, (3:4)]) #add estimated DON to plot data

# box plot DON ####
plots$habitat <- ordered(plots$habitat, 
                         levels=c("Open shrub", "Grassland","Open forest", 
                                  "Bracken", "Dense shrub", "Forest"))

ggplot(data = plots, aes(x = habitat, y =n.transect), color="black")+
  geom_jitter(width = 0.25, colour = "grey50", shape=1) + theme_bw(base_size = 12) +
  theme(panel.grid.major = element_blank())+
  stat_summary(fun.data=data_summary, color="black")+
  ylab(expression(bold(paste("Density of questing nymphs (DON), ind./100 ", m^2)))) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle=45, hjust=1, size=12, color="black"), 
        axis.text.y = element_text(size=12, color="black"))

rm(est.DON, m1, sim.res)

###############################################################################
# FIELD DATA ANALYSIS
###############################################################################

# Add community competence per habitat to plot data
per.habitat$habitat <- gsub('_', ' ', per.habitat$habitat)

plots <- plots %>% 
  left_join(per.habitat, by="habitat")

p <- plots %>% #subset to include only plots with ticks
  filter(successf.analyzed >=1)

# total analysed and positives for manuscript
p %>% summarize(total.analysed = sum(successf.analyzed),
                total.positive.b = sum(positives.b),
                total.positive.a = sum(positives.a))

tmp <- p %>% group_by(habitat) %>% # calculate prevalence per habitat type
  summarize(total.analysed = sum(successf.analyzed),
            total.positive.b = sum(positives.b),
            total.positive.a = sum(positives.a),
            prev.b = mean(prev.b),
            prev.a = mean(prev.a)) %>% 
  mutate(prev.sum.b = total.positive.b / total.analysed,
         prev.sum.a = total.positive.a / total.analysed) %>% 
  select(habitat, prev.sum.b, prev.sum.a)

per.habitat <- per.habitat %>% #add to per habitat data
  left_join(tmp, by = "habitat")

p <- p %>% 
  left_join(tmp, by = "habitat") #add to plot data
p <- p %>% 
  mutate(DIN.b.all = DON * prev.sum.b, #DIN per habitat
         DIN.a.all = DON * prev.sum.a)

rm(tmp)

p %>% group_by(habitat) %>% #summarize DIN as mean of all plots
  summarize(din.b = round(mean(DIN.b.all), digits = 3),
            din.a = round(mean(DIN.a.all), digits = 3))


# Barplots B. burgdorferi s.l. ####
p$habitat <- ordered(p$habitat, #keep same order
                     levels=c("Open shrub", "Grassland", "Open forest", 
                              "Bracken", "Dense shrub", "Forest"))

per.habitat$habitat <- ordered(per.habitat$habitat, #keep same order
                               levels=c("Open shrub", "Grassland", "Open forest", 
                                        "Bracken", "Dense shrub", "Forest"))

per.habitat$CC = "CC"

borrelia <- grobTree(textGrob("B. burgdorferi s.l.", x=0.025,  y=0.94, hjust=0, #create text
                              gp=gpar(col="black", fontsize=12, fontface="bold.italic")))

# Barplot NIP ####
nip.b <- ggplot(data = per.habitat, aes(x=habitat, y=prev.sum.b), color="black")+
  geom_bar(stat="identity", fill="slategray2") + theme_linedraw(base_size=12) +
  theme(panel.grid.major = element_blank())+
  geom_point(data=p, aes(x=habitat, y=prev.b, colour=successf.analyzed), size=1, shape=1,
             position = position_jitter(width = 0.25)) +
  scale_color_gradient(low="grey70", high="grey5", name="N")+
  #ylab("Nymphal infection prevalence (NIP)") +
  ylim(0,0.4)+
  annotation_custom(borrelia) +
  geom_point(data=per.habitat, aes(x=habitat, y=cc.bsl, fill=CC),
             color='red2', size=8, shape=42)+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), 
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=12),
        legend.position = "none")

# Barplot DIN #####
df1 <- p %>% group_by(habitat) %>% #Borrelia
  summarize(din.b = round(mean(DIN.b.all), digits = 2),
            absolute.cc = round(mean(absolute.cc.bsl), digits=2))
df2 <- melt(df1, id.vars='habitat')
df2$variable <- ordered(df2$variable, 
                        levels=c("absolute.cc", "din.b"))
df1$RCC = "RCC"
df1$habitat <- ordered(df1$habitat, #keep same order
                       levels=c("Open shrub", "Grassland", "Open forest", 
                                "Bracken", "Dense shrub", "Forest"))

din.b <- ggplot(data = df1, aes(x=habitat, y=din.b), color="black")+
  geom_bar(stat="identity", fill="slategray2") + theme_linedraw(base_size=12) +
  theme(panel.grid.major = element_blank())+
  #ylab(expression(bold(paste("Density of infected nymphs, ind./100 ", m^2)))) +
  ylim(0,2)+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1, size=12), 
        axis.text.y = element_text(size=12),
        legend.position = "none")+
  annotation_custom(borrelia) +
  geom_point(data=df1, aes(x=habitat, y=absolute.cc, fill=RCC),
             color='red2', size=8, shape=42)


# Barplots Anaplasma ####
anaplasma <- grobTree(textGrob("A. phagocytophilum", x=0.025,  y=0.94, hjust=0,
                               gp=gpar(col="black", fontsize=12, #create text
                                       fontface="bold.italic")))

# Barplot NIP Anaplasma
nip.a <- ggplot(data = per.habitat, aes(x = habitat, y = prev.sum.a), color="black") +
  geom_bar(stat="identity", fill="slategray2") + theme_linedraw(base_size = 12) +
  theme(panel.grid.major = element_blank()) +
  geom_point(data=p, aes(x=habitat, y=prev.a, colour=successf.analyzed), 
             size=1, shape=1, position = position_jitter(width = 0.25)) +
  scale_color_gradient(low = "grey70", high = "grey5", name="N") +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  ylim(0,0.4) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size=12),
        legend.position = "none") +
  annotation_custom(anaplasma) +
  geom_point(data=per.habitat, aes(x=habitat, y=cc.ap, fill=CC),
             color='red2', size=8, shape=42)

# Barplot DIN Anaplasma
df.ap <- p %>% group_by(habitat) %>% #Anaplasma
  summarize(din.a = round(mean(DIN.a.all), digits = 3),
            absolute.cc = round(mean(absolute.cc.ap), digits=3))
df.ap2 <- melt(df.ap, id.vars='habitat')
df.ap2$variable <- ordered(df.ap2$variable, 
                           levels=c("absolute.cc", "din.a"))
df.ap$RCC = "RCC"
df.ap$habitat <- ordered(df.ap$habitat, #keep same order
                         levels=c("Open shrub", "Grassland", "Open forest", 
                                  "Bracken", "Dense shrub", "Forest"))

din.a <- ggplot(data = df.ap, aes(x=habitat, y=din.a), color="black")+
  geom_bar(stat="identity", fill="slategray2") + theme_linedraw(base_size=12) +
  theme(panel.grid.major = element_blank())+
  ylim(0,2)+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1, size=12), 
        axis.text.y = element_text(size=12),
        legend.position = "none")+
  annotation_custom(anaplasma) +
  geom_point(data=df.ap, aes(x=habitat, y=absolute.cc, fill=RCC),
             color='red2', size=8, shape=42)

# add up plots
fig3 <- ggarrange(nip.b, nip.a, din.b, din.a, 
                  labels = c("A", "B", "C", "D"), label.x=0.9, label.y=0.96,
                  ncol = 2, nrow = 2, align = "hv")
print(fig3)

rm(nip.b, nip.a, din.b, din.a, anaplasma, borrelia,df1,df2,fig3, df.ap, df.ap2)

# GLM NIP B. burgdorferi s.l. ####
p <- p %>% # create binary variables
  mutate(negatives.b = successf.analyzed - positives.b)

p %>% group_by(habitat) %>% 
  summarize(all = sum(successf.analyzed),
            positives = sum(positives.b)) %>% 
  mutate(prev = positives / all)

# Binomial GLM B. burgdorferi s.l. ####

summary(m1 <- glmer(cbind(positives.b, negatives.b) ~ cc.bsl + (1|volgorde), 
                    data=p, family=binomial))

#model including CC of surrounding habitat
summary(m2 <- glmer(cbind(positives.b, negatives.b) ~ cc.bsl + sur.cc.bsl + (1|volgorde), 
                    data=p, family=binomial))

#check model assumptions
testDispersion(m2)

sim.res <- simulateResiduals(fittedModel = m2, plot = T)
testDispersion(sim.res)

exp(7.8231)-1

rm(m1, m2, sim.res)

# GLM NIP Anaplasma ####
p <- p %>% mutate(negatives.a = successf.analyzed - positives.a)

summary(m3 <- glmer(cbind(positives.a, negatives.a) ~ cc.ap + (1|volgorde), 
                    data=p, family=binomial))

summary(m4 <- glmer(cbind(positives.a, negatives.a) ~ cc.ap + sur.cc.ap + (1|volgorde), 
                    data=p, family=binomial))

# daynr
summary(m5 <- glmer(cbind(positives.a, negatives.a) ~ cc.ap + daynr + (1|volgorde), 
                    data=p, family=binomial)) #anaplasma

summary(m6 <- glmer(cbind(positives.b, negatives.b) ~ cc.bsl + daynr + (1|volgorde), 
                    data=p, family=binomial)) #borrelia

#check model assumptions
testDispersion(m5)
sim.res <- simulateResiduals(fittedModel = m5, plot = T)
testDispersion(sim.res)

rm(m3, m4, m5, m6, p, sim.res)

