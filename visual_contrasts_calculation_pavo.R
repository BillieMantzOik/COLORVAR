#Visual modelling with pavo
#scripts by Ariel Rodriguez and Vasiliki Mantzana Oikonomaki
#Data from O. granulifera, O. pumilio collected in Panama and COsta Rica, May-June 2022
setwd('C:/visualmodelling_Vasiliki')
library(pavo)
# First import the frog dorsal reflectances
pumilio.data <- as.rspec(read.csv("./Frog_reflectance_spectrometer_Panama.csv"),
                         lim = c(300, 700))


granu.data <- as.rspec(read.csv
                       ("./Spectrometer_data_O_granulifera_new.csv"),
                       lim = c(300, 700))


############CONNECT FROG DATA###########
all.data <- cbind(granu.data, pumilio.data[-1])
localities <- colnames(all.data[-1])
locality <- gsub("\\.[0-9]+","",localities)

# Check that there are 4 x 8 =32 meassurements by locality
table(locality)
head(all.data)
frog.spectra <- procspec(all.data, 
                         fixneg = "zero", 
                         opt = 'smooth', 
                         span = 0.2)

mspectra <- aggspec(frog.spectra, 
                    by = locality, 
                    FUN = mean)

pdf("frogs_spectra_by_id.pdf")
explorespec(frog.spectra, by=4, ylim=c(0,60))
dev.off()

head(mspectra)
write.csv(mspectra, 'allfrogsmean.csv')
plot(mspectra)

pdf("frogs_spectra_by_locality.pdf")
explorespec(mspectra, by=8, ylim=c(0,60))
dev.off()
############SUBSTRATE DATA##########
# Now Import the substrate reflectances
data.panama <- read.csv("./PanamaSub_OK.csv")
data.costarica <- read.csv("./Substrate_CostaRica_OK.csv")


subpanama <- as.rspec(data.panama,
                      lim = c(300, 700))


subcosta <- as.rspec(data.costarica,
                     lim = c(300, 700))

all.substrates <- cbind(subcosta,subpanama[-1])
all.substrates <- procspec(all.substrates, 
                           fixneg = "zero", 
                           opt = 'smooth', 
                           span = 0.2)


head(all.substrates)
sum.substrates<-summary(all.substrates)
write.csv(sum.substrates,"./summary_stats_substrates.csv")


# localities and substrates list
subst.desc <- colnames(all.substrates[-1])
head(subst.desc)
locality.types <- gsub("\\.[0-9]+","",subst.desc)
table(locality.types)
loc.mean <- aggspec(all.substrates, 
                    by = locality.types, 
                    FUN = mean)
pdf("substrate_spectra_by_locality.pdf")					
explorespec(loc.mean, by=3)
dev.off()
GreenLeaf.means<-as.rspec(cbind(loc.mean[1],loc.mean[,grepl("GreenLeaf", colnames(loc.mean))]))
LeafLitter.means<-as.rspec(cbind(loc.mean[1],loc.mean[,grepl("LeafLitter", colnames(loc.mean))]))
Trunk.means<-as.rspec(cbind(loc.mean[1],loc.mean[,grepl("trunk", colnames(loc.mean))]))

pdf("./GreenLeaf.means.pdf")
explorespec(GreenLeaf.means,by=4) # SAM & DAM are not OK
dev.off()

pdf("./LeafLitter.means.pdf")
explorespec(LeafLitter.means,by=4) # is OK
dev.off()

pdf("./Trunk.means.pdf")
explorespec(Trunk.means,by=4) # ALM is not OK
dev.off()

# Hitoy did not fit the plots, doing it now:
HIT.means<-as.rspec(cbind(loc.mean[1],loc.mean[,grepl("HIT", colnames(loc.mean))]))
pdf("./HIT.means.pdf")
explorespec(HIT.means) # is OK
dev.off()

# It was impossible to get a good representation of all three substrate types by locality.
# We decided to average all the messuremnts of each substrate type across all localities
substrate.types <- gsub("^[A-z]+\\.","",locality.types)
table(substrate.types)
# double check correspondence:
cbind(names(all.substrates[-1]), substrate.types)

subst.mean <- aggspec(all.substrates, by = substrate.types, 
                      FUN = mean)
pdf("average_substrates_used.pdf")
plot(subst.mean)			
dev.off()

###########BIRD VISUAL MODEL AND CONTRASTS#######################
#Visual Modelling
# Bird visual model
#all frogs against all background
vismodel.data <- cbind(mspectra,subst.mean[-1])

vismod.bird <- vismodel(vismodel.data, 
                        visual = "avg.uv", 
                        achromatic = "bt.dc",
                        qcatch = "Qi",
                        illum = "forestshade",
                        relative = FALSE)
bird.contrasts <- coldist(vismod.bird, noise = "neural",
                          achro = TRUE, 
                          n = c (1, 2, 3, 3), 
                          #Hart et al.,2000c
                          weber = 0.1, 
                          weber.achro = 0.05)
#(Siddiqi et al. 2004 J Exp Biol)
###########ORGANIZE CONTRAST RESULTS FOR BIRD BY LOCALITY#############################################
frogs <- colnames(mspectra[-1])
bird.contrast.subst<-bird.contrasts[bird.contrasts$patch1 %in% frogs & 
                                      bird.contrasts$patch2 %in% unique(substrate.types),]
bird.contrast.subst$predator<-"bird"
#write.csv(bird.contrast.subst, "frog-subst.contrasts_bird.csv")
frog.contrast.bird<-bird.contrasts[bird.contrasts$patch1 %in% frogs & 
                                     bird.contrasts$patch2 %in% frogs,]
bird <- jnd2xyz(frog.contrast.bird,
                ref1 = "l",
                axis1 = c(0,1,1), axis2 = c(1,0,0),
                ref2 = NULL,
                center = TRUE,
                rotate = TRUE,
                rotcenter = c("mean"))

# Loading the metadata info from the sampling
metadata<-read.csv("frogs_metadata.csv")
morph <- metadata$morph
morph <- gsub('red', 'indianred2', morph)
morph <- gsub('green', 'forestgreen', morph)
species <- metadata$species
species <- as.factor(species)
jndplot(btit, 
        arrow = "relative", 
        arrow.p=0.8,
        achro = FALSE,
        arrow.labels = TRUE,
        labels.cex = 1.5,
        arrow.col = "gray27",
        pch = c(24, 21)[species],
        bg = morph,
        col = 'white',
        cex = 3.5, 
        cex.axis = 1.5,
        cex.lab = 1.5)
legend(x = "topleft",
       inset = c(-0.1, -0.04),
       box.col = "black",
       bg = "white",
       cex=1.5,
       text.font=3,
       title="Species", 
       legend=c("O. pumilio", "O. granulifera"),
       col = "black",
       pch = c(21, 24))
###########LIZARD VISUAL MODEL AND CONTRASTS#######################################################
# Lizard visual model - all frogs - all background
lizard <- sensmodel(c(370, 495, 550, 590), lambdacut = c(330, 371, 463, 507), oiltype = c("T", "C", "Y", "R"), om = TRUE, beta = 0.05) 
#Loew et al. 2002; Fleishman et al. 2011
vismod.lizard <- vismodel(vismodel.data, 
                          visual = lizard, 
                          achromatic = "ml",
                          qcatch = "Qi",
                          illum = "forestshade",
                          relative = FALSE)

lizard.contrasts <- coldist(vismod.lizard, noise = "neural",
                            achro = TRUE, 
                            n = c (1, 1, 1, 3), 
                            #Loew et al. 2002; Fleishman et al. 2020
                            weber = 0.1,  #(Fleishman, Perez et al., 2016).
                            # Loew et al. 2002; Fleishman et al. 2020
                            weber.achro = 0.1)

lizard.contrast.subst<-lizard.contrasts[lizard.contrasts$patch1 %in% frogs & 
                                          lizard.contrasts$patch2 %in% unique(substrate.types),]
lizard.contrast.subst$predator<-"lizard"
#write.csv(lizard.contrast.subst, "frog-subst.contrasts_lizard.csv")
frog.contrast.lizard<-lizard.contrasts[lizard.contrasts$patch1 %in% frogs & 
                                         lizard.contrasts$patch2 %in% frogs,]

###########FROGS IN LIZARD COLOUR SPACE MODEL#####################################################
#plot lizard - frogs against frogs
lizard <- jnd2xyz(frog.contrast.lizard,
                  ref1 = "lmax590",
                  ref2 = NULL, 
                  axis1 = c(0,1,1), 
                  axis2 = c(1,0,0),
                  center = TRUE, 
                  rotate = TRUE,
                  rotcenter = c("mean"))
jndplot(lizard, 
        arrow = "relative", 
        arrow.p=0.8,
        achro = FALSE,
        arrow.labels = TRUE,
        labels.cex = 1.5,
        arrow.col = "gray27", 
        pch = c(24, 21)[species],
        bg = morph,
        col = 'white',
        cex = 3.5, 
        cex.axis = 1.5,
        cex.lab = 1.5)
legend(x = "topleft",
       inset = c(0.04, 0),
       box.col = "black",
       bg = "white",
       cex=1.5,
       text.font=3,
       title="Species", 
       legend=c("O. pumilio", "O. granulifera"),
       col = "black",
       pch = c(21, 24))

###########CRAB VISUAL MODEL AND CONTRASTS##########################################
# Crab visual model - all frogs - all background
crab <- sensmodel(c(430, 590)) #. Horch et al. (2002)
vismod.crab <- vismodel(vismodel.data,
                        visual = crab, 
                        achromatic = "ml",
                        qcatch = "Qi",
                        illum = "forestshade",
                        relative = FALSE) 

crab.contrasts <- coldist(vismod.crab, noise = "neural",
                          achro = TRUE, 
                          n = c (1, 2),
                          #Horch et al.2002
                          weber = 0.12,
                          weber.achro = 0.1)
#Peitsch 1992; Hempel De Ibarra et al. 2000

crab.contrast.subst<-crab.contrasts[crab.contrasts$patch1 %in% frogs & 
                                      crab.contrasts$patch2 %in% unique(substrate.types),]
crab.contrast.subst$predator<-"crab"
#write.csv(crab.contrast.subst, "frog-subst.contrasts_crab.csv")
frog.contrast.crab<-crab.contrasts[crab.contrasts$patch1 %in% frogs & 
                                     crab.contrasts$patch2 %in% frogs,]

#write.csv(crab.contrast.sfr, "crabfrogs.csv")
###########CRAB colour space#######

crab <- jnd2xyz(frog.contrast.crab,
                ref1 = "lmax590",                  
                ref2 = NULL, 
                axis1 = c(0,1,1), 
                axis2 = c(1,0,0),
                center = TRUE, 
                rotate = TRUE,
                rotcenter = c("mean"))
jndplot(crab, arrow = "none",#"relative",
        achro = TRUE,
        arrow.labels = TRUE,
        arrow.col = "gray27",
        pch = c(24, 21)[species],
        bg = morph, 
        col="white",
        cex = 3,
        cex.axis = 2)
legend(x = "topleft",
       inset = c(0.01, 0.01),
       box.col = "black",
       bg = "white",
       cex=1.4,
       text.font=3,
       title="Species", 
       legend=c("O. pumilio", "O. granulifera"),
       col = "black",
       pch = c(21, 24))

###########COMBINE EVERYTHING TOGETHER#############################
contrast.subst.all <- rbind(bird.contrast.subst, lizard.contrast.subst, crab.contrast.subst)
names(contrast.subst.all)<-c("individual","substrate","dS","dL", "predator")
dd<-merge(contrast.subst.all,metadata, all.x=TRUE)
contrast.subst.all.metadata<-dd[order(dd$species, dd$morph, dd$locality, dd$predator, dd$substrate),]
write.csv(contrast.subst.all.metadata, 'visual_contrast_predators_on_substrates.csv', row.names=FALSE)