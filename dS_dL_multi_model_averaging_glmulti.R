# Script created for analyzing the visual contrasts of Oophaga frogs 
# contrasts were meassured against a leaf litter substrate
# Script by A. Rodríguez, data by V. Mantzana-Oikonomaki
setwd('./visualmodelling_Vasiliki')

library(lme4)
library(sjPlot)
library(MuMIn)

data<-read.csv("visual_contrast_predators_on_substrates_new_b.csv", stringsAsFactors=TRUE)
str(data)
head(data)

####### models #####################################################
# Modeling ds and dl from the effect of predator/ phenotype/ substrate with no interaction
#species fixed effect
# We first check for departures from normal distribution
hist(data$dS)
hist(data$dL)
# Both histograms show strong departures from normality
shapiro.test(data$dS)
shapiro.test(data$dL)
# Tests are highly significant
# As the data have structure we also perform the test by groups:
library(data.table)
DT <- data.table(data)
DT[,.(W = shapiro.test(dS)$statistic, P.value = shapiro.test(dS)$p.value),
   by = .(species, morph, predator)]
#	As per the previous table the normality assumption is not meet even after 
#	accounting for the structure in the data
#  We now perform the test on log-transformed data
shapiro.test(log(data$dS))
shapiro.test(log(data$dL))
DT[,.(W = shapiro.test(log(dS))$statistic, P.value = shapiro.test(log(dS))$p.value),
   by = .(species, morph, predator)]
DT[,.(W = shapiro.test(log(dL))$statistic, P.value = shapiro.test(log(dL))$p.value),
   by = .(species, morph, predator)]
# We now inspect the dato for outliers 
library(ggplot2)
p1<-ggplot(data, 
           aes(x = predator, y = dS, fill=substrate)) + 
  xlab("predator") +  
  geom_boxplot() + 
  theme_classic(base_size = 12) + 
  scale_fill_manual(values = alpha(c("darkolivegreen","tan3","sienna4"), .6)) +
  labs(fill = "substrate") + facet_wrap(~species + locality, ncol = 3) #scale="free", 
pdf("boxplots_dS_species.pdf", width = 11.7, height = 8.25)
p1
dev.off()

p2<-ggplot(data, 
           aes(x = predator, y = dL, fill=substrate)) + 
  xlab("predator") +  
  geom_boxplot() + 
  theme_classic(base_size = 12) + 
  scale_fill_manual(values = alpha(c("darkolivegreen","tan3","sienna4"), .6)) +
  labs(fill = "substrate") + facet_wrap(~species + locality, ncol = 3)
pdf("boxplots_dL_species.pdf", width = 11.7, height = 8.25)
p2
dev.off()

library(lme4)
mod1.1<- lmer(dS ~ predator + morph + substrate + species + 
                (1|individual), data=data)
mod1.1.std<- lmer(scale(dS) ~ predator + morph + substrate + species + 
                    (1|individual), data=data)
# Plot of model diagnostics	
plot_model(mod1.1, type="diag")
# Test for Overdispersion in the model
# We tested for overdispersion using parametric overdispersion tests implemented in the DHARMa package (https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html)
# This approach was followed after reading: https://royalsocietypublishing.org/doi/full/10.1098/rstb.2017.0281
library(DHARMa)
simulationOutput <- simulateResiduals(fittedModel = mod1.1)
plot(simulationOutput)
testDispersion(simulationOutput)
# Now the real model we need, with interactions
#interaction
#all interactions and species
mod1.2 <- lmer(dS ~predator*morph *substrate * species +
                 (1|individual), data=data, REML=FALSE) 
# control=lmerControl(optimizer="optimx", optCtrl=list(method='nlminb')))
plot_model(mod1.2, type="diag", grid=TRUE)
# Departure from normality of residuals does not seem so dramatic
library(glmulti)
library(lme4)
# Using wrapper function update from Ben Bolker : https://gist.github.com/bbolker/4ae3496c0ddf99ea2009a22b94aecbe5
## The wrapper function for linear mixed-models
mm_glmulti <- function(formula, data, random="", FUN = lme4::lmer,
                       extra_args = NULL, ...) {
  ## FIXME: check/throw error if REs not parenthesized?
  ff <- glmmTMB::addForm(formula, substitute(random))
  arglist <- c(list(ff, data = data), list(...), extra_args)
  do.call(FUN, arglist)
}
lmer.glmulti <- function(...) mm_glmulti(..., extra_args = list(REML = FALSE))
glmer.glmulti <- function(...) mm_glmulti(..., FUN = lme4::glmer)
glmmTMB.glmulti <- function(...) mm_glmulti(..., FUN = glmmTMB::glmmTMB,
                                            extra_args = list(REML = FALSE))
setMethod('getfit', 'merMod',
          function(object, ...) {
            summ <- coef(summary(object))
            summ1 <- summ[, c("Estimate", "Std. Error"), drop = FALSE]
            ## FIXME: update for lmerTest?
            cbind(summ1, df=rep(10000, length(fixef(object))))
          }
)

setMethod('getfit', 'glmmTMB',
          function(object, ...) {
            ## only handles conditional
            ## warn if disp, zi detected??
            summ <- coef(summary(object))$cond
            summ1 <- summ[, c("Estimate", "Std. Error"), drop = FALSE]
            cbind(summ1, df=rep(10000, length(fixef(object)$cond)))
          }
)

individual<- as.numeric(data$individual)
locality<-as.numeric(data$locality)
predator<-as.numeric(data$predator)
substrate<-as.numeric(data$substrate)
data.frame(individual=data$individual,id=individual,locality=locality,predator=predator,substrate=substrate)
mixed_model.dS <- glmulti(
  dS ~ predator* morph * substrate * species,
  random  = (1|individual),
  crit    = aicc,
  data    = data,
  method  = "h",
  fitfunc = lmer.glmulti,
  marginality = TRUE,
  level   = 2)
pdf("AICc_plots_dS_revised.pdf")
par(mfcol=c(2,1), mar=c(5,6,4,3))
#plot(mixed_model.dS, type = "p")
plot(mixed_model.dS, type = "s")
plot(mixed_model.dS, type = "w")
dev.off()
dS.coef<-coef(mixed_model.dS, select=0.95,
              varweighting="Johnson",
              icmethod="Burnham")
model.table.dS<-weightable(mixed_model.dS)
write.csv(model.table.dS, "glmulti_model_table_dS_new.csv")
write.csv(dS.coef, "model-averaged_estimates.95_dS_new.csv")
mixed_model.dS
# Plots using the best model
best.mixed_model.dS<-mixed_model.dS@objects[[1]]
par(mfcol=c(1,3))
theme_set(theme_sjplot(sjplot_pal(pal = "viridis")))
tab_model(best.mixed_model.dS)
plot_model(best.mixed_model.dS, type = "std", 
           show.values = TRUE, value.offset = .3, sort.est=TRUE)
plot_model(best.mixed_model.dS, type="int", 
           connect.lines = TRUE)
plot_model(best.mixed_model.dS, type="eff",  
           terms=c("predator", "substrate", "morph"))
mdrt.values = "meansd"

# Diagnostic plot to check for assumptions of linear models
plot_model(best.mixed_model.dS, type="diag")
# Interaction plots variant two
library(interactions)
cat_plot(best.mixed_model.dS, pred = substrate, modx = predator, plot.points = TRUE)
cat_plot(best.mixed_model.dS, pred = predator, 
         modx =morph, mod2=substrate, plot.points = TRUE, colors=c("seagreen", "orangered4"), alpha=0.6)
# Statistics for brightness contrasts
mixed_model.dL <- glmulti(
  dL ~ predator* morph * substrate * species,
  random  = (1|individual),
  crit    = aicc,
  data    = data,
  method  = "h",
  fitfunc = lmer.glmulti,
  marginality = TRUE,
  intercept = TRUE,
  level   = 2)
pdf("AICc_plots_dL_new.pdf")
par(mfcol=c(2,1), mar=c(5,6,4,3))
#plot(mixed_model.dL, type = "p")
plot(mixed_model.dL, type = "s")
plot(mixed_model.dL, type = "w")
dev.off()
dL.coef<-coef(mixed_model.dL, select=0.95,
              varweighting="Johnson",
              icmethod="Burnham")
model.table.dL<-weightable(mixed_model.dL)
write.csv(model.table.dL, "glmulti_model_table_dL_new.csv")
# model averaged estimates using models that sum up to 95% evidence weight
write.csv(dL.coef, "model-averaged_estimates.95_dL_new.csv")

# Plots using the best model
best.mixed_model.dL<-mixed_model.dL@objects[[1]]
par(mfcol=c(1,3))
theme_set(theme_sjplot())
tab_model(best.mixed_model.dL)
plot_model(best.mixed_model.dL, type = "std", 
           show.values = TRUE, value.offset = .3, sort.est=TRUE)
plot_model(best.mixed_model.dL, type="int")
plot_model(best.mixed_model.dL, type="pred",  
           terms=c("predator", "substrate", "morph"))
mdrt.values = "meansd"
# Diagnostic plot to check for assumptions of linear models
plot_model(best.mixed_model.dL, type="diag")

# Plot dS and dl
p1<-ggplot(data[data$species=="granulifera",], aes(x = dL, y = dS, color=morph)) + 
  xlab("dL") +  geom_point(size = 3) + 
  theme_classic(base_size = 12) %+replace% 
