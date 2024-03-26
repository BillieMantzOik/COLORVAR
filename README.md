# COLORVAR
Color variance analysis between red and green frogs of Oophaga species.
Spectral data on dorsal side of frogs and substrate collected from frog localities.
To measure conspicuousness we implemented viewer dependent and taxon-specific
detection models for three potential viewers, including a UV-sensitive bird,
an Anolis lizard and a dichromatic crab. ). In all models, quantum catches at 
each photoreceptor were calculated from the spectra, and the conspicuousness of 
frogs was evaluated in terms of color (ΔS) and luminance (ΔL) contrast relative 
to an estimated average background substrate of three different substrate types, 
commonly found in the sampled localities (leaf litter, trunks, green leaves). 

We first implemented a tetrachromatic UV-sensitive visual system representing 
an average avian UV system using “avg.uv” in the vismodel command of pavo (Maia et al., 2019). 
The avian average UV system (“avg.uv” in pavo) was selected to simulate the 
visual system of a potential bird predator. For the second detection model,
we used a tetrachromatic lizard  visual system as described for Anolis cristatellus 
(Fleishman et al., 1997). The visual system  of this species is similar in physiology
and anatomy to those of other Anoline species studied 
(Fleishman et al., 1995; Fleishman et al., 1997; Persons et al., 1999; Loew et al., 2002). 
Anoline visual systems are adapted for high-acuity diurnal vision and the retina contains
four classes of cones modified by oil droplets (Fleishman et al., 1995). Spectral sensitivity 
for UVS of A. cristatellus peaks at 370 nm and SWS peaks at 495 nm. The peak sensitivity of 
the MWS is at 550 nm and for LWS at 590 nm (Fleishman et al., 2020). Oil droplet cut off was
estimated at 330 for UVS, 371 for SWS, 463 for MWS and 507 for LWS (Loew et al., 2002). 
Luminance receptor stimulation was estimated based on the longest wavelength photoreceptor 
(Loew et al., 2002; Fleishman et al., 2020). LWS photoreceptors densities is three times 
as high as the other photoreceptors and photoreceptor densities were set at 1:1:1:3 
(Fleishman et al., 2020). 

Finally, for the third detection model, we implemented a dichromatic visual system 
representing a crab predator, based on long-wavelength-sensitive (LWS) cone 
absorbance spectra and electrophysiological measures of a short-wavelength-sensitive 
(SWS) cone response for Uca thayeri (Horch et al. 2002).

Spectral analysis, Color perception analysis done using pavo (Maia et al., 2019)
Estmate of color and luminance ocntrasts between frogs and background were
estimated for three potential viewers. To investigate the variance between
color (dS) and luminance (dL) contrast, we tested the relationship between the
predictors and the estimated contrasts as responce variables. 

In separated analyses we fitted linear mixed effects models with individual 
as a random factor (to account for the replicated measurements per individual) 
and frog species, morph, substrate and predator (and their interactions) as
predictors of ΔS and ΔL variation. The models included individual as a random 
factor to account for the repeated measurements and were fitted with the 
function lmer from the lme4-package in R (R software team, 2023). 

With glmulti a multi-averaging model selction approach was used to select for
the best fitted models based on AICCc criteria.
