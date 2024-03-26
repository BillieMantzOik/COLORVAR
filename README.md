# COLORVAR

contents:
visual_contrasts_calculation_pavo.R --> script for spectral analysis, color perception analysis and calculating visual contrasts between frogs and background for three different viewers.
dS_dL_multi_model_averaging_glmulti.R --> script for statistical analysis using linear models and model averaging for investigating the relationship between predictors and responce variables.
data --> subfolder containing input and output data


Color variance analysis between red and green frogs of Oophaga species.

Spectral data collected on the dorsal side of frogs collected from different localities
and spectral data on substrate collected from these localities were used.

To measure conspicuousness we implemented viewer dependent and taxon-specific
detection models for three potential viewers, including a UV-sensitive bird,
an Anolis lizard and a dichromatic crab, using pavo package (Maia et al., 2019).

In all models, quantum catches at each photoreceptor were calculated from the 
spectra, and the conspicuousness of  frogs was evaluated in terms of color (ΔS) 
and luminance (ΔL) contrast relative to an estimated average background substrate 
of three different substrate types, commonly found in the sampled localities 
(leaf litter, trunks, green leaves). 

We first implemented a tetrachromatic UV-sensitive visual system representing 
an average avian UV system using “avg.uv” in the vismodel command of pavo. 
The avian average UV system (“avg.uv” in pavo) was selected to simulate the 
visual system of a potential bird predator. For the second detection model,
we used a tetrachromatic lizard  visual system as described for Anolis cristatellus 
(Fleishman et al., 1997). 

Finally, for the third detection model, we implemented a dichromatic visual system 
representing a crab predator, based on long-wavelength-sensitive (LWS) cone 
absorbance spectra and electrophysiological measures of a short-wavelength-sensitive 
(SWS) cone response for Uca thayeri (Horch et al. 2002).

To investigate the variance of spectral contrasts between frogs, we tested the 
relationship between predictors and the estimated color and luminance contrasts 
as responce variables. 

In separated analyses we fitted linear mixed effects models with individual 
as a random factor (to account for the replicated measurements per individual) 
and frog species, morph, substrate and predator (and their interactions) as
predictors of ΔS and ΔL variation. The models were fitted with the 
function lmer from the lme4-package in R (R software team, 2023). 

With glmulti a multi-averaging model selction approach was used to select for
the best fitted models based on AICCc criteria.
