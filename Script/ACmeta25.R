# ACmeta25.R
# Meta-analysis script for testing consensus of Arctic cod population
# structure from the studies that have attempted to estimate it over the years.
# The data are compiled into a single dataset except for Fevolden et al. 1999
# and Palsson et al. 2009, which do not report all Fst values specifically.
# For these latter two studies, the pairwise Fst values are estimated using the 
# values they do report.

# The mapping part of the script uses ggOceanmap to create a visual 
# rendering of the where studies have been conducted. The initial map image 
# look pretty crowded on the image viewing window, but it ends up looking 
# much better as a jpeg saved at 1500 X 800 px from the "save as Image" 
# option in the plots window.

# Most of the analyses are performed using the metafor package. The graphs and figures
# are generated using both metafor and, in most cases, ggplot2. 

# written by DR and ML June 2025

############################# 1. PRELIMINAIRIES ###########################
# clearing instances and resetting R
rm(list=ls())

# Installing packages that might be needed (only run once)

#install.packages("dplyr")
#install.packages("vctrs")

# loading all libraries

{
  library(dplyr)
  library(metafor)
  library(vctrs)
  library(tidyr)
  library(tidyverse)
  library(ggplot2)
  library(emmeans)
  library(stringr)
  library(ggplot2)
  library(ggnewscale)
  library(circlize)
  library(scales)
}

# Loading libraries for Mapping
{
  library(ggOceanMaps)
  library(ggspatial)    # for data plotting
}

# set working directory
setwd("path-to-file-on-computer")

##################### 2. PALSSON et al, 2009 ###########################
# Dealing with Palsson study which uses mtDNA to estimate genetic differences,
# but only reports 1 or 2 out 136 pairwise comparisons. All others "vary around
# 0, with a few (one in particular) being significant after correcting for multiple
# tests. 

# Set sample size
ss <- 136

# Use beta distribution with left skew to focus Fsts around 0.000
# with a maximal value near 0.171 as reported in the paper.
# in rbeta shape 1 controls accumulation near 1 larger -> more near 1
# shape2 controls accumulation near 0 larger -> more near 0 
paldat <- rbeta(ss, shape1 = 0.2, shape2 = 5)
scpaldat <- paldat * 0.171

# Set the current plot area parameters to 'op' 
op <- par(cex = 1.0, font = 2, mgp = c(3, 1, 0))

# Plot the histogram (quickly!) to see the skew and the overall data
par(mar=c(8,8,2,1), mgp = c(6, 1, 0))
hist(scpaldat, col = "cadetblue", cex.lab = 2, cex.axis = 2, xlab = "FST", main = "")

par(op)

## She also specified that 12 vs 17 has fst = 0.171, corresponding to line 119
scpaldat[119] <- 0.171
mean(scpaldat)

# Data created for Palsson will be inserted in real data where currently
# values are placeholders. Real and only reported comparison value between 
# site 12(NE-Greenland) and 17(Spitsbergen) = 0.171 is also inserted in data. 
# scpaldat will replace placeholder data.


###################### 3. FEVOLDEN et al, 1999 ###########################
# Need to manipulate Fevolden data who report only AMOVA. We know the total
# variance, among-locations variance, and thus the global 
# Fst. To approximate the 3 pairwise Fst values (since there are 3 sites = 3 
# pairwise), we randomly partition among-locations variance into 3 parts 
# and compute Fst_ij = var_ij / var_total

# Data provided by Fevolden et al. 1999
varam <- 0.086
varwi <- 4.952
varto <- 5.039

# Mean can be recovered as among/overall variance
# over 3 populations
mu <- varam / varto
npairs <- 3

fevfst <- function(k) {
  x <- rgamma(npairs, shape = k, scale = 1)
  fstfev <- mu * npairs * x / sum(x)
  return(fstfev)
}
fevfst(3)


################# 4. READ AND MANIPULATE DATA ##################
# Import data - looking for metaACdata.csv
mcoord <- read.csv("metaACdata.csv", header = T, stringsAsFactors = T)

# Confirming number of studies listed
length(unique(mcoord$study))
nrow(mcoord)

# The data were entered as text?! so convert to standardised lat longs 
# with 5 significant digits
mcoord[c("lat1", "lon1", "lat2", "lon2", "Fst")] <-
  lapply(mcoord[c("lat1", "lon1", "lat2", "lon2", "Fst")],
         function(x) round(as.numeric(x), 5))

# Set empty data frame to store data (trying to retrieve
# locations of each sampling site from different studies)
mapcoord <- data.frame()

# Loop runs through studies and finds unique lats 
# longs to extract individual sampling points for each.
for (curstu in unique(mcoord$study)) {
  
  # Because Wilson has two markers types but uses same sampling
  # sites - we only need to use those once
  if (curstu == "Wilson et al., 2019b") curstu <- "Wilson et al., 2019a"
  
  # "test" object retrieves current study's mcoord and returns 
  # the ones that are not repeated
  test <- mcoord[mcoord$study==curstu,]
  unilat1 <- which(!duplicated(test[, c("lat1", "lon1")]))
  
  # fills a tmp df that has study, site, and lats/lons
  tmp1 <- cbind.data.frame(study = test$study[unilat1],site = test$site1[unilat1],
                           lat = test$lat1[unilat1], lon = test$lon1[unilat1])
  
  # Retrieves unique locations for lat2, and puts study 
  # unique location in tmp2
  unilat2 <- which(!duplicated(test[, c("lat2", "lon2")]))
  tmp2 <- cbind.data.frame(study = test$study[unilat2], site = test$site2[unilat2],
                           lat = test$lat2[unilat2], lon = test$lon2[unilat2])
  
  # Next, combine tmp dfs and look for duplicates.
  # End with list of unique locations for each study in an
  # object mapcoord
  tmpcomb <- rbind(tmp1, tmp2)
  tmpcomb <- tmpcomb[!duplicated(tmpcomb[, c("lat", "lon")]), ]
  
  mapcoord <- rbind(mapcoord, tmpcomb)
}

######################### 5. R-MAP #############################
# Save above mapcoord and look at them on a plotted map to see what 
# coordinates of the larger 7 groups are. Make a file I can also read
# into R called grpcoord
ggrps <- read.csv("grpcoord.csv", header = T, stringsAsFactors = T)

# Assign label data.frame containing information for major land 
# masses
landlabs <- data.frame(
  Long = c(-111, 50, -2),          # adjust X positions
  Lat  = c(59, 59, 50),            # adjust Y positions
  name = c("Canada", "Russia", "Europe")  # text on the map
)

# Since Greenland is smaller - needs its own
greenl <- data.frame(
  Long = c(-44),          # adjust X positions
  Lat  = c(76),           # adjust Y positions
  name = c("Greenland")   # text on the map
)

# Finally, also need to name groups we'll use in regional analyses
regnam <- data.frame(
  Long = c(-149, -90, -48, -27, 19, 65, 161),     # adjust X positions
  Lat  = c(65, 62, 53, 61, 72, 77, 81 ),          # adjust Y positions
  name = c("Western\n Arctic", "Baffin Bay \n West Greenland",
           "Northwest\n Atlantic", "East Greenland\n Iceland",
           "Svalbard", "Novaya\n Zemlya", "Eastern\n Arctic")  # text on the map
)

# Assign the studies a specific colour here:
grpcol <- c(
  "Bringloe et al. 2024" = "dodgerblue",
  "Emilianova et al. 2023" = "palegreen1",
  "Fevolden  et al. 1999" = "hotpink",
  "Goordeva and Mishin 2019" = "red",
  "Madsen et al. 2016" = "firebrick4",
  "Maes et al. 2021" = "palegoldenrod",
  "Maes et al. 2025" = "darkorchid",
  "Nelson et al. 2020" = "cyan",
  "Palsson et al. 2009" = "darkorange",
  "Quintela et al. 2021" = "gold",
  "Wilson et al. 2019a" = "darkolivegreen",
  "Wilson et al. 2019b" = "darkolivegreen"
)

# Due to some limitations of ggOceanMaps package, the Goordeva and Mishin
# coordinates get plotted as Nelson. So need to replot at end of ggplot to
# make them appear.
GMfix <- filter(mapcoord, study == "Goordeva and Mishin 2019")


# Make the basemap as per direction in ggOceanMaps outlining circumpolar 
# North. Add sample sites and labels, and major land features. Also 
# include North pointing Arrow, scale, and legend for the depth profile
basemap(limits = c(-165, 65, 48, 80), shapefiles = "Arctic", rotate = T, 
        bathymetry = T, bathy.style = "rbb", land.col = "darkgray")+
  ggnewscale::new_scale_fill() +
  geom_spatial_polygon(data = ggrps, aes(x = lon, y = lat, group = area),
                       fill = "black", color = "#ff4d4d", alpha = 0.2, linewidth = 1.0, 
                       linetype = "solid")+
  geom_spatial_point(data = mapcoord, aes(x = lon, y = lat, colour = study), size = 5, stroke = 1, alpha = 0.7)+
  scale_colour_manual(values = grpcol) +
  geom_spatial_text(data = landlabs, aes(x = Long, y = Lat, label = name),
                    size = 18, fontface = "italic", colour = "black") +
  geom_spatial_text(data = greenl, aes(x = Long, y = Lat, label = name, angle = -35),
                    size = 11, fontface = "italic", colour = "black") +
  geom_spatial_text(data = regnam, aes(x = Long, y = Lat, label = name),
                    size = 9, fontface = "italic", colour = "black") +
  labs(fill = "Depth (m)", x = "Longitude (° DD)", y = "Latitude (° DD)") +
  annotation_scale(location = "tr", text_cex = 1.9, text_col = "ghostwhite") + 
  annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0.8, "cm"), 
                         pad_y = unit(0.8, "cm"), style = north_arrow_fancy_orienteering,
                         height = unit(2.5, "cm"), width = unit(2.5, "cm"))+
  theme(axis.text.x = element_text(size = 17),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        axis.text.y = element_text(size = 17),
        axis.title = element_text(face = "bold"),
        axis.line = element_line(linewidth=1.5),
        axis.ticks = element_line(linewidth = 1.5),
        axis.ticks.length = unit(.3,"cm")) +
  theme(legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.background = element_blank())+
  geom_spatial_point(data = GMfix, aes(x = lon, y = lat), fill = grpcol["Goordeva and Mishin 2019"], 
                     colour = "red", size = 5, stroke = 1, alpha = 0.7, inherit.aes = F)
  

####################### 6. LOOK AT THE DATA #####################
# look at the actual data especially the FST values
head(mcoord)

# Insert simulated Fst for Fevolden and Palsson into mcoord dataframe
mcoord$Fst[which(mcoord$study=="Fevolden  et al., 1999")] <- fevfst(3)
mcoord$Fst[which(mcoord$study=="Palsson et al., 2009")] <- scpaldat

# For each FST extracted estimate variance on that value using formula:
# Vfst <- FST * (1-FST)/(N1 +N2) from Quinn and Keough 2002; 
mcoord$V <- mcoord$Fst*(1-mcoord$Fst)/(mcoord$n1+mcoord$n2)

# Above can generate - FST variances which are not usable
# in the remaining analyses (e.g., -ve values). So need to
# set a minimum value of the Fst variance to always be +ve.
minV <- 0.000001
mcoord$V[mcoord$V < minV] <- minV

# Calculate sample sizes used to estimate individual pairwise 
# Fsts in each study
mcoord$np <- mcoord$n1 + mcoord$n2

# verify that Fsts and nps are numeric
mcoord$Fst <- round(as.numeric(mcoord$Fst), digits = 5)
mcoord$V <- round(as.numeric(mcoord$V), digits = 6)
mcoord$np  <- as.numeric(mcoord$np)

# R doesn't like reading the "µsat" as listed in csv file. So, it 
# is converted to character - fixed - then returned to a factor.
mcoord$markers <- as.character(mcoord$markers)
mcoord$markers[which(as.character(mcoord$markers) == "\xb5sat")] <- "µsat"
mcoord$markers <- as.factor(mcoord$markers)

# Calculate the summary stats for the Fsts by study
# using the dplyr functions
sumstat <- mcoord %>%
  group_by(study) %>%
  summarise(
    mfst = mean(Fst),
    minfst = min(Fst),
    maxfst = max(Fst),
    vfst = var(Fst),
    sefst = sd(Fst)/sqrt(n()),
    pairwise = n(),
    m_np = round(mean(np)/2, digits = 0),
    marker = first(markers)
  )

# Calculate the weights or inverse standard errors for each study to 
# outline weight it has in impacting overall effect size (Fst).
# weights <- 1 / sqrt(sumstat$vfst)   
# sumstat$weights <- weights / max(weights)   # normalize so max = 1

####################### 7. INITIAL META ANALYSES ##################
# We first run random-effects meta-analysis with yi = mfst and sei equal to
# standard error of study-level mean FST, using restricted maximum
# likelihood (REML) to estimate between-study heterogeneity. This model includes
# random study effect representing deviations of each study’s true FST from
# the overall mean (i.e., ui in the model yi = μ + ui + εi, where μ is the
# consensus FST, ui is the study-specific random effect, and εi is sampling
# error). 

# Differences among studies in the number of pairwise population comparisons 
# influences precision of study-level mean estimates and are therefore 
# reflected in the reported standard errors and inverse-variance
# weighting, rather than being modeled as a separate random factor.
metares1 <- rma.uni(yi = mfst, sei = sefst, method = "REML", data = sumstat) 

# Look at metares1 results 
summary(metares1)

# Model Results are straight forward with consensus FST ~ 0.0062 and etc.
# indicating statistically significant but very weak population structure 
# on average.

# Heterogeneity: Q test shows clearly statistically significant 
# overall population structure across studies (consensus FST = 0.0062).
# I² = 97.5%; substantial heterogeneity present among studies, with most 
# variation in reported FSTs reflecting real biological differences rather 
# than sampling error.

# While estimated between-study variance (τ² = 0.0001) is small in
# absolute terms due to bounded and near-zero scale of FST, study-level
# estimates differ considerably in sampling precision. As a
# result, majority of observed variability among studies attributable
# to true heterogeneity rather than sampling error (I² = 97.5%).

# Metafor uses strange settings to organise forest plots so setting 
# it up is useful. Here setting margin limits
par(mgp = c(3.5, 1, 0))

# Visualise the metares1 result with a forest plot 
forest(metares1, slab = sumstat$study, header = F, addpred = F, showweights = T, 
       xlim=c(-0.18, 0.10), shade=T, at=seq(-0.04, 0.04, by = 0.005),
       digits = 3, ilab = cbind(as.character(sumstat$marker), sumstat$m_np, 
       sumstat$pairwise), ilab.xpos = c(-0.11, -0.08, -0.05), cex = 1.4, 
       xlab = c(expression(italic(F[ST]))), cex.lab = 1.8, refline = c(0, metares1$b),
       colout = "black", colshade = "gray87", mlab = expression(Consensus ~ italic(F)[ST]),
       col = "dodgerblue4")

# Forest plots in themselves aren't that great unless you annotate the 
# produced figure. The above creates the base and below cmds to annotate
# the figure better. best resolution here is also 1350 width 818 height.
op <- par(cex = 1.0, font = 2, mgp = c(3, 1, 0))

text(0.074, length(sumstat$study) + 2, expression("Mean " * italic(F)[ST] * " \u00B1 95% CI"), 
     cex = 1.4, font = 2)

text(0.035, length(sumstat$study) + 2, "Weight", cex = 1.4, font = 2)
text(-0.164, length(sumstat$study) + 2, "Study", cex = 1.6)
text(-0.11, length(sumstat$study) + 2, "Marker", cex = 1.4)
text(-0.08, length(sumstat$study) + 2, "Mean N\nfor pairs", cex = 1.4)
text(-0.05, length(sumstat$study) + 2, "# Comp", cex = 1.4)
par(op)

# Adding the lines for the overall effect size and its 95% CI.
abline(v = as.numeric(metares1$b), col = "dodgerblue4", lty = 1, lwd = 2)
abline(v = as.numeric(metares1$ci.lb), col = "dodgerblue4", lty = 2, lwd = 1)
abline(v = as.numeric(metares1$ci.ub), col = "dodgerblue4", lty = 2, lwd = 1)

# Looking at Funnel plot for data to identify missingness of studies. 
# Large gaps in the efforts to estimate the population structure using Fst.

# Calculate inverse-variance weights from the random-effects 
# meta-analysis, reflecting each study’s contribution to the 
# pooled FST estimate.
weights1 <- weights(metares1)   
nweights1 <- weights1 / sum(weights1)   # normalise so max = 1

# Scale the point sizes (ps) so points reflect study weight in 
# meta-analysis or how much it influences the consensus FST 
# (inverse of the sampling variance + random effect if included)
ps <- 40 * nweights1

# Save the current plotting parameters as pp.
pp <- par(no.readonly = T)   # save current settings

# Change the plotting parameters to make the funnel plots cleaner
# and easier to see.
par(mar=c(8,8,2,1), mgp = c(6, 1, 0))

# Funnel plot for the first model metares1.
funnel(metares1, yaxis = "sei", cex.axis = 1.6, xlim = c(-0.03, 0.04), 
       ylim = c(0, 0.01), steps = 10, las = 1, xlab = "", cex.lab = 2.0, 
       ylab = "Standard Error", addtau2 = F, back = "gray87",
       pch = 19, col = adjustcolor("dodgerblue4", alpha.f = 0.7), cex = ps)

mtext(c(expression(italic(F[ST]))), side = 1, line = 3, cex = 2)
abline(v = 0, col = "black", lty = 3, lwd = 3)
abline(v = as.numeric(metares1$b), col = "dodgerblue4", lty = 1, lwd = 2)
abline(v = as.numeric(metares1$ci.lb), col = "dodgerblue4", lty = 2, lwd = 2)
abline(v = as.numeric(metares1$ci.ub), col = "dodgerblue4", lty = 2, lwd = 2)
par(pp)

# Metares1 funnel plot shows that on average, Fst calculated from studies 
# included in model are > 0 (black dotted line), but not by much ~ 0.0062. 
# The blue dots (mean Fst/study) should be "close" to overall consensus (blue
# line), but are scattered. So, studies show variation around consensus.
# Most studies have low SE and so are relatively precise (lower SE values and 
# higher on the plot). Size of the each study's point reflects the study
# weight (related to their SE and # of comparisons).


# Use metafor to estimate the number of studies missing, or that would be 
# required to estimate an even smaller overall FST (0.001 here), that is 
# essentially 0.
fsn(x = metares1, target = 0.001, type = "General")

# So, there would have to be 60 + studies with FSTs lower than 0.0062 to 
# nullify the significant FST we've observed here.


####################### 8. BY MARKER TYPE ##################
# Previous versions of analyses tried to re-run analyses by marker type
# and included a summary step aggregating data by marker rather than
# by study. This approach incorrect for meta-analysis because it
# collapses effect sizes across studies, removing between-study
# variability and obscuring important differences in study design,
# sampling, and precision.

# To ensure studies (or individual effect sizes) remain the unit of
# analysis, we instead analyse the full dataset and include marker type
# as moderator. This preserves original study-level information while
# formally testing whether marker type explains heterogeneity among
# FST estimates, as recommended in metafor package documentation.

# Run analyses using marker as moderator to see if it helps 
# explain heterogeneity in the included studies FST indices.
metares2 <- rma.mv(yi = Fst, V = mcoord$V, mods = ~ markers, random = ~ 1|study,
                    data = mcoord, method = "REML")

# Summarise model
summary(metares2)

# Assess results of marker specifically (moderator effects) or
# Omnibus test.
anova(metares2)

# rma.mv is not supported by emmeans <-  so all calculations
# that would be provided by emmeans must be done by hand.

# Extracting regression coefficients from model for each 
# of the moderator levels
regco <- coef(metares2)

# Also extract the variance–covariance matrix for the 
# different moderator values
mr2vcv <- vcov(metares2)

# Make L matrix specifying which factor levels to compare
# for example the first characterises mtDNA, and the next one 
# looks at RAPD compared with mtDNA, etc,.. Sort of like 
# model.matrix cmd in base R.
L <- rbind(
  mtDNA  = c(1, 0, 0, 0),
  RAPD  = c(1, 1, 0, 0),
  SNP   = c(1, 0, 1, 0),
  µsat  = c(1, 0, 0, 1)
)

# Setting column names for L matrix to those of the extracted
# coefficients
colnames(L) <- names(regco)

# Estimate predicted mean FST per marker type by multiplying L 
# matrix with that of the regco (regression coefficients) -
# essentially using mtDNA (intercept) to convert other values.
meanFST <- as.numeric(L %*% regco)

# Calculating standard error (FULL covariance)
seFST <- sqrt(diag(L %*% mr2vcv %*% t(L)))

# Normal value for 95% CI (97.5 % quantile)
z <- 1.96 

# Calculating 95% CI for meanFST based on seFST from the 
# model
ci.lb <- meanFST - z * seFST
ci.ub <- meanFST + z * seFST

# Gather marker data into a small data.frame for plotting
# using ggplot
markest <- data.frame(marker = rownames(L), n = as.numeric(table(mcoord$markers)),
                      meanFST = meanFST, se = seFST, ci.lb = ci.lb, ci.ub = ci.ub)

# Set colours for the modeled FSTs by marker type
markcol <- c("mtDNA" = "firebrick4", "RAPD" = "springgreen4", "SNP" = "darkcyan", "µsat" = "darkorange")


# Plot using ggplot 
ggplot(markest, aes(x = marker, y = meanFST)) +
  geom_jitter(data = mcoord, aes(x = markers, y = Fst, colour = markers), width = 0.15,
    alpha = 0.3, size = 5, stroke = 1, show.legend = F) +
  geom_point(size = 5, colour = 'black', alpha = 0.7, stroke = 1) +
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), 
                width = 0.1, colour = "black", linewidth = 1.5) +
  geom_text(aes(y = 0.14, label = paste0("italic(n) == ", n)), parse = T, ,size = 10) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black", linewidth = 0.7) +
  scale_colour_manual(values = markcol)+
  scale_y_continuous(limits = c(-0.02, 0.18), breaks = seq(-0.02, 0.18, by = 0.04)) +
  labs(x = "Marker type", y = expression("Estimated mean " * italic(F)[ST])) +
  theme_classic() +
  theme(
    panel.grid.major = element_line(color = "grey96"),
    axis.title = element_text(size = 30, colour = "black"),
    axis.text.x = element_text(size = 26), 
    axis.text.y = element_text(size = 26),
    axis.line = element_line(linewidth = 2),
    axis.ticks = element_line(linewidth = 2),
    axis.ticks.length = unit(.4, "cm"))


####################### 9. WITHIN vs. AMONG/BETWEEN ##################
# Similar to the above, we want to test whether variation in FSTs in 
# different studies (so, again taking study into account) is higher when 
# calculated within vs. among regions. Our hypothesis is that overall 
# FST heterogeneity would be better explained by among region than within. 

# Follow similar analysis as for markers (section 8). To ensure studies 
# (or individual effect sizes) remain unit of analysis, we analyse 
# full dataset and include comptype as moderator. This preserves original 
# study-level information while formally testing whether comptype type 
# explains heterogeneity among FST estimates, as recommended in metafor.

# Sorting the comptypes
mcoord$comptype <- factor(mcoord$comptype,
                          levels = c("within", "among"))

# Run analyses using comptype type as moderator to see if it helps 
# explain heterogeneity in the included studies FST indices.
metares3 <- rma.mv(yi = Fst, V = mcoord$V, mods = ~ comptype, random = ~ 1|study,
                    data = mcoord, method = "REML")

# Summarise model
summary(metares3)

# Create a small df with the two levels of interest that we can use to dump
# information into.
sortdat <- data.frame(comptype = factor(c("within", "among"), 
                                         levels = levels(mcoord$comptype)))

# Get the model estimated means and CI for the different comptypes
compred <- predict(metares3, newmods = model.matrix(~ comptype, sortdat)[,-1])

# Set names of predicted values to comptype as re-sorted above to make within 
# appear before among
compred$slab <- as.character(sortdat$comptype)

# Make results into a dataframe easier to plot in ggplot.
compest <- cbind.data.frame(comptype = compred$slab, meanFSTC = compred$pred,
                            seFSTC = compred$se, ciC.lb = compred$ci.lb,
                            ciC.ub = compred$ci.ub, n = as.numeric(table(mcoord$comptype)))

# Make sure comptype appear as levels in the factor
compest$comptype <- factor(compest$comptype, levels = c("within", "among"))

# Set colours for the modeled FSTs by marker type
compcol <- c("within" = "peachpuff", "among" = "powderblue")

# Plot using ggplot 
ggplot(compest, aes(x = comptype, y = meanFSTC)) +
  geom_jitter(data = mcoord, aes(x = comptype, y = Fst, colour = comptype), width = 0.15,
              alpha = 0.6, size = 5, stroke = 1, show.legend = F) +
  geom_point(size = 5, colour = 'black', alpha = 0.7, stroke = 1) +
  geom_errorbar(aes(ymin = ciC.lb, ymax = ciC.ub), 
                width = 0.1, colour = 'black', linewidth = 1.5) +
  annotate("text", x = "within", y = 0.18, label = "X", color = "black", size = 12, family = "Courier") +
  annotate("text", x = "among", y = 0.18, label = "Y", color = "black", size = 12, family = "Courier") +
  geom_text(aes(y = 0.14, label = paste0("italic(n) == ", n)), parse = T, ,size = 10) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black", linewidth = 0.5) +
  scale_colour_manual(values = compcol)+
  scale_y_continuous(limits = c(-0.02, 0.18), breaks = seq(-0.02, 0.18, by = 0.04)) +
  labs(x = "Comparisons", y = expression("Estimated mean " * italic(F)[ST])) +
  theme_classic() +
  theme(
    panel.grid.major = element_line(color = "grey96"),
    axis.title = element_text(size = 30, colour = "black"),
    axis.text.x = element_text(size = 26), 
    axis.text.y = element_text(size = 26),
    axis.line = element_line(linewidth = 2),
    axis.ticks = element_line(linewidth = 2),
    axis.ticks.length = unit(.4, "cm"))


####################### 10. WITHIN GRPS ##################
# Here, interested in asking in which region is population Structure 
# highest, taking into account the different studies and variation 
# within them. 

# Can do similar to marker type analyses above. 

# Filter to consider only the within groups comparisons 
# getting rid of among groups.  
windat <- mcoord[!is.na(mcoord$wingrp),]

# Look at the regions listed in the within group data. 
reglist <- unique(windat$wingrp)

# Set an empty variable to eventually store regional code 
# data as a factor.
regcod <- as.character()

# There is a typo in data with 2 codes for same region "Nortwest_Atlantic" 
# and "Northwest Atlantic" - R considers these separate.
# To fix collapsed into single region - "NWAT".
idx <- which(windat$wingrp == "Nortwest_Atlantic" | windat$wingrp == "Northwest Atlantic")
regcod[idx] <- "NWAT"

# Since the above were levels one and two in the reglist, we can remove these.
reglist <- reglist[-c(1,2)]

# Set a new variable with new regional codes to use. These are ordered 
# to the levels in reglist and so if we change reglist - we would need to 
# also change rc too.
rc <- c("WARC", "NOZE", "EARC", "EGRE", "SVAL", "BBWG")

# Short loop populating the regcod variable according to 
# the reglist levels.
for (a in seq_along(reglist)) {
    idx <- which(windat$wingrp == as.character(reglist[a]))
    regcod[idx] <- rc[a]
}

# Once finished can set the regcod as a new variable in windat data.
windat$regcode <- regcod

# Analyses of windat data using regcode as moderator but also taking study as 
# random factor to see if it helps explain heterogeneity in included studies FSTs.
metares4 <- rma.mv(yi = Fst, V = windat$V, mods = ~ regcode, random = ~ 1|study,
                   data = windat, method = "REML")

# Summarise model
summary(metares4)

# Assess results of region specifically (moderator effects) or
# Omnibus test.
anova(metares4)

# rma.mv is not supported by emmeans <-  so all calculations
# that would be provided by emmeans must be done by hand.

# Extract regression coefficients from model for each moderator level
regco4 <- coef(metares4)

# Also extract variance–covariance matrix for moderator levels
mr4vcv <- vcov(metares4)

# Making an LR matrix specifying which factor levels to compare;
# for example the first characterises BBWG, and the next one 
# looks at EARC compared with BBWG, etc,.. Sort of like 
# model.matrix cmd in base R.
LR <- rbind(
  BBWG = c(1, 0, 0, 0, 0, 0, 0),
  EARC = c(1, 1, 0, 0, 0, 0, 0),
  EGRE = c(1, 0, 1, 0, 0, 0, 0),
  NOZE = c(1, 0, 0, 1, 0, 0, 0),
  NWAT = c(1, 0, 0, 0, 1, 0, 0), 
  SVAL = c(1, 0, 0, 0, 0, 1, 0),
  WARC = c(1, 0, 0, 0, 0, 0, 1)
)

# Setting column names for LR matrix to those of the extracted
# coefficients
colnames(LR) <- names(regco4)

# Estimate predicted mean FST per region by multiplying LR 
# matrix with that of the regco (regression coefficients) -
# essentially using BBWG (intercept) to convert other values.
mFSTR <- as.numeric(LR %*% regco4)

# Calculating the standard error (FULL covariance)
sFSTR <- sqrt(diag(LR %*% mr4vcv %*% t(LR)))

# Calculating 95% CI for mFSTR based on sFSTR from the 
# model
ci.lbR <- mFSTR - z * sFSTR
ci.ubR <- mFSTR + z * sFSTR

# Gather region data into a small dataframe for plotting
# using ggplot
regest <- data.frame(region = rownames(LR), n = as.numeric(table(windat$regcode)),
                      mFSTR = mFSTR, sFSTR = sFSTR, ci.lbR = ci.lbR, ci.ubR = ci.ubR)

# Sort the data from west to east
regest$region <- factor(regest$region,
                        levels = c("EARC", "NOZE", "SVAL", "EGRE", "NWAT", "BBWG", "WARC"))

# Use regcol to set plotting colours
regcol <- c("BBWG" = "skyblue3", "EARC" = "goldenrod4", "EGRE" = "springgreen4", "NOZE" = "darkorange",
             "NWAT" = "navy", "SVAL" = "salmon", "WARC" = "mediumpurple")

# Plot using ggplot 
ggplot(regest, aes(x = region, y = mFSTR, colour = region)) +
  geom_jitter(data = windat, aes(x = regcode, y = Fst, colour = regcode), width = 0.15,
              alpha = 0.3, size = 5, stroke = 1, show.legend = F) +
  geom_point(size = 5, alpha = 0.7, colour = "black", stroke = 1) +
  geom_errorbar(aes(ymin = ci.lbR, ymax = ci.ubR), 
                width = 0.1, colour = "black", linewidth = 1.5) +
  geom_text(aes(y = 0.10, label = paste0("italic(n) == ", n)), 
            colour = "black", parse = T, size = 7.5) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black", linewidth = 0.5) +
  scale_colour_manual(values = regcol, guide = "none")+
  scale_y_continuous(limits = c(-0.02, 0.12), breaks = seq(-0.02, 0.12, by = 0.02)) +
  labs(x = "Geographic region", y = expression("Estimated mean " * italic(F)[ST])) +
  theme_classic() +
  theme(
    panel.grid.major = element_line(color = "grey96"),
    axis.title = element_text(size = 30, colour = "black"),
    axis.text.x = element_text(size = 26), 
    axis.text.y = element_text(size = 26),
    axis.line = element_line(linewidth = 2),
    axis.ticks = element_line(linewidth = 2),
    axis.ticks.length = unit(.4, "cm"))


####################### 11. SCOPE ##########################
# In this analysis, interested in asking if the mean FST differ systematically 
# depending on how many geographic regions a study includes (1, 2, 3, or 7)?
# Some studies include just 1 while others are circumpolar. Do those that
# include more have higher mean FSTs? Here we treat scope as a category 
# with 1 = 1 region, 2 = 2 regions etc,.. except that 7 = 4 or more.

# the prediction is that studies with greater scope likely will have 
# higher mean FSTs than those that are smaller.

# We can use scope variable added to the dataset to include this in a meta-
# regression as we have done for marker, comparison, and geographic 
# region

# First, change scope to a factor (not continuous)
mcoord$scope <- as.factor(mcoord$scope)

# Run a model that includes scope as moderator (fixed) and study as random. 
# Tests for differences in mean FST across studies
metaresS <- rma.mv(yi = Fst, V = V, mods = ~ scope, random = ~ 1|study, data = mcoord)

# Look at the summary
summary(metaresS)

# look at the levels and table them to get a feel for the numbers
levels(mcoord$scope)
table(mcoord$scope)

# Because the metares8b results are relative to the intercept group
# need model matrix to convert to estimated values.
Xs <- model.matrix(~ scope, data = data.frame(scope = levels(mcoord$scope)))[, -1]

# Predict the mean FST for the data and correct them with model matrix.
scpred <- predict(metaresS, newmods = Xs)

# Can also list the number of comparisons used in each study scope.
nsc <- as.numeric(table(mcoord$scope))

# Make results into a dataframe easier to plot in ggplot.
scest <- data.frame(scope = levels(mcoord$scope), meanFSTS = scpred$pred,
                    seFSTS = scpred$se, ciS.lb = scpred$ci.lb,
                    ciS.ub = scpred$ci.ub, n = nsc)

# Set colours for the modeled FSTs by scope level
sccol <- c("1" = "dodgerblue4", "2" = "goldenrod4", "3" = "lightcoral", "7" = "mediumpurple3")

# Plot using ggplot as previous
ggplot(scest, aes(x = scope, y = meanFSTS)) +
  geom_jitter(data = mcoord, aes(x = scope, y = Fst, colour = scope), width = 0.15,
              alpha = 0.4, size = 5, stroke = 1, show.legend = F) +
  geom_point(size = 5, alpha = 0.7, colour = "black", stroke = 1) +
  geom_errorbar(aes(ymin = ciS.lb, ymax = ciS.ub), 
                width = 0.10, linewidth = 1.5, colour = "black") +
  geom_text(aes(y = 0.14, label = paste0("italic(n) == ", n)), 
            colour = "black", parse = T, size = 10, show.legend = F) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black", linewidth = 0.5) +
  scale_colour_manual(values = sccol, guide = "none") +
  scale_y_continuous(limits = c(-0.02, 0.18), breaks = seq(-0.02, 0.18, by = 0.04)) +
  labs(x = "# of regions included", y = expression("Estimated mean " * italic(F)[ST])) +
  theme_classic() +
  theme(
    panel.grid.major = element_line(color = "grey96"),
    axis.title = element_text(size = 30, colour = "black"),
    axis.text.x = element_text(size = 26), 
    axis.text.y = element_text(size = 26),
    axis.line = element_line(linewidth = 2),
    axis.ticks = element_line(linewidth = 2),
    axis.ticks.length = unit(.4, "cm"))


####################### 12. AMONG GRPS ##################
# Here, interested in asking among which regions is population structure 
# most different, taking into account the different studies and the 
# variation within them.

# different from the other analyses above because units are pairwise
# comparisons spanning different regions and so cannot be easily
# clumped into single regions. 

# First, filter to consider only among groups comparisons 
# getting rid of within groups.  
amgdat <- mcoord[is.na(mcoord$wingrp),]

# Initialise new character vectors to store regional codes.
arg1 <- as.character()
arg2 <- as.character()

# Sanity check to look at new amgdat structure
amgdat

# Make factor vector renaming regions as within analyses 
# above
rgc <- c("NWAT", "BBWG", "NOZE", "SVAL", "EGRE", "EARC", "WARC")

# Extract region list from the group data in amgdat
rgl <- unique(amgdat$group1) 

# Short loop populating the regcode variable according to 
# the reglist levels.
for (a in seq_along(rgl)) {
  idx1 <- which(amgdat$group1 == as.character(rgl[a]))
  idx2 <- which(amgdat$group2 == as.character(rgl[a]))
  arg1[idx1] <- rgc[a]
  arg2[idx2] <- rgc[a]
}

# Tagging the new indices to the end of amgdat
amgdat$arg1 <- as.factor(arg1)
amgdat$arg2 <- as.factor(arg2)

# Structuring data to include "regpair" variable listing regions being 
# compared for each estimated Fst in data generated from among-region 
# comparisons. regpair should be standardised with no directional value. 
amgdat <- amgdat |>
  mutate(r1 = pmin(as.character(arg1), as.character(arg2)),
         r2 = pmax(as.character(arg1), as.character(arg2)),
         regpair = factor(paste(r1,r2, sep = "_")))

# Analyses of amgdat data using regpair as moderator but taking study as 
# random factor to assess differences among among-region FSTs.
metares5 <- rma.mv(yi = Fst, V = V, mods = ~ regpair, random = ~ 1 | study, 
                   data = amgdat, method = "REML")

# Summerise model into object that will be used later
sm <- summary(metares5)
sm

# Create new df listing the among region pairs 
newdat <- data.frame(regpair = levels(amgdat$regpair))

# Use newdat to define model matrix for mean FST estimates.
# Should make 18 x 18 matrix
X <- model.matrix(~ regpair, newdat)
dim(X)

# Reform the model matrix to exclude first column as this is reference 
X <- X[, -1, drop = F]

# Convert coefficients derived in models to actual estimates using refined 
# model.matrix
preds <- predict(metares5, newmods = X)

# Make new df accumulating modeled coefficients and their
# error and tracking among region comparisons. 
plotdf <- data.frame(regpair = newdat$regpair, pfst = preds$pred, 
                     ci.lbp = preds$ci.lb, ci.ubp = preds$ci.ub)

# Test whether the confidence interval of each estimated FST includes 
# 0 or not. If it doesn't - star, if it does - no star. Tests whether 
# FST is significantly different from 0.
plotdf <- plotdf |>
  mutate(
    psig = !(ci.lbp <= 0 & ci.ubp >= 0),
    star = ifelse(psig, "*", "")
  )

# Estimated among-region FST for each region pair, accounting 
# for study-level random effects.
plotdf

# Split the regpair into r1 and r2 for input into circlise
# which needs a to: and a from:
chorddf <- plotdf |>
  separate(regpair, into = c("r1", "r2"), sep = "_")

# Sanity check to verify all comparisons in amgdat are included in 
# chorddf, and whether any chorddfs are duplicated. First should be
# T, and second should be F.
nrow(chorddf) == length(levels(amgdat$regpair))
any(duplicated(chorddf$regpair))

# Generate n, the sampling effort / representation of that region 
# across all among-region FSTs. Extract region appearances 
# from amgdat
regcts <- bind_rows(
  amgdat |> select(region = arg1),
  amgdat |> select(region = arg2)
)

# This should show how many times each region appears in the 432 
# among region analyses. It should count to 432*2 = 864
regcts

# Count how many times each region appears in regcts. This is the 
# value needed in the sectors of our circle figure.
secsize <- regcts |>
  count(region, name = "n")

# Sanity check <- each row in amgdat contributes 2 regions
# so with 2 *432 comparisons there should be 864 total ns
# in secsize
sum(secsize$n) == 2 * nrow(amgdat)


# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
## Circlise plot of among region pairwise FSTs ##

# Use regcol from above 

# Reset plotting par and canvas for circle
par(mar = c(0, 0, 0, 0))
circos.clear()

# Reorder sectors from west to east.
regord <- c("EARC", "NOZE", "SVAL", "EGRE", "NWAT", "BBWG", "WARC")
chorddf$r1 <- factor(chorddf$r1, levels = regord)
chorddf$r2 <- factor(chorddf$r2, levels = regord)

# Count how many significant comparisons touch each region
star_counts <- table(c(
  chorddf$r1[chorddf$star == "*"],
  chorddf$r2[chorddf$star == "*"]
))

# Set the circlise parameters for figure
circos.par(start.degree = 90, gap.degree = 4, track.margin = c(0.02, 0.02),
           cell.padding = c(0, 0, 0, 0), canvas.xlim = c(-1, 1),
           canvas.ylim = c(-1, 1))

# Extract the frequency of use for each region in pairwise comparisons  
secsize2 <- secsize[match(regord, secsize$region), ]

# Initialise the values to be used as sectors in the circle 
circos.initialize(factors = secsize2$region, xlim = cbind(rep(0, nrow(secsize2)), secsize2$n))

# Fsts can be small (especially here). So, good idea to make them 
# more visible by scaling up using rescale from scales package.              
fstcw <- scales::rescale(chorddf$pfst, to = c(0.8, 3))

# Set colours of the links in the circle to be alternating to make them 
# easier to see <- but colours do not mean anything here.
alt_cols <- c("powderblue", "lightpink", "palegreen")

# Assigns colour to each chord, cycling through alt_cols, and 
# making them semi-transparent.
chordcol <- adjustcolor(alt_cols[(seq_len(nrow(chorddf)) - 1) %% length(alt_cols) + 1],
                        alpha.f = 0.6)

# Plotting the chords in the circle figure which represent the magnitude of the 
# FST between the regions (r1, r2), and plotting them in the order of regord, 
# and with chordcol colours, with a thin border to make them easier to see.
# The last few cmds tell it to not plot tracks, sectors, but allocates space 
# for them outside circle.
chordDiagram(
  x = chorddf[, c("r1", "r2", "pfst")],
  order = regord,
  col = chordcol,
  link.lwd = fstcw,
  link.border = "azure",
  directional = 0,
  annotationTrack = NULL,   # do NOT redraw sectors
  grid.col = NULL,
  symmetric = T,
  preAllocateTracks = list(track.height = 0.1, track.height = 0.09))

# Plots track sectors (regions) by regcol and whose size reflects 
# frequency of use in pairwise comparisons. 
circos.trackPlotRegion(track.index = 1, ylim = c(0, 1), bg.col = regcol, 
                       bg.border = NA, #track.height = 0.09, 
                       panel.fun = function(x, y) {
                         
                         sector = CELL_META$sector.index
                         
                         circos.rect(CELL_META$xlim[1], 0,
                                     CELL_META$xlim[2], 1,
                                     col = regcol[sector],
                                     border = NA)
                         
                         circos.text(x = CELL_META$xcenter, y = CELL_META$ylim[2] + mm_y(5),
                                     labels = sector,
                                     facing = "bending",
                                     niceFacing = T,
                                     adj = c(0.5, 0.5),
                                     cex = 1.5)
                         
                         # number of comparisons for this sector
                         nc <- secsize2$n[secsize2$region == sector]
                         
                         # added inside colour band
                         circos.text(x = CELL_META$xcenter, y = 0.45, labels = nc,
                                     col = "white", cex = 1.5, font = 2)
                       }
)

# $$$$$$$$$$$$$$$$$$ SETTING STARS $$$$$$$$$$$$$$$$$$$$
# Set the stars in the plot to indicate which fsts are sign. 
# different from 0. 

# Set object storing the sectors that need stars from chorddf
starsec <- rep(c("BBWG", "WARC", "NWAT"), times = c(3, 1, 3))

# Set x value/sector that needs stars. x values are sector 
# specific. The commented code below inserted after the Chord.diagram
# will plot grid to help place stars (points). y values are constant 
# @ -0.8.

# Derived from plotted grid 
starx <- c(0.0062, 0.0158, 0.025, 0.02, 0.0031, 0.0115, 0.027)
stardf <- cbind.data.frame(sec = starsec, x = starx)

# Loop through stardf and adding x,y positions of the stars 
for (d in 1:length(starx)) {
  circos.points(track.index  = 1,
                sector.index = stardf$sec[d],
                x = stardf$x[d],
                y = -0.8, 
                pch = 8,
                cex = 1.4)
}

# $$$$$$$$$$$$$$$$$$ ADD LEGEND $$$$$$$$$$$$$$$$$$$$
# Adding legend manually using segment lengths
fstref <- c(0.001, 0.005, 0.01)

# scaling the fstref using similar scale as for actual fst.
fstreflen <- scales::rescale(fstref, to = c(0.04, 0.12))

# Allows plotting outside circle, going to left and bottom
par(xpd = NA)  
x0 <- -1.4     
y0 <- -0.75

# Size of separation between legend items
ystep <- 0.1

# loop placing segments and their estimated fst values
# below and to the left of the circle.
for (b in seq_along(fstref)) {
  segments(x0, y0 - (b - 1) * ystep, x0 + fstreflen[b], y0 - (b - 1) * ystep,
           lwd = 2, col = "grey30")
  
  text(x0 + fstreflen[b] + 0.01, y0 - (b - 1) * ystep, 
       labels = bquote(italic(F)[ST] == .(fstref[b])), adj = 0, cex = 1.5)
}

# Adding the legend title (needed here)
text(x0, y0 + ystep, "Chord width", adj = 0, font = 2, cex = 1.5)

## END circlise.
circos.clear()


####################### 13. OVERALL IBD ##########################
# In this part of the analyses, we are asking whether there is evidence
# throughout the included studies of isolation by distance, i.e., is FST 
# impacted by how far samples are from one another. As before, unit of 
# analysis should be FST comparison, but keeping track of the 
# non-independence of values within study. 

# Also important is to linearise FSTs as FST/(1-FST) values, and then to also 
# consider the distances in log forms.   

# First - linearise fst as per most IBD studies.
mcoord$lFst <- mcoord$Fst/(1-mcoord$Fst)

# However probably need to adjust FST variance to match linearised 
# version, so as to take variance into account in meta-analysis 
# framework. Turns out the way to linearise FST variance is to 
# derive it using the formula V = V/(1-FST)^4 (derivative).
mcoord$Vl <- mcoord$V / (1-mcoord$Fst)^4

# Distances should also be log scaled to account for the fact 
# that some values will be very large relative to others.
mcoord$lwwd <- log(mcoord$waterdist)
mcoord$lcud <- log(mcoord$curdist)

# With converted data, can run the IBDs using meta-analysis framework
# with lwwd/lcud as moderator(s), and considering study as a random factor
# (account for non-independence in values within study).

# First compute the ibdnull model for comparisons which tests the Ho that
# Genetic differentiation (lFst) does not depend on geographic distance, or 
# after accounting for study level differences, pairwise Fsts are drawn from 
# the same distribution regardless of distance. 
ibdnull <- rma.mv(yi = lFst, V = Vl, random = ~ 1|study,
                  data = mcoord, method = "REML")
summary(ibdnull)

# Then with waterway distance as moderator - Tests whether adding distance
# (here log(distance)) explain some systematic variation in lFst beyond 
# study level differences.
metares6 <- rma.mv(yi = lFst, V = Vl, mods = ~ lwwd, random = ~ 1|study,
                   data = mcoord, method = "REML")

# Save IBD lwwd summary to object
smIBD1 <- summary(metares6)
smIBD1

# Testing if including distance is better than not for lwwd. 
anova(ibdnull, metares6, refit = T)


# Second, we try it again with the lcud
# Then with waterway distance as moderator
metares7 <- rma.mv(yi = lFst, V = Vl, mods = ~ lcud, random = ~ 1|study,
                   data = mcoord, method = "REML")

# Save IBD lcud summary to object
smIBD2 <- summary(metares7)
smIBD2

# Testing if including distance is better than not for lwwd. 
anova(ibdnull, metares7, refit = T)

# Compute pseudo-R²s by comparing the sigma2s from both models to null.
taunull <- sum(ibdnull$sigma2)
taum6 <- sum(metares6$sigma2)
taum7 <- sum(metares7$sigma2)

# Taus calculated and normalised for lwwd and lcud models
r2wd <- (taunull - taum6)/taunull
r2wd

r2cd <- (taunull - taum7)/taunull
r2cd

# Generate new df to plot the data and fit the meta-regression line.
# First for lwwd.

# Generate a ordered sequence of 200 lwwds 
prdconf1 <- data.frame(lwwd = seq(min(mcoord$lwwd), 
                                 max(mcoord$lwwd), length.out = 200))

# Use these to predict what the linearised FSTs should be | model.
lpred1 <- predict(metares6, newmods = prdconf1$lwwd)

# Store predictions and CIs into predconf df
prdconf1$fit <- lpred1$pred
prdconf1$ci.lb <- lpred1$ci.lb
prdconf1$ci.ub <- lpred1$ci.ub

# Plot using ggplot 
ggplot(mcoord, aes(x = lwwd, y = lFst, colour = study)) +
  geom_point(size = 5, alpha = 0.6, stroke = 1) +
  geom_ribbon(data = prdconf1, aes(x = lwwd, ymin = ci.lb, ymax = ci.ub),
              inherit.aes = F, fill = "gray70", alpha = 0.4) +
  geom_line(data = prdconf1, aes(x = lwwd, y = fit), inherit.aes = F, 
                                colour = "firebrick4", linewidth = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black", linewidth = 0.5) +
  annotate("text", x = 2.4, y = 0.20, label = paste0("pseudo-", "italic(R)^2 == ", 
  round(r2wd, 3)), parse = T, size = 12) +
  annotate("text", x = 2.4, y = 0.18, label = paste0("slope = ", 
  formatC(as.numeric(coef(metares6)[2]), format = "f", digits = 4)), size = 10)+
  annotate("text", x = 2.4, y = 0.16, label = paste0("italic(p) == '", 
  formatC(as.numeric(metares6$pval[2]), format = "f", digits = 3),"'"), parse = T, size = 10) +
  scale_y_continuous(limits = c(-0.02, 0.22), breaks = seq(-0.02, 0.22, by = 0.04)) +
  scale_x_continuous(limits = c(-0, 9), breaks = seq(0, 9, by = 2)) +
  scale_colour_manual(values = grpcol, guide = "none" )+
  #scale_colour_manual(values = grpcol) + comment above and decomment here to see legend
  labs(x = "Log(distance)", y = expression("Linearised " * italic(F)[ST])) +
  theme_classic() +
  theme(
    panel.grid.major = element_line(color = "grey96"),
    axis.title = element_text(size = 30, colour = "black"),
    axis.text.x = element_text(size = 26), 
    axis.text.y = element_text(size = 26),
    axis.line = element_line(linewidth = 1.3),
    axis.ticks = element_line(linewidth = 1),
    axis.ticks.length = unit(.2, "cm")) # + uncomment this and below to see legend
#  theme(legend.text = element_text(size = 15),
#        legend.title = element_blank(),
#        legend.position = "bottom",
#        legend.background = element_blank())
  

# Then do the same for the lcud model.
# New df to plot data and fit meta-regression line for lcud.

# Generate a ordered sequence of 200 lcuds 
prdconf2 <- data.frame(lcud = seq(min(mcoord$lcud), 
                                  max(mcoord$lcud), length.out = 200))

# Use these to predict what the linearised FSTs should be | model.
lpred2 <- predict(metares7, newmods = prdconf2$lcud)

# Store predictions and CIs into predconf df
prdconf2$fit <- lpred2$pred
prdconf2$ci.lb <- lpred2$ci.lb
prdconf2$ci.ub <- lpred2$ci.ub

# Plot using ggplot 
ggplot(mcoord, aes(x = lcud, y = lFst, colour = study)) +
  geom_point(size = 5, alpha = 0.6, stroke = 1) +
  geom_ribbon(data = prdconf2, aes(x = lcud, ymin = ci.lb, ymax = ci.ub),
              inherit.aes = F, fill = "gray70", alpha = 0.4) +
  geom_line(data = prdconf2, aes(x = lcud, y = fit), inherit.aes = F, 
            colour = "firebrick4", linewidth = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black", linewidth = 0.5) +
  annotate("text", x = 3, y = 0.20, label = paste0("pseudo-", "italic(R)^2 == ", 
     round(r2cd, 3)), parse = T, size = 12) +
  annotate("text", x = 3, y = 0.18, label = paste0("slope = ", 
     formatC(as.numeric(coef(metares7)[2]), format = "f", digits = 4)), size = 10)+
  annotate("text", x = 3, y = 0.16, label = paste0("italic(p) == '", 
     formatC(as.numeric(metares7$pval[2]), format = "f", digits = 3),"'"), parse = T, size = 10) +
  scale_y_continuous(limits = c(-0.02, 0.22), breaks = seq(-0.02, 0.22, by = 0.04)) +
  scale_x_continuous(limits = c(-0, 13), breaks = seq(0, 13, by = 2)) +
  scale_colour_manual(values = grpcol, guide = "none" )+
  labs(x = "Log(distance)", y = expression("Linearised " * italic(F)[ST])) +
  theme_classic() +
  theme(
    panel.grid.major = element_line(color = "grey96"),
    axis.title = element_text(size = 30, colour = "black"),
    axis.text.x = element_text(size = 26), 
    axis.text.y = element_text(size = 26),
    axis.line = element_line(linewidth = 1.3),
    axis.ticks = element_line(linewidth = 1),
    axis.ticks.length = unit(.2, "cm")) 


####################### 14. IBD SEPARATE SLOPES LWWD ##########################
# While the above found that over the entire set of studies IBD was small (LWWD)
# to negligible (lcud), we want to know whether different regions have different 
# IBD within them. This asks whether the IBD is more or less constrained in different
# geographical regions based on the studies gathered.

# To do this we need to filter the data to include only within region FSTs as 
# among regions would create bias upward by including among region comparisons.

# We already have windat (regional comparisons) but need to calculate the 
# linearised Fst, lFst variance, and log distances for this data too.
windat$lFst <- windat$Fst/(1-windat$Fst)
windat$Vl <- windat$V / (1-windat$Fst)^4

# Distances should also be log scaled to account for the fact 
# that some values will be very large relative to others.
windat$lwwd <- log(windat$waterdist)
windat$lcud <- log(windat$curdist)

# Reordering the results in the way that we ordered regions above.
windat$regcode <- factor(windat$regcode, levels = regord)

# This new analysis will look specifically at within region fsts and test for IBD 
# within the regions considering only the Fsts within those regions. 
metares8 <- rma.mv(yi = lFst, V = Vl, mods = ~ lwwd * regcode, random = ~ 1 | study,
                   data = windat, method = "REML")

# Summerise the model 
summary(metares8)

# When considering only within-region pairwise FSTs, we found no evidence of IBD, 
# and no indication that the relationship between pairwise FST and geographic distance
# differed among regions.

# Again, the summary shows that there is no difference in the IBD relationships 
# estimated within the different regions (in this case using the EARC
# as the initial group). The full moderator set explains basically nothing.

# Overall IBD within regions: Is the baseline slope (reference region) different from 0?
# Regional differences in IBD strength: Do other regions have slopes different from 
# the reference region?

# No evidence of IBD in the reference region. lwwd est. ~ 0.00 p = 0.98
# All interaction terms: lwwd:regcodeXXXX p >> 0.05

# Generating the prediction distances to plot the modeled relationship
# that considers both the lwwd and region, but also the non-independence 
# of pairwise fsts estimated within studies.
winreshp <- windat |>
  group_by(regcode) |>
  summarise(minlwwd = min(lwwd, na.rm = T),
            maxlwwd = max(lwwd, na.rm = T),
            .groups = "drop") |>
  rowwise() |>
  reframe(regcode, lwwd = seq(minlwwd, maxlwwd, length.out = 100)) 

# Again reordering regions for plotting
winreshp$regcode <- factor(winreshp$regcode, levels = regord)

# Using the model predictions to generate the best fit linearised
# relationship for each region. First extract the model matrix and 
# remove the intercept
wrpred <- predict(metares8, newmods = model.matrix(~ lwwd * regcode, winreshp)[,-1])

# Setting the linearised data into a dataframe for plotting
winreshp$pred <- wrpred$pred
winreshp$ci.lb <- wrpred$ci.lb
winreshp$ci.ub <- wrpred$ci.ub

# To add slopes and p-values to each facet in the data means these
# must be calculated from the metares8 results

# Get coefficients and the vcov matrix of the model 
mc <- coef(metares8)
mvcv <- vcov(metares8)

# Again listing the regions in the order we want them in.
regions <- c("EARC","NOZE","SVAL","EGRE","NWAT","BBWG","WARC")

# Create an empty df calculating IBD slopes, vars, and ps for each region.
slopedf <- data.frame(regcode = character(), 
                      slope = numeric(), mse = numeric(), p = numeric())

# Use a loop to calculate slopes and predicted values for the 100 lwwd 
# values for each region generated above. Here, need to add coef from each 
# region to the "lwwd" estimate which is the reference region (EARC). 
# vars (for 95%CI) and to estimate uncertainty in the coef needs to use 
# the formula: Var(A + B) = Var(A) + Var(B) + 2Cov(A,B)
# This is beacuse the coef is part of the interaction of lwwd and region
# and so the var in both region and lwwd must be taken together to track
# regional coef confidence (and assessments such as p values). 
for(e in regions){
  
  # Define coefficient names
  base  <- "lwwd"
  inter <- paste0("lwwd:regcode", e)
  
  # Calculate slope
  if(inter %in% names(mc)) { slope <- mc[base] + mc[inter]
  } else {
    slope <- mc[base]
  }
  
  # Calculate variance of slope
  vars <- mvcv[base, base]
  
  if(inter %in% colnames(mvcv)){
    vars <- vars + mvcv[inter, inter] + 2 * mvcv[base, inter]
  }
  
  # Standard error
  mse <- sqrt(vars)
  
  # Test statistic
  z <- slope/mse
  
  # p value 
  p <- 2 * pnorm(abs(z), lower.tail = FALSE)
  
  # Add row to slopedf for storage
  slopedf <- rbind(slopedf, data.frame(regcode = e, slope = slope, 
                                       mse = mse,p = p))
}

# Add a label to the df for plotting 
slopedf$lab <- sprintf("slope = %.5f\np = %.4f", slopedf$slope, slopedf$p)

# Make sure the slopedf is in the same order as plotting
slopedf$regcode <- factor(slopedf$regcode, levels = regord)

# Combine the plot of the raw within region pairwise Fst, with the modeled 
# predictions - putting these in a facet plot allows all to be displayed at once
ggplot(windat, aes(lwwd, lFst)) +
  geom_point(aes(colour = regcode), size = 5, alpha = 0.5, stroke = 1)+
  geom_line(data = winreshp, aes(x = lwwd, y = pred), colour = "firebrick4", inherit.aes = F) +
  geom_ribbon(data = winreshp, aes(x = lwwd, ymin = ci.lb, ymax = ci.ub), alpha = 0.2, inherit.aes = F) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black", linewidth = 0.5) +
  geom_text(data = slopedf, aes(x = -Inf, y = Inf, label = lab), 
           hjust = -0.1, vjust = 2, size = 7)+
  labs(x = "Log(distance)", y = expression("Linearised " * italic(F)[ST])) +
  scale_colour_manual(values = regcol) +
  facet_wrap(~ regcode, strip.position = "top") +
  theme_classic()+   
  theme(panel.grid.major = element_line(color = "grey96"),
    axis.title = element_text(size = 26, colour = "black"),
    axis.text.x = element_text(size = 22), 
    axis.text.y = element_text(size = 22),
    axis.line = element_line(linewidth = 1.2),
    axis.ticks = element_line(linewidth = 1.2),
    axis.ticks.length = unit(.2, "cm"),
    legend.position = "none", 
    plot.title = element_text(family = "Helvetica", face = "bold.italic", size = (20)),
    strip.background = element_blank(),
    strip.text = element_text(size = 20)) 



####################### 15. IBD SEPARATE SLOPES LCUD ##########################
# Can do the same for IBD with along current distances / region
# Same reordering the results in the way that we ordered regions above.
windat$regcode <- factor(windat$regcode, levels = regord)

# This analysis looks at within region fsts and test for IBD using lcud 
# within the regions considering only the Fsts within those regions. 
metares9 <- rma.mv(yi = lFst, V = Vl, mods = ~ lcud * regcode, random = ~ 1 | study,
                   data = windat, method = "REML")

# Summerise the model 
summary(metares9)

# When considering only within-region pairwise FSTs, we found no evidence of IBD, 
# and no indication that the relationship between pairwise FST and geographic distance
# differed among regions.

# Again, the summary shows that there is no difference in the IBD relationships 
# estimated within the different regions (in this case using the EARC
# as the initial group). The full moderator set explains basically nothing.

# Overall IBD within regions: Is the baseline slope (reference region) different from 0?
# Regional differences in IBD strength: Do other regions have slopes different from 
# the reference region?

# No evidence of IBD in the reference region. lcud est. ~ 0.00 p = 0.98
# All interaction terms: lcud:regcodeXXXX p >> 0.05

# Generating the prediction distances to plot the modeled relationship
# that considers both the lcud and region, but also the non-independence 
# of pairwise fsts estimated within studies.
winreshp2 <- windat |>
  group_by(regcode) |>
  summarise(minlcud = min(lcud, na.rm = T),
            maxlcud = max(lcud, na.rm = T),
            .groups = "drop") |>
  rowwise() |>
  reframe(regcode, lcud = seq(minlcud, maxlcud, length.out = 100)) 

# Reordering regions for plotting
winreshp2$regcode <- factor(winreshp2$regcode, levels = regord)

# Using the model predictions to generate the best fit linearised
# relationship for each region. First extract the model matrix and 
# remove the intercept
wrpred2 <- predict(metares9, newmods = model.matrix(~ lcud * regcode, winreshp2)[,-1])

# Setting the linearised data into a dataframe for plotting
winreshp2$pred <- wrpred2$pred
winreshp2$ci.lb <- wrpred2$ci.lb
winreshp2$ci.ub <- wrpred2$ci.ub

# To add slopes and p-values to each facet in the data means these
# must be calculated from the metares9 results

# Get coefficients and the vcov matrix of the model 
mc2 <- coef(metares9)
mvcv2 <- vcov(metares9)

# Again listing the regions in the order we want them in.
regions <- c("EARC","NOZE","SVAL","EGRE","NWAT","BBWG","WARC")

# Create an empty df calculating IBD slopes, vars, and ps for each region.
slopedf2 <- data.frame(regcode = character(), 
                      slope = numeric(), mse = numeric(), p = numeric())

# Use a loop calculate the slopes and the predicted values for the 
# 100 lcud values for each region generated above. Here, need to add
# coef from each region to the "lcud" estimate which is the reference
# region (EARC). vars (for 95%CI) and to estimate uncertainty in the 
# coef needs to use the formula: Var(A + B) = Var(A) + Var(B) + 2Cov(A,B)
# This is beacuse the coef is part of the interaction of lwwd and region
# and so the var in both region and lwwd must be taken together to track
# regional coef confidence (and assessments such as p values). 
for(f in regions){
  
  # Define coefficient names
  base  <- "lcud"
  inter <- paste0("lcud:regcode", f)
  
  # Calculate slope
  if(inter %in% names(mc2)) { slope <- mc2[base] + mc2[inter]
  } else {
    slope <- mc2[base]
  }
  
  # Calculate variance of slope
  vars <- mvcv2[base, base]
  
  if(inter %in% colnames(mvcv2)){
    vars <- vars + mvcv2[inter, inter] + 2 * mvcv2[base, inter]
  }
  
  # Standard error
  mse <- sqrt(vars)
  
  # Test statistic
  z <- slope/mse
  
  # p value 
  p <- 2 * pnorm(abs(z), lower.tail = F)

  # Add row to slopedf for storage
  slopedf2 <- rbind(slopedf2, data.frame(regcode = f, slope = slope, 
                                       mse = mse, p = p))
}

# Add a label to the df for plotting 
slopedf2$lab <- sprintf("slope = %.5f\np = %.4f", slopedf2$slope, slopedf2$p)

# Make sure the slopedf is in the same order as plotting
slopedf2$regcode <- factor(slopedf2$regcode, levels = regord)

# Combine the plot of the raw within region pairwise Fst, with the modeled 
# predictions - putting these in a facet plot allows all to be displayed at once
ggplot(windat, aes(lcud, lFst)) +
  geom_point(aes(colour = regcode), size = 5, alpha = 0.5, stroke = 1)+
  geom_line(data = winreshp2, aes(x = lcud, y = pred), colour = "firebrick4", inherit.aes = F) +
  geom_ribbon(data = winreshp2, aes(x = lcud, ymin = ci.lb, ymax = ci.ub), alpha = 0.2, inherit.aes = F) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black", linewidth = 0.5) +
  geom_text(data = slopedf2, aes(x = -Inf, y = Inf, label = lab), 
            hjust = -0.1, vjust = 2, size = 7)+
  labs(x = "Log(distance)", y = expression("Linearised " * italic(F)[ST])) +
  scale_colour_manual(values = regcol) +
  facet_wrap(~ regcode, strip.position = "top") +
  theme_classic()+   theme(
    panel.grid.major = element_line(color = "grey96"),
    axis.title = element_text(size = 26, colour = "black"),
    axis.text.x = element_text(size = 22), 
    axis.text.y = element_text(size = 22),
    axis.line = element_line(linewidth = 1.2),
    axis.ticks = element_line(linewidth = 1.2),
    axis.ticks.length = unit(.2, "cm"),
    legend.position = "none", 
    plot.title = element_text(family = "Helvetica", face = "bold.italic", size = (20)),
    strip.background = element_blank(),
    strip.text = element_text(size = 20)) 


################################## 16.NEARSHORE OFFSHORE ############################
# Create Nearshore - Offshore indices based on distance to shore from sample location. 
# This was done in using a short script below but it was written for manipulating
# excel files. This is because many studies have tried to show that there may be
# fine-scale structure between nearshore/fjord individuals and offshore ones. But, 
# this is not convincing so far. Here we try with meta-data.

# Next few cmds run once. Can then be commented out as the manipulation was permanently
# added to working data.

#  mutate(
#    distance_to_shore1 = as.numeric(distance_to_shore1),
#    distance_to_shore2 = as.numeric(distance_to_shore2),
#    shore_type1 = ifelse(distance_to_shore1 < 50000, "nearshore", "offshore"),
#    shore_type2 = ifelse(distance_to_shore2 < 50000, "nearshore", "offshore")
#  )

#write_xlsx(data, "workingdata_with_shoretype.xlsx")

# Create "Dist_type" variable based on pairwise comparisons
# using dplyr cmds
mcoord <- mcoord %>%
  mutate(
    Dist_type = case_when(
      shore_type1 == "nearshore" & shore_type2 == "nearshore" ~ "NN",
      shore_type1 == "offshore" & shore_type2 == "offshore" ~ "OO",
      T ~ "NO"
    )
  )

# First run model on the entire dataset 
metares10 <- rma.mv(yi = Fst, V = mcoord$V, mods = ~ Dist_type, random = ~ 1|study,
                   data = mcoord, method = "REML")

# Summarise model
summary(metares10)

# Assess results of marker specifically (moderator effects) or
# Omnibus test.
anova(metares10)

# rma.mv is not supported by emmeans <-  so all calculations
# that would be provided by emmeans must be done by hand.

# Extracting regression coefficients from model for each 
# of the moderator levels
regco <- coef(metares10)

# Also extract the variance–covariance matrix for the 
# different moderator values
mr10vcv <- vcov(metares10)

# Make L matrix specifying which factor levels to compare
# for example the first characterises mtDNA, and the next one 
# looks at RAPD compared with mtDNA, etc,.. Sort of like 
# model.matrix cmd in base R.
L <- rbind(
  NN  = c(1, 0, 0),
  NO  = c(1, 1, 0),
  OO  = c(1, 0, 1)
)

# Setting column names for L matrix to those of the extracted
# coefficients
colnames(L) <- names(regco)

# Estimate predicted mean FST per marker type by multiplying L 
# matrix with that of the regco (regression coefficients) -
# essentially using mtDNA (intercept) to convert other values.
meanFST <- as.numeric(L %*% regco)

# Calculating the standard error (FULL covariance)
seFST <- sqrt(diag(L %*% mr10vcv %*% t(L)))

# Normal value for 95% CI (97.5 % quantile)
z <- 1.96 

# Calculating 95% CI for meanFST based on seFST from the 
# model
ci.lb <- meanFST - z * seFST
ci.ub <- meanFST + z * seFST

# Gather marker data into a small dataframe for plotting
# using ggplot
shore <- data.frame(shore = rownames(L), n = as.numeric(table(mcoord$Dist_type)),
                      meanFST = meanFST, se = seFST, ci.lb = ci.lb, ci.ub = ci.ub)

# Set colours for the modelled FSTs by Dist to shore
shorecol <- c("NN" = "firebrick4", "OO" = "springgreen4", "NO" = "darkcyan")

# Plot using ggplot 
ggplot(shore, aes(x = shore, y = meanFST)) +
  geom_jitter(data = mcoord, aes(x = Dist_type, y = Fst, colour = Dist_type), width = 0.15,
              alpha = 0.3, size = 5, stroke = 1, show.legend = F) +
  geom_point(size = 5, colour = 'black', alpha = 0.7, stroke = 1) +
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), 
                width = 0.1, colour = "black", linewidth = 1.5) +
  geom_text(aes(y = 0.14, label = paste0("italic(n) == ", n)), parse = T, ,size = 10) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black", linewidth = 0.7) +
  scale_colour_manual(values = shorecol)+
  scale_y_continuous(limits = c(-0.02, 0.18), breaks = seq(-0.02, 0.18, by = 0.04)) +
  labs(x = "Shore type", y = expression("Estimated mean " * italic(F)[ST])) +
  theme_classic() +
  theme(
    panel.grid.major = element_line(color = "grey96"),
    axis.title = element_text(size = 30, colour = "black"),
    axis.text.x = element_text(size = 26), 
    axis.text.y = element_text(size = 26),
    axis.line = element_line(linewidth = 2),
    axis.ticks = element_line(linewidth = 2),
    axis.ticks.length = unit(.4, "cm"))


# Run same model but on only Among dataset 
# Keep only among 
mcoord_among <- mcoord %>%
  filter(comptype == "among")

# Run the model
metares11 <- rma.mv(yi = Fst,V = V, mods = ~ Dist_type, random = ~ 1 | study, data = mcoord_among, method = "REML")

# Summarise model
summary(metares11)

# assign coefficients
regco11 <- coef(metares11)

# assign vcv
mr11vcv <- vcov(metares11)

# Make L matrix specifying which factor levels to compare
# for example the first characterises mtDNA, and the next one 
# looks at RAPD compared with mtDNA, etc,.. Sort of like 
# model.matrix cmd in base R.
L <- rbind(
  NN  = c(1, 0, 0),
  NO  = c(1, 1, 0),
  OO  = c(1, 0, 1)
)

# Setting column names for L matrix to those of the extracted
# coefficients
colnames(L) <- names(regco11)

# Estimate predicted mean FST per marker type by multiplying L 
# matrix with that of the regco (regression coefficients) -
# essentially using mtDNA (intercept) to convert other values.
meanFST11 <- as.numeric(L %*% regco11)

# Calculating the standard error (FULL covariance)
seFST11 <- sqrt(diag(L %*% mr11vcv %*% t(L)))

# Normal value for 95% CI (97.5 % quantile)
z <- 1.96 

# Calculating 95% CI for meanFST based on seFST from the 
# model
ci.lb <- meanFST11 - z * seFST11
ci.ub <- meanFST11 + z * seFST11

# Gather marker data into a small dataframe for plotting
# using ggplot
shore <- data.frame(shore = rownames(L), n = as.numeric(table(mcoord_among$Dist_type)),
                    meanFST11 = meanFST11, se = seFST11, ci.lb = ci.lb, ci.ub = ci.ub)

# Set colours for the modelled FSTs by Dist to shore
shorecol <- c("NN" = "firebrick4", "OO" = "springgreen4", "NO" = "darkcyan")


# Plot using ggplot 
ggplot(shore, aes(x = shore, y = meanFST11)) +
  geom_jitter(data = mcoord_among, aes(x = Dist_type, y = Fst, colour = Dist_type), width = 0.15,
              alpha = 0.3, size = 5, stroke = 1, show.legend = F) +
  geom_point(size = 5, colour = 'black', alpha = 0.7, stroke = 1) +
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), 
                width = 0.1, colour = "black", linewidth = 1.5) +
  geom_text(aes(y = 0.14, label = paste0("italic(n) == ", n)), parse = T, ,size = 10) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black", linewidth = 0.7) +
  scale_colour_manual(values = shorecol)+
  scale_y_continuous(limits = c(-0.02, 0.18), breaks = seq(-0.02, 0.18, by = 0.04)) +
  labs(x = "Shore type", y = expression("Estimated mean " * italic(F)[ST])) +
  theme_classic() +
  theme(
    panel.grid.major = element_line(color = "grey96"),
    axis.title = element_text(size = 30, colour = "black"),
    axis.text.x = element_text(size = 26), 
    axis.text.y = element_text(size = 26),
    axis.line = element_line(linewidth = 2),
    axis.ticks = element_line(linewidth = 2),
    axis.ticks.length = unit(.4, "cm"))


#### Within specific Region ? #####
# Remove lines with NA that correspond to among group
mcoord_clean <- mcoord %>%
  filter(!is.na(wingrp))

# Identify the wingrp that have the 3 Dist_type
wingrp_valid <- mcoord_clean %>%
  group_by(wingrp) %>%
  summarise(types = n_distinct(Dist_type)) %>%
  filter(types == 3) %>%
  pull(wingrp)

# Filter that dataset to only keep those 
mcoord_filtered <- mcoord_clean %>%
  filter(wingrp %in% wingrp_valid)

# Check to see which regions are included
table(mcoord_filtered$Dist_type, mcoord_filtered$wingrp)

# Run the new model
metares12 <- rma.mv(yi = Fst, V = mcoord_filtered$V, mods = ~ Dist_type * wingrp, random = ~ 1 | study,
  data = mcoord_filtered, method = "REML")

# Summarise
summary(metares12)

# Reordering and filtering levels that have enough comparisons to be useful
mcoord_filtered$Dist_type <- factor(mcoord_filtered$Dist_type, levels = c("NN", "NO", "OO"))
mcoord_filtered$wingrp    <- factor(mcoord_filtered$wingrp, levels = unique(mcoord_filtered$wingrp))

# Extracting and listing existing combinations 
combos_obs <- unique(mcoord_filtered[, c("Dist_type", "wingrp")])

# Create L matrix for calculating summaries of comparisons.
L <- model.matrix(~ Dist_type * wingrp, data = combos_obs)
coef_names <- names(coef(metares12))

# Assign coefficients 
regco12 <- coef(metares12)
mr12vcv <- vcov(metares12)

# Recalculating the adjusted means, SE and CIs
meanFST12 <- as.numeric(as.matrix(L) %*% regco12)
seFST12   <- sqrt(rowSums((as.matrix(L) %*% mr12vcv) * as.matrix(L)))
preds_df <- combos_obs %>%
  mutate(
    meanFST11 = meanFST12,
    ci.lb   = meanFST12 - 1.96 * seFST12,
    ci.ub   = meanFST12 + 1.96 * seFST12
  )

n_df <- mcoord_filtered %>%
  group_by(wingrp, Dist_type) %>%
  summarise(n = n(), .groups = "drop")

# Plotting with ggplot
shorecol <- c("NN" = "firebrick4", "OO" = "springgreen4", "NO" = "darkcyan")

ggplot(mcoord_filtered, aes(x = Dist_type, y = Fst, colour = Dist_type)) +
  geom_jitter(width = 0.2, size = 5, alpha = 0.5, stroke = 1) +                   
  geom_point(data = preds_df, aes(x = Dist_type, y = meanFST12),
             size = 5, shape = 18, colour = "black", inherit.aes = F) +
  geom_errorbar(data = preds_df, aes(x = Dist_type, ymin = ci.lb, ymax = ci.ub),
                width = 0.15, linewidth = 1, colour = "black", inherit.aes = F) +                     
  geom_text(data = n_df, aes(x = Dist_type, y = 0.14, label = paste0("italic(n) == ", n)),
            parse = T, size = 10, inherit.aes = F) + 
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black", linewidth = 0.5) +
  scale_y_continuous(limits = c(-0.02, 0.16), breaks = seq(-0.02, 0.16, by = 0.04)) +
  facet_wrap(~ wingrp, strip.position = "top") +
  scale_colour_manual(values = shorecol) +
  labs(x = "Shore type", y = expression("Estimated mean F" [ST])) +
  theme_classic() +
  theme(
    panel.grid.major = element_line(color = "grey96"),
    axis.title = element_text(size = 28, colour = "black"),
    axis.text.x = element_text(size = 22), 
    axis.text.y = element_text(size = 22),
    axis.line = element_line(linewidth = 1.2),
    axis.ticks = element_line(linewidth = 1.2),
    axis.ticks.length = unit(.2, "cm"),
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(size = 24)
  )


#################################### END #########################################

# TEMPORARY GRID FOR CIRCLISE 

# circos.trackPlotRegion(
#  track.index = 1,
#  bg.border = "grey80",
#  panel.fun = function(x, y) {
#
#    xlim <- get.cell.meta.data("xlim")
#    ylim <- get.cell.meta.data("ylim")
#
#    # X grid
#    xs <- seq(xlim[1], xlim[2], length.out = 6)
#    for (xx in xs) {
#      circos.lines(c(xx, xx), ylim, col = "grey85", lty = 3)
#      circos.text(xx, ylim[1],
#                  labels = round(xx, 3),
#                  cex = 0.5,
#                  adj = c(0.5, 1.2),
#                  col = "grey40")
#    }
#
#    # Y grid
#    ys <- seq(ylim[1], ylim[2], by = 0.25)
#    for (yy in ys) {
#      circos.lines(xlim, c(yy, yy), col = "grey85", lty = 3)
#      circos.text(xlim[1], yy,
#                  labels = yy,
#                  cex = 0.5,
#                  adj = c(1.2, 0.5),
#                  col = "grey40")
#    }
#  }
#)
#

# END TEMPORARY GRID


