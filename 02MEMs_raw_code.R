#Load libraries
library(ape)
library(spdep)
library(ade4)
library(vegan) # contains the pcnm function
library(packfor)
library(MASS)

# Import the data
setwd("C:\\Users\\adminuser\\Desktop\\EVERYTHING THOMAS\\gitprojects\\Moran_Eigenvector_Maps_MOECC")
a <- read.csv("chrysophyte-dorset.csv")
# b <- environmental variables to come
c <- read.csv("chrysophyte-dorset-xy.csv")

# Hellinger-transform species data
algae.h <- decostand(a, "hellinger")
# Check for significant trends in response data
anova(rda(algae.h, c)) 
# Detrend the transformed species data
algae.h.det <- resid(lm(as.matrix(algae.h) ~ ., data = c))

# Construction of spatial variables at all relevant scales
# PCNM analysis of data using Vegan's "pcnm" code
algae.pcnm.vegan <- pcnm(dist(c))
algae.pcnm <- as.data.frame(algae.pcnm.vegan$vectors)
nb.ev <- length(which(algae.pcnm.vegan$values > 0.0000001))

# Run a PCNM on the detrended species data and determine if analysis is significant
algae.pcnm.rda <- rda(algae.h.det, algae.pcnm)
anova.cca(algae.pcnm.rda) # it is significant

# Compute Adjusted R2 value and run forward selection of the PCNM variables on the species data.  You do this because the analysis above was significant
(algae.r2a <- RsquareAdj(algae.pcnm.rda)$adj.r.squared)

#Run forward selection to determine which MEMs are the most important predictors
(algae.pcnm.fwd <- forward.sel(algae.h.det, as.matrix(algae.pcnm), adjR2thresh = algae.r2a))

### Note: In the above selection procedure, sometimes 4 PCNMs come out as important because the fourth is either just under or just over the significance threshold (0.05)

## Gives the number of SIGNIFICANT PCMNs
(nb.sig.pcnm <- nrow(algae.pcnm.fwd))

# The below plots are purely exploratory. They show how increasing PCNMs decreases the scale of the period from broad to fine
par(mfrow=c(1,3))
plot(algae.pcnm$PCNM2, type="o", col="black", ylab = "PCNM 2")
plot(algae.pcnm$PCNM4, type="o", col="black", ylab = "PCNM 4")
plot(algae.pcnm$PCNM8, type="o", col="black", ylab = "PCNM 8")

# Identify the significant PCNMs in increasing order of significance
(pcnm.sign <- sort(algae.pcnm.fwd[,2]))
# Write the 4 significant PCNMs to a new object
pcnm.significant <- algae.pcnm[,c(pcnm.sign)]
pcnm.significant <- pcnm.significant[,c(1:3)]
write.csv(pcnm.significant, "MEM variables.csv", row.names = FALSE)

# New PCNM analysis with 3 significant PCNM predictors
algae.pcnm.rda2 <- rda(algae.h.det ~ ., data = pcnm.significant)
# Adjusted R2 after forward selection: R2adj = 0.2215
(algae.fwd.R2a <- RsquareAdj(algae.pcnm.rda2)$adj.r.squared) 
anova.cca(algae.pcnm.rda2)
axes.test <- anova.cca(algae.pcnm.rda2, by="axis")
(nb.ax <- length(which(axes.test[,4] <= 0.05))) # Number of significant axes
