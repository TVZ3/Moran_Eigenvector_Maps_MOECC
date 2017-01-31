#Load libraries
library(ape)
library(spdep)
library(ade4)
library(vegan) # contains the pcnm function
library(packfor)
library(MASS)
library(dplyr)
library(tidyr)

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
algae.h.det <- as.data.frame(algae.h.det)

colnames(algae.h.det)

#Select a subset of variables by column name
Diatoms <- select(algae.h.det, contains("Diat"))
Chlorophytes <- select(algae.h.det, contains("Chlo"))
Chrysophytes <- select(algae.h.det, contains("Chry"))
Cryptophytes <- select(algae.h.det, contains("Cryp"))
Cyanobactera <- select(algae.h.det, contains("Cyan"))
Dinoflagellates <- select(algae.h.det, contains("Dino"))

# Construction of spatial variables at all relevant scales
# PCNM analysis of data using Vegan's "pcnm" code
algae.pcnm.vegan <- pcnm(dist(c))
algae.pcnm <- as.data.frame(algae.pcnm.vegan$vectors)
nb.ev <- length(which(algae.pcnm.vegan$values > 0.0000001))

# Run a PCNM on the detrended species data and determine if analysis is significant
Diatoms.pcnm <- rda(Diatoms, algae.pcnm)
Chlorophytes.pcnm <- rda(Chlorophytes, algae.pcnm)
Chrysophytes.pcnm <- rda(Chrysophytes, algae.pcnm)
Cryptophytes.pcnm <- rda(Cryptophytes, algae.pcnm)
Cyanobactera.pcnm <- rda(Cyanobactera, algae.pcnm)
Dinoflagellates.pcnm <- rda(Dinoflagellates, algae.pcnm)

anova.cca(Diatoms.pcnm) # Not significant
anova.cca(Chlorophytes.pcnm) # Significant
anova.cca(Chrysophytes.pcnm) # Significant
anova.cca(Cryptophytes.pcnm) # Significant
anova.cca(Cyanobactera.pcnm) # Not significant
anova.cca(Dinoflagellates.pcnm) # Significant

# Compute Adjusted R2 value and run forward selection of the PCNM variables on the species data.  You do this because the analysis above was significant
(Chlo.r2a <- RsquareAdj(Chlorophytes.pcnm)$adj.r.squared)
(Chry.r2a <- RsquareAdj(Chrysophytes.pcnm)$adj.r.squared)
(Cryp.r2a <- RsquareAdj(Cryptophytes.pcnm)$adj.r.squared)
(Dino.r2a <- RsquareAdj(Dinoflagellates.pcnm)$adj.r.squared)

#Run forward selection to determine which MEMs are the most important predictors
(Chlorophytes.pcnm.fwd <- forward.sel(Chlorophytes, as.matrix(algae.pcnm), adjR2thresh = Chlo.r2a))

(Chrysophytes.pcnm.fwd <- forward.sel(Chrysophytes, as.matrix(algae.pcnm), adjR2thresh = Chry.r2a))

(Cryptophytes.pcnm.fwd <- forward.sel(Cryptophytes, as.matrix(algae.pcnm), adjR2thresh = Cryp.r2a))

(Dinoflagellates.pcnm.fwd <- forward.sel(Dinoflagellates, as.matrix(algae.pcnm), adjR2thresh = Dino.r2a))

### Note: In the above selection procedure, sometimes 4 PCNMs come out as important because the fourth is either just under or just over the significance threshold (0.05)

## Gives the number of SIGNIFICANT PCMNs
(Chlo.sig.pcnm <- nrow(Chlorophytes.pcnm.fwd)) # Chlorophytes
(Chry.sig.pcnm <- nrow(Chrysophytes.pcnm.fwd)) # Chrysophytes
(Cryp.sig.pcnm <- nrow(Cryptophytes.pcnm.fwd)) # Cryptophytes
(Dino.sig.pcnm <- nrow(Dinoflagellates.pcnm.fwd)) # Dinoflagellates

# Plots of significant PCNMs over our thirty year time period
par(mfrow=c(3,3))
plot(algae.pcnm$PCNM2, type="o", col="black", ylab = "PCNM 2") # Chlorophytes, Chrysophytes, Dinoflagellates
plot(algae.pcnm$PCNM3, type="o", col="black", ylab = "PCNM 3") # Cryptophytes
plot(algae.pcnm$PCNM4, type="o", col="black", ylab = "PCNM 4") # Chlorophytes, Chrysophytes, Cryptophytes
plot(algae.pcnm$PCNM5, type="o", col="black", ylab = "PCNM 5") # Chlorophytes
plot(algae.pcnm$PCNM8, type="o", col="black", ylab = "PCNM 8") # Chrysophytes, Cryptophytes
plot(algae.pcnm$PCNM9, type="o", col="black", ylab = "PCNM 9") # Chrysophytes
plot(algae.pcnm$PCNM17, type="o", col="black", ylab = "PCNM 17") # Chlorophytes


# Write the significant PCNMs for each species to a new object and save the file
colnames(algae.pcnm)
pcnm.significant <- algae.pcnm[,c(2, 3, 4, 5, 8, 9, 17)]
write.csv(pcnm.significant, "MEMs_Individual_species.csv", row.names = FALSE)
