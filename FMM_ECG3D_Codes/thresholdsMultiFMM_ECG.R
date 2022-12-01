relevantLeads<-c("I","II","V2","V5")

#### Noisy Set ####
noisyAlpha<-c(pi-0.25,pi+0.25)
noisyMinOmega<-0.008; noisyMaxOmega<-0.5
minExtraVar<-0.15

#### R Wave ####
rLeftDistanceToAlphaN<-c(0.15,0.2,0.3); rRightDistanceToAlphaN<-0.2
rMaxOmega<-c(0.25,0.5); rMinVarMaxOmega<-0.25

## R Confusions
lowVarR_P<-0.15
lowVarR_Q<-c(0.15,0.25); minVarQ_to_R<-0.12
lowVarR_S<-0.25
highVarQ<-0.15
highVarS<-0.15
minVarU_to_R<-0.05

#### P Wave ####
pMaxOmega<-c(0.35, 0.40); pMinVarMaxOmega<-0.10
minDistPR<-c(0.065, 0.25); maxDistPR<-c(1.5, 1.8)
minDistPN<-0.065

## Typical P Wave
typMinDistPR<-c(0.25, 0.1)
typMinOmegaP<-0.05; typMaxOmegaP<-0.25

## P Comprobation
minExtraVarP<-0.15; maxVarDecrease<-0.10; maxRelativeVarDecrease<-0.75

#### Q Wave ####
qMaxOmega<-c(0.25,0.30); qMinVarMaxOmega<-0.15
maxDistQR<-c(0.25,0.40)

#### S Wave ####
sMaxOmega<-c(0.30,0.40); sMinVarMaxOmega<-0.25
maxDistRS<-0.5; maxDistSN<-0.5
omegaS_or_T<-0.08

#### T Wave ####
minDistRT<-c(0.30,0.40)
minDistNT<-0.3