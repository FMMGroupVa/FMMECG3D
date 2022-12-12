
path <- getwd()
source(paste0(path,"/FMM_ECG3D_Codes/auxMultiFMM_ECG.R"), chdir = TRUE)
source(paste0(path,"/requiredFunctionsPreprocessing_v4.1.R"))
source(paste0(path,"/runPreprocessing_v4.1.R"))

# Patients #1 and #2 from PTB-XL ECG Database:
load("Examples12leads.RData")

preprocessedOutput <- givePreprocessing_app(dataIn = patient1, freqHz = 500)
# preprocessedOutput[[1]] # Beat ID | R location | beat start | beat end | status of leads {1=OK, 0=bad}
# preprocessedOutput[[2]] # Preprocessed ECG

selectedBeat <- 5
beatRefsPatient1 <- preprocessedOutput[[1]][preprocessedOutput[[1]][,"beat_id"] == selectedBeat,
                                            c("iniRef", "finRef", "annoRef")]
ecgDataPatient1 <- preprocessedOutput[[2]]
selectedBeatData <- ecgDataPatient1[as.numeric(beatRefsPatient1[1]):as.numeric(beatRefsPatient1[2]),]

paramsPerLeadPre <- fitMultiFMM_ECG(vDataMatrix = selectedBeatData, 
                                    annotation = as.numeric(beatRefsPatient1[3] - beatRefsPatient1[1]), 
                                    maxIter = 10, parallelize = TRUE)



