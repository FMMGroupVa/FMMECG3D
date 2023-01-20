
path <- getwd()
source(paste0(path,"/requiredFunctionsPreprocessing_v4.1.R")) # Load preprocessing functions
source(paste0(path,"/FMM_ECG3D_Codes/auxMultiFMM_ECG.R"), chdir = TRUE) # Load analysis functions

# Patients #1 from PTB-XL Database:
load("Patient1_PTBXL.RData") # 10 seconds of ECG records across leads (5000X12) 


###############################
#	Preprocessing		#
###############################


preprocessedOutput <- givePreprocessing_git(dataIn = Patient1_PTBXL, freqHz = 500) # For sampling frequency see database info
# preprocessedOutput[[1]] # Beat ID | QRS location | beat start | beat end | logical matrix indicating beat processing by lead {1=Yes, 0=No}
# preprocessedOutput[[2]] # Preprocessed ECG data (5000X12)


###############################
#	Data Analysis		#
###############################

# Beat selection
selectedBeat <- 5

# Data heartbeat adquisition for from processed data  
beatRefsPatient1 <- preprocessedOutput[[1]][preprocessedOutput[[1]][,"beat_id"] == selectedBeat,
                                            c("iniRef", "finRef", "annoRef")]
selectedBeatData <- preprocessedOutput[[2]][as.numeric(beatRefsPatient1[1]):as.numeric(beatRefsPatient1[2]),]

# 3D FMMM_ecg
paramsPerLeadPre <- fitMultiFMM_ECG(vDataMatrix = selectedBeatData,
                                    annotation = as.numeric(beatRefsPatient1[3] - beatRefsPatient1[1] + 1))
# paramsPerLeadPre # List named by leads with FMM parameter estimates
# Note: Plot with 3D FMM_ecg prediction is shown by screen. Use extra arguments of paramsPerLeadPre to save locally (advanced)



