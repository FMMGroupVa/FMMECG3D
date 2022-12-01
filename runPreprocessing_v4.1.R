source("requiredFunctionsPreprocessing_v4.1.R")
require(stringr)
require(rsleep)

##############################################
#		INPUTS
##############################################
#dataIn: Matrix nObs x 12. Each column a lead in the following order: I,II,III,AVL,AVR,AVF,V1,V2,V3,V4,V5,V6
#	   Colnames must be empty or be THE VECTOR: c("I","II","III","AVL","AVR","AVF","V1","V2","V3","V4","V5","V6")
#freqHz: is a numeric value indicating the sampling frequency in Hz


##############################################
#		OUTPUTS
##############################################
#List of two elements.
#1.- The first element is a matrix. Each row corresponds to an analyzed heartbeart. The columns are described as follows: 
#	beat_id: numeric value indicating the position of the corresponding t_QRS with regard to the median annotation list.
#	annoRef: numeric value indicating the median t_QRS position of the beat
#	iniRef: numeric value indicating the position at which the heartbeat starts 
#	finRef: numeric value indicating the position at which the heartbeat finishes
#	I: logical value (1/0) indicating if the beartbeat is processed for lead I 
#	II: logical value (1/0) indicating if the beartbeat is processed for lead II
#	V1: logical value (1/0) indicating if the beartbeat is processed for lead V1
#	V2: logical value (1/0) indicating if the beartbeat is processed for lead V2
#	V3: logical value (1/0) indicating if the beartbeat is processed for lead V3
#	V4: logical value (1/0) indicating if the beartbeat is processed for lead V4
#	V5: logical value (1/0) indicating if the beartbeat is processed for lead V5
#	V6: logical value (1/0) indicating if the beartbeat is processed for lead V6
#2.- The second element is  a matrix nObs x 12 contained baseline corrected data

givePreprocessing_app<-function(dataIn=dataIn,freqHz=freqHz){#
	dataIn<-t(dataIn)
	leadName<-c("I","II","III","AVL","AVR","AVF","V1","V2","V3","V4","V5","V6")
	less4<-FALSE
	uniqueLead<-FALSE
	# Order leads and correct lead names
	dataIn<-dataIn[match(tolower(rownames(dataIn)), tolower(leadName)),]
	rownames(dataIn)<-leadName
	
	numNonNA<-sum((apply(dataIn,1,sum,na.rm=TRUE)!=0))
	#if(sum(is.na(dataIn))==dim(dataIn)[1]*dim(dataIn)[2]) stop("Empty data matrix")
	# if(sum(is.na(dataIn[1,]))==dim(dataIn)[2] & sum(is.na(dataIn[2,]))==dim(dataIn)[2] & 
	# 		sum(is.na(dataIn[8,]))==dim(dataIn)[2]  & sum(is.na(dataIn[11,]))==dim(dataIn)[2] ) stop("Please, include data for leads I, II, V2 or V5")
	if(sum((apply(dataIn,1,sum,na.rm=TRUE)!=0))<4) less4<-TRUE
	if((sum(!is.na(dataIn[1,]))==dim(dataIn)[2] & sum(is.na(dataIn[-1,]))==11*dim(dataIn)[2]) |
		(sum(!is.na(dataIn[2,]))==dim(dataIn)[2] & sum(is.na(dataIn[-2,]))==11*dim(dataIn)[2]) |
		(sum(!is.na(dataIn[8,]))==dim(dataIn)[2] & sum(is.na(dataIn[-8,]))==11*dim(dataIn)[2]) |
		(sum(!is.na(dataIn[11,]))==dim(dataIn)[2] & sum(is.na(dataIn[-11,]))==11*dim(dataIn)[2])) uniqueLead<-TRUE

	pre_I<-leadPrePanTom_app(dataIn[1,],freqHz)#leadPrePanTomYolanda POR leadPrePanTom y quito la dependdencia de paci
	pre_II<-leadPrePanTom_app(dataIn[2,],freqHz)
	pre_III<-leadPrePanTom_app(dataIn[3,],freqHz)
	pre_AVL<-leadPrePanTom_app(dataIn[4,],freqHz)
	pre_AVR<-leadPrePanTom_app(dataIn[5,],freqHz)
	pre_AVF<-leadPrePanTom_app(dataIn[6,],freqHz)
	pre_V1<-leadPrePanTom_app(dataIn[7,],freqHz)
	pre_V2<-leadPrePanTom_app(dataIn[8,],freqHz)
	pre_V3<-leadPrePanTom_app(dataIn[9,],freqHz)
	pre_V4<-leadPrePanTom_app(dataIn[10,],freqHz)
	pre_V5<-leadPrePanTom_app(dataIn[11,],freqHz)
	pre_V6<-leadPrePanTom_app(dataIn[12,],freqHz)





	############################
	#	Multi
	############################

	
	lead8<-1:nrow(dataIn)
	
		
	annos<-c(pre_I[[1]],pre_II[[1]],pre_III[[1]],pre_AVL[[1]],pre_AVR[[1]],pre_AVF[[1]],
				pre_V1[[1]],pre_V2[[1]],pre_V3[[1]],pre_V4[[1]],pre_V5[[1]],pre_V6[[1]])
	leads<-rep(1:nrow(dataIn),times=c(length(pre_I[[1]]),length(pre_II[[1]]),length(pre_III[[1]]),
					length(pre_AVL[[1]]),length(pre_AVR[[1]]),length(pre_AVF[[1]]),
					length(pre_V1[[1]]),length(pre_V2[[1]]),length(pre_V3[[1]]),
					length(pre_V4[[1]]),length(pre_V5[[1]]),length(pre_V6[[1]])))
	annoLeads<-rbind(annos,leads)
	annoLeads<-annoLeads[,order(annoLeads[1,])]
	#aplicamos multi a las anotaciones		
	
	if(!less4){
		getAnnosMedian<-giveAnnoMed_app(annoLeads,ncol(dataIn),freqHz)
		annosPac<-getAnnosMedian[[1]]
		leadsPac<-getAnnosMedian[[2]]
		segPac<-getAnnosMedian[[3]]
	}else{
		if(!uniqueLead){
			getAnnosMedian<-giveAnnoMed2_app(annoLeads,ncol(dataIn),freqHz)
			annosPac<-getAnnosMedian[[1]]
			leadsPac<-getAnnosMedian[[2]]
			segPac<-getAnnosMedian[[3]]
		}else{
			if(sum(!is.na(dataIn[1,]))==dim(dataIn)[2]){
				annosPac<-pre_I[[1]]
				leadsPac<-rep(1,length(annosPac))
				segPac<-pre_I[[2]]
			}
			if(sum(!is.na(dataIn[2,]))==dim(dataIn)[2]){
				annosPac<-pre_II[[1]]
				leadsPac<-rep(2,length(annosPac))
				segPac<-pre_II[[2]]
			}
			if(sum(!is.na(dataIn[8,]))==dim(dataIn)[2]){
				annosPac<-pre_V2[[1]]
				leadsPac<-rep(8,length(annosPac))
				segPac<-pre_V2[[2]]
			}
			if(sum(!is.na(dataIn[11,]))==dim(dataIn)[2]){
				annosPac<-pre_V5[[1]]
				leadsPac<-rep(11,length(annosPac))
				segPac<-pre_V5[[2]]
			}
		}
	}


	#eliminar primero y último
	#eliminar artefactos
	#eliminar 35-45
	#contar tendnecia

	pos_I<-leadPreMulti_app(pre_I,annosPac,segPac,freqHz)
	pos_II<-leadPreMulti_app(pre_II,annosPac,segPac,freqHz)
	pos_III<-leadPreMulti_app(pre_III,annosPac,segPac,freqHz)
	pos_AVL<-leadPreMulti_app(pre_AVL,annosPac,segPac,freqHz)
	pos_AVR<-leadPreMulti_app(pre_AVR,annosPac,segPac,freqHz)
	pos_AVF<-leadPreMulti_app(pre_AVF,annosPac,segPac,freqHz)
	pos_V1<-leadPreMulti_app(pre_V1,annosPac,segPac,freqHz)
	pos_V2<-leadPreMulti_app(pre_V2,annosPac,segPac,freqHz)
	pos_V3<-leadPreMulti_app(pre_V3,annosPac,segPac,freqHz)
	pos_V4<-leadPreMulti_app(pre_V4,annosPac,segPac,freqHz)
	pos_V5<-leadPreMulti_app(pre_V5,annosPac,segPac,freqHz)
	pos_V6<-leadPreMulti_app(pre_V6,annosPac,segPac,freqHz)


	#Cambiamos a 8 derivaciones
	
	finales<-list(pos_I,pos_II,pos_V1,pos_V2,pos_V3,pos_V4,pos_V5,pos_V6)
	leadChar<-c("I","II","V1","V2","V3","V4","V5","V6")



	lead8<-c(1,2,7:12)
	

	#CAMBIO: addLength, sera TRUE si hay que eliminarlo 
	##################################################
	#	Elimino si hay menos de tres
	##################################################
	addLength<-FALSE
	lRef<-length(annosPac)
	if(lRef<3)addLength<-TRUE


	
	#ANOTACIONES PARA LAS LEADS I, II
	if((!uniqueLead & (sum(!is.na(dataIn[1,]))==dim(dataIn)[2] | sum(!is.na(dataIn[2,]))==dim(dataIn)[2])) | 
			(uniqueLead & (sum(!is.na(dataIn[1,]))==dim(dataIn)[2] | sum(!is.na(dataIn[2,]))==dim(dataIn)[2]))){
		annosA<-c(pre_I[[1]],pre_II[[1]])
		leadsA<-rep(1:2,times=c(length(pre_I[[1]]),length(pre_II[[1]])))
		annoLeadsA<-rbind(annosA,leadsA)
		annoLeadsA<-annoLeadsA[,order(annoLeadsA[1,])]
	}else{
		if(!uniqueLead){
			annosA<-c(pre_V1[[1]],pre_V2[[1]],pre_V3[[1]],pre_V4[[1]],pre_V5[[1]],pre_V6[[1]])
			leadsA<-rep(7:12,times=c(length(pre_V1[[1]]),length(pre_V2[[1]]),length(pre_V3[[1]]),length(pre_V4[[1]]),length(pre_V5[[1]]),length(pre_V6[[1]])))
			annoLeadsA<-rbind(annosA,leadsA)
			annoLeadsA<-annoLeadsA[,order(annoLeadsA[1,])]
		}else{
			if(sum(!is.na(dataIn[8,]))==dim(dataIn)[2]){
				annosA<-c(pre_V2[[1]])
				leadsA<-rep(8,times=c(length(pre_V2[[1]])))
				annoLeadsA<-rbind(annosA,leadsA)
				annoLeadsA<-annoLeadsA[,order(annoLeadsA[1,])]
			}else{
				annosA<-c(pre_V5[[1]])
				leadsA<-rep(11,times=c(length(pre_V5[[1]])))
				annoLeadsA<-rbind(annosA,leadsA)
				annoLeadsA<-annoLeadsA[,order(annoLeadsA[1,])]
			}
		}
	}


	#ANOTACIONES PARA LAS LEADS V1,V2,V3,V4,V5,V6
	if((!uniqueLead & (sum(!is.na(dataIn[8,]))==dim(dataIn)[2] | sum(!is.na(dataIn[11,]))==dim(dataIn)[2])) | 
			(uniqueLead & (sum(!is.na(dataIn[8,]))==dim(dataIn)[2] | sum(!is.na(dataIn[11,]))==dim(dataIn)[2]))){
		annosB<-c(pre_V1[[1]],pre_V2[[1]],pre_V3[[1]],pre_V4[[1]],pre_V5[[1]])
		leadsB<-rep(7:11,times=c(length(pre_V1[[1]]),length(pre_V2[[1]]),length(pre_V3[[1]]),length(pre_V4[[1]]),length(pre_V5[[1]])))
		annoLeadsB<-rbind(annosB,leadsB)
		annoLeadsB<-annoLeadsB[,order(annoLeadsB[1,])]
	}else{
		if(!uniqueLead){
			annosB<-c(pre_I[[1]],pre_I[[1]])
			leadsB<-rep(1:2,times=c(length(pre_I[[1]]),length(pre_II[[1]])))
			annoLeadsB<-rbind(annosB,leadsB)
			annoLeadsB<-annoLeadsB[,order(annoLeadsB[1,])]
		}else{
			if(sum(!is.na(dataIn[1,]))==dim(dataIn)[2]){
				annosB<-c(pre_I[[1]])
				leadsB<-rep(1,times=c(length(pre_I[[1]])))
				annoLeadsB<-rbind(annosB,leadsB)
				annoLeadsB<-annoLeadsB[,order(annoLeadsB[1,])]
			}else{
				annosB<-c(pre_II[[1]])
				leadsB<-rep(2,times=c(length(pre_II[[1]])))
				annoLeadsB<-rbind(annosB,leadsB)
				annoLeadsB<-annoLeadsB[,order(annoLeadsB[1,])]
			}
		}
	}

	#BEFORE V3
	#annosB<-c(pre_V1[[1]],pre_V2[[1]],pre_V3[[1]],pre_V4[[1]])
	#leadsB<-rep(7:10,times=c(length(pre_V1[[1]]),length(pre_V2[[1]]),length(pre_V3[[1]]),length(pre_V4[[1]])))
	#annoLeadsB<-rbind(annosB,leadsB)
	#annoLeadsB<-annoLeadsB[,order(annoLeadsB[1,])]
	
	medI_II<-giveAnnoMed2_app(annoLeadsA,ncol(dataIn),freqHz)
	medV1_V6<-giveAnnoMed2_app(annoLeadsB,ncol(dataIn),freqHz)

	#CAMBIO: addPace, sera TRUE si hay que eliminarlo 
	##################################################
	#	Elimino si hay menos de tres
	##################################################
	addPace<-FALSE
	vAux<-sort(c(medI_II[[1]],medV1_V6[[1]]))
	vAux2<-c(vAux,vAux)
	if(length(vAux)>0){
		diffe<-vAux2[2:(length(vAux))]-vAux[1:(length(vAux)-1)]
		limit1<-103*freqHz/1000
		if(sum(diffe>limit1)>0.5*length(vAux) & min(diffe,na.rm=TRUE)>=limit1)addPace<-TRUE
	}

	#CAMBIO: elimina será TRUE si lo es addLength o addPace 
	elimina<-FALSE
	if(addLength | addPace)elimina<-TRUE

	#############################################################################
	#	Actualizacion de finales para los de menos 3 latidos o paces
	#############################################################################
	if(elimina){
		for(j in 1:length(lead8)){
			finales[[j]][[1]]<-NA
			finales[[j]][[2]]<-c(NA,NA)
			finales[[j]][[3]]<-NA
			ff<-which(elimina[i]==finales[[j]][[4]][,1],arr.ind=TRUE)
			finales[[j]][[4]][ff,3]<-rep(10,length(ff))
		}
	}


	#CAMBIO: paciOut prescindo de ello y lo cambio por !elimina
	#PARA ELIMINAR LOS DE AMPLITUD MUY VARIABLE 			
	addAnnoAll<-c()
	cvDivision<-c()
	addAnnoWrong<-c()
	if(!elimina & sum(c(is.null(finales[[1]][[5]]),is.null(finales[[2]][[5]]),
					is.null(finales[[3]][[5]]),is.null(finales[[4]][[5]]),
					is.null(finales[[5]][[5]]),is.null(finales[[6]][[5]]),
					is.null(finales[[7]][[5]]),is.null(finales[[8]][[5]])))==0   ){
		ampDivisionMat<-matrix(NA,length(lead8),5)
		addAnnoWrong<-c()
		for(j in 1:length(lead8)){
			divisiones<-seq(1,length(finales[[j]][[5]])+1,by=ncol(dataIn)/5)#Theshold
			ampDivision<-c()
			dataDivi<-finales[[j]][[5]]
			for(k in 1:(length(divisiones)-1)){
				dataDivision<-dataDivi[divisiones[k]:(divisiones[k+1]-1)]
				ampDivision<-c(ampDivision,quantile(dataDivision,0.9)-quantile(dataDivision,0.1))
			}
			cvDivision[j]<-sd(ampDivision)/abs(mean(ampDivision))
			ampDivisionMat[j,]<-ampDivision
		}
		medVal<-quantile(cvDivision,probs=c(0.4))#Threshold
		if(sum(cvDivision>1.75*medVal)>0 & max(cvDivision)>0.4){#Threshold
			sele<-which(cvDivision>1.75*medVal & (cvDivision)>0.4,arr.ind=TRUE)#Threshold
			for(kk in 1:length(sele)){
				col<-which(ampDivisionMat[sele[kk],]>2*median(ampDivisionMat[sele[kk],]),arr.ind=TRUE)
				if(length(col)>0){
					for(kkk in 1:length(col)){
						addAnnoWrong<-rbind(addAnnoWrong,c(sele[kk],divisiones[col[kkk]],(divisiones[col[kkk]+1]-1)))
					}
				}else{
					col<-NA
				}
			}
		}
		if(length(addAnnoWrong)>0)if(length(col)>0)addAnnoAll<-rbind(addAnnoAll,rbind(cbind(rep(1234,nrow(addAnnoWrong)),addAnnoWrong)))
	}

	#############################################################################
	#	Actualizacion de finales para los de gran amplitud:77
	#############################################################################

	if(length(addAnnoWrong)>0){
		sss<-0
		#recorreria addAnnoAll y pondria NA si se encuentra alguna anotacion ahi
		for(i in 1:nrow(addAnnoAll)){
			derivation<-addAnnoAll[i,2]
			segOutIni<-addAnnoAll[i,3]
			segOutFin<-addAnnoAll[i,4]
			if(sum(!is.na(match(finales[[derivation]][[1]],segOutIni:segOutFin)))>0){
				aux<-match(finales[[derivation]][[1]],segOutIni:segOutFin)
			
				annosFound<-c(segOutIni:segOutFin)[aux[!is.na(aux)]]
				posOut<-c()
				for(k in 1:length(annosFound)){
					posOut<-c(posOut,which(annosFound[k]==finales[[derivation]][[1]]))
				}
				for(k in 1:length(posOut)){
					sss<-sss+1
					finales[[derivation]][[1]][posOut[k]]<-NA
					finales[[derivation]][[2]][posOut[k],]<-c(NA,NA)
					finales[[derivation]][[3]][posOut[k]]<-NA
					ff<-which(posOut[k]==finales[[derivation]][[4]][,2])
					finales[[derivation]][[4]][ff,3]<-77
				}
			}
		}
	}



	#para eliminar los de rango o problemas de conduccion 
	#CAMBIO: addStrange, addStrange3, addCeros, addRangos lo cambio a TRUE/FALSE
	#CAMBIO: strange, strange2, ceros, rango lo cambio a vector NA de longitud las leads que sera TRUE en aquella que haya cambio
	addStrange<-FALSE;addStrange3<-FALSE;addCeros<-FALSE;addRangos<-FALSE;
	strange<-rep(NA,length(lead8));strange3<-rep(NA,length(lead8));ceros<-rep(NA,length(lead8));rangos<-rep(NA,length(lead8))
	options(warn=-1)
	for(j in 1:length(lead8)){
		#COMO EN EL 4410
		ss7<-0
		if(sum(!is.na(finales[[j]][[1]]))==1){##NEW##
			otherLeads<-(1:8)[-j]
			for(k in 1:length(otherLeads)){
				if(sum(!is.na(finales[[otherLeads[k]]][[1]]))>=3)ss7<-ss7+1
			}
		}
		if(ss7==7){
			addStrange<-TRUE
			strange[j]<-TRUE
		}
			
		ss11<-0;ss22<-0
		cond1<-FALSE
		if(sum(!is.na(finales[[j]][[1]]))>0){
			if(!is.na(finales[[j]][[1]][1]))cond1<-finales[[j]][[1]][1]==annosPac[1]
		}
		if(!is.na(cond1)){
			if(cond1 ){
				ss11<-ss11+1
				otherLeads<-(1:8)[-j]
				sigue<-TRUE;
				for(k in 1:length(otherLeads)){
					sigue2<-TRUE
					if(sum(!is.na(finales[[otherLeads[k]]][[1]]))>0 ){
						if( !is.na(finales[[otherLeads[k]]][[1]][1])){
							if(finales[[otherLeads[k]]][[1]][1]==annosPac[1] ){#&
								ss11<-ss11+1
								sigue2<-FALSE
							}
							if(length(finales[[otherLeads[k]]][[1]])>1){
								if(sum(annosPac[2]==finales[[otherLeads[k]]][[1]],na.rm=TRUE)>0){
										ss22<-ss22+1
										sigue2<-FALSE
								}
									
								if(sum(annosPac[2]==finales[[j]][[1]],na.rm=TRUE)>0){
									ss22<-ss22+1
									sigue<-FALSE
									sigue2<-FALSE
								}
									
							}
						}
					}
				}
			}
		}
		if(ss11==1 & ss22==0){
			addStrange3<-TRUE
			strange3[j]<-TRUE
		}
		if(sum(!is.na(finales[[j]][[1]]))==0){#cuando todos son NA, es redundante?
			addCero<-TRUE
			ceros[j]<-TRUE
		}
		rangos[j]<-max(finales[[j]][[5]])-min(finales[[j]][[5]])##NEW##10565
	}
	if(max(rangos)>15000)addRangos<-TRUE#Threshold
	options(warn=0)

	#############################################################################
	#	Actualizacion de finales para los de gran amplitud:77
	#############################################################################

	for(j in 1:length(lead8)){
		if(!is.na(ceros[j])){
			if(ceros[j]){
				finales[[j]][[1]]<-NA
				finales[[j]][[2]]<-c(NA,NA)
				finales[[j]][[3]]<-NA
		
				#recodificacion#para ceros codigo 11
				finales[[j]][[4]][,3]<-rep(11,length(finales[[j]][[4]][,3]))

			}
		}else{
			if(rangos[j]>15000 | !is.na(strange[j])   | !is.na(strange3[j])){
				finales[[j]][[1]]<-NA
				finales[[j]][[2]]<-c(NA,NA)
				finales[[j]][[3]]<-NA
	
				#recodificacion#
				if(rangos[j]>15000){#codigo 22
					finales[[j]][[4]][,3]<-rep(22,length(finales[[j]][[4]][,3]))
				}
				if(!is.na(strange[j])){#codigo 44
					finales[[j]][[4]][,3]<-rep(44,length(finales[[j]][[4]][,3]))
				}
				if(!is.na(strange3[j])){#codigo 333
					finales[[j]][[4]][,3]<-rep(333,length(finales[[j]][[4]][,3]))
				}
			}
		} 
	}


	###############################################################
	#	Añado nueva condicion apra eliminar latidos mas cortos de2/3 de la longitu mediana de los que quedan y los actualizo en finales
	##############################################################


	for(j in 1:length(lead8)){
		if(sum(finales[[j]][[2]],na.rm=TRUE)>2){
			if(length(finales[[j]][[2]])>2){
				rrs<-finales[[j]][[2]][,2]-finales[[j]][[2]][,1]
			}else{
				rrs<-finales[[j]][[2]][2]-finales[[j]][[2]][1]
			}
			medS<-median(rrs,na.rm=TRUE)
			if(sum(rrs<2/3*medS,na.rm=TRUE)>0){
				posi<-which(!is.na(rrs) & rrs<2/3*medS,arr.ind=TRUE)
				finales[[j]][[1]][posi]<-rep(NA,length(posi))
				for(k in 1:length(posi)){
					finales[[j]][[2]][posi[k],]<-rep(NA,2)
					ff<-which(posi[k]==finales[[j]][[4]][,2],arr.ind=TRUE)
					finales[[j]][[4]][ff,3]<-100
				}
				finales[[j]][[3]][posi]<-rep(NA,length(posi))
			}
		}
	}

	
	##################################################################################
	#		Cuantas derivaciones con menos de tres latidos 
	##################################################################################

	if( sum(length(finales[[1]][[1]]))){
		s<-0
		wrongDeri<-0
		addDeri<-c()
		for(j in 1:length(lead8)){
			if(sum(finales[[j]][[4]][,3]!=11)>0){#distinto del codigo 11 
				
				if(( sum(finales[[j]][[4]][,3]==0))<3){#que los latidos de ese pac en esa der sin codigo 2,5,3,77,22,44,88,333
					s<-s+1
					addDeri<-c(addDeri,j)
					if(j==1 | j==2 | j==8 | j==11)wrongDeri<-wrongDeri+1
				}
			}
		}
		menos3<-s
		menos3Der<-addDeri
		menos3DerImp<-wrongDeri
	}else{
		menos3<-0
		menos3Der<-0
		menos3DerImp<-0
	}


	#Si hay 4 o más derivaciones con menos de 3 latidos addPaci es TRUE 
	addPaci<-FALSE																	
	if(menos3>3)addPaci<-TRUE

	#Numero de latidos analizdos por paciente
	suma<-0;sumAnno<-0
	for(j in 1:length(lead8)){
		suma<-suma+sum(finales[[j]][[4]][,3]==0)
	}
	sumAnno<-suma

	#Eliminamos o 4 o mas derivaciones con menos de 3 latidos o menos de 24 latidos en total analizados
	eliminaPac<-FALSE
	#if((sumAnno>0 & (sumAnno<24 & !less4)) | addPaci) eliminaPac<-TRUE
	if((sumAnno>0 & (sumAnno<0.2*ncol(dataIn)/freqHz*numNonNA & !less4)) | addPaci) eliminaPac<-TRUE


	#¿se elimina el paciente?
	eliminaPacAll<-eliminaPac


	
	#voy a sacr la forma en leer los latidos
	uno<-c();dos<-c();tres<-c();cuatro<-c()
	if(length(annosPac)>1){
		uno<-cbind(rep(1234,length(annosPac)),1:length(annosPac),annosPac,segPac)
	}else{
		uno<-c(rep(1234,length(annosPac)),1:length(annosPac),annosPac,segPac)
	}
	dos<-matrix(0,length(annosPac),length(lead8))
	for(j in 1:length(annosPac)){
		for(k in 1:length(lead8)){
			if(!is.na(match(annosPac[j],finales[[k]][[1]])))dos[j,k]<-1
		}
	}
	if(length(annosPac)>1){
		tres<-rbind(tres,cbind(uno,dos))
	}else{
		tres<-rbind(tres,c(uno,dos))
	}

	colnames(tres)<-c("ecg_id","beat_id","annoRef","iniRef","finRef",leadChar)
	if(nrow(tres)>1){
		if(sum(apply(tres[,leadChar],1,sum)==0)>0){
			cuatro<-tres[-which(apply(tres[,leadChar],1,sum)==0,arr.ind=TRUE),]#elimino los que suman cero's
		}else{
			cuatro<-tres
		}
	}else{
		if(sum(tres[,leadChar])==0){
			cuatro<-c()
		}else{
			cuatro<-tres
		}

	}
	if(length(cuatro)==13){
		cuatro<-t(as.matrix(cuatro))
		final<-cbind(cuatro,cuatro[,"I"]+cuatro[,"II"]+cuatro[,"V2"]+cuatro[,"V5"])
	}else{
		if(length(cuatro)==0){
			final<-c()
		}else{
			final<-cbind(cuatro,apply(cuatro[,c("I","II","V2","V5")],1,sum))
		}
	}
	

	# Dos grupos de anotaciones comparandolas con las iniciales
	iniciales<-list(pre_I,pre_II,pre_V1,pre_V2,pre_V3,pre_V4,pre_V5,pre_V6)
	tabCheck<-c()
	filasPacDelete<-c()
	addFilasDelete<-c()

	if(length(final)>0){
	for(k in 1:nrow(final)){
		annosF<-final[k,"annoRef"]
		annosI<-c()
		for(j in 1:length(leadChar)){
			if(sum((final[k,"annoRef"]+limit1)>=iniciales[[j]][[1]] & 
			(final[k,"annoRef"]-limit1)<=iniciales[[j]][[1]])==1){
				pos<-which((final[k,"annoRef"]+limit1)>=iniciales[[j]][[1]] & 
					(final[k,"annoRef"]-limit1)<=iniciales[[j]][[1]])
				if(final[k,leadChar[j]]==1){
					annosI<-c(annosI,iniciales[[j]][[1]][pos])
				}else{
					annosI<-c(annosI,NA)
				}
			}else{
				annosI<-c(annosI,NA)
			}

		}
		tabCheck<-rbind(tabCheck,c(final[k,"ecg_id"],final[k,"beat_id"],
			annosF,annosI,median(annosI,na.rm=TRUE),abs(annosF-median(annosI,na.rm=TRUE))))
	}
	colnames(tabCheck)<-c("ecg_id","beat_id","annoFinal",paste0("annoIni_",leadChar),"MedianIni","Dif")
	#lo decido con los percentiles 99.9 y para que no se eliminen pacientes
	toDelete<-tabCheck[tabCheck[,"Dif"]>12 & !is.na(tabCheck[,"Dif"]),c("ecg_id","beat_id")]#Thershold

	if(length(toDelete)>0){

		#en pacientes cuantos se eliminarian 
		cuenta12<-tabCheck[tabCheck[,"Dif"]>12 & !is.na(tabCheck[,"Dif"]),c("ecg_id")]
		pac12<-as.numeric(names(table(cuenta12)))
		nLat12<-c()
		filas<-which(pac12==final[,"ecg_id"],arr.ind=TRUE)
		resta<-(table(cuenta12))[match(as.character(pac12),names(table(cuenta12)))]
		nLat12<-c(nLat12,length(filas)-resta)
		tabla12<-cbind(pac12,nLat12)
	

		filasPacDelete<-c()
		#si hubiera que eliminar paciente completo lo haria asi 
		if(min(nLat12)<3){
			pacDelete<-tabla12[which(nLat12<3,arr.ind=TRUE),1]
			for(j in 1:length(pacDelete)){
				filasPacDelete<-c(filasPacDelete,which(pacDelete[j]==final[,"ecg_id"]))
			}
		}
	}

	
	#modificamos los códigos de los >12 y les eliminamos de finales
	if(length(filasPacDelete)==0){
		addFilasDelete<-c()
	}else{
		addFilasDelete<-filasPacDelete
	}

	if(length(toDelete)>0){
		if(length(toDelete)<=2){
			nn<-1
		}else{
			nn<-nrow(toDelete)
		}
		if(nn>1){
			for(i in 1:nn){
				addFilasDelete<-unique(c(addFilasDelete,
					which(toDelete[i,1]==final[,"ecg_id"] & toDelete[i,2]==final[,"beat_id"],arr.ind=TRUE)))
				for(j in 1:length(lead8)){
					finales[[j]][[4]][which(toDelete[i,1]==finales[[j]][[4]][,1] & 
							toDelete[i,2]==finales[[j]][[4]][,2],arr.ind=TRUE),3]<-6
					auxF<-final[final[,"ecg_id"]==toDelete[i,1] & final[,"beat_id"]==toDelete[i,2],5+j]
					if(auxF==1){
						auxF2<-final[final[,"ecg_id"]==toDelete[i,1] & final[,"beat_id"]==toDelete[i,2],"annoRef"]
						posi<-which(auxF2==finales[[j]][[1]])
						finales[[j]][[1]]<-(finales[[j]][[1]])[-posi]
						if(length(finales[[j]][[2]])>2){
							finales[[j]][[2]]<-(finales[[j]][[2]])[-posi,]
						}else{
							finales[[j]][[2]]<-c(NA,NA)
						}
						finales[[j]][[3]]<-(finales[[j]][[3]])[-posi]
					}
				}
			}
		}else{
			addFilasDelete<-unique(c(addFilasDelete,
			which(toDelete[1]==final[,"ecg_id"] & toDelete[2]==final[,"beat_id"],arr.ind=TRUE)))
			for(j in 1:length(lead8)){
				finales[[j]][[4]][which(toDelete[1]==finales[[j]][[4]][,1] & 
						toDelete[2]==finales[[j]][[4]][,2],arr.ind=TRUE),3]<-6
				auxF<-final[final[,"ecg_id"]==toDelete[1] & final[,"beat_id"]==toDelete[2],5+j]
				if(auxF==1){
					auxF2<-final[final[,"ecg_id"]==toDelete[1] & final[,"beat_id"]==toDelete[2],"annoRef"]
					posi<-which(auxF2==finales[[j]][[1]])
					finales[[j]][[1]]<-(finales[[j]][[1]])[-posi]
					if(length(finales[[j]][[2]])>2){
						finales[[j]][[2]]<-(finales[[j]][[2]])[-posi,]
					}else{
						finales[[j]][[2]]<-c(NA,NA)
					}
					finales[[j]][[3]]<-(finales[[j]][[3]])[-posi]
				}
			}

		}
	}
	}
	#actualizamos la tabla final y la guardamos
	finalDef<-c()
	if(length(addFilasDelete)>0){
		finalDef<-final[-addFilasDelete,]
	}else{
		finalDef<-final
	}

	############################
	#	Output 	Funcion 
	######################
	finalDefAll<-c()
	if(eliminaPacAll){
		#10/01/2021
	  errorPreprocessing<-TRUE
		#finalDefAll<-c()
	}else{
	  errorPreprocessing<-FALSE
		if(length(final)>0)finalDefAll<-finalDef[,-c(1,14)]
	}

	
	
	detrend_I<-pre_I[[5]]
	if(is.null(detrend_I))detrend_I<-rep(NA,ncol(dataIn))
	detrend_II<-pre_II[[5]]
	if(is.null(detrend_II))detrend_II<-rep(NA,ncol(dataIn))
	detrend_III<-pre_III[[5]]
	if(is.null(detrend_III))detrend_III<-rep(NA,ncol(dataIn))
	detrend_AVL<-pre_AVL[[5]]	
	if(is.null(detrend_AVL))detrend_AVL<-rep(NA,ncol(dataIn))
	detrend_AVR<-pre_AVR[[5]]
	if(is.null(detrend_AVR))detrend_AVR<-rep(NA,ncol(dataIn))
	detrend_AVF<-pre_AVF[[5]]
	if(is.null(detrend_AVF))detrend_AVF<-rep(NA,ncol(dataIn))
	detrend_V1<-pre_V1[[5]]
	if(is.null(detrend_V1))detrend_V1<-rep(NA,ncol(dataIn))
	detrend_V2<-pre_V2[[5]]
	if(is.null(detrend_V2))detrend_V2<-rep(NA,ncol(dataIn))
	detrend_V3<-pre_V3[[5]]
	if(is.null(detrend_V3))detrend_V3<-rep(NA,ncol(dataIn))
	detrend_V4<-pre_V4[[5]]
	if(is.null(detrend_V4))detrend_V4<-rep(NA,ncol(dataIn))
	detrend_V5<-pre_V5[[5]]
	if(is.null(detrend_V5))detrend_V5<-rep(NA,ncol(dataIn))
	detrend_V6<-pre_V6[[5]]
	if(is.null(detrend_V6))detrend_V6<-rep(NA,ncol(dataIn))

	mDetrend<-cbind(detrend_I,detrend_II,detrend_III,detrend_AVL,detrend_AVR,detrend_AVF,
				detrend_V1,detrend_V2,detrend_V3,detrend_V4,detrend_V5,detrend_V6)
	colnames(mDetrend)<-leadName

	return(list(finalDefAll,mDetrend,errorPreprocessing))
}



