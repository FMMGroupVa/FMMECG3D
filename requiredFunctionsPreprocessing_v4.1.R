
######################################################################################################################################
#
#
#									FUNCTIONS
#
#
#######################################################################################################################################
#data=dataIn[1,]
#freq=freqHz
require(stringr)
require(rsleep)

givePreprocessing_git<-function(dataIn=dataIn,freqHz=freqHz){#
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



leadPrePanTom_app<-function(data,freq){

data[is.na(data)]<-rep(0,length(data))
if( sum(data!=0)>0.75*length(data)){

	#ANALITICALLY
	dd<-data
	loessFilter_II<- loess(dd ~ c(1:length(dd)), span = 0.2) ## control span to define the amount of smoothing
	loessFitted_II<- loessFilter_II$fitted
	loessData_II<-dd-loessFitted_II
	loessRPeaks_II<-detect_rpeaks(loessData_II,freq)*freq

	#en los que hay uno muy outliaer y pan tompking falla
	while(length(loessRPeaks_II)<=2 ){ 
		q1<-quantile(dd,probs=0.9999)
		if(max(dd)>=q1){
			dd[which(dd>=q1,arr.ind=TRUE)]<-rep(median(dd),length(which(dd>=q1,arr.ind=TRUE)))
		}
			q2<-quantile(dd,probs=0.0001)
			if(min(dd)<=q2){
				dd[which(dd<=q2,arr.ind=TRUE)]<-rep(median(dd),length(which(dd<=q2,arr.ind=TRUE)))
			}
		
		
		loessFilter_II<- loess(dd ~ c(1:length(dd)), span = 0.2) ## control span to define the amount of smoothing
		loessFitted_II<- loessFilter_II$fitted
		loessData_II<-dd-loessFitted_II
		loessRPeaks_II<-detect_rpeaks(loessData_II,freq)*freq
	}
	allECG<-loessData_II
	allRPeaks<-loessRPeaks_II
	loessSeg_II<-giveSegmentsECG(loessData_II,loessRPeaks_II,freq)
	rrAux<-c()
	for(j in 2:length(loessRPeaks_II)){
		rrAux<-c(rrAux,loessRPeaks_II[j]-loessRPeaks_II[j-1])
	}	
	loessRR_II<-rrAux
	
	#tabla de todos pacientes y latidos y codigo
	#codigo 0 es ok
	#codigo asociado a problema

	tabCodes_II<-cbind(rep(1234,length(loessRPeaks_II)),
		1:length(loessRPeaks_II),
		rep(0,length(loessRPeaks_II)))


	#preprocesado
	paso1_II<-dropAnno1_app(loessData_II,loessRPeaks_II,loessRR_II,loessSeg_II,freq)
	#marco TRUE/FALSE y CODIGO
	tabCodes_II[(1:nrow(tabCodes_II))[-paso1_II[[2]]],3]<-rep(1,nrow(tabCodes_II)-length(paso1_II[[2]]))
	
	#REAJUSTO 
	seg1<-giveSegmentsECG(loessData_II,paso1_II[[1]],freq)
	rrAux<-c()
	for(j in 2:length(paso1_II[[1]])){
		rrAux<-c(rrAux,paso1_II[[1]][j]-paso1_II[[1]][j-1])
	}	
	rr1<-rrAux



	
	#finales
	if(length(seg1)>2)if(seg1[1,1]<0)seg1[1,1]<-1
	if(length(seg1)>2)if(seg1[nrow(seg1),2]>length(data))seg1[nrow(seg1),2]<-length(data)
	addR<-c();addSeg<-c();addRR<-c()
	if(length(paso1_II[[1]])>0){
		union<-paso1_II[[1]]
	}else{
		union<-c()
	}
	if(length(union)>0){
		union<-unique(union)
		c1<-FALSE
		if(length(union)==length(paso1_II[[1]]))c1<-sum(union==paso1_II[[1]])!=length(paso1_II[[1]])
		if(length(union)!=length(paso1_II[[1]]) | c1 ){
			for(j in 1:length(union)){
				if(!is.na(match(union[j],paso1_II[[1]]))){
					addR<-c(addR,paso1_II[[1]][which(union[j]==paso1_II[[1]])])
					addSeg<-rbind(addSeg,seg1[which(union[j]==paso1_II[[1]]),])
					addRR<-c(addRR,rr1[which(union[j]==paso1_II[[1]])-1])
				}
			}
			loessRPeaksEnd_II<-addR
			loessSegEnd_II<-addSeg
			loessRREnd_II<-addRR
		}else{
			loessRPeaksEnd_II<-paso1_II[[1]]
			loessSegEnd_II<-seg1
			loessRREnd_II<-rr1
		}
	}else{
		loessRPeaksEnd_II<-paso1_II[[1]][1]		#################################
		loessSegEnd_II<-seg1[1,]				#		decidir que hacer aqui
		loessRREnd_II<-rr1[1]					##############################
	}	


tabFinal_II<-tabCodes_II

}else{
	loessRPeaksEnd_II<-c()
	loessSegEnd_II<-c()
	loessRREnd_II<-c()
	
	loessData_II<-c()
	loessFitted_II<-c()
	addCount_II<-c()
	paso1_II<-c()
	tabCodes_II<-cbind(rep(1234,1),
		1:1,
		rep(11,1))
	tabFinal_II<-tabCodes_II

}
#EVALUEVMOS EN CUANTOS HAY PROBLEMAS CON DISTINTOS VALORES EN LOS EXRTEMOS
addCount_II<-c()
contPaci_II<-0
contador_II<-0
#tabFinal_II<-c()
	

	
return(list(loessRPeaksEnd_II,loessSegEnd_II,loessRREnd_II,tabFinal_II,loessData_II,loessFitted_II,addCount_II,
			paso1_II))
}




dropAnno1_app<-function(data,anno,rr,seg,freq){
	thrDNew<-2/3*median(rr)
	annoDrop<-c()
	rrNew<-rr
	annoNew<-anno
	segNew<-seg
	#la distancia ente dos anno es mayor que 2/mediana RR
	distAnnoNew<-c()
	for(j in 2:length(annoNew)){
		distAnnoNew[j-1]<-annoNew[j]-annoNew[j-1]
	}
	seguir<-(sum(rrNew<thrDNew)>0 | sum(distAnnoNew>(4/3*median(rrNew)))>=1) | (nrow(segNew)>(4/3*length(data)/freq) & (sd(rrNew)/mean(rrNew))>0.2)

	time<-1
	while(seguir){
		
		aqui<-FALSE
		cond<-rrNew<thrDNew
		if(sum(distAnnoNew>(4/3*median(rrNew)))>=1 & sum(rrNew<thrDNew)==0 & time==1 ){
			cond<-distAnnoNew>(4/3*median(rrNew))
			time<-2
			aqui<-TRUE
		}else{
			if(nrow(segNew)>(4/3*length(data)/freq)  & (sd(rrNew)/mean(rrNew))>0.2)cond<-c(rep(TRUE,nrow(segNew)-2),FALSE)# & sum(cond)==0
		}
		
		#if(sum(rrNew<thrDNew)>0)time<-2
		a1<-segNew[which(cond,arr.ind=TRUE),1]
		a2<-segNew[which(cond,arr.ind=TRUE),2]
		b5000<-rep(length(data),length(which(cond,arr.ind=TRUE)))
		b1<-rep(1,length(which(cond,arr.ind=TRUE)))
		a1_1<-segNew[which(cond,arr.ind=TRUE)+1,1]
		a2_1<-segNew[which(cond,arr.ind=TRUE)+1,2]
		b5000_1<-rep(length(data),length(which(cond,arr.ind=TRUE)+1))
		b1_1<-rep(1,length(which(cond,arr.ind=TRUE)+1))
		
		extMismo1<-(annoNew[which(cond,arr.ind=TRUE)]-ceiling(0.02*freq))
		extMismo2<-(annoNew[which(cond,arr.ind=TRUE)]+ceiling(0.02*freq))
		extMasUno1<-(annoNew[which(cond,arr.ind=TRUE)+1]-ceiling(0.02*freq))
		extMasUno2<-(annoNew[which(cond,arr.ind=TRUE)+1]+ceiling(0.02*freq))
		extMenosUno1<-(annoNew[which(cond,arr.ind=TRUE)]-ceiling(0.02*freq))
		extMenosUno2<-(annoNew[which(cond,arr.ind=TRUE)]+ceiling(0.02*freq))

		mismo<-c()
		masUno<-c()
		menosUno<-c()
		for(j in 1:length(extMismo1)){
			mismo[j]<-max(data[extMismo1[j]:extMismo2[j]])-min(data[extMismo1[j]:extMismo2[j]])
			masUno[j]<-max(data[extMasUno1[j]:extMasUno2[j]])-min(data[extMasUno1[j]:extMasUno2[j]])
			menosUno[j]<-max(data[extMenosUno1[j]:extMenosUno2[j]])-min(data[extMenosUno1[j]:extMenosUno2[j]])
		}
	
		une<-c()
		for(j in 1:length(which(cond,arr.ind=TRUE))){
			dd<-data[max(1,segNew[which(cond,arr.ind=TRUE)[j],1]):min(length(data),segNew[which(cond,arr.ind=TRUE)[j],2])]
			cond1<-abs(mismo[j]-masUno[j])>1/4*(max(dd)-min(dd))
			cond2<-(max(dd)-min(dd))>2*abs(mismo[j]) & (max(dd)-min(dd))>2*abs(masUno[j])
			cond3<-abs(mismo[j]-masUno[j])<1/4*(max(dd)-min(dd))
			if( (cond1 | (cond2 & cond3)) & 
				 (segNew[j+1,1]-segNew[j,2])<1/3*median(rrNew) ){
				if(mismo[j]>masUno[j] ){
					une<-c(une,(which(cond,arr.ind=TRUE)+1)[j])
				}else{
					une<-c(une,(which(cond,arr.ind=TRUE))[j])
				}
			}
			
		}
		if(length(une)>0)une<-unique(une)
		if(length(une)>0){
			annoDrop<-annoNew[-une]
			annoNew<-annoDrop
		}
		if(length(une)>0)segNew<-segNew[-une,]
		rrAux<-c()
		for(j in 2:length(annoNew)){
			rrAux<-c(rrAux,annoNew[j]-annoNew[j-1])
		}	
		rrNew<-rrAux
		elimina<-c()
		#revisar si entre los que nos han quedado todavia laguno tiene en un entrorno de la anota la amplitud mucho mas baja que ant o post
		if(length(annoDrop)>2 & time==1){
			elimina<-c()
			for(j in 2:(length(annoDrop)-1)){
				inf<-annoDrop[j]-ceiling(0.02*freq)
				sup<-annoDrop[j]+ceiling(0.02*freq)
				infSig<-annoDrop[j+1]-ceiling(0.02*freq)
				supSig<-annoDrop[j+1]+ceiling(0.02*freq)
				infAnt<-annoDrop[j-1]-ceiling(0.02*freq)
				supAnt<-annoDrop[j-1]+ceiling(0.02*freq)
				condAnt<-(2*max(data[infAnt:supAnt])<max(data[inf:sup]))
				condSup<-(2*max(data[infSig:supSig])<max(data[inf:sup]))
				if( (condAnt | condSup) & (min(length(data),segNew[,2][j])-max(1,segNew[,1][j]))<2/3*freq & sd(rrNew)/mean(rrNew)<0.2){
					if(condAnt){
						elimina<-c(elimina,j-1)
					}else{
						elimina<-c(elimina,j+1)
					}
				}
			}
			if(length(elimina)>0){
				elimina<-unique(elimina)
				annoDrop<-annoDrop[-elimina]
				segNew<-segNew[-elimina,]
			}
		}
		time<-2


		if(length(elimina)>0)annoNew<-annoDrop
		rrAux<-c()
		for(j in 2:length(annoNew)){
			rrAux<-c(rrAux,annoNew[j]-annoNew[j-1])
		}	
		rrNew<-rrAux
		thrDNew<-2/3*median(rrNew)
		if(length(annoNew)<=1){
			seguir<-FALSE
		}else{	
			distAnnoNew<-c()
			for(j in 2:length(annoNew)){
				distAnnoNew[j-1]<-annoNew[j]-annoNew[j-1]
			}
			seguir<-(sum(rrNew<thrDNew)>0 )#| sum(distAnnoNew>(4/3*median(rrNew)))>=2)# & time==1
			if(length(une)==0 & length(elimina)==0)seguir<-FALSE
		}
	}
	posAnnoNew<-c()
	for(j in 1:length(annoNew)){
		posAnnoNew<-c(posAnnoNew,which(annoNew[j]==anno))
	}
	return(list(annoNew,posAnnoNew))
}

giveSegmentsECG<-function(allECG,allRPeaks,fr){
  sucBeats<-giveSucesiveECG(allRPeaks,fr)[[1]]
  medRR<-giveSucesiveECG(allRPeaks,fr)[[2]]
  if(!is.na(medRR)){
    datpac<-allECG
    boundBeats<-matrix(0,length(allRPeaks),2)
    conti<-c()
    ya<-FALSE
    for(j in 1:length(allRPeaks)){

      if(sum(sucBeats[j,]==c(1,0,1))==3){
        conti[j]<-"R&L Cont"
      }else{
        if(sum(sucBeats[j,]==c(0,0,1))==3){
          conti[j]<-"R Cont"
        }else{
          if(sum(sucBeats[j,]==c(1,0,0))==3){
            conti[j]<-"L Cont"
          }else{
            if(sum(sucBeats[j,]==c(0,1,0))==3){
              conti[j]<-"Alone"
            }
          }
        }
      }

      if(sum(sucBeats[j,]==c(0,0,1))==3){
        RRbeat<-allRPeaks[j+1]-allRPeaks[j]
        if(sum(sucBeats[j+1,]==c(0,1,0))!=3){
          if(j==1 & floor(allRPeaks[j]-0.4*RRbeat)<(1:length(datpac))[1]){
            boundBeats[j,]<-c((1:length(datpac))[1],ceiling(allRPeaks[j]+0.6*RRbeat))
          }else{
            boundBeats[j,]<-c(floor(allRPeaks[j]-0.4*RRbeat),ceiling(allRPeaks[j]+0.6*RRbeat))
          }
        }else{
          boundBeats[j,]<-c(floor(allRPeaks[j]-0.4*medRR),ceiling(allRPeaks[j]+0.6*medRR))
        }
      }
      if(sum(sucBeats[j,]==c(1,0,0))==3 | (sum(sucBeats[j,]==c(0,1,0))==3 & j==length(allRPeaks))){
        RRbeat<-allRPeaks[j]-allRPeaks[j-1]
        if(sum(sucBeats[j-1,]==c(0,1,0))!=3 & (sum(sucBeats[j,]==c(0,1,0))!=3)){
          if(j==length(allRPeaks) & ceiling(allRPeaks[j]+0.6*RRbeat)>(1:length(datpac))[length((1:length(datpac)))]){
            boundBeats[j,]<-c(floor(allRPeaks[j]-0.4*RRbeat),(1:length(datpac))[length((1:length(datpac)))])
            ya<-TRUE
          }else{
            boundBeats[j,]<-c(floor(allRPeaks[j]-0.4*RRbeat),ceiling(allRPeaks[j]+0.6*RRbeat))
          }
        }else{
          boundBeats[j,]<-c(floor(allRPeaks[j]-0.4*medRR),ceiling(allRPeaks[j]+0.6*medRR))
        }
      }
      if(sum(sucBeats[j,]==c(1,0,1))==3){
        RRbeatL<-allRPeaks[j]-allRPeaks[j-1]
        RRbeatR<-allRPeaks[j+1]-allRPeaks[j]
        boundBeats[j,]<-c(floor(allRPeaks[j]-0.4*RRbeatL),
                          ceiling(allRPeaks[j]+0.6*RRbeatR))
      }
      if(sum(sucBeats[j,]==c(0,1,0))==3 & ya==FALSE){
        RRbeat<-medRR
        boundBeats[j,]<-c(floor(allRPeaks[j]-0.4*medRR),ceiling(allRPeaks[j]+0.6*medRR))
      }
    }
    #corregimos para no ajustar dos veces los mismos puntos
    for(j in 2:nrow(boundBeats)){
      if(boundBeats[j,1]<=boundBeats[j-1,2])boundBeats[j,1]<-boundBeats[j-1,2]+1
    }
  }else{
    boundBeats<-c()
  }

  return(boundBeats)
}

#vAnnotN=allRPeaks[[1]]
#frec=fr
giveSucesiveECG<-function(vAnnotN,frec){
  if(length(vAnnotN)>1){
    lRR<-c()
    for(i in 2:length(vAnnotN)){
      lRR[i-1]<-vAnnotN[i]-vAnnotN[i-1]
    }
#si en al menos el 75% de ellos la longitud es mayor que la frequencia y tenemos mas de 3 latinos	
    if(sum(lRR>frec)>0.75*length(lRR) & length(lRR)>3){
      lRR<-sort(lRR)[-c(length(lRR):(length(lRR)-2))]
    }else{
      #si tenemos mas de tres latidos y en al menos uno la long es mayor que la frequencia 
      if(sum(lRR>frec)>0 & length(lRR)>3){
        lRR<-(lRR)[-which(lRR>frec)]##
      }else{
		#si tenemos tres latidos y  en alguno la longitud mayor a la frequencia
		if(sum(lRR>frec)>0 & length(lRR)==3){
        		lRR<-sort(lRR)[-length(lRR)]
      	}
	}
    }

    M<-median(lRR)
    #print(paste("es m",M))
    fin<-c()
    for(i in 1:length(vAnnotN)){
      if(i==1 &
         (vAnnotN[2]-vAnnotN[1])<1.5*M)idSuc<-c(0,0,1)#podr?a se problema si RR_12 != RR_23, asumimos que no ocurre
      if(i==1 &
         (vAnnotN[2]-vAnnotN[1])>=1.5*M)idSuc<-c(0,1,0)
      if(i==length(vAnnotN) &
         (vAnnotN[length(vAnnotN)]-vAnnotN[length(vAnnotN)-1])<1.5*M)idSuc<-c(1,0,0)
      if(i==length(vAnnotN) &
         (vAnnotN[length(vAnnotN)]-vAnnotN[length(vAnnotN)-1])>=1.5*M)idSuc<-c(0,1,0)
      if( i >=2 & i<=(length(vAnnotN)-1)){
        RRl<-vAnnotN[i]-vAnnotN[i-1]
        RRd<-vAnnotN[i+1]-vAnnotN[i]
        if(RRl<1.5*M & RRd<1.5*M)idSuc<-c(1,0,1)
        if(RRl<1.5*M & RRd>1.5*M)idSuc<-c(1,0,0)
        if(RRl>1.5*M & RRd<1.5*M)idSuc<-c(0,0,1)
        if(RRl>1.5*M & RRd>1.5*M)idSuc<-c(0,1,0)
      }
      fin<-rbind(fin,idSuc)
    }
  }else{
    M<-NA
    fin<-c()
  }

  return(list(fin,M))
}



#annoLeads2=annoLeads
#nObsSignal=ncol(dataIn)
#freq=freqHz

giveAnnoMed_app<-function(annoLeads2,nObsSignal,freq){

		auxBeats<-annoLeads2[1,]
		auxLeads<-annoLeads2[2,]
		Vbeat<-c()
		Vlead<-c()
		Vbeat[1]<-auxBeats[1]
		Vlead[1]<-auxLeads[1]
		beats<-list()
		leads<-list()
		finalBeats<-list()
		finalLeads<-list()
		b<-1
		ii<-2
		j<-2
		limit1<-103*freq/1000#103ms
		limit2<-322*freq/1000#322ms

		#primera parte: eliminar los que nos osn
		while(j<=length(auxBeats)){
			if(auxBeats[j]>(auxBeats[j-1]+limit1)){
				#Diferent beat, evaluate previous
				if(length(Vbeat)>3){
					#hacer la mediana
					beats[[b]]<-ceiling(median(Vbeat))
					leads[[b]]<-Vlead
					b<-b+1	
					ii<-2
					Vbeat<-c()
					Vlead<-c()
					Vbeat[1]<-auxBeats[j]
					Vlead[1]<-auxLeads[j]
					j<-j+1
				}else{
					#reject
					ii<-2
					Vbeat<-c()
					Vlead<-c()
					Vbeat[1]<-auxBeats[j]
					Vlead[1]<-auxLeads[j]
					j<-j+1
				}
			}else{
				if(j!=length(auxBeats)){
					Vbeat[ii]<-auxBeats[j]
					Vlead[ii]<-auxLeads[j]
					ii<-ii+1
					j<-j+1
					#print(j)
				}else{
					#lo he añadido yo para que terminase
					if(length(Vbeat)>3){
						#hacer la mediana
						beats[[b]]<-ceiling(median(Vbeat))
						leads[[b]]<-Vlead
						j<-j+1
					}else{
						j<-j+1
					}
				}
			}
		}
	

		#segunda parte: posibles falsos positivos
		ii<-2
		k<-2
		finalBeats[[1]]<-beats[[1]]
		finalLeads[[1]]<-leads[[1]]
	
		#pongo <= en lugar de <
		while(ii<=length(beats)){
			if(beats[[ii]]>(beats[[ii-1]]+limit2)){
				finalBeats[[k]]<-beats[[ii]]
				finalLeads[[k]]<-leads[[ii]]
				ii<-ii+1
				k<-k+1
			}else{
				#falso positivo: esto lo añado yo, porque sino no acaba
				ii<-ii+1
			}
		}


		##############################3
		#	Lead de referencia
		###############################
		annosPac<-unlist(finalBeats)
		leadsPac<-unlist(finalLeads)
		#lo he cambiado, para no tener q	ue poner los datos,
		#sino la longitud, que es lo unico necesario
		if(length(giveSegmentsECG_multi(lengthData=nObsSignal,allRPeaks=annosPac,fr=freq))>0){
			segPac<-giveSegmentsECG_multi(lengthData=nObsSignal,allRPeaks=annosPac,fr=freq)
			if(segPac[1,1]<1)segPac[1,1]<-1
			if(segPac[length(annosPac),2]>nObsSignal )segPac[length(annosPac),2]<-nObsSignal
		}else{
			segPac<-c(NA,NA)
		}

		
return(list(annosPac,leadsPac,segPac))
}


giveSegmentsECG_multi<-function(lengthData,allRPeaks,fr){
  sucBeats<-giveSucesiveECG(allRPeaks,fr)[[1]]
  medRR<-giveSucesiveECG(allRPeaks,fr)[[2]]
  if(!is.na(medRR)){
 
    boundBeats<-matrix(0,length(allRPeaks),2)
    conti<-c()
    ya<-FALSE
    for(j in 1:length(allRPeaks)){

      if(sum(sucBeats[j,]==c(1,0,1))==3){
        conti[j]<-"R&L Cont"
      }else{
        if(sum(sucBeats[j,]==c(0,0,1))==3){
          conti[j]<-"R Cont"
        }else{
          if(sum(sucBeats[j,]==c(1,0,0))==3){
            conti[j]<-"L Cont"
          }else{
            if(sum(sucBeats[j,]==c(0,1,0))==3){
              conti[j]<-"Alone"
            }
          }
        }
      }

      if(sum(sucBeats[j,]==c(0,0,1))==3){
        RRbeat<-allRPeaks[j+1]-allRPeaks[j]
        if(sum(sucBeats[j+1,]==c(0,1,0))!=3){
          if(j==1 & floor(allRPeaks[j]-0.4*RRbeat)<(1:lengthData)[1]){
            boundBeats[j,]<-c((1:lengthData)[1],ceiling(allRPeaks[j]+0.6*RRbeat))
          }else{
            boundBeats[j,]<-c(floor(allRPeaks[j]-0.4*RRbeat),ceiling(allRPeaks[j]+0.6*RRbeat))
          }
        }else{
          boundBeats[j,]<-c(floor(allRPeaks[j]-0.4*medRR),ceiling(allRPeaks[j]+0.6*medRR))
        }
      }
      if(sum(sucBeats[j,]==c(1,0,0))==3 | (sum(sucBeats[j,]==c(0,1,0))==3 & j==length(allRPeaks))){
        RRbeat<-allRPeaks[j]-allRPeaks[j-1]
        if(sum(sucBeats[j-1,]==c(0,1,0))!=3 & (sum(sucBeats[j,]==c(0,1,0))!=3)){
          if(j==length(allRPeaks) & ceiling(allRPeaks[j]+0.6*RRbeat)>(1:lengthData)[length((1:lengthData))]){
            boundBeats[j,]<-c(floor(allRPeaks[j]-0.4*RRbeat),(1:lengthData)[length((1:lengthData))])
            ya<-TRUE
          }else{
            boundBeats[j,]<-c(floor(allRPeaks[j]-0.4*RRbeat),ceiling(allRPeaks[j]+0.6*RRbeat))
          }
        }else{
          boundBeats[j,]<-c(floor(allRPeaks[j]-0.4*medRR),ceiling(allRPeaks[j]+0.6*medRR))
        }
      }
      if(sum(sucBeats[j,]==c(1,0,1))==3){
        RRbeatL<-allRPeaks[j]-allRPeaks[j-1]
        RRbeatR<-allRPeaks[j+1]-allRPeaks[j]
        boundBeats[j,]<-c(floor(allRPeaks[j]-0.4*RRbeatL),
                          ceiling(allRPeaks[j]+0.6*RRbeatR))
      }
      if(sum(sucBeats[j,]==c(0,1,0))==3 & ya==FALSE){
        RRbeat<-medRR
        boundBeats[j,]<-c(floor(allRPeaks[j]-0.4*medRR),ceiling(allRPeaks[j]+0.6*medRR))
      }
    }
    #corregimos para no ajustar dos veces los mismos puntos
    for(j in 2:nrow(boundBeats)){
      if(boundBeats[j,1]<=boundBeats[j-1,2])boundBeats[j,1]<-boundBeats[j-1,2]+1
    }
  }else{
    boundBeats<-c()
  }

  return(boundBeats)
}


#objPrePanTom<-pre_I
#annosRef<-annosPac
#segRef<-segPac
#freq<-freqHz

leadPreMulti_app<-function(objPrePanTom,annosRef,segRef,freq){
paso2_II<-c();paso3_II<-c();paso5_II<-c()
contTend<-c()
totalLat<-c()
if( length(annosRef)>=3 & sum(objPrePanTom[[5]]!=0)>0.75*length(objPrePanTom[[5]])){
	#print(i)
	#ANALITICALLY
	loessFitted_II<- objPrePanTom[[6]]
	loessData_II<- objPrePanTom[[5]]
	loessRPeaks_II<-annosRef
	loessSeg_II<-segRef
	rrAux<-c()
	for(j in 2:length(loessRPeaks_II)){
		rrAux<-c(rrAux,loessRPeaks_II[j]-loessRPeaks_II[j-1])
	}	
	loessRR_II<-rrAux
	
	#tabla de todos pacientes y latidos y codigo
	#codigo 0 es ok
	#codigo asociado a problema

	tabCodes_II<-cbind(rep(1234,length(loessRPeaks_II)),
		1:length(loessRPeaks_II),
		rep(0,length(loessRPeaks_II)))


	
	#preprocesado#eliomar primero o ultimo valido si son menos del 855 d ela RR mediana
	paso2_II<-dropAnno2_app(loessData_II,loessRPeaks_II,loessRR_II,loessSeg_II,loessRPeaks_II,freq)
	#marco TRUE/FALSE y CODIGO
	tabCodes_II[(1:nrow(tabCodes_II))[paso2_II[[2]]],3]<-rep(2,length(paso2_II[[2]]))
	
	#preprocesado#eliminar latidos con anotaciones mal colocadas fuera del 35-45 de la longitud del latido 
	paso3_II<-dropAnno3_app(loessRPeaks_II,loessRR_II,loessSeg_II,loessRPeaks_II)
	#marco TRUE/FALSE y CODIGO
	tabCodes_II[(1:nrow(tabCodes_II))[paso3_II[[2]]],3]<-rep(3,length(paso3_II[[2]]))


		

	#preprocesado#extemos distintos 
	paso5_II<-dropAnno5_app(loessData_II,loessRPeaks_II,loessRR_II,loessSeg_II,loessRPeaks_II)
	#marco TRUE/FALSE y CODIGO
	tabCodes_II[(1:nrow(tabCodes_II))[paso5_II[[2]]],3]<-rep(5,length(paso5_II[[2]]))
	sumaT<-length(paso5_II[[2]])
	




	
	#finales
	if(length(loessRPeaks_II)>2)if(loessSeg_II[1,1]<0)loessSeg_II[1,1]<-1
	if(length(loessRPeaks_II)>2)if(loessSeg_II[nrow(loessSeg_II),2]>length(objPrePanTom[[5]]))loessSeg_II[nrow(loessSeg_II),2]<-objPrePanTom[[5]]
	addR<-c();addSeg<-c();addRR<-c()
	if(length(paso2_II[[1]])>0){
		if(length(paso3_II[[1]])>0){
			contTend<-sumaT
			totalLat<-length(intersect(paso2_II[[1]],paso3_II[[1]]))
			if(length(paso5_II[[1]])>0){
				union<-intersect(paso2_II[[1]],intersect(paso3_II[[1]],paso5_II[[1]]))
			}else{
				union<-c()#intersect(paso2_II[[1]],paso3_II[[1]])
			}
		}else{
			contTend<-sumaT
			totalLat<-length(paso2_II[[1]])
			if(length(paso5_II[[1]])>0){
				union<-c()#intersect(paso2_II[[1]],paso5_II[[1]])
			}else{
				union<-c()#paso2_II[[1]]
			}
		}
	}else{
		if(length(paso3_II[[1]])>0){
			contTend<-sumaT
			totalLat<-length(paso3_II[[1]])
			if(length(paso5_II[[1]])>0){
				union<-c()#intersect(paso3_II[[1]],paso5_II[[1]])
			}else{
				union<-c()#paso3_II[[1]]
			}
		}else{
			contTend<-sumaT
			totalLat<-length(loessRPeaks_II)
			if(length(paso5_II[[1]])>0){
				union<-c()#paso5_II[[1]]
			}else{
				union<-c()
			}
		}
	}
	if(length(union)>0){
		union<-unique(union)
		
		c1<-FALSE
		if(length(union)==length(loessRPeaks_II))c1<-sum(union==loessRPeaks_II)!=length(loessRPeaks_II)
		if(length(union)!=length(loessRPeaks_II) | c1 ){
			for(j in 1:length(union)){
				if(!is.na(match(union[j],loessRPeaks_II))){
					addR<-c(addR,loessRPeaks_II[which(union[j]==loessRPeaks_II)])
					addSeg<-rbind(addSeg,loessSeg_II[which(union[j]==loessRPeaks_II),])
					addRR<-c(addRR,loessRR_II[which(union[j]==loessRPeaks_II)-1])
				}
			}
			loessRPeaksEnd_II<-addR
			loessSegEnd_II<-addSeg
			loessRREnd_II<-addRR
		}else{
			loessRPeaksEnd_II<-loessRPeaks_II
			loessSegEnd_II<-loessSeg_II
			loessRREnd_II<-loessRPeaks_II
		}
	}else{
		loessRPeaksEnd_II<-c()#loessRPeaks_II[1]		#################################
		loessSegEnd_II<-c()#loessSeg_II[1,]			#		decidir que hacer aqui
		loessRREnd_II<-c()#loessRPeaks_II[1]					##############################
	}	

}else{
		loessFitted_II<- objPrePanTom[[6]]
	loessData_II<- objPrePanTom[[5]]
	loessRPeaks_II<-annosRef
	loessSeg_II<-segRef

		tabCodes_II<-cbind(rep(1234,length(loessRPeaks_II)),
		1:length(loessRPeaks_II),
		rep(10,length(loessRPeaks_II)))
		loessRPeaksEnd_II<-c()#loessRPeaks_II[1]		#################################
		loessSegEnd_II<-c()#loessSeg_II[1,]			#		decidir que hacer aqui
		loessRREnd_II<-c()#loessRPeaks_II[1]					##############################

}	




#EVALUEVMOS EN CUANTOS HAY PROBLEMAS CON DISTINTOS VALORES EN LOS EXRTEMOS
addCount_II<-c()
contPaci_II<-0
contador_II<-0
tabFinal_II<-c()
	#if(length(loessSegEnd_II)>6){
	#	for(j in 1:nrow(loessSegEnd_II)){
	#		ini1<-1
	#		ini2<-ceiling(0.02*(loessSegEnd_II[j,2]-loessSegEnd_II[j,1]))
	#		fin1<-loessSegEnd_II[j,2]-loessSegEnd_II[j,1]+1-ceiling(0.02*(loessSegEnd_II[j,2]-loessSegEnd_II[j,1]))
	#		fin2<-loessSegEnd_II[j,2]-loessSegEnd_II[j,1]
	#		dataNorm<-minimax(loessData_II[loessSegEnd_II[j,1]:loessSegEnd_II[j,2]])
	#		valIni<-median(dataNorm[ini1:ini2])
	#		valFin<-median(dataNorm[fin1:fin2])
	#		contador_II<-contador_II+1
	#		if(abs(valFin-valIni)>0.1){
	#			addCount_II<-rbind(addCount_II,c(1234,loessRPeaksEnd_II[j],abs(valFin-valIni)))
	#		}
	#		
	#	}
	#}else{
	#	contPaci_II<-contPaci_II+1
	#}
	tabFinal_II<-rbind(tabFinal_II,tabCodes_II)

return(list(loessRPeaksEnd_II,loessSegEnd_II,loessRREnd_II,tabFinal_II,loessData_II,loessFitted_II,addCount_II,
			paso2_II,paso3_II,paso5_II,contTend,totalLat))
}



dropAnno2_app<-function(data,anno,rr,seg,annoOrig,freq){
	annoNew<-anno
	rrNew<-rr
	segNew<-seg
	primero<-FALSE
	segundo<-FALSE
	a<-data[max(1,segNew[1,1]):max((segNew[1,1]+ceiling(0.02*(min(length(data),segNew[1,2])-max(1,segNew[1,1])))),1+ceiling(0.02*(min(length(data),segNew[1,2])-max(1,segNew[1,1]))))]
	b<-data[(min(length(data)-ceiling(0.02*(min(length(data),segNew[length(annoNew),2])-max(1,segNew[length(annoNew),1]))),segNew[length(annoNew),2])-ceiling(0.02*(min(length(data),segNew[length(annoNew),2])-max(1,segNew[length(annoNew),1])))):min(length(data),segNew[length(annoNew),2])]




	##NEW##9306, 11342 En el primer latdo hay una caida/subida muy pronunciada antes 5%
			ini1<-max(1,segNew[1,1])##NEW##
			ini2<-ini1-1+ceiling(0.15*(min(length(data),segNew[1,2])-max(1,segNew[1,1])))##NEW##
			fin1<-min(length(data),segNew[1,2])+1-ceiling(0.15*(min(length(data),segNew[1,2])-max(1,segNew[1,1])))##NEW##
			fin2<-min(length(data),segNew[1,2])##NEW##
			dataNorm<-(data[ini1:fin2])##NEW##
			vMin<-which.min(dataNorm[(ini1-ini1+1):(fin2-ini1+1)])##NEW##
			vMax<-which.max(dataNorm[(ini1-ini1+1):(fin2-ini1+1)])##NEW##
			valMin<-min(dataNorm[(ini1-ini1+1):(fin2-ini1+1)])##NEW##
			valMax<-max(dataNorm[(ini1-ini1+1):(fin2-ini1+1)])##NEW##

			diff1<-abs(valMin-valMax)##NEW##
			if(!is.na(match(vMin,(ini1-ini1+1):(ini2-ini1+1))) | !is.na(match(vMax,(ini1-ini1+1):(ini2-ini1+1))) )primero<-TRUE##NEW##
			

			diff<-c()
			for(j in 2:nrow(segNew)){
				ini1<-max(1,segNew[j,1])##NEW##
				ini2<-ini1-1+ceiling(0.15*(min(length(data),segNew[j,2])-max(1,segNew[j,1])))##NEW##
				fin1<-min(length(data),segNew[j,2])+1-ceiling(0.15*(min(length(data),segNew[j,2])-max(1,segNew[j,1])))##NEW##
				fin2<-min(length(data),segNew[j,2])##NEW##
				dataNorm<-(data[ini1:fin2])##NEW##
				vMin<-which.min(dataNorm[(ini1-ini1+1):(fin2-ini1+1)])##NEW##
				vMax<-which.max(dataNorm[(ini1-ini1+1):(fin2-ini1+1)])##NEW##
				valMin<-min(dataNorm[(ini1-ini1+1):(fin2-ini1+1)])##NEW##
				valMax<-max(dataNorm[(ini1-ini1+1):(fin2-ini1+1)])##NEW##
				diff<-c(diff,abs(valMin-valMax))

			}
			if(diff1>2/3*median(diff) | diff1<2/3*median(diff)){
				#primero<-TRUE##NEW##
				if(primero)segundo<-TRUE
			}

	seguir<-((segNew[1,2]-segNew[1,1])<(0.9*median(segNew[,2]-segNew[,1])) & length(annoNew)>1)|#NEW
			(sd(a)/abs(mean(a))<0.01 )|
			(sd(a)/abs(mean(a))<0.02 & length(annoNew)>=0.04*freq) | 
			segundo##NEW##
			


	while(seguir){
		annoDrop<-annoNew[-1]
		annoNew<-annoDrop
		rrAux<-c()
		for(j in 2:length(annoNew)){
			rrAux<-c(rrAux,annoNew[j]-annoNew[j-1])
		}	
		rrNew<-rrAux
		segNew<-segNew[-1,]
		if(length(annoNew)<=1){
			seguir<-FALSE
		}else{
			seguir<-(segNew[1,2]-segNew[1,1])<(0.85*median(segNew[,2]-segNew[,1])) & length(annoNew)>1
		}
		
	}
	if(length(annoNew)>1)seguir<-((segNew[nrow(segNew),2]-segNew[nrow(segNew),1])<(0.9*median(segNew[,2]-segNew[,1]))) | (sd(b)/abs(mean(b))<0.01)|(sd(b)/abs(mean(b))<0.02 & length(annoNew)>=0.04*freq)#NEW
	while(seguir){
		annoDrop<-annoNew[-nrow(segNew)]
		annoNew<-annoDrop
		rrAux<-c()
		for(j in 2:length(annoNew)){
			rrAux<-c(rrAux,annoNew[j]-annoNew[j-1])
		}	
		rrNew<-rrAux
		segNew<-segNew[-nrow(segNew),]
		if(length(annoNew)<=1){
			seguir<-FALSE
		}else{
			seguir<-(segNew[nrow(segNew),2]-segNew[nrow(segNew),1])<(0.85*median(segNew[,2]-segNew[,1])) & length(annoNew)>1
		}
	}
	posAnnoNew<-c()
	diferencia<-setdiff(anno,annoNew)
	for(j in 1:length(diferencia)){
		posAnnoNew<-c(posAnnoNew,which(diferencia[j]==annoOrig))
	}

	return(list(annoNew,posAnnoNew))
}



dropAnno3_app<-function(anno,rr,seg,annoOrig){
	annoNew<-anno
	rrNew<-rr
	segNew<-seg
	if(length(annoNew)>1){
		if(length(segNew)>0)seguir<-sum((((segNew[,2]-segNew[,1])*0.35)>(annoNew-segNew[,1]+1)) | (annoNew-segNew[,1]+1)>(0.45*(segNew[,2]-segNew[,1])))>=1
		if(seguir){
			donde<-which(((segNew[,2]-segNew[,1])*0.35)>(annoNew-segNew[,1]+1) | (annoNew-segNew[,1]+1)>(0.45*(segNew[,2]-segNew[,1])),arr.ind=TRUE)
			annoDrop<-annoNew[-donde]
			annoNew<-annoDrop
			rrAux<-c()
			for(j in 2:length(annoNew)){
				rrAux<-c(rrAux,annoNew[j]-annoNew[j-1])
			}	
			rrNew<-rrAux
			segNew<-segNew[-donde,]
		}
	}
	posAnnoNew<-c()
	diferencia<-setdiff(anno,annoNew)
	for(j in 1:length(diferencia)){
		posAnnoNew<-c(posAnnoNew,which(diferencia[j]==annoOrig))
	}

	return(list(annoNew,posAnnoNew))
}

#data=loessData_II
#anno=loessRPeaks_II
#rr=loessRR_II
#seg=loessSeg_II
#annoOrig=loessRPeaks_II
dropAnno5_app<-function(data,anno,rr,seg,annoOrig){
	annoNew<-anno
	rrNew<-rr
	segNew<-seg
	seguir<-FALSE
	if(length(annoNew)>1){
		contador_II<-1
		addCount_II<-c()
		dife<-c()
		amp<-c()
		condiKey<-c()
		revisarKey<-c()
		for(j in 1:length(annoNew)){
			ini1<-max(1,segNew[j,1])
			ini2<-ini1-1+ceiling(0.05*(min(length(data),segNew[j,2])-max(1,segNew[j,1])))
			fin1<-min(length(data),segNew[j,2])+1-ceiling(0.05*(min(length(data),segNew[j,2])-max(1,segNew[j,1])))
			fin2<-min(length(data),segNew[j,2])
			dataNorm<-(data[ini1:fin2])
			#dataNorm<-(data[ini1:fin2])
			valIni<-median(dataNorm[(ini1-ini1+1):(ini2-ini1+1)])
			valFin<-median(dataNorm[(fin1-ini1+1):(fin2-ini1+1)])
			contador_II<-contador_II+1
			dife[j]<-abs(valFin-valIni)
			amp[j]<-abs(max(dataNorm)-min(dataNorm))
			condiKey[j]<-(dife[j])>0.15*(max(dataNorm)-min(dataNorm))##NEW##
			if(condiKey[j])revisarKey<-c(revisarKey,j)
		}
		if(sum(condiKey,na.rm=TRUE)>=4){
			condiKey<-rep(TRUE,length(annoNew))
			dife<-rep(999,length(annoNew))
			amp<-rep(1,length(annoNew))
		}
		if(length(revisarKey)>0)seguir<-TRUE
		
		if(seguir){
			if(length(revisarKey)>0)donde<-unique(c(which(dife>amp,arr.ind=TRUE),revisarKey))
			annoDrop<-annoNew[-donde]
			annoNew<-annoDrop
			rrAux<-c()
			for(j in 2:length(annoNew)){
				rrAux<-c(rrAux,annoNew[j]-annoNew[j-1])
			}	
			rrNew<-rrAux
			segNew<-segNew[-donde,]
		}
	}
	posAnnoNew<-c()
	diferencia<-setdiff(anno,annoNew)
	if(length(diferencia)>0){
		for(j in 1:length(diferencia)){
			posAnnoNew<-c(posAnnoNew,which(diferencia[j]==annoOrig))
		}
	}
	return(list(annoNew,posAnnoNew))
}

minimax<-function(y){
        return(-1+2*(y-min(y))/(max(y)-min(y)))
}







giveAnnoMed2_app<-function(annoLeads2,nObsSignal,freq){


		auxBeats<-annoLeads2[1,]
		auxLeads<-annoLeads2[2,]
		Vbeat<-c()
		Vlead<-c()
		Vbeat[1]<-auxBeats[1]
		Vlead[1]<-auxLeads[1]
		beats<-list()
		leads<-list()
		finalBeats<-list()
		finalLeads<-list()
		b<-1
		ii<-2
		j<-2
		limit1<-103*freq/1000#103ms
		limit2<-322*freq/1000#322ms

		#primera parte: eliminar los que nos osn
		while(j<=length(auxBeats)){
			if(auxBeats[j]>(auxBeats[j-1]+limit1)){
				#Diferent beat, evaluate previous
				if(length(Vbeat)>0){#OJO QUE ESTO EN EL CASO DE LAS 12 LEADS ES >3
					#hacer la mediana
					beats[[b]]<-ceiling(median(Vbeat))
					leads[[b]]<-Vlead
					b<-b+1	
					ii<-2
					Vbeat<-c()
					Vlead<-c()
					Vbeat[1]<-auxBeats[j]
					Vlead[1]<-auxLeads[j]
					j<-j+1
				}else{
					#reject
					ii<-2
					Vbeat<-c()
					Vlead<-c()
					Vbeat[1]<-auxBeats[j]
					Vlead[1]<-auxLeads[j]
					j<-j+1
				}
			}else{
				if(j!=length(auxBeats)){
					Vbeat[ii]<-auxBeats[j]
					Vlead[ii]<-auxLeads[j]
					ii<-ii+1
					j<-j+1
					#print(j)
				}else{
					#lo he añadido yo para que terminase
					if(length(Vbeat)>0){#OJO QUE ESTO EN EL CASO DE LAS 12 LEADS ES >3
						#hacer la mediana
						beats[[b]]<-ceiling(median(Vbeat))
						leads[[b]]<-Vlead
						j<-j+1
					}else{
						j<-j+1
					}
				}
			}
		}
	

		#segunda parte: posibles falsos positivos
		ii<-2
		k<-2
		finalBeats[[1]]<-beats[[1]]
		finalLeads[[1]]<-leads[[1]]
	
		#pongo <= en lugar de <
		while(ii<=length(beats)){
			if(beats[[ii]]>(beats[[ii-1]]+limit2)){
				finalBeats[[k]]<-beats[[ii]]
				finalLeads[[k]]<-leads[[ii]]
				ii<-ii+1
				k<-k+1
			}else{
				#falso positivo: esto lo añado yo, porque sino no acaba
				ii<-ii+1
			}
		}


		##############################3
		#	Lead de referencia
		###############################
		annosPac<-unlist(finalBeats)
		leadsPac<-unlist(finalLeads)
		#lo he cambiado, para no tener q	ue poner los datos,
		#sino la longitud, que es lo unico necesario
		if(length(giveSegmentsECG_multi(lengthData=nObsSignal,allRPeaks=annosPac,fr=freq))>0){
			segPac<-giveSegmentsECG_multi(lengthData=nObsSignal,allRPeaks=annosPac,fr=freq)
			if(segPac[1,1]<1)segPac[1,1]<-1
			if(segPac[length(annosPac),2]>nObsSignal)segPac[length(annosPac),2]<-nObsSignal
		}else{
			segPac<-NA
		}

	return(list(annosPac,leadsPac,segPac))
}

