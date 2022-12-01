
######################################################################################################################################
#
#
#									FUNCTIONS
#
#
#######################################################################################################################################
#data=dataIn[1,]
#freq=freqHz

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

