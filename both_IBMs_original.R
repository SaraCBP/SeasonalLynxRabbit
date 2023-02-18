#Both

setwd("~/uni/master/TFM/R")
library(truncnorm)
library(dbplyr)
library(ggplot2)
library(tidyverse)

#Inputs 
data.fine_gcm <- read.csv("data_fine_gcm.csv")  #climate data
nucleus_and_their_areas <- read.csv("nucleus_and_their_areas_lynx_size.csv", sep=";",header=T)
nucleus.vec <- as.factor(nucleus_and_their_areas$Nucleus)
lynx.sizes <- nucleus_and_their_areas$lynx_size   #lynx number in 2021 in each nuclei, used as initial population size

# Ages: I sampled an age structure for each nucleus and kept them if they were realistic. 
edades <- read.csv("edades.csv")

##Vital Rates Rabit

##Survival
PSA <- function(n,Fm){
  S1=exp(AA-(ds*(n/K)))
  S2=1-exp(-v/Fm)
  S=(S1/(1+S1))*S2
  return(S^30.4)
}
AA <-  10   #age effect on survival adults
#AA.values <- 5:10
ds <- 1.5   #density dependence in survival
#ds.values <- 1.5:5
v <- 40   #food availability factor
#V.values <- 14:45

PSJ <- function(n,Fm){
  S1=exp(AJ-(ds*(n/K)))
  S2=1-exp(-v/Fm)
  S=(S1/(1+S1))*S2
  return(S^30.4)
}
AJ <- 5    #age effect on survival juvenils
#AJ.values <- 2:4

#Babies
PSB <- function(MB){
  1-MB
}
MB <- 0.4 #newborn probability of mortaliy
#MB.values <- 0.4:0.8

#Breeding season  #if >= 0.5 --> reproduction
PB <- function(Ta,D,DI,W){
  B1=4.542-(0.605*Ta)+(0.029*(Ta^2))-(0.006*D)-(0.017*DI)-W
  return(1/(1+exp(B1)))
}
#Reproduction
PR <- function(n,rA){
  R1=rA-dr*(n/K)
  return(exp(R1)/(1+exp(R1)))
}
#dr.values <- 1.5:5
dr <- 3 #density-dependence in fecundity
r4 <- 0   #age parameter on reproduction months 4-6 
r6 <- 3   #age parameter on reproduction months 6-9
r46 <- r4+r6
r9 <- 6   #age parameter on reproduction months >9
r49 <- r4+r9

###VITAL RATES LYNX
#MONTHLY SURV

S0=(0.8406)^(1/12)        #newborns survival
SF1 = (0.87)^(1/12)       #(female survival from 12 to 24 months)
SM1 = (0.56)^(1/12)
SF2 = (0.99)^(1/12)     
SM2 = (1)^(1/12)        #( survival of males, from 24 months to 9 years)
SF9 = (0.8)^(1/12)       #( survival of females, from 10 years)
SM9 = (0.8)^(1/12) 

#Added mortality without a territory

AddMort=0.09

# function to calculate the monthly carrying capacity depending on rabbit adundance
klynx_function <- function(A){
  klynx=A/1000
  return(klynx)
}




months=12 # months of simulations 
years=5
scen=c("MPI-ESM-MR","GISS-E2-H","ACCESS1-0",
       "bcc-csm1-1-m","IPSL-CM5A-LR","CNRM-CM5","MRI-CGCM3","MIROC-ESM-CHEM","GISS-E2-R",
       "EC-EARTH") #cmip5 models that work fine

#hold rabbit
ibm.data1=NULL
dens.rabbit=NULL

# hold lynx
ibm.datal1=NULL
dens.lynx=NULL



for (z in 1:1){   # 10 runs of each nucleus x climate scenario
  for(nuc in 1:2){  
    for(s in 1:2){
      
      #initial rabbit pop
      #REFERENCIA --> 5.6 conejos/ha
      sizer0=90 #unidad espacial 0,1km2
      K <- sizer0+10 #carrying capacity
      {temp.conejos= data.frame(ID=paste(1:sizer0,"0",sep="_"),
                                gcm=scen[s],
                                edad=round(rnorm(sizer0,17,5),0), #?Est? bien? tendr?a que ser gaussiana y truncada
                                estado=factor(x=(sample(c("A","J","B"),size=sizer0,replace = TRUE)),levels = c("A","J","B")),
                                sexo=sample(c("M","F"),size=sizer0,replace = TRUE),
                                stringsAsFactors = T)
        temp.conejos$estado[temp.conejos$edad%in%(0)]="B"
        temp.conejos$estado[temp.conejos$edad%in%(1:3)]="J"
        temp.conejos$estado[temp.conejos$edad>3]="A"
      }
      
      #initial lynx pop
      
      sizel0=lynx.sizes[nuc]
      area_pob=nucleus_and_their_areas$Total.area..km2.[nuc] #area in km2
      K0=round(nucleus_and_their_areas$Total.area..km2.[nuc]/10) 
      
      
      temp.lynx <- data.frame(ID=paste(1:sizel0,"0",sep="_"),
                              gcm=scen[s],
                              edad=edades$edad[edades$nucleus%in%nucleus.vec[nuc]],
                              estado=sample(c("C","PD","D"),size=sizel0,replace = TRUE),
                              sexo=sample(c("M","F"),size=sizel0,replace = TRUE),
                              territory=sample(c("Yes","No","NotYet"),size=sizel0,replace = TRUE),
                              motherID=NA,
                              stringsAsFactors = TRUE
      )
      temp.lynx$estado[temp.lynx$edad%in%(0:12)]="C" #cubs
      temp.lynx$estado[temp.lynx$edad%in%(13:24)]="PD" #pre dispersal
      temp.lynx$estado[temp.lynx$edad>24]="D" #dispersed
      temp.lynx$territory[temp.lynx$estado%in%"C"]="NotYet"
      temp.lynx$territory[temp.lynx$estado%in%"PD"]="NotYet"
      temp.lynx$territory[temp.lynx$estado%in%"D"]="No"  #they will be assigned later 
      temp.lynx$age=trunc(temp.lynx$edad/12) #"edad" is the age in months, "age" is the age in years
      
      for(y in 1:years){
        for(m in 1:months){
          
          #rabbit
          temp.conejos$month=m
          temp.conejos$year=y
          temp.conejos$run=z
          temp.conejos$surv=NA
          temp.conejos$offspring=NA
          temp.conejos$densi=NA
          temp.conejos$densf=NA
          temp.conejos$bs=NA
          
          temp.conejos$Kits=NA
          temp.conejos$Juveniles=NA
          temp.conejos$Adults=NA
          
          temp.conejos$nucleus=nucleus.vec[nuc]
          date=data.fine_gcm$date[data.fine_gcm$yearl%in%y&data.fine_gcm$month%in%m&
                                    data.fine_gcm$gcm%in%scen[s]&data.fine_gcm$nucleus%in%nucleus.vec[nuc]] #data.fine_gcm$nucleus%in%coord_points$Nucleus[point],data.fine_gcm$point%in%coord_points$Point[point]
          
          temp.conejos$date=date
          
          FmL=data.fine_gcm$Fm[data.fine_gcm$yearl%in%y&data.fine_gcm$month%in%m&data.fine_gcm$gcm%in%scen[s]&data.fine_gcm$nucleus%in%nucleus.vec[nuc]]  #data.gcm$nucleus%in%coord_points$Nucleus[point],data.gcm$point%in%coord_points$Point[point]
          dens=as.numeric(nrow(temp.conejos))
          TaL=data.fine_gcm$ta[data.fine_gcm$yearl%in%y&data.fine_gcm$month%in%m&data.fine_gcm$gcm%in%scen[s]&data.fine_gcm$nucleus%in%nucleus.vec[nuc]] #data.gcm$nucleus%in%coord_points$Nucleus[point],data.gcm$point%in%coord_points$Point[point]
          DL=data.fine_gcm$mdl[data.fine_gcm$yearl%in%y&data.fine_gcm$month%in%m&data.fine_gcm$gcm%in%scen[s]&data.fine_gcm$nucleus%in%nucleus.vec[nuc]] #data.gcm$nucleus%in%coord_points$Nucleus[point],data.gcm$point%in%coord_points$Point[point]
          DIL=data.fine_gcm$DI[data.fine_gcm$yearl%in%y&data.fine_gcm$month%in%m&data.fine_gcm$gcm%in%scen[s]&data.fine_gcm$nucleus%in%nucleus.vec[nuc]] #data.gcm$nucleus%in%coord_points$Nucleus[point],data.gcm$point%in%coord_points$Point[point]
          WL=data.fine_gcm$W[data.fine_gcm$yearl%in%y&data.fine_gcm$month%in%m&data.fine_gcm$gcm%in%scen[s]&data.fine_gcm$nucleus%in%nucleus.vec[nuc]] #data.gcm$nucleus%in%coord_points$Nucleus[point],data.gcm$point%in%coord_points$Point[point]
          
          #initial density for this month
          temp.conejos$densi=dens
          
          #survival  
          temp.conejos$surv[temp.conejos$estado%in%"A"]=rbinom(nrow(temp.conejos[temp.conejos$estado%in%"A", ]),1,PSA(dens,FmL))
          temp.conejos$surv[temp.conejos$estado%in%"J"]=rbinom(nrow(temp.conejos[temp.conejos$estado%in%"J", ]),1,PSJ(dens,FmL))
          
          #offspring
          #desde el repr previo
          temp.conejos$offspring[temp.conejos$repr%in%1&temp.conejos$surv%in%1]=rpois(nrow(temp.conejos[temp.conejos$repr%in%1&temp.conejos$surv%in%1, ]),4) 
          born.babies=as.numeric(sum(temp.conejos$offspring,na.rm = T)) #son los q han nacido y su madre ha sobrevivido a este mes (si ha muerto ya no le calculo offspring)
          
          #para acabar el mes añadir los conejos recien nacidos
          if(born.babies>0){
            data.babies=data.frame(ID=paste(1:(born.babies),m,y,sep = "_"),gcm=scen[s],edad=0,
                                   estado="B",sexo=sample(c("M","F"),size=born.babies,replace = TRUE),
                                   month=m,run=z,surv=1,offspring=NA,repr=NA,year=y,densi=dens,densf=NA,
                                   bs=NA,date=date, stringsAsFactors = T,nucleus=nucleus.vec[nuc],
                                   Kits=Kits,Juveniles=Juveniles,Adults=Adults) #nucleus=coord_points$Nucleus[point],point=coord_points$Point[point]
            temp.conejos <- rbind(temp.conejos,data.babies,stringsAsFactors = T)
          }
          
          
          #breeding season
          #calcular si ese mes hay breeding season para determinar reproduction
          is.BS=PB(TaL,DL,DIL,WL)
          ifelse(is.BS<0.5,(temp.conejos$bs=0),(temp.conejos$bs=1))
          
          #reproduction
          temp.conejos$repr=NA 
          if(is.BS<0.5) {temp.conejos$repr=NA
          }else{
            temp.conejos$repr[temp.conejos$sexo%in%"F"&temp.conejos$estado%in%"A"&
                                temp.conejos$edad<6&temp.conejos$surv==1]=
              rbinom(nrow(temp.conejos[temp.conejos$sexo%in%"F"
                                       &temp.conejos$estado%in%"A"&temp.conejos$edad<6
                                       &temp.conejos$surv%in%1, ]),1,PR(dens,r4))
            temp.conejos$repr[temp.conejos$sexo%in%"F"&temp.conejos$estado%in%"A"&
                                temp.conejos$edad%in%6:9&temp.conejos$surv==1]=
              rbinom(nrow(temp.conejos[temp.conejos$sexo%in%"F"
                                       &temp.conejos$estado%in%"A"&temp.conejos$edad%in%6:9
                                       &temp.conejos$surv%in%1, ]),1,PR(dens,r46))
            temp.conejos$repr[temp.conejos$sexo%in%"F"&temp.conejos$estado%in%"A"&
                                temp.conejos$edad>9&temp.conejos$surv==1]=
              rbinom(nrow(temp.conejos[temp.conejos$sexo%in%"F"
                                       &temp.conejos$estado%in%"A"&temp.conejos$edad>9
                                       &temp.conejos$surv%in%1, ]),1,PR(dens,r49))
          }
          
          #final density
          densf <- nrow(temp.conejos[temp.conejos$surv%in%1, ])
          temp.conejos$densf= densf
          
          Kits= as.numeric(nrow(temp.conejos[temp.conejos$estado%in%"B",]))
          temp.conejos$Kits=Kits
          
          Juveniles= as.numeric(nrow(temp.conejos[temp.conejos$estado%in%"J",]))
          temp.conejos$Juveniles=Juveniles
          
          Adults= as.numeric(nrow(temp.conejos[temp.conejos$estado%in%"A",]))
          temp.conejos$Adults=Adults
          
          dens.temp=data.frame(scen=scen[s],year=y,month=m,densi=dens,densf=densf,
                               Kits=Kits,Juveniles=Juveniles,Adults=Adults,date=date,
                               nucleus=nucleus.vec[nuc],run=z)
          dens.rabbit<-rbind(dens.rabbit,dens.temp,stringsAsFactors=T)
          
          ibm.data1=rbind(ibm.data1,temp.conejos,stringsAsFactors = T)
          
          ## break si la no quedan conejos
          if((sum(temp.conejos$densf))<1) break
          
          #delete surv=0
          temp.conejos <- temp.conejos[temp.conejos$surv%in%1, ] 
          
          #update age and stage
          temp.conejos$edad <- temp.conejos$edad + 1
          temp.conejos$estado[temp.conejos$edad%in%(1:3)]="J"
          temp.conejos$estado[temp.conejos$edad>3]="A"
          temp.conejos$estado <- as.factor(temp.conejos$estado)
          
          ##############################################
          ######################
          ######################
          ###############################################
          
          #lynx
          temp.lynx$month=m
          temp.lynx$year=y
          temp.lynx$run=z
          temp.lynx$surv=NA
          temp.lynx$offspring=NA
          temp.lynx$densi=NA
          temp.lynx$densf=NA
          temp.lynx$K=NA
          
          date.lynx=data.fine_gcm$date[data.fine_gcm$yearl%in%(y)&data.fine_gcm$month%in%m&
                                         data.fine_gcm$gcm%in%scen[s]&data.fine_gcm$nucleus%in%nucleus.vec[nuc]] #data.fine_gcm$nucleus%in%nucleus_and_their_areas$Nucleus[nuc],data.gcm$point%in%1
          temp.lynx$date=date.lynx
          
          temp.lynx$nucleus=nucleus.vec[nuc]
          
          temp.lynx$Cubs=NA
          temp.lynx$PreDispersal=NA
          temp.lynx$Dispersal=NA
          
          
          #initial density for this month
          densi=as.numeric(nrow(temp.lynx))
          temp.lynx$densi=densi
          
          CSize=as.numeric(nrow(temp.lynx[temp.lynx$estado%in%"C",]))
          temp.lynx$Cubs=CSize
          
          PDSize=as.numeric(nrow(temp.lynx[temp.lynx$estado%in%"PD",]))
          temp.lynx$PreDispersal=PDSize
          
          DSize=as.numeric(nrow(temp.lynx[temp.lynx$estado%in%"D",]))
          temp.lynx$Dispersal=DSize
          
          
          #Carrying capacity
          rsizem=dens.rabbit$densi[dens.rabbit$scen%in%scen[s]&dens.rabbit$year%in%y&
                                   dens.rabbit$month%in%m&dens.rabbit$nucleus%in%nucleus.vec[nuc]&dens.rabbit$run%in%z] ##rabbit initial density for this month
          if(length(rsizem)> 0){rsizem=rsizem}else{rsizem=0} #avoid numeric(0)
          r_toar <- area_pob*10*rsizem #rabbit number for the total area
          r_per_l <- round(r_toar/K0) #rabbit per potential lynx area
          
          
          if (r_per_l>1000){
            K=K0
          }else{
            K=round(klynx_function(r_toar))
          }
          
          temp.lynx$K=K
          
          terr.sub=temp.lynx[temp.lynx$territory%in%"Yes",]
          
          ideal.age1=c(4,5,6,7,3,2,8,9) #this is the age ranking for territorial competition
          
          if(m%in%c(10,11,12,1)){ 
            
            new.terr=terr.sub$ID[terr.sub$age==(ideal.age1[1])]
            
            if(length(new.terr)<K){
              
              new.terr = c (new.terr,temp.lynx$ID[temp.lynx$age==(ideal.age1[1])&temp.lynx$territory%in%c("NotYet","No")])
            }
            
            track.terr=new.terr
            
            for(x in 2:length(ideal.age1)){ 
              
              
              if(length(track.terr)<K){ 
                
                new.terr=terr.sub$ID[terr.sub$age==(ideal.age1[x])]
                
                if(length(new.terr)<K){
                  if(x==6){     #this age is only aplied to females, 2 yo males apparently don't compete
                    new.terr = c(new.terr,temp.lynx$ID[temp.lynx$age==(ideal.age1[x])&
                                                         temp.lynx$territory%in%c("NotYet","No")&temp.lynx$sexo%in%"F"]) 
                    
                  }else{new.terr = c(new.terr,temp.lynx$ID[temp.lynx$age==(ideal.age1[x])&temp.lynx$territory%in%c("NotYet","No")])
                  }
                }
                
                track.terr=c(track.terr,new.terr) 
                
              }
              
              
            }
            
          }else{ 
            
            track.terr=terr.sub$ID[order(terr.sub$age)][1:K]
            
            
          }
          
          # reassign territories
          temp.lynx$territory[temp.lynx$territory%in%"Yes"]="No"
          temp.lynx$territory[temp.lynx$territory%in%"NotYet"&temp.lynx$estado%in%"D"]="No"
          temp.lynx$territory[temp.lynx$ID%in%track.terr]="Yes"
          
          
          #survival  
          
          temp.lynx$surv[temp.lynx$estado%in%"C"]=rbinom(nrow(temp.lynx[temp.lynx$estado%in%"C", ]),1,S0)  
          
          temp.lynx$surv[temp.lynx$estado%in%"PD"&temp.lynx$sexo%in%"F"]=rbinom(nrow(temp.lynx[temp.lynx$estado%in%"PD"&temp.lynx$sexo %in%"F", ]),1,SF1)  
          temp.lynx$surv[temp.lynx$estado%in%"PD"&temp.lynx$sexo%in%"M"]=rbinom(nrow(temp.lynx[temp.lynx$estado%in%"PD"&temp.lynx$sexo %in%"M", ]),1,SM1)  
          
          temp.lynx$surv[temp.lynx$estado%in%"D"&temp.lynx$sexo%in%"F"&temp.lynx$territory%in%"Yes"&temp.lynx$age%in%2:9]=rbinom(nrow(temp.lynx[temp.lynx$estado%in%"D"&temp.lynx$sexo%in%"F"&temp.lynx$territory%in%"Yes"&temp.lynx$age%in%2:9, ]),1,SF2)  
          temp.lynx$surv[temp.lynx$estado%in%"D"&temp.lynx$sexo%in%"M"&temp.lynx$territory%in%"Yes"&temp.lynx$age%in%2:9]=rbinom(nrow(temp.lynx[temp.lynx$estado%in%"D"&temp.lynx$sexo%in%"M"&temp.lynx$territory%in%"Yes"&temp.lynx$age%in%2:9, ]),1,SM2)  
          
          temp.lynx$surv[temp.lynx$estado%in%"D"&temp.lynx$sexo%in%"F"&temp.lynx$territory%in%"No"&temp.lynx$age%in%2:9]=rbinom(nrow(temp.lynx[temp.lynx$estado%in%"D"&temp.lynx$sexo%in%"F"&temp.lynx$territory%in%"No"&temp.lynx$age%in%2:9, ]),1,(SF2-AddMort))  
          temp.lynx$surv[temp.lynx$estado%in%"D"&temp.lynx$sexo%in%"M"&temp.lynx$territory%in%"No"&temp.lynx$age%in%2:9]=rbinom(nrow(temp.lynx[temp.lynx$estado%in%"D"&temp.lynx$sexo%in%"M"&temp.lynx$territory%in%"No"&temp.lynx$age%in%2:9, ]),1,(SM2-AddMort))  
          
          temp.lynx$surv[temp.lynx$estado%in%"D"&temp.lynx$sexo%in%"F"&temp.lynx$territory%in%"Yes"&temp.lynx$age>9]=rbinom(nrow(temp.lynx[temp.lynx$estado%in%"D"&temp.lynx$sexo%in%"F"&temp.lynx$territory%in%"Yes"&temp.lynx$age>9, ]),1,SF9)  
          temp.lynx$surv[temp.lynx$estado%in%"D"&temp.lynx$sexo%in%"M"&temp.lynx$territory%in%"Yes"&temp.lynx$age>9]=rbinom(nrow(temp.lynx[temp.lynx$estado%in%"D"&temp.lynx$sexo%in%"M"&temp.lynx$territory%in%"Yes"&temp.lynx$age>9, ]),1,SM9)  
          
          temp.lynx$surv[temp.lynx$estado%in%"D"&temp.lynx$sexo%in%"F"&temp.lynx$territory%in%"No"&temp.lynx$age>9]=rbinom(nrow(temp.lynx[temp.lynx$estado%in%"D"&temp.lynx$sexo%in%"F"&temp.lynx$territory%in%"No"&temp.lynx$age>9, ]),1,(SF9-AddMort))  
          temp.lynx$surv[temp.lynx$estado%in%"D"&temp.lynx$sexo%in%"M"&temp.lynx$territory%in%"No"&temp.lynx$age>9]=rbinom(nrow(temp.lynx[temp.lynx$estado%in%"D"&temp.lynx$sexo%in%"M"&temp.lynx$territory%in%"No"&temp.lynx$age>9, ]),1,(SM9-AddMort))  
          
          temp.lynx$surv[temp.lynx$edad%in%0:6&temp.lynx$motherID%in%temp.lynx$ID[temp.lynx$surv%in%0]]=0 #cubs younger than 6 months die if their mother does
          
          
          #reproduction
          
          temp.lynx$repr=NA 
          if(m%in%3){ temp.lynx$repr[temp.lynx$sexo%in%"F"&temp.lynx$age%in%3:9&
                                       temp.lynx$surv%in%1&temp.lynx$territory%in%"Yes"]=1
          
          temp.lynx$offspring[temp.lynx$repr%in%1]=rpois(nrow(temp.lynx[temp.lynx$repr%in%1, ]),3.5) 
          motherID=rep.int(temp.lynx$ID[temp.lynx$repr%in%1],times=temp.lynx$offspring[temp.lynx$repr%in%1])
          
          born.cubs=as.numeric(sum(temp.lynx$offspring,na.rm = T)) 
          
          }else{temp.lynx$repr=NA
          born.cubs=0
          }
          
          
          #newborns are added to the population
          if(born.cubs>0){
            data.cubs=data.frame(ID=paste(1:(born.cubs),m,y,sep = "_"),gcm=scen[s],edad=0,estado="C",
                                 sexo=sample(c("M","F"),size=1,replace = TRUE),territory="NotYet",
                                 month=m,run=z,surv=rbinom(born.cubs,1,S0),offspring=NA,repr=NA,year=y,densi=densi,densf=NA,
                                 K=K,Cubs=CSize,PreDispersal=PDSize,Dispersal=DSize,
                                 date=date.lynx,nucleus=nucleus.vec[nuc],age=0,motherID=motherID)   #nucleus=nucleus_and_their_areas$Nucleus[nuc]
            temp.lynx <- rbind(temp.lynx,data.cubs)
          }
          
          
          
          #final density
          densf=as.numeric(nrow(temp.lynx[temp.lynx$surv%in%1, ]))
          temp.lynx$densf=densf
          
          ibm.datal1=rbind(ibm.datal1,temp.lynx)
          dens.lynx.temp=data.frame(scen=scen[s],year=y,month=m,densi=densi,densf=densf,Cubs=CSize,
                                    PreDispersal=PDSize,Dispersal=DSize,date=date.lynx,
                                    nucleus=nucleus_and_their_areas$Nucleus[nuc],run=z)
          dens.lynx=rbind(dens.lynx,dens.lynx.temp)
          
          
          if(densf<1) break
          
          #new temp
          
          #delete surv=0
          temp.lynx <- temp.lynx[temp.lynx$surv%in%1, ]
          
          #update ages and stage
          temp.lynx$edad= temp.lynx$edad + 1
          
          temp.lynx$estado[temp.lynx$edad%in%(0:12)]="C" #cubs
          temp.lynx$estado[temp.lynx$edad%in%(13:24)]="PD" #pre dispersal
          temp.lynx$estado[temp.lynx$edad>24]="D" #dispersal
          temp.lynx$age=trunc(temp.lynx$edad/12)
        }
        if((sum(temp.lynx$densf))<1) break
        if((sum(temp.conejos$densf))<1) break
      }
      if((sum(temp.lynx$densf))<1) next
      if((sum(temp.conejos$densf))<1) next
    }
    if((sum(temp.lynx$densf))<1) next
    if((sum(temp.conejos$densf))<1) next
  }
}


+*write.csv(ibm.data1,"ibm_rabbit.csv")
write.csv(dens.rabbit,"dens_rabbit.csv")
write.csv(ibm.datal1,"ibm_lynx_trial_OK_df_2.csv",row.names = F)
write.csv(dens.lynx,"dens_lynx_trial_OK_df_2.csv",row.names = F)


  ibm.data1$run=z   #mal
  dens.rabbit$run=z


