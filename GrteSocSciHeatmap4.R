setwd("d:/Rprojects/GrteSocSci/Rcode")
library(here)
library(readxl) # data input
library(data.table) # data input
library(lubridate) # text to POSIXct
library(mokken) # NIRT functions
library(Cairo) # plotting to files
library(ggplot2) # plotting framework
library(patchwork)
library(twosamples)
library(qualpalr)
library(party)
library(future.apply)
library(progressr)
library(beepr)
library(purrr)
# skip the IF statement to run all processing in a single thread
# KF: reset Rprofile.site to be simple, not Intel oneAPI
# could use parallely::supportsMulticore()
# may need to expand the memory allocated to global variables
# options(future.globals.maxSize= 500*1024^2)
if ("Darwin" == Sys.info()["sysname"]){ # parallel ops using fork
  plan(multicore, workers=7)
}else{ # WinDoze
plan(multisession, workers=14, gc=TRUE)
}
cfctrl <- cforest_control(ntree=1000)
custom.col <- c("#FFDB6D", "#C4961A", "#F4EDCA", "#D16103", "#C3D7A4", "#52854C", "#4E84C4", "#293352")
# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# string constants for graphics
quest1 <- paste0("For the current lighting conditions created by\n",
                 "the streetlights in Colter Bay, please indicate how\n",
                 "true the following statements are of your experience\n",
                 "(Select one number for each item.)\n",
                 "The current lighting conditions created by the\n",
                 "streetlights in Colter Bay...")
quest2 <- paste0("Please indicate the degree to which you oppose or support\n",
                 "the following management actions designed to protect\n",
                 "the quality of stargazing/viewing the night sky at this park.")
questParts1 <- c(
  "... makes nighttime recreation activities\n      more pleasurable.",
  "... makes it easy for my eyes to transition\n     from light to dark places.",
  "... makes it more difficult for me to do other\n     nighttime activities in the park.",
  "... makes it less safe.",
  "... makes it easier to navigate while doing\n     nighttime activities.",
  "... promotes more natural wildlife behaviors.",
  "... provides benefits for wildlife and insects.",
  "... reduces human impacts on wildlife.",
  "... makes it seem inappropriately bright.",
  "... makes it seem inappropriately dark.",
  paste0("... makes it difficult to find areas and points",
         "\n     of interest."),
  "... affects the Greater Yellowstone Ecosystem.",
  "... makes it easy for me to see in unlit areas\n     (in the shadows).")
questParts2 <- c(
  "Setting lights to the minimum necessary brightness.",
  "Reducing the number of park lights.",
  "Restricting the number of lights visitors can\nuse at night.",
  "Creating a shield on lights that direct light\nonly to intended areas.",
  "Adjusting hues of lights to be wildlife friendly.",
  "Adjusting hues of lights to preserve\nhumans' night vision."
  )
questVals1 <- NULL
answerVals1 <- c("Not at all true (1)", "Slightly True (2)",
                 "Moderately True(3)", "Very True (4)", "Completely True (5)")
answerVals2 <- c("Completely Oppose", "Oppose", "Neither Oppose nor Support",
                 "Support", "Completely Support")
revCodeH <- function(likert, vChk=F, naDrop=T){
  if(nullNames <- is.null(names(likert)))
    lkNames <- dimnames(likert)[[2]]
  else
    lkNames <- names(likert)
  nix <- apply(is.na(likert),1,any)
  likertNA <- likert[!nix,]
  if(F == vChk) vChk <- seq(along = likertNA[1,])
  hIV <- coefH(as.data.frame(likertNA), se = F, results = F)
  hBest <- hIV$H # overall rating of the current scale
  # look for variables that merit reverse coding
  vList <- order(apply(hIV$Hij * (hIV$Hij<0),1,sum))
  vList <- vList[vList %in% vChk]
  while(length(vList)>0){
    iv <- vList[1]; vChk <- vChk[iv != vChk]
    likertNA[,iv] <- 6 - likertNA[,iv]
    hIV <- coefH(as.data.frame(likertNA), se = F, results = F)
    if(hBest > hIV$H) # undo, not better
      likertNA[,iv] <- 6 - likertNA[,iv]
    else{ #better, update best
      hBest <- hIV$H
      vList <- order(apply(hIV$Hij * (hIV$Hij<0),1,sum))
      lkNames[iv] <- paste(lkNames[iv],"-",sep="")
    }
    vList <- vList[vList %in% vChk]
  }
  if(!naDrop){ # apply answer converse to all data
    revIX <- sapply(lkNames, function(x)
      substring(x,nchar(x),nchar(x)))=="-"
    likertNA <- likert
    likertNA[,revIX] <- 6 - likertNA[,revIX]
  }
  if(nullNames)
    dimnames(likertNA)[[2]] <- lkNames
  else
    names(likertNA) <- lkNames
  return(list(Likert = likertNA, dropped = nix))
} # end revCodeH()
itemReplicates <- function(itmXsimMat,reGroup=T){
  # compresses the columns of the matrix argument, keeping track of the number of
  # duplicates with the column names
  itmXsimMat <- matrix(as.numeric(itmXsimMat),nrow=dim(itmXsimMat)[1],
                       dimnames = dimnames(itmXsimMat))
  itmXsimMat[is.na(itmXsimMat)] <- -999
  res <- itmXsimMat # to be populated with updated group numbers
  # first, optionally remove the artifact of group order (convert to shared basis)
  if(reGroup){
    res <- res*0
    allClmnGrps <- apply(itmXsimMat,2,split,x=dimnames(itmXsimMat)[[1]])
    allGrps <- lapply(unlist(allClmnGrps,recursive = F),sort)
    ikeep <- names(allGrps) != "1.0" # do not care about unscaled items
    names(allGrps) <- c(1,cumsum(diff(as.numeric(names(allGrps)))<0) + 1) # names become column indices
    allGrps <- allGrps[ikeep]
    uniqGrps <- unique(allGrps) # all unique scales
    uniqGrps <- uniqGrps[order(sapply(uniqGrps,length),sapply(uniqGrps,paste,collapse=""),
                               decreasing=T)]
    mjunk <- match(allGrps,uniqGrps)
    smjunk <- split(mjunk,names(allGrps))
    smjunk <- smjunk[order(as.numeric(names(smjunk)))]
    for (iCol in seq(along=smjunk))
      res[unlist(uniqGrps[smjunk[[iCol]]]),iCol] <- 
      rep(order(smjunk[[iCol]]),times=sapply(uniqGrps[smjunk[[iCol]]],length))
  } # if reGroup
  colIx <- seq(along=res[1,])[colDup <- duplicated(res,MARGIN=2)]
  colRmvd <- NULL
  for(Cl in colIx){
    cMatch <- seq(along=res[1,])[!apply(res[,Cl] != res,2,any)]
    cMatch <- cMatch[!(cMatch %in% colRmvd)]
    if(length(cMatch)>1){
      dimnames(res)[[2]][min(cMatch)] <-
        sum(as.numeric(dimnames(res)[[2]][cMatch]))
      colRmvd <- c(colRmvd,cMatch[cMatch != min(cMatch)])
    }
  }
  res <- res[,!colDup]
  res <- res[,order(as.numeric(dimnames(res)[[2]]),
                    decreasing = T)]
  res[res==-999] <- NA
  return(res)
} #itemReplicates()
getIllum <- function(fname){
  mtemp <- fread(file=fname) #EDIT FILE TO REMOVE DASHED LINE
  names(mtemp) <- gsub(patt=" ",replacement = "",names(mtemp))
  if(is.character(mtemp[,DateTime])) mtemp[,DateTime:=ymd_hms(DateTime)]
  return(mtemp)
} #getIllum()
source(here("Rcode", "heatMkf6.R"))
natLum <- getIllum(here("Data", "SunMoon2019.txt"))
natLum[,DateTime := with_tz(DateTime,tzone="America/Denver")] # MDT
setkey(natLum,DateTime)
hntrUV <- fread(here("Data", "06_AggBatDataWuv.csv"))
nix <- seq(along=names(hntrUV))[!(names(hntrUV) %in%
                                    c("SampleDate","Color","Intensity"))]
hntrUV[,(nix) := NULL]
hntrUV[, ":="(SampleDate = ymd(SampleDate),
              Color = as.factor(Color))]
setnames(hntrUV, "Intensity", "N10city")
hntrUV <- unique(hntrUV)
setkey(hntrUV,SampleDate)
{ # prepare red data matrix
  red <- read_excel(here("Data", "GRTE lighting.xlsx"))
  red[red==999] <- NA
  red[,grep(patt="awe",names(red),value = T,ignore.case = T)] <- NULL
  red[,grep(patt="Biodiver",names(red),value=T,ignore.case=T)] <- NULL
  red[,grep(patt="Non_resp",names(red),value=T,ignore.case=T)] <- NULL
  red[,grep(patt="BEL",names(red),value=T)] <- NULL
  red[,grep(patt="center",names(red),value=T,ignore.case=T)] <- NULL
  red$SurveyMinutes <- as.numeric(red$EndDate-red$StartDate)
  red$SawBats <- !is.na(red$NumberBats)
  red$NumberBats[is.na(red$NumberBats)] <- 0
  varKeep <- c("StartDate","SurveyMinutes","Cloud_Cover","Darkness",
               "StreetlightCondition","Hue","Intensity","Location",
               "NightSkyCondition","FirstTimeVisitor",
               "YearFirstTripGRTE","PreviousVisitCB",
               "PreviousEveningCB","CampCB_ThisTrip",
               "NumberNightsCB","PreviousCampCB",
               "YearFirstCampCB","CBActivity_ranger",
               "CBActivity_walk","CBActivity_stargze",
               "CBActivity_Other",
               "NumberBats","Streetlights_activitiesmorepleasurable",
               "Streetlights_eyetransition",
               "Streetlights_activitiesdifficult",
               "Streetlights_lesssafe","Streetlights_easiernavigate",
               "Streetlights_wildlifebehavior",
               "Streetlights_wildlifebenefits",
               "Streetlights_reducehumanimpacts",
               "Streetlights_bright","Streetlights_dark",
               "Streetlights_pointsofinterest","Streetlights_affectGYE",
               "Streetlights_unlitareas","Management_minimumbright",
               "Management_reducelights",
               "Management_restrictvisitorlights","Management_shield",
               "Management_hueswildlife","Management_huespeople",
               "Yearborn","Gender",
               "Permanentresident","Personalgroup","Language",
               "Illuminance_Vertical_0","Illuminance_Vertical_90",
               "Illuminance_Vertical_180","Illuminance_Vertical_270",
               "Illuminance_Vertical_Sun",
               "Illuminance_Vertical_BrightestLight",
               "Illuminance_Horizontal","SawBats")
  red <- red[,varKeep]
  factorVars <- c("StreetlightCondition","Hue","Location",
                  "FirstTimeVisitor","PreviousVisitCB",
                  "PreviousEveningCB","CampCB_ThisTrip",
                  "PreviousCampCB","CBActivity_ranger",
                  "CBActivity_walk","CBActivity_stargze",
                  "CBActivity_Other","Gender","SawBats",
                  "Permanentresident", "Language")
  red[, factorVars][is.na(red[, factorVars])] <- 0
  red[,factorVars] <- lapply(red[,factorVars],as.factor)
  orderVars <- c("Darkness","NightSkyCondition")
  red[, orderVars][is.na(red[, orderVars])] <- 0
  red[, orderVars] <- lapply(red[, orderVars],as.ordered)
  red[,grep(patt="YearFirst", names(red), value=T)][
	is.na(red[,grep(patt="YearFirst", names(red), value=T)])] <- 2020
  red$StartDate <- force_tz(red$StartDate,tzone="America/Denver")
  red$Illuminance_Vertical_Integrated <-
    as.matrix(red[,grep(patt="Illuminance",names(red))][,1:4]) %*% rep(1,4)
  red$NumberNightsCB <- as.numeric(red$NumberNightsCB)
  red$NumberNightsCB[is.na(red$NumberNightsCB)] <- -1
} # prepare red data matrix
# read brightness weighting from:
# Spectral considerations for outdoor lighting:
# Designing for perceived scene brightness
brtCoef <- read.csv(here("Radiometry/", "Vb2BrightnessCoefficients.csv"), header=F)
brtCoef2 <- t(matrix(unlist(brtCoef), nrow=5,
	dimnames=list(c("nm", "pRelHps", "pRelLed", "Vlambda", "V2lambda"), NULL)))
brtCoef <- as.data.table(brtCoef2)
load(here("Radiometry", "AsphaltReflectance.RDat"))
# eval.parent, quote, and substitute to access the name of the list element?
HntrSpec <- fread(here("Radiometry", "LightMeasurementsLong.csv"))
names(HntrSpec)[3] <- "Erg"
HntrSpec[, Erg:=Erg/max(Erg)]
HntrSpec[, LightSet := sub(patt="Red(..)$", repl="Red0\\1", HntrSpec[, LightSet])]
HntrSpec[, LightSet := sub(patt="White(..)$", repl="White0\\1", HntrSpec[, LightSet])]
setkey(HntrSpec, wl, LightSet)
if(0){
	ggplot(LightsLong, aes(wl, `Watts/sr/m2`, color = LightSet)) +
  geom_line() +
  scale_color_manual(breaks = c("Red10", "Red20", "Red30", "Red40", "Red50", "Red60", "Red70",
                                "Red80", "Red90", "Red100", "White10", "White20", "White30",
                                "White40", "White50", "White60", "White70", "White80",
                                "White90", "White100"),
                     labels = c("Red 10%", "Red 20%", "Red 30%", "Red 40%", "Red 50%", "Red 60%", "Red 70%",
                                "Red 80%", "Red 90%", "Red 100%", "White 10%", "White 20%", "White 30%",
                                "White 40%", "White 50%", "White 60%", "White 70%", "White 80%",
                                "White 90%", "White 100%"),
                     values = c("#ff5a5a", "#ff4242", "#ff2b2b", "#ff1414", "#e51212",
                                "#cc1010", "#b20e0e", "#990c0c", "#7f0a0a", "#660808",
                                "#14b4ff", "#12a2e5", "#1090cc", "#0e7db2", "#0c6c99",
                                "#0a5a7f", "#084866", "#06364c", "#042433", "#021219")) +
  labs(x = "Wavelength (nm)", y = bquote("W *"~sr^-1~m^-2),
       color = "Light Settings") +
  theme_PREM(MajorGrid = "grey90", MinorGrid = "grey95")
	}
kolrs <- c(rgb(red=(seq(from=32, to=248, by=24))/255, green=(seq(from=32, to=248, by=24))/1023,
	blue=rep(0,10)),
	rgb(red=(seq(from=32, to=248, by=24))/1023, green=(seq(from=32, to=248, by=24))/384,
		blue=(seq(from=32, to=248, by=24))/255))
specPlot <- ggplot(data=HntrSpec, aes(x=wl, y=Erg, color=LightSet))+
	geom_line(linewidth=1.1) + scale_color_manual(values=kolrs) +
	theme(axis.title.x = element_text(size=20),
		axis.text.x = element_text(size=16),
		axis.title.y = element_text(size=20),
		axis.text.y = element_text(size=16),
		legend.title = element_text(size=20),
		legend.text = element_text(size=16)) +
	labs(x = "Wavelength (nm)", y = "Relative spectral intensity") +
	guides(color=guide_legend(title="Color/\nIntensity"))
# lumens(red) = 5040.45917063
# calculate lumens(white)
unlist(HntrSpec[LightSet=="White100", Erg])%*%Hspec2[, "Vphot"]*
	5040.45917063/unlist(HntrSpec[LightSet=="Red100",Erg])%*%Hspec2[, "Vphot"]
HntrSpec[, Color:=NULL]
Hspec2 <- t(matrix(unlist(HntrSpec[, Erg]), nrow=20,
	dimnames=list(HntrSpec[, unique(LightSet)], HntrSpec[, unique(wl)])))
waveLn <- sort(unique(HntrSpec$wl))
tst <- approx(x=brtCoef[, .(nm, Vlambda)], xout=waveLn, rule=1)$y
tst[is.na(tst)] <- 0
Hspec2 <- cbind(Hspec2, Vphot=tst)
tst <- approx(x=brtCoef[, .(nm, V2lambda)], xout=waveLn, rule=1)$y
tst[is.na(tst)] <- 0
Hspec2 <- cbind(Hspec2, Vb2=tst)
asphMlt <- approx(x=as.matrix(asphDat[["B"]]), xout=as.integer(dimnames(Hspec2)[[1]]), rule=2)[[2]]
asphMlt <- asphMlt/max(asphMlt)
Hspec2[, "Vphot"] <- Hspec2[, "Vphot"] * asphMlt
Hspec2[, "Vb2"] <- Hspec2[, "Vb2"] * asphMlt
# add the photometric luminance row
Hspec2 <- rbind(Hspec2, Photopic=0)
for(jcol in 1:10){
	Hspec2[402, jcol] <- Hspec2[,jcol] %*% Hspec2[, "Vphot"]
	Hspec2[402, jcol+10] <- Hspec2[,jcol+10] %*% Hspec2[, "Vphot"]
}
# add the Vb2 luminance row
Hspec2 <- rbind(Hspec2, B2=0)
for(jcol in 1:10){
	Hspec2[403, jcol] <- Hspec2[,jcol] %*% Hspec2[, "Vb2"]
	Hspec2[403, jcol+10] <- Hspec2[,jcol+10] %*% Hspec2[, "Vb2"]
}
write.csv(cbind((Hspec2["Photopic", 1:10]), (Hspec2["Photopic", 11:20]),
	(Hspec2["B2", 1:10]), (Hspec2["B2", 11:20])),
	file=here("Output", "BrightnessComparison.csv"))
zQix <- # ordered as in the questionnaire (and initial data table)
  grep(patt="^Streetlights.+", names(red))
dIx <- findInterval(red$StartDate,natLum$DateTime) + 1 # aligning red .. natLum
# dIx marks a sol/luna datum that is within five minutes of StartDate
red$dTime <- lubridate::hour(red$StartDate)+
  lubridate::minute(red$StartDate)/60+lubridate::second(red$StartDate)/3600
red$yDay <- yday(red$StartDate)
red <- cbind(red,natLum[dIx,-"DateTime"])
dIx <- match(date(red$StartDate),date(hntrUV$SampleDate))
red<- cbind(red,hntrUV[dIx,.(Color,N10city)]) #AFTER NEGQIX, POSQIX
# add subjective brightness metrics
dimVals <- as.integer(sub(patt="Red(...)", repl="\\1",
	names(Hspec2["Photopic", grep(patt="Red", names(Hspec2["Photopic",]))])))
PhotFun <- approxfun(x=dimVals,
	y=Hspec2["Photopic", grep(patt="Red", names(Hspec2["Photopic",]))],
	rule=2)
red$PhotBrite[red$Color=="R"] <- PhotFun(red$N10city[red$Color=="R"])
PhotFun <- approxfun(x=dimVals,
	y=Hspec2["Photopic", grep(patt="White", names(Hspec2["Photopic",]))],
	rule=2)
red$PhotBrite[red$Color=="W"] <- PhotFun(red$N10city[red$Color=="W"])
PhotFun <- approxfun(x=dimVals,
	y=Hspec2["B2", grep(patt="Red", names(Hspec2["B2",]))],
	rule=2)
red$B2Brite[red$Color=="R"] <- PhotFun(red$N10city[red$Color=="R"])
PhotFun <- approxfun(x=dimVals,
	y=Hspec2["B2", grep(patt="White", names(Hspec2["B2",]))],
	rule=2)
red$B2Brite[red$Color=="W"] <- PhotFun(red$N10city[red$Color=="W"])

tst <- ymd_hms(red$StartDate)
unique(date(tst)) # Surveys administered on how many dates?
# data set for aisp calculations
red2 <- revCodeH(red[ ,zQix], naDrop=F)
negQix <- grep(patt="-",names(red2$Likert))
posQix <- seq(along=names(red2$Likert)); posQix <- posQix[!(posQix %in% negQix)]
All3S <- apply(red2$Likert==3,1,all) # also get rid of all "3"
# Odds of all 3s by chance
prod(apply(red2$Likert[!red2$dropped & !All3S,]==3,2,sum)/
	dim(red2$Likert[!red2$dropped & !All3S,])[1])
# Odds of any replicated surveys
redTmp <- data.frame(lapply(red[, grep(patt="Streetlights_", names(red))], as.factor))
streetFreqs <- sapply(redTmp, summary)/dim(redTmp)[1]
1 - (1 - prod(apply(streetFreqs^2, 2, sum)))^dim(redTmp)[1]
###### add lighting and weather to red2
dateIX <- match(date(red$StartDate),date(hntrUV$SampleDate)) # adding lighting
All3S[is.na(All3S)] <- F # get rid of all "3" for complexPlot
wthr <- fread(here("Data", "MooseWeatherHr.csv"))
wthr[,ytime :=as.POSIXct(
  julianHr*3600,origin=as.POSIXct("1970-01-01",tz="GMT"),tz="America/Denver")]
setkey(wthr,ytime)
# look at duplicates by color
negQix <- negQix[-6] # get rid of affctGE-
sum(apply(red2$Likert[,posQix]>3,1,all),na.rm=T) # how many surveys had all >3 for "+Q"
sum(apply(red2$Likert[,negQix]==5,1,all),na.rm=T) # how many surveys had all 1 for "-Q"
safety5Ix <- red2$Likert$Streetlights_lesssafe==1 
safety5Ix[is.na(safety5Ix)] <- F
naCount <- apply(is.na(red2$Likert), 1, sum, na.rm=T)
red5 <- t(red2$Likert)
dimnames(red5)[[2]] <- rep(1,dim(red5)[2])
dupSurveys <- itemReplicates(red5, reGroup = F)
write.csv(dupSurveys[,as.numeric(dimnames(dupSurveys)[[2]])>1],
          file=here("Output", "StreetlightDuplications.csv"))
stlt <- red2$Likert
names(stlt) <- sapply(names(stlt), sub, patt="-", repl="")
skales <- list(
  "1"=names(stlt)[6:7],
  "1a"=names(stlt)[6:8],
  "2"=names(stlt)[1:2],
  "1b"=names(stlt)[c(1,2,6:8)],
  "1c"=names(stlt)[c(1,2,6:9)],
  "3"=names(stlt)[c(5,13)],
  "4"=names(stlt)[c(3,4)],
  "5"=names(stlt)[c(10,11)]
)
dateIX <- match(date(red$StartDate),date(hntrUV$SampleDate)) # adding lighting
stlt <- cbind(stlt, hntrUV[dateIX,.(Color,N10city)])
dateIX2 <- apply(red2$Likert[,posQix]>3,1,all); dateIX2[is.na(dateIX2)] <- F
summary(stlt[dateIX2,"Color"])
dateIX2 <- apply(red2$Likert[,negQix]==5,1,all); dateIX2[is.na(dateIX2)] <- F
summary(stlt[dateIX2,"Color"])
strtNirt <- complexPlot(rmat=stlt, xCensor = All3S,
	scales=skales, Quest = quest1, qVals=questParts1, aVals=answerVals1)
write.csv(strtNirt, file=here("Output", "StreetlightScalesTable.csv"))
# test for differences in Safety for incomplete surveys
incTst <- cbind(
	summary(as.factor(stlt$Streetlights_lesssafe[stlt$Color=="R" &
		apply(is.na(stlt[, c(1:3, 5:13)]), 1, any)])),
	c(0, summary(as.factor(stlt$Streetlights_lesssafe[stlt$Color=="W" &
		apply(is.na(stlt[, c(1:3, 5:13)]), 1, any)]))))
chisq.test(incTst[-6,]) #, simulate=T)
############# NOW REPEAT FOR MANAGEMENT QUESTIONS
mQix <- grep(patt="Management_",names(red))
mgmt <- red[,mQix]
mgmTr <- t(mgmt); dimnames(mgmTr)[[2]] <- rep(1,dim(mgmTr)[2])
dupSurveys2 <- itemReplicates(mgmTr,reGroup = F)
repProb <- function(tval, datmat){
  replicIx <- apply(mgmt, 1, sd)==0
  replicIx[is.na(replicIx)] <- F
  tvProb <- apply(datmat[!replicIx,]==tval, 2, sum, na.rm=T)/dim(datmat)[1]
  return(tvProb)
}
sum(apply(is.na(mgmt[, 1:6]), 1, any))
print("expected number of uniform responses, after omitting them (5, 4, 3)")
prod(repProb(5, mgmt[, 1:6]))*dim(mgmt)[1]
prod(repProb(4, mgmt[, 1:6]))*dim(mgmt)[1]
prod(repProb(3, mgmt[, 1:6]))*dim(mgmt)[1]
prod(repProb(2, mgmt[, 1:6]))*dim(mgmt)[1]
prod(repProb(1, mgmt[, 1:6]))*dim(mgmt)[1]
AllSame <- apply(mgmt[, 1:6], 1, sd) == 0
summary(as.factor(apply(mgmt[AllSame, 1:6], 1, paste, collapse="")))
AllSame[is.na(AllSame)] <- F
write.csv(dupSurveys2[,as.numeric(dimnames(dupSurveys2)[[2]])>3],
          file=here("Output","ManagementDuplications.csv"))
mgmt2 <- revCodeH(mgmt, naDrop=F)
dateIX <- match(date(red$StartDate),date(hntrUV$SampleDate)) # adding lighting
mgmt <- cbind(mgmt2$Likert,hntrUV[dateIX,.(Color,N10city)])
complexPlot(rmat=mgmt, xCensor=AllSame, Quest = quest2,
            qVals=questParts2, aVals=answerVals2)
skalesM <- list(
  "1"=names(mgmt)[4:5],
  "2"=names(mgmt)[1:2],
  "1a"=names(mgmt)[c(1,2,4,5)],
  "1b"=names(mgmt)[1:5],
  "1c"=names(mgmt)[1:6]
  )
mgmtNirt <- complexPlot(rmat=mgmt, xCensor = AllSame, scales=skalesM, Quest = quest2,
            qVals=questParts2, aVals=answerVals2)
write.csv(mgmtNirt, here("Output", "ManagementScalesTable.csv"))
############# Response scales by color
predVar <- c("SurveyMinutes", "Cloud_Cover", "Darkness", "Location",
			"FirstTimeVisitor", "YearFirstTripGRTE",
			grep(patt="CB", names(red), value=T),
			grep(patt="Bats", names(red), value=T),
			"Yearborn", "Gender", "Permanentresident",
			"Personalgroup", "Language",
			grep(patt="Illuminance.*", names(red), value=T),
			 "SAz", "SAlt", "MAz", "MAlt", "Phase", "MSD", "SDir",
			"SSky", "SIllm", "MDir", "MSky", "MIllm", "TIllm", "Color",
			"N10city", grep(patt="Brite", names(red), value=T))
sclMake <- function(xDframe, skls, vnames, sknames){
  scls <- matrix(data=0, nrow=dim(xDframe)[2], ncol=length(skls),
                 dimnames=list(dimnames(xDframe)[[2]], sknames))
  for(sc in seq(along=skls))
    scls[skls[[sc]], sc] <- 1
  scls <- sweep(scls,2,apply(scls,2,sum),"/")
  sclDframe <- as.data.frame(as.matrix(xDframe) %*% scls)
  sclDframe <- lapply(sclDframe, as.ordered)
  return(as.data.frame(sclDframe))
}
redSc <- sclMake(stlt[, !(names(stlt) %in% c("Color", "N10city"))],
	skales, names(stlt),
	c("Wild2", "Wild3", "Visi2", "Integ5", "Integ6",
		"Navi2", "Saft2", "Dark2"))
tmp <- !(names(stlt) %in% c("Color", "N10city"))
rdIx <- apply(is.na(stlt[,tmp]), 1, any) | apply(red[,tmp]==3, 1, all)
redStrt <- red[!rdIx, ]
redStrt <- cbind(redStrt, redSc[!rdIx, ]) # after "fixing" NA factors
tmp <- !(names(mgmt) %in% c("Color", "N10city"))
rdIx <- apply(is.na(mgmt[,tmp]), 1, any) | AllSame
mgmtSc <- sclMake(mgmt[, !(names(mgmt) %in% c("Color", "N10city"))],
	skalesM, names(mgmt2),
	c("mWild2", "mLess2", "mWild4", "mWild5", "mWild6"))
redMgm <- red[!rdIx, ]
redMgm <- cbind(redMgm, mgmtSc[!rdIx,])
redStrt$Color <- factor(3-as.integer(redStrt$Color),
	labels=rev(levels(redStrt$Color)))
redMgm$Color <- factor(3-as.integer(redMgm$Color),
	labels=rev(levels(redMgm$Color)))
######### scale CDF plot
ksSize <- 10
scaleFUN <- function(x) sprintf("%.2f", x)
scaleEcdfPlot <- function(redDat, SclNam, LablNam, binaryVar="Color", axSize=30){
	gplt <- ggplot(data=redDat, mapping=aes(x=get(SclNam))) +
		stat_ecdf(mapping=aes(linewidth=get(binaryVar),
			colour=get(binaryVar), alpha=0.3),
			geom="step", pad=T) +
		labs(x=LablNam, colour=binaryVar, linewidth=binaryVar) +
		coord_flip() + scale_y_continuous(0:5) +
		geom_text(size=ksSize, mapping=aes(
			x=3,y=0.9,
			label=sprintf("DTS test:\np=%9.2e",
				with(redDat,
# rev(levels()) used to omit NA that were translated to zeroes
					dts_test(
						get(SclNam)[get(binaryVar)==rev(levels(get(binaryVar)))[1]],
						get(SclNam)[get(binaryVar)==rev(levels(get(binaryVar)))[2]],
						nboots=10000)["P-Value"])))) +
		guides(alpha="none")
	if("Color" == binaryVar)
		gplt <- gplt +
			scale_color_manual(values=c("lightblue", "red")) +
			scale_linewidth_manual(breaks=c("W", "R"), values=c(8,1))
	if("Navigate" == LablNam | "Integrated Management" == LablNam){
		gplt <- gplt + theme(axis.text = element_text(size=16),
			axis.title.x = element_blank(),
			axis.title.y = element_text(size=axSize),
			legend.position = c(0.07, 0.85), legend.title=element_blank(),
			legend.text = element_text(size=axSize))
			} else {
					gplt <- gplt + theme(axis.text = element_text(size=16),
						axis.title.x = element_blank(),
						axis.title.y = element_text(size=axSize),
						legend.position="none")
						}
	return(gplt)
}
redScales <- c("Wild2", "Wild3", "Visi2", "Integ5", "Integ6",
					"Navi2", "Saft2", "Dark2")
if(is.ordered(redStrt$Wild2))
	redStrt[, redScales] <- sapply(redScales,
		function(x) as.numeric(as.character((redStrt[, x]))))
pWld <- scaleEcdfPlot(redStrt, "Wild3", "Wildlife")
pVis <- scaleEcdfPlot(redStrt, "Visi2", "Vision")
pNav <- scaleEcdfPlot(redStrt, "Navi2", "Navigate")
pSaf <- scaleEcdfPlot(redStrt, "Saft2", "Safety")
pDrk <- scaleEcdfPlot(redStrt, "Dark2", "Dark")
redScales <- c("mWild2", "mLess2", "mWild4", "mWild5", "mWild6")
if(is.ordered(redMgm$mWild2))
	redMgm[, redScales] <- sapply(redScales,
		function(x) as.numeric(as.character((redMgm[, x]))))
pMwld <- scaleEcdfPlot(redMgm, "mWild2", "Manage for Wildlife", axSize=24)
pMlss <- scaleEcdfPlot(redMgm, "mLess2", "Manage for Less Light", axSize=24)
pMint <- scaleEcdfPlot(redMgm, "mWild5", "Integrated Management", axSize=24)
pLang <- scaleEcdfPlot(redMgm, "mWild2", "Manage for Wildlife",
	binaryVar="Language", axSize=24)
pLang <- pLang + theme(legend.position = c(0.5, 0.3),
			legend.text = element_text(size=24))
pLang <- pLang + guides(colour=guide_legend(title="Primary\nLanguage"),
		linewidth=guide_legend(title="Primary\nLanguage")) +
	theme(legend.title=element_text(size=24))
#CairoSVG
#pWld + pVis + pNav + pSaf + pDrk +
#	pMwld + pMlss + pMint + plot_layout(ncol=2)
CairoPNG(filename = here("Output", "ScaleByColor.png"),
         width=1200,height=1500)
pWld + pVis + pNav + pSaf + plot_layout(ncol=1)
dev.off()
CairoPNG(filename = here("Output", "ManageScaleByColor.png"),
         width=1200,height=1500)
pMwld + pMlss + pMint + pLang + plot_layout(ncol=1)
dev.off()
red$Color <- factor(3-as.integer(red$Color),
	labels=rev(levels(red$Color)))
hplt <- scaleEcdfPlot(red[red$SAlt < -6, ], "Illuminance_Horizontal",
	"Horizontal Illuminance (lux)")
hplt + scale_x_continuous(trans="log10") +
	labs(caption="Fraction of survey locations after civil twilight (N=168)") +
	theme(plot.caption=element_text(size=30, hjust=0.5),
		legend.position = c(0.07, 0.85), legend.title=element_blank(),
		legend.text = element_text(size=24))
tmp <- split(red, red$Color)
tmp2 <- lapply(tmp, function(x) x[x$SAlt < -6, ])
sapply(tmp2, function(x) summary(x$Illuminance_Horizontal))
######### variable importance analyses
plllVimp <- function(tgt, redDat, predVar, parReps, cfctrl){
  handlers("progress")
  frmla <- as.formula(paste0(tgt, " ~ ", paste(predVar, collapse=" + ")))
  dta <- redDat[, c(tgt, predVar)]
  tstrt <- Sys.time()
  with_progress({
    p <- progressor(along = parReps)
    varIMP <- future_lapply(parReps, function(x, ...) {
      cfRed <- cforest(frmla, control=cfctrl,
                       data=dta)
      target <- redDat[[tgt]]
      if(is.factor(target)){
        tPred <- predict(cfRed, OOB=T, type="prob")
        tmp <- list(
          vImp = varimp(cfRed, nperm=1, conditional=T),
# largest positive difference from the mean conditional probability
          predMiss = (summary(as.factor(
            abs(apply(sweep(jnk<-sapply(tPred, function(x)
              x[1,]), 1, apply(jnk,1,mean), "-"), 2, which.max) -
                as.numeric(target))))/length(target))[
                  1:as.numeric(substr(tgt, nchar(tgt), nchar(tgt)))])
      } else if(is.numeric(target)){
        tmp <- list(
          vImp = varimp(cfRed, nperm=1, conditional=T),
          varExpl = list(
            madFrac = mad(target - predict(cfRed, OOB=T))/
              mad(target),
            varFrac = var(target - predict(cfRed, OOB=T))/
              var(target)))
      }
      p(sprintf("x=%g, tdiff=%g", x, Sys.time()-tstrt))
      return(tmp)
    }, future.seed=T)
  })
  print(Sys.time() - tstrt)
  return(varIMP)
} # plllVimp
sapply(redStrt[, sapply(redStrt, is.numeric)], function(v)
  range(v[v>0], na.rm=T))
nWorkers <- as.numeric(as.character(attr(plan(), "call"))[3])
parReps <- 1: (nWorkers*ceiling(100/nWorkers)) # redStrt 40, redMgm 80 min
redStrt <- redStrt[!apply(is.na(redStrt), 1, any), ]
redMgm <- redMgm[!apply(is.na(redMgm), 1, any),
	-grep(patt="Streetlights", names(redMgm))]
write.csv(predVar, here("Output", "Independent Variables.csv"))
tstrtAll <- Sys.time()
Wild3 <- plllVimp("Wild3", redStrt, predVar, parReps, cfctrl)
save(Wild3, file=here("Output","Wild3b.RData"))
Visi2 <- plllVimp("Visi2", redStrt, predVar, parReps, cfctrl)
save(Visi2, file=here("Output","Visi2b.RData"))
Integ6 <- plllVimp("Integ6", redStrt, predVar, parReps, cfctrl)
save(Integ6, file=here("Output","Integ6b.RData"))
Navg2 <- plllVimp("Navi2", redStrt, predVar, parReps, cfctrl)
save(Navg2, file=here("Output","Navg2b.RData"))
Saft2 <- plllVimp("Saft2", redStrt, predVar, parReps, cfctrl)
save(Saft2, file=here("Output","Saft2b.RData"))
Dark2 <- plllVimp("Dark2", redStrt, predVar, parReps, cfctrl)
save(Dark2, file=here("Output","Dark2b.RData"))
Mwild2 <- plllVimp("mWild2", redMgm, predVar, parReps, cfctrl)
save(Mwild2, file=here("Output","Mwild2b.RData"))
Mless2 <- plllVimp("mLess2", redMgm, predVar, parReps, cfctrl)
save(Mless2, file=here("Output","Mless2b.RData"))
Mwild4 <- plllVimp("mWild4", redMgm, predVar, parReps, cfctrl)
save(Mwild4, file=here("Output","Mwild4b.RData"))
Mwild5 <- plllVimp("mWild5", redMgm, predVar, parReps, cfctrl)
save(Mwild5, file=here("Output","Mwild5b.RData"))
Mwild6 <- plllVimp("mWild6", redMgm, predVar, parReps, cfctrl)
save(Mwild6, file=here("Output","Mwild6b.RData"))
print(Sys.time()-tstrtAll)
########## Investigating directional illuminance results
with(redStrt, plot(log10(Illuminance_Vertical_90),
  as.numeric(as.character(Wild3))))
jnk <- with(redStrt, supsmu(log10(Illuminance_Vertical_90),
          as.numeric(as.character(Wild3))))
lines(jnk$x, jnk$y)
## visitors at the store and roving locations were less supportive of mWild2
## when north illuminance was high (daytime)
locPlots <- function(redDat, vars=c("SAlt", "Saft2"), splt="Location", keep=2:4){
	tst <- split(redDat[, vars], redDat[, splt])
	par(mfrow=c(3,1))
	lapply(tst[keep], function(x){
		plot(x[,1], x[,2], ylim=c(1,5), xlab=vars[1], ylab=vars[2])
		jnk <- supsmu(x[, 1], x[, 2], span=0.3)
		lines(jnk$x, jnk$y)
		})
	par(mfrow=c(1,1))
}
locPlots(redMgm, vars=c("SAlt", "mWild2"))
with(redMgm, plot(Illuminance_Vertical_0, SAlt, log="x"))
# First time visitor and Navigation
tmp <-with(redStrt,table(FirstTimeVisitor, Navi2))
sprintf("%0.2f", apply(sweep(tmp, 1, apply(tmp, 1, sum), "/"), 2, diff))[
	seq(from=2, to=18, by=2)]
varSumry <- function(fname){
  vnam <- load(here("Output", fname))
  resList <- get(vnam)
  vImp <- reduce(lapply(resList, "[[", "vImp"), rbind)
  zVimp <- apply(vImp,2,mean)/apply(vImp,2,sd)
  zVimp <- sort(zVimp, decreasing = T)
  zVimp <- c(
    apply(reduce(lapply(resList, "[[", "predMiss"), rbind),
                 2, mean), zVimp)
  return(zVimp)
}
Flist <- c("Integ6b.RData", "Wild3b.RData", "Visi2b.RData", "Navg2b.RData",
           "Saft2b.RData", "Dark2b.RData", "Mwild2b.RData",
           "Mless2b.RData", "Mwild4b.RData", "Mwild5b.RData",
           "Mwild6b.RData")
for(fl in Flist){
  if(match(fl, Flist)==1){
    tmp <- varSumry(fl)
    zVarList <- matrix(NA, nrow=length(tmp), ncol=length(Flist),
                       dimnames = list(names(tmp), Flist))
    zVarList[names(tmp), fl] <- tmp
  } else{
      tmp <- varSumry(fl)
      zVarList[names(tmp), fl] <- tmp
  }
}
zVarList[zVarList<0] <- 0
# first six rows are fractions of predictions that missed the mark
# by I-1 steps
ixSorted <- c(c(1:6),
  order(apply(zVarList[-(1:6), 1:6], 1, max), decreasing = T) + 6)
write.csv(zVarList[ixSorted, 1:6],
          file=here("Output","VarImpZStreetlightC.csv"))
ixSorted <- c(c(1:6),
              order(apply(zVarList[-(1:6), 7:11], 1, max), decreasing = T) + 6)
write.csv(zVarList[ixSorted, 7:11],
          file=here("Output","VarImpZManagementB.csv"))
zVarList <-
  zVarList[c(c(1:6),
             order(apply(zVarList[-(1:6),], 1, max), decreasing = T) + 6), ]
#ix <- apply(zVarList>2, 1, any); ix[1:6] <- T
#zVarList <- zVarList[ix, ]
dimnames(zVarList)[[2]] <- c("Integrated(1,2,6,7,8,9)", "Wildlife(6,7,8)", 
                             "Vision(1,2)", "Navigate(5,13)", "Safety(3,4)", "Dark(10,11)",
                             "mWildlife(4,5)", "mLessLight(1,2)", "mCombined4(1,2,4,5)",
                             "mCombined5(1,2,3,4,5)", "mCombined6(1,2,3,4,5,6)")
write.csv(zVarList,
          file=here("Output","VarImpZscoresC.csv"))
partDependProbs <- function(tgt, redDat, xVar){
# if tgt is numeric, returns the mean difference
# if tgt is ordered, returns a vector of differences for each level
  frmla <- as.formula(paste0(tgt, " ~ ", paste(predVar, collapse=" + ")))
  dta <- redDat[, c(tgt, predVar)]
  cfRed <- cforest(frmla, control=cfctrl,
                   data=dta)
  xvLvls <- sort(unique(redDat[[xVar]]))
  xVals <- rev(xvLvls)[1:2]
#selects the last two options for xVar
  redDat[["xVar"]] <- xVals[2]
  tPred0 <- predict(cfRed, OOB=T, type="prob", newdata=redDat)
  tPred0 <- reduce(tPred0, rbind)
  redDat[[xVar]] <- xVals[1]
  tPred1 <- predict(cfRed, OOB=T, type="prob", newdata=redDat)
  tPred1 <- reduce(tPred1, rbind)
  xiTion <- tPred1 - tPred0
  return(apply(xiTion, 2, mean))
}
if(is.numeric(redStrt$Wild2))
	redStrt[dim(redStrt)[2]:(dim(redStrt)[2]-8)] <- lapply(rev(redStrt)[1:8],
		function(x) as.ordered(x))
if(is.ordered(redStrt$Wild2))
	redStrt[dim(redStrt)[2]:(dim(redStrt)[2]-8)] <- lapply(rev(redStrt)[1:8],
		function(x) as.numeric(as.character(x)))
nIter <- 40
cfctrl <- cforest_control(ntree=1000)
tstrt <- Sys.time()
Int6Lst <- future_lapply(1:nIter, function(x)
	partDependProbs("Integ6", redStrt, "B2Brite"), future.seed=T)
print(Int6Avg <- apply(reduce(Int6Lst, rbind),2, mean))
## the next calculation is slightly different than the average real number predictions
if(length(Int6Avg)>1) print(
	sum(as.numeric(sub(patt="Integ6\\.", repl="", names(Int6Avg))) * Int6Avg))
print(Sys.time()-tstrt)
tstrt <- Sys.time()
NaviLst <- future_lapply(1:nIter, function(x)
	partDependProbs("Navi2", redStrt, "FirstTimeVisitor"), future.seed=T)
print(NaviAvg <- apply(reduce(NaviLst, rbind),2, mean))
FirstLessLst <- future_lapply(1:nIter, function(x)
	partDependProbs("mLess2", redMgm, "FirstTimeVisitor"), future.seed=T)
print(FirstLessAvg <- apply(reduce(FirstLessLst, rbind),2, mean))
VisiLst <- future_lapply(1:nIter, function(x)
	partDependProbs("Visi2", redStrt, "CBActivity_walk"), future.seed=T)
print(NaviAvg <- apply(reduce(VisiLst, rbind),2, mean))
NaviLstB <- future_lapply(1:nIter, function(x)
	partDependProbs("Navi2", redStrt, "SawBats"), future.seed=T)
print(NaviAvgB <- apply(reduce(NaviLstB, rbind),2, mean))
SaftLst <- future_lapply(1:nIter, function(x)
	partDependProbs("Saft2", redStrt, "SawBats"), future.seed=T)
print(SaftAvg <- apply(reduce(SaftLst, rbind),2, mean))
SolarLst <- future_lapply(1:nIter, function(x)
	partDependProbs("mWild2", redMgm, "SAlt"), future.seed=T)
print(SolarAvg <- apply(reduce(SolarLst, rbind),2, mean))
LangLst <- future_lapply(1:nIter, function(x)
	partDependProbs("mWild2", redMgm, "Language"), future.seed=T)
print(LangAvg <- apply(reduce(LangLst, rbind),2, mean))
## the next calculation is slightly different than the average real number predictions
if(length(Int6Avg)>1) print(
	sum(as.numeric(sub(patt="Integ6\\.", repl="", names(Int6Avg))) * Int6Avg))
print(Sys.time()-tstrt)
## graphical explorations
ggplot(data=redStrt, mapping=aes(x=as.factor(Saft2))) +
  geom_bar(mapping=aes(fill=as.factor(B2Brite)), position="dodge")
# what is going on with mWild2 and Language?
pLang <- scaleEcdfPlot(redMgm, "mWild2", "Manage for Wildlife", "Language")
pLang <- pLang + theme(axis.text = element_text(size=16),
			axis.title.x = element_blank(),
			axis.title.y = element_text(size=axSize),
			legend.position = c(0.85, 0.2), legend.title=element_blank(),
			legend.text = element_text(size=axSize))
CairoPNG(filename = here("Output", "mWildLang.png"),
         width=1200,height=1500)
pLang
dev.off()
ggplot(data=redMgm, mapping=aes(x=as.factor(mWild2))) +
  geom_bar(mapping=aes(fill=as.factor(Language)), position="dodge")
ggplot(data=redMgm, mapping=aes(x=Language, y=mWild2)) + geom_boxplot()
lapply(split(as.numeric(as.character(redMgm$mWild2)), redMgm$Language), summary)
tst <- with(redMgm, split(mWild2, Language))
dts_test(tst[[1]], tst[[2]])
dts_test(tst[[1]], tst[[3]])
dts_test(tst[[2]], tst[[3]])
# Illuminance North?
tst <- split(redMgm[, c("Illuminance_Vertical_0", "mWild2")], redMgm$Location)
par(mfrow=c(3,1))
lapply(tst[-1], function(x) plot(x[,1], x[,2])
# lower astronomical darkness increases Vision scores
tst <- sapply(split(redStrt$Visi2, redStrt$Darkness), summary)
apply(sweep(tst, 2, apply(tst, 2, sum), "/"), 2, cumsum)
