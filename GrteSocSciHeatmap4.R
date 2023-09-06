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
  red <- read_excel(here("Data", "GRTE lighting.xlsx"));
  red[red==999] <- NA
  red[,grep(patt="awe",names(red),value = T,ignore.case = T)] <- NULL
  red[,grep(patt="Biodiver",names(red),value=T,ignore.case=T)] <- NULL
  red[,grep(patt="Non_resp",names(red),value=T,ignore.case=T)] <- NULL
  red[,grep(patt="BEL",names(red),value=T)] <- NULL
  red[,grep(patt="center",names(red),value=T,ignore.case=T)] <- NULL
  red$SurveyMinutes <- as.numeric(red$EndDate-red$StartDate)
  red$SawBats <- !is.na(red$NumberBats)
  varKeep <- c("StartDate","SurveyMinutes","Cloud_Cover","Darkness",
               "StreetlightCondition","Hue","Intensity","Location",
               "NightSkyCondition","FirstTimeVisitor",
               "YearFirstTripGRTE","PreviousVisitCB",
               "PreviousEveningCB","CampCB_ThisTrip",
               "NumberNightsCB","PreviousCampCB",
               "YearFirstCampCB","CBActivity_ranger",
               "CBActivity_walk","CBActivity_stargze",
               "CBActivity_Other",
               "NumberBats","NightSky_important","NightSky_expect",
               "NightSky_experience","NightSky_pristine",
               "NightSky_reputation","NightSky_futuregenerations",
               "NightSky_resources","NightSky_modifylight",
               "Streetlights_activitiesmorepleasurable",
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
               "Acceptable_Sky_Conditions","Yearborn","Gender",
               "Permanentresident","PrimaryZip",
               "Personalgroup","Language",
               "Illuminance_Vertical_0","Illuminance_Vertical_90",
               "Illuminance_Vertical_180","Illuminance_Vertical_270",
               "Illuminance_Vertical_Sun",
               "Illuminance_Vertical_BrightestLight",
               "Illuminance_Horizontal","SawBats")
  red <- red[,varKeep]
  factorVars <- c("Darkness","NightSkyCondition",
                  "StreetlightCondition","Hue","Location",
                  "FirstTimeVisitor","PreviousVisitCB",
                  "PreviousEveningCB","CampCB_ThisTrip",
                  "PreviousCampCB","CBActivity_ranger",
                  "CBActivity_walk","CBActivity_stargze",
                  "CBActivity_Other","Gender","SawBats",
                  "Permanentresident")
  red[,factorVars] <- lapply(red[,factorVars],as.factor)
  red$StartDate <- force_tz(red$StartDate,tzone="America/Denver")
  red$Illuminance_Vertical_Integrated <-
    as.matrix(red[,grep(patt="Illuminance",names(red))][,1:4]) %*% rep(1,4)
  red$PrimaryZip <- as.factor(substring(red$PrimaryZip,1,1))
  red$NumberNightsCB <- as.numeric(red$NumberNightsCB)
} # prepare red data matrix
# now check the alternative data matrix
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
tst <- ymd_hms(red$StartDate)
unique(date(tst)) # Surveys administered on how many dates?
# data set for aisp calculations
red2 <- data.frame(red[,zQix])
rcAns <- revCodeH(data.frame(red2))
red2r <- rcAns$Likert
negQix <- grep(patt="-",names(red2r))
posQix <- seq(along=names(red2r)); posQix <- posQix[!(posQix %in% negQix)]
red2[, negQix] <- 6 - red2[, negQix]
All3S <- apply(red2r==3,1,all) # also get rid of all "3"
# Odds of all 3s by chance
prod(apply(red2r[!All3S,]==3,2,sum)/dim(red2r[!All3S,])[1])
# Odds of any replicated surveys
redTmp <- data.frame(lapply(red[, grep(patt="Streetlights_", names(red))], as.factor))
streetFreqs <- sapply(redTmp, summary)/dim(redTmp)[1]
1 - (1 - prod(apply(streetFreqs^2, 2, sum)))^dim(redTmp)[1]
###### add lighting and weather to red2
dateIX <- match(date(red$StartDate),date(hntrUV$SampleDate)) # adding lighting
All3S <- apply(red2==3,1,all)
All3S[is.na(All3S)] <- F # get rid of all "3" for complexPlot
red2 <- cbind(red2,hntrUV[dateIX,.(Color,N10city)])
# add weather and natural illumination
wthr <- fread(here("Data", "MooseWeatherHr.csv"))
wthr[,ytime :=as.POSIXct(
  julianHr*3600,origin=as.POSIXct("1970-01-01",tz="GMT"),tz="America/Denver")]
setkey(wthr,ytime)
# look at duplicates by color
negQix <- negQix[-6] # get rid of affctGE-
sum(apply(red2[,posQix]>3,1,all),na.rm=T) # how many surveys had all >3 for "+Q"
sum(apply(red2[,negQix]==5,1,all),na.rm=T) # how many surveys had all 1 for "-Q"
dateIX2 <- apply(red2[,posQix]>3,1,all); dateIX2[is.na(dateIX2)] <- F
summary(red2[dateIX2,"Color"])
dateIX2 <- apply(red2[,negQix]==5,1,all); dateIX2[is.na(dateIX2)] <- F
summary(red2[dateIX2,"Color"])
safety5Ix <- red2$Streetlights_lesssafe==1 
safety5Ix[is.na(safety5Ix)] <- F
naCount <- apply(is.na(red2), 1, sum, na.rm=T)
red5 <- t(red2)
dimnames(red5)[[2]] <- rep(1,dim(red5)[2])
dupSurveys <- itemReplicates(
  red5[grep(patt="Streetlights", dimnames(red5)[[1]]), ], reGroup = F)
write.csv(dupSurveys[,as.numeric(dimnames(dupSurveys)[[2]])>1],
          file=here("Output", "StreetlightDuplications.csv"))
skales <- list(
  names(red2)[6:7],
  names(red2)[6:8],
  names(red2)[1:2],
  names(red2)[c(1,2,6:8)],
  names(red2)[c(1,2,6:9)],
  names(red2)[c(5,13)],
  names(red2)[c(3,4)],
  names(red2)[c(10,11)]
)
names(skales) <- c("1", "1a", "2", "1b", "1c", "3", "4", "5")
complexPlot(rmat=red2, xCensor = apply(red2==3,1,all),
  scales=skales, Quest = quest1, qVals=questParts1, aVals=answerVals1)
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
mgmt2 <- revCodeH(mgmt, naDrop=F);
dateIX <- match(date(red$StartDate),date(hntrUV$SampleDate)) # adding lighting
mgmt <- cbind(mgmt,hntrUV[dateIX,.(Color,N10city)])
complexPlot(rmat=mgmt, xCensor=AllSame, Quest = quest2,
            qVals=questParts2, aVals=answerVals2)
skalesM <- list(
  names(mgmt)[4:5],
  names(mgmt)[1:2],
  names(mgmt)[c(1,2,4,5)],
  names(mgmt)[1:5],
  names(mgmt)[1:6]
  )
names(skalesM) <- c("1", "2", "1a", "1b", "1c")
complexPlot(rmat=mgmt, xCensor = AllSame, scales=skalesM, Quest = quest2,
            qVals=questParts2, aVals=answerVals2)
############# Response scales by color
predVar <- c("SurveyMinutes", "Cloud_Cover", "Darkness",
             "NightSkyCondition", "FirstTimeVisitor", "YearFirstTripGRTE",
             "PreviousVisitCB", "PreviousEveningCB", "CampCB_ThisTrip",
             "NumberNightsCB", "PreviousCampCB", "YearFirstCampCB",
             "CBActivity_ranger", "CBActivity_walk", "CBActivity_stargze",
             "Yearborn", "Gender", "Permanentresident", "PrimaryZip",
             "Personalgroup", "Language",
             names(red)[grep(patt="Illuminance.*", names(red))],
             "SAz", "SAlt", "MAz", "MAlt", "Phase", "MSD", "SDir",
             "SSky", "SIllm", "MDir", "MSky", "MIllm", "TIllm", "Color",
             "N10city")
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
redSc <- sclMake(red2[, !(names(red2) %in% c("Color", "N10city"))],
                skales, names(red2),
                c("Wild2", "Wild3", "Visi2", "Integ5", "Integ6",
                  "Navi2", "Saft2", "Dark2"))
rdIx <- apply(is.na(red[,zQix]), 1, any) | apply(red[,zQix]==3, 1, all)
redStrt <- red[!rdIx, ]
##########################
# REVISIT THIS CONVERSION OF FACTORS TO INTEGERS, AND CHANGING NA TO -1
##########################
redStrt[sapply(redStrt, is.factor)] <- lapply(
  redStrt[sapply(redStrt, is.factor)], as.numeric)
redStrt[is.na(redStrt)] <- -1
redStrt <- cbind(redStrt, redSc[ !rdIx, ]) # after "fixing" NA factors
rdIx <- apply(is.na(red[,mQix]), 1, any) | AllSame
mgmtSc <- sclMake(mgmt2$Likert[!rdIx, ], skalesM, names(mgmt2),
                  c("mWild2", "mLess2", "mWild4", "mWild5", "mWild6"))
redMgm <- red[!rdIx, ]
redMgm[sapply(redMgm, is.factor)] <- lapply(
  redMgm[sapply(redMgm, is.factor)], as.numeric)
redMgm[is.na(redMgm)] <- -1
redMgm <- cbind(redMgm, mgmtSc)
######### variable importance analyses
plllVimp <- function(tgt, redDat, predVar, parReps, cfctrl, nWorkers){
  handlers("progress", "beepr")
  plan(multisession, workers=nWorkers)
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
    })
  })
  print(Sys.time() - tstrt)
  plan(sequential)
  return(varIMP)
} # plllVimp
sapply(redStrt[, sapply(redStrt, is.numeric)], function(v)
  range(v[v>0]))
cfctrl <- cforest_control(ntree=1000)
nWorkers <- 14
parReps <- 1: (nWorkers*round(100/nWorkers)) # redStrt 40, redMgm 80 min
Wild3 <- plllVimp("Wild3", redStrt, predVar, parReps, cfctrl, nWorkers)
save(Wild3, file=here("Output","Wild3b.RData")); plan(.skip=T, .cleanup=T)
Visi2 <- plllVimp("Visi2", redStrt, predVar, parReps, cfctrl, nWorkers)
save(Visi2, file=here("Output","Visi2b.RData")); plan(.skip=T, .cleanup=T)
Integ6 <- plllVimp("Integ6", redStrt, predVar, parReps, cfctrl, nWorkers)
save(Integ6, file=here("Output","Integ6b.RData")); plan(.skip=T, .cleanup=T)
Navg2 <- plllVimp("Navi2", redStrt, predVar, parReps, cfctrl, nWorkers)
save(Navg2, file=here("Output","Navg2b.RData")); plan(.skip=T, .cleanup=T)
Saft2 <- plllVimp("Saft2", redStrt, predVar, parReps, cfctrl, nWorkers)
save(Saft2, file=here("Output","Saft2b.RData")); plan(.skip=T, .cleanup=T)
Dark2 <- plllVimp("Dark2", redStrt, predVar, parReps, cfctrl, nWorkers)
save(Dark2, file=here("Output","Dark2b.RData")); plan(.skip=T, .cleanup=T)
Mwild2 <- plllVimp("mWild2", redMgm, predVar, parReps, cfctrl, nWorkers)
save(Mwild2, file=here("Output","Mwild2b.RData")); plan(.skip=T, .cleanup=T)
Mless2 <- plllVimp("mLess2", redMgm, predVar, parReps, cfctrl, nWorkers)
save(Mless2, file=here("Output","Mless2b.RData")); plan(.skip=T, .cleanup=T)
Mwild4 <- plllVimp("mWild4", redMgm, predVar, parReps, cfctrl, nWorkers)
save(Mwild4, file=here("Output","Mwild4b.RData")); plan(.skip=T, .cleanup=T)
Mwild5 <- plllVimp("mWild5", redMgm, predVar, parReps, cfctrl, nWorkers)
save(Mwild5, file=here("Output","Mwild5b.RData")); plan(.skip=T, .cleanup=T)
Mwild6 <- plllVimp("mWild6", redMgm, predVar, parReps, cfctrl, nWorkers)
save(Mwild6, file=here("Output","Mwild6b.RData")); plan(.skip=T, .cleanup=T)
##########
with(redStrt, plot(log10(Illuminance_Vertical_90),
  as.numeric(as.character(Wild3))))
jnk <- with(redStrt, supsmu(log10(Illuminance_Vertical_90),
          as.numeric(as.character(Wild3))))
lines(jnk$x, jnk$y)
with(redMgm, plot(log10(Illuminance_Vertical_0),
  as.numeric(as.character(mWild2))))
jnk <- with(redMgm, supsmu(log10(Illuminance_Vertical_0),
          as.numeric(as.character(mWild2))))
lines(jnk$x, jnk$y)

varSumry <- function(fname){
  vnam <- load(here("Output", fname))
  resList <- get(vnam)
  vImp <- reduce(lapply(resList, "[[", "vImp"), rbind)
  zVimp <- apply(vImp,2,mean)/apply(vImp,2,sd)
  zVimp <- zVimp[order(zVimp, decreasing = T)]
  zVimp <- c(
    apply(reduce(lapply(resList, "[[", "predMiss"), rbind),
                 2, mean), zVimp)
  return(zVimp)
}
Flist <- c("Integ6.RData", "Wild3.RData", "Visi2.RData", "Navg2.RData",
           "Saft2.RData", "Dark2.RData", "Mwild2b.RData",
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
          file=here("Output","VarImpZStreetlightB.csv"))
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
                             "Vision(1,2)", "Navigate(5, 13)", "Safety(3,4)",
                             "mWildlife(4,5)", "mLessLight(1,2)", "mCombined4(1,2,4,5)",
                             "mCombined5(1,2,3,4,5)", "mCombined6(1,2,3,4,5,6)")
write.csv(zVarList,
          file=here("Output","VarImpZscoresB.csv"))
partDependProbs <- function(tgt, redDat, xVar){
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
nIter <- 40
Navi2Avg <- partDependProbs("Navi2", redStrt, "FirstTimeVisitor")
for(itr in 2:nIter)
  Navi2Avg <- Navi2Avg + partDependProbs("Navi2", redStrt, "FirstTimeVisitor")
Navi2Avg <- Navi2Avg/nIter
sum(as.numeric(sub(patt="Navi2\\.", repl="", names(Navi2Avg))) * 
  Navi2Avg)
sum(c(-0.00040, -0.00014, 0.0015, 0.0019, 0.000022, -0.0034, -0.0042, 0.00075, 0.0070))
Less2Avg <- partDependProbs("mLess2", redMgm, "Color")
for(itr in 2:nIter)
  Less2Avg <- Less2Avg + partDependProbs("mLess2", redMgm, "Color")
Less2Avg <- Less2Avg/nIter
sum(as.numeric(sub(patt="mLess2\\.", repl="", names(Less2Avg))) * 
  Less2Avg)
ggplot(data=redMgm, mapping=aes(x=as.factor(mLess2))) +
  geom_bar(mapping=aes(fill=as.factor(Color)), position="dodge")

