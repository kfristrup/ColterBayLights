coefHdist <- function(xmat){
  # generates a 1-D ordering based on MDS of 1-Hij
  if(diff(dim(xmat)) > 0) xmat <- t(xmat)
  nvars <- dim(xmat)[2]
  xmCh <- coefH(na.exclude(xmat), se = F, results = F)
  oneDlocs <- cmdscale(1-xmCh$Hij, 1)
  return(oneDlocs[ , 1])
}
lastGrp <- function(firstGrp, mrgs){
  # finds the last merge that contains variables to be merged
  tst <- firstGrp
  while(!is.na(tst)){
    tst <- match(tst, mrgs) %% dim(mrgs)[1]
    if(!is.na(tst))
      firstGrp <- tst
  }
  return(firstGrp)
}
findNeg <- function(mrgs){
  # returns a vector of all cases connected to the last entry in the
  # merge matrix (Nx2, conforming with hclust) <negative values>
  mDims <- dim(mrgs)
  if(is.null(mDims)){
    lastRow <- 1
    lastVals <- mrgs
  } else{
    lastRow <- dim(mrgs)[1]
    lastVals <- mrgs[lastRow,]
  }
  posIx <- (0 < lastVals) # T/F, length=2
  retVal <- lastVals[!posIx]
  if(any(posIx))
    retVal <- c(retVal, unlist(
      sapply(lastVals[posIx], function(rw) findNeg(mrgs[1:rw,]))))
  return(retVal)
}
findRows <- function(mrgs){
  # returns a vector of all rows connected to the last entry in the
  # merge matrix (Nx2, conforming with hclust)
  mDims <- dim(mrgs)
  if(is.null(mDims)){
    lastRow <- 1
    lastVals <- mrgs
  } else {
    lastRow <- dim(mrgs)[1]
    lastVals <- mrgs[lastRow,]
  }
  retVal <- lastRow
  posIx <- (0 < lastVals) # T/F, length=2
  if(any(posIx))
    retVal <- c(retVal, unlist(
      sapply(lastVals[posIx], function(rw) findRows(mrgs[1:rw,]))))
  return(retVal)
}
hclGrps <- function(mrgs){
  grps <- vector(mode="list", length=0)
  # mrgs is Nx2 as in hclust()
  nextGrp <- 1
  rwsLeft <- 1:sum(!apply(is.na(mrgs),1,any))
  while(length(rwsLeft)){
    rws <- findRows(mrgs[1:max(rwsLeft), ])
    grps[[nextGrp]] <- mrgs[rws,][mrgs[rws,]<0]
    names(grps)[nextGrp] <- max(rws)
    nextGrp <- nextGrp + 1
    rwsLeft <- rwsLeft[!(rwsLeft %in% rws)]
  }
  return(grps)
}
maxH <- function(mList, rMatR){
  mList <- lapply(mList, abs) # convert merged negative numbers to positive
  nGrps <- sum(sapply(mList, length) > 1)
  lLgth <- length(mList)
  lPairs <- cbind(
    lFrom <- unlist(sapply(1:nGrps, function(x) rep(x, lLgth - x))),
    unlist(sapply(1:max(lFrom), function(x) (1+x):lLgth)))
  xH <- apply(lPairs, 1, function(x){
    tst <- coefH(rMatR[, c(mList[[x[1]]], mList[[x[2]]])], se = F, results = F)
    return(c(x, tst$H))
  })
  imx <- which.max(xH[3, ])
  mrgd2 <- sapply(xH[1:2, imx],
                  function(x) ifelse(1==length(mList[[x]]),
                                     -mList[[x]], names(mList)[x]))
  return(list(hMax=1-xH[3, imx], mErge=as.numeric(mrgd2)))
}
removeHij <- function(vrs, lowTri){
  if(1 > length(lowTri))
    return()
  # removes elements with pairs containing one of vrs
  # assumes a vector structures like mat[low.tri[mat]]
  vrs <- abs(vrs)
  lTlen <- length(lowTri)
  nVars <- (1 + sqrt(1 + 8*lTlen)) / 2
  nix <- NULL
  for(rw in vrs){
    rowMxs <- c(0, cumsum(seq(from=nVars - 1, to=1)))
    if(rw > 1)
      nix <- c(nix, (rw-1):1 + rowMxs[1:(rw-1)])
    if(rw < nVars)
      nix <- c(nix, (1 + rowMxs[rw]):rowMxs[rw + 1])
  }
  return(lowTri[-unique(nix)])
}
coefHclust <- function(rmat){
  allCols <- seq(along=rmat[1, ])
  rmatr <- rmat[!apply(is.na(rmat), 1, any), ]
  hcCH <- vector(mode="list", length = 7) # class hclust is list of 7
  names(hcCH) <- c("merge", "height", "order", "labels",
                   "method", "call", "dist.method")
  attr(hcCH, "class") <- "hclust"
  hcCH[["merge"]] <- array(NA, dim=c(dim(rmatr)[2] - 1, 2))
  hcCH[["height"]] <- rep(0, dim(rmatr)[2] - 1)
  hcCH[["order"]] <- rep(0, dim(rmatr)[1])
  hcCH[["labels"]] <- sub(patt="Streetlights_", repl="",
                          dimnames(rmatr)[[2]][1:13])
  hcCH[["call"]] <- "MrKurt"
  hcCH[["dist.method"]] <- "1MinusCoefH"
  mergNum <- 1
  vHij <- coefH(rmatr, se = F, results = F)$Hij
  vHij <- vHij[lower.tri(vHij)] # n*(n-1)/2
  for(cl in 1:length(hcCH[["height"]])){
    mergdCols <- NA
    if(!all(is.na(hcCH[["merge"]]))){
      mergList <- hclGrps(hcCH[["merge"]])
      mergdCols <- unlist(mergList);
      mergdCols <- -mergdCols[mergdCols<0]
      availCols <- allCols[!(allCols %in% mergdCols)]
      mergList <- c(mergList, availCols)
      bestMerg <- maxH(mergList, rmatr)
    } else {
      availCols <- allCols
      bestMerg <- list(hMax = 1, mErg <- c(NA, NA))
    }
    if(hijLen <- length(vHij)){
      # find the best new cluster pair in the unclustered instances
      bestPair <- which.max(vHij)
      # N*(N-1)/2 = k
      # N^2 - N - 2k = 0
      # N = (1 + sqrt(1 + 8k))/2
      nVars <- (1 + sqrt(1 + 8*hijLen)) / 2
      rowMxs <- cumsum(seq(from=nVars - 1, to=1))
      bclmn <- which.min(rowMxs < bestPair) # first rowMxs > bestPair
      brw <- ifelse(bclmn > 1,
                    bestPair - rowMxs[bclmn - 1] + bclmn,
                    bestPair + 1)
      mergeTmp <- -c(availCols[brw], availCols[bclmn])
      heightTmp <- 1 - vHij[bestPair]
    }
    if(heightTmp < bestMerg$hMax){ # pair forms new group
      hcCH[["merge"]][mergNum, ] <- mergeTmp
      hcCH[["height"]][mergNum] <- heightTmp
      vHij <- removeHij(c(brw, bclmn), vHij)
    } else { # merge one group and one var or two groups
      tst <- bestMerg$mErge < 0
      if(any(tst)){
        rmIxs <- match(abs(bestMerg$mErge[tst]), availCols)
        vHij <- removeHij(rmIxs, vHij)
      }
      hcCH[["merge"]][mergNum, ] <- bestMerg$mErge
      hcCH[["height"]][mergNum] <- bestMerg$hMax
    }
    # remove the merged variable from rmatr
    mergNum <- mergNum + 1
  }
  return(hcCH)
}
clustGramKf <- function (xmat, rmRows, Q1, Qparts, ansVals) {
  if(diff(dim(xmat)) > 0) xmat <- t(xmat) #rows are surveys
  cexAll <- 1
  numQ <- length(Qparts)
  colDrop <- match(c("Color", "N10city"), names(xmat))
  if(is.null(rmRows))
      xmatr <- revCodeH(xmat[, -colDrop], naDrop=F)
  else
      xmatr <- revCodeH(xmat[!rmRows, -colDrop], naDrop=F)
# Create the survey dendrogram
# Add Color*50 for main separation, NA,_8 for secondary
  xmatTmp <- cbind(xmat,
    srvSep = as.numeric(xmat$Color) * 50)
  xmatTmp[is.na(xmatTmp)] <- 8
  rowMns <- apply(xmatTmp[,-colDrop], 1, mean, na.rm=T)
  rowHc <- hclust(dist(xmatTmp[, -colDrop], method = "manhattan"))
  hccDgram <- as.dendrogram(rowHc)
  hccDgram <- reorder(hccDgram, rowMns)
# tweak to improve graphical similarity of red v. white subsets for Management
  if(all((1:6) == grep(patt="Management_", dimnames(xmat)[[2]]))){
    tmp <- hccDgram[[c(2,2,2,2,1)]]
    hccDgram[[c(2,2,2,2,1)]] <- hccDgram[[c(2,2,2,2,2)]]
    hccDgram[[c(2,2,2,2,2)]] <- tmp
  }
  colInd <- order.dendrogram(hccDgram)
  # omit censored surveys for item clustering
  colHc <- coefHclust(xmatr$Likert)
  hclDgram <- as.dendrogram(colHc)
  lH <- coefH(na.omit(xmatr$Likert), se = F, results = F)
  hclDgram <- reorder(hclDgram, lH$Hi)
  rowInd <- order.dendrogram(hclDgram)
  par(oma=c(rep(0.4, 3), 8.4), mar=rep(0,4))
  layout(matrix(c(4, 5, 6, 3, 1, 2), ncol=2), heights=c(1,2, 0.35),
         widths = c(2,2))
  klrs <- gray(seq(from=0.5, to=4.5, by=1)/5)
  # 1. Survey response image
  image(as.matrix(xmat[colInd, rowInd]), col=klrs, axes=F)
  mtext("questions", side=4, line=3.5, outer=T, cex=cexAll*3,
        at=sum(par()$usr[1:2] %*% c(0.6, 0.4)))
  mtext("aggregate difference\n in visitor responses",
        side=4, line=5.75, outer=T, cex=cexAll*1.75,
        at=sum(par()$usr[1:2] %*% c(1/6, 5/6)))
  # 2. Answer legend panel
  plot(1, 1, type="n", axes=F, xlab="", ylab="")
  legend(x="top", legend=c(ansVals, "No Answer"), fill=c(klrs, "white"),
         cex=cexAll*2, ncol=2,
         x.intersp = 0.25, y.intersp = 0.8,
         title="visitor survey responses", bty="n")
  # 3. Dendrogram panel
  dendYmax <- 1.08*max(attr(hccDgram[[1]],"height"),
                       attr(hccDgram[[2]],"height"))
  pres <- plot(hccDgram, leaflab = "none", axes=F, xaxs = "i",
               ylim=c(-0.25, dendYmax), yaxs="i")
  dendUsr <- par()$usr
  for(lns in 1:3){
    abline(a=lns*numQ+0.1, b=0, col=gray(0.7), lty=2)
    mtext(text=lns*numQ, side = 4, line=0.5,
         at = lns*numQ, cex=cexAll, las = 1)
#    text(lab=lns*numQ, x=dendUsr[1:2] %*% c(0.92, 0.08),
#         y=lns*numQ+0.1, adj=c(1,-0.2), cex=cexAll + 0.25)
  }
  text(x=attr(hccDgram[[1]], "midpoint"),
       y=0.97*par()$usr[4], cex=cexAll + 0.5, adj=1.1,
       lab="RED")
  text(x=attr(hccDgram[[1]], "midpoint"),
       y=0.97*par()$usr[4], cex=cexAll + 0.5, adj=-0.2,
       lab="LIGHT")
  text(x=attr(hccDgram[[2]], "midpoint") + attr(hccDgram[[1]], "members"),
       y=0.97*par()$usr[4], cex=cexAll + 0.5, adj=1.1,
       lab="WHITE")
  text(x=attr(hccDgram[[2]], "midpoint") + attr(hccDgram[[1]], "members"),
       y=0.97*par()$usr[4], cex=cexAll + 0.5, adj=-0.2,
       lab="LIGHT")
  # 4. Question panel
  plot(1, 1, type="n", axes=F, xlab="", ylab="")
  text(lab=Q1, x=sum(par()$usr[1:2] %*% c(0.98, 0.02)), y=1.65*par()$usr[3],
       adj=0, cex=cexAll*2.5)
  # 5. Question parts and scales
  plot(hclDgram, leaflab = "none", axes=F, xaxs = "i", yaxs = "i", horiz = T,
       xlim=c(1.7, 0))
  print(imgUsr <- par()$usr)
  ySpace <- diff(imgUsr[3:4])/numQ
  vAdj <- -1/16
  for(qq in seq(along=Qparts)){
    qText <- Qparts[rowInd][qq]
    text(lab = qText,
          x=imgUsr[1]-0.05, y = (qq-vAdj)*ySpace, adj=c(0, 0.5), cex=cexAll*1.75)
  }
  Llvls <- c(0.8, 0.6, 0.4)
  for(lns in seq(along=Llvls)){
    lines(x=1 - rep(Llvls[lns], 2), y=imgUsr[3:4], col=gray(0.7), lty=2)
    text(lab=Llvls[lns], x=1 - Llvls[lns],
         y=imgUsr[3:4] %*% c(0.97, 0.03), adj=c(-0.2, 1), cex=cexAll + 0.5)
  }
  # 6. Loevinger's panel
  image(as.matrix(xmat[, -colDrop]), col="white", axes=F)
  imgUsr <- par()$usr
  text(x=imgUsr[1:2] %*% c(0.15, 0.85), y=imgUsr[3:4] %*% c(0.3, 0.7),
       lab="Loevinger's H", cex=cexAll*2)
  return(invisible(list(survOrder=colInd, varOrder=rowInd)))
}
complexPlot <- function(rmat, xCensor=NULL, scales=NULL, Quest, qVals, aVals){
  # expects global odir for write.csv()
  if(is.null(xCensor))
     redMat <- rmat[, -match(c("Color", "N10city"), names(rmat))] else
       redMat <- rmat[!xCensor, -match(c("Color", "N10city"), names(rmat))]
  aispHrange <- seq(from=0.05, to=0.9, by=0.05)
  # AISP results exclude the censored data and any NA (via revCodeH)
  redGps <- aisp(revCodeH(redMat)$Likert, lowerbound = aispHrange)
  redGps <- t(unique(t(redGps)))
  redGps <- redGps[,rev(dimnames(redGps)[[2]])]
  redGps <- redGps[,!apply(redGps==0, 2, all)]
  write.csv(redGps, file=paste0(odir,
                                sub(patt="_.*", repl="", names(redMat)[1]),
                                "Aisp.csv"))
  if(is.null(scales)) return(redGps)
  # This section creates a table that is no longer plotted
  LH <- vector(mode="list",length=length(scales))
  ############# psych::alpha seems to throw warnings for 2-item scales
  redMatr <- revCodeH(redMat)$Likert
  for (sc in seq(along=scales)){
    tmp <- coefH(redMatr[,scales[[sc]]], se = T, results = F)
    LH[[sc]] <- tmp
  } # get Loevinger's H and Coefficient Alpha
  ##### create LoevH rows for complex figure
  lHrows <- sapply(sapply(unlist(lapply(LH,"[[", "H")),
                          gsub, patt="[)( ]", repl=""), as.numeric)
  lHrows <- matrix(sprintf("%2.2f", lHrows), nrow=2)
  lHrows <- cbind(rep(NA, 2), lHrows)
  vPolars <- sapply(names(redMatr), function(s)
    substr(s, nchar(s), nchar(s)))=="-"
  vPolars <- sapply(vPolars, ifelse, "-", "+")
  vPolars[grep(patt="GYE", names(vPolars))] <- "?"
  scls <- matrix(data=0, nrow=dim(redMat)[2], ncol=length(scales),
                 dimnames=list(dimnames(redMatr)[[2]], NULL))
  for(sc in seq(along=scales))
    scls[scales[[sc]], sc] <- names(scales)[sc]
  scls <- cbind("Polarity" = vPolars, scls)
  scls <- rbind(scls, lHrows)
# end of unused table
  ####### Create Streetlight compound figure
  kolr <- gray(seq(from=0, to=5, length.out=5)/5)
#  CairoPNG(filename = paste0(odir,
#                             sub(patt="_.*", repl="", names(redMat)[1]),
#                             "Composite.png"),
#           width=2800,height=1750)
  hmRes <- clustGramKf(rmat, xCensor, Q1=Quest, Qparts = qVals, ansVals=aVals)
#  dev.off()
  return(scls)
}
