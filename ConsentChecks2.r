library(here)
library(memisc)
library(party)
fact2logical <- function(x){
	levels(x) <- c("F", "T")
	return(as.logical(x))
}
fl <- "MASTER_GRTE_NightSkies_Cleaned_NRB2.sav"
chk <- as.data.set(spss.system.file(
    here("Data",fl)))
chk2 <- as.data.frame(chk[grep("NRB", names(chk))])
chk2[2:5] <- lapply(chk2[2:5], fact2logical)
chk2 <- chk2[chk2$NRB_binary %in% c("REFUSAL", "AGREE"),]
chk2$NRB_binary <- droplevels(chk2$NRB_binary)
cfctrl <- cforest_control(ntree = 1000, replace = F, fraction = 0.632, mtry=2)
nrbRf <- cforest(NRB_binary ~ ., data=chk2,
	controls = cfctrl)
sum(diag(table(chk2$NRB_binary,predict(nrbRf, OOB=T))))/dim(chk2)[1]
varimp(nrbRf)
chk3 <- split(chk2, chk2$NRB_binary)
ynam <- names(chk2)[-1]
for(nm in ynam){
	tabl <- table(chk2[["NRB_binary"]], chk2[[nm]])
	print(nm)
	print(fisher.test(tabl))
}
