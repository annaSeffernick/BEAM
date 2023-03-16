###################################
# Create Example Data for
# BEAM R Package
# AES
# 2023-03-16
######################################


#rm(list=ls())

ddir <- "C:/Users/aseffern/Box/BookChapterApril2022"
source(paste0(ddir, "/Xueyuan-BEAM/prep-beam-data.R"))
library(survival)

lesion<-readRDS(paste0(ddir, "/Abdel-GRIN/TALL-binary-lsn-mtx-atleast5patients.RDS"))
#lesion <- cbind(lesion.ID=rownames(lesion), lesion)
lesion.sm <- lesion[1:20,]

rdat<-"C:/Users/aseffern/Box/BookChapterApril2022/DataSet/v1/TALL-dataset.Rdata"
load(rdat)
geneann<-read.csv(paste0(ddir, "/DataSet/v1/ensembl.annotation.csv"))
RNA<-RNA[!is.na(RNA$ensembl.ID), ]
rownames(RNA)<-RNA$ensembl.ID
RNA.sm <- RNA[1:20,]
geneann.sm <- geneann[which(geneann$ensembl.ID %in% rownames(RNA.sm)),]

#RNA<-RNA[rownames(RNA)%in%substring(rownames(lesion), 1, 15), ]
#lesion<-lesion[substring(rownames(lesion), 1, 15)%in%rownames(RNA), ]

omicdat<-list(Lesion=as.matrix(lesion.sm),
              RNA=as.matrix(RNA.sm[, -1]))
omicann<-list(Lesion=cbind.data.frame(prob.id=rownames(lesion.sm),
                                      gene.id=substring(rownames(lesion.sm), 1, 15)),
              RNA=geneann.sm)


clin$EFScensor[clin$First_Event%in%c("None", "Censored")]<-0
clin$EFScensor[clin$First_Event%in%c("Death","Progression", "Relapse" , "Second Malignant Neoplasm")]<-1

clin$OScensor[clin$Vital_Status%in%"Alive"]<-0
clin$OScensor[clin$Vital_Status%in%"Dead"]<-1

rownames(clin)<-clin$ID
clin$"RNA.clm" <- clin$ID
clin$"Lesion.clm" <- clin$ID
clin$RNA.id<-clin$Lesion.id<-clin$ID
clin$EFS<-Surv(clin$Event_Days, clin$EFScensor)
clin$OS<-Surv(clin$OS_Days, clin$OScensor)

clinf<-clin[, c(1, 8, 15:ncol(clin))]

omicann<-list(Lesion=cbind.data.frame(id=rownames(lesion.sm),
                                      gene=substring(rownames(lesion.sm), 1, 15)),
              RNA=cbind.data.frame(id=rownames(RNA.sm),
                                   gene=rownames(RNA.sm)))
setdat<-cbind.data.frame(set.id=c(rownames(RNA.sm), substring(rownames(lesion.sm), 1, 15)),
                         mtx.id=c(rep("RNA", nrow(RNA.sm)), rep("Lesion", nrow(lesion.sm))),
                         row.id=c(rownames(RNA.sm), rownames(lesion.sm)))


usethis::use_data(omicdat, overwrite = TRUE)
usethis::use_data(clinf, overwrite = TRUE)
usethis::use_data(omicann, overwrite = TRUE)
usethis::use_data(setdat, overwrite = TRUE)
