#load libraries
library("seqinr")
library("data.table")
library("stringr")
library("dplyr")

#set options
options(stringsAsFactors = FALSE)

#set wd
setwd("~/Desktop/Project_Suillus_comp_genomics/R")

#read in the input files
Suivar1_TMHMM<- data.frame(read.csv("Suivar1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Suitom1_TMHMM<- data.frame(read.csv("Suitom1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Suisub1_TMHMM<- data.frame(read.csv("Suisub1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Suisu1_TMHMM<- data.frame(read.csv("Suisu1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Suipla1_TMHMM<- data.frame(read.csv("Suipla1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Suipic1_TMHMM<- data.frame(read.csv("Suipic1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Suipal1_TMHMM<- data.frame(read.csv("Suipal1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Suiocc1_TMHMM<- data.frame(read.csv("Suiocc1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Suilu4_TMHMM<- data.frame(read.csv("Suilu4_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Suilak1_TMHMM<- data.frame(read.csv("Suilak1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Suihi1_TMHMM<- data.frame(read.csv("Suihi1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Suigr1_TMHMM<- data.frame(read.csv("Suigr1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Suidec1_TMHMM<- data.frame(read.csv("Suidec1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Suicot1_TMHMM<- data.frame(read.csv("Suicot1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Suicli1_TMHMM<- data.frame(read.csv("Suicli1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Suibr2_TMHMM<- data.frame(read.csv("Suibr2_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Suibov1_TMHMM<- data.frame(read.csv("Suibov1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Suiamp1_TMHMM<- data.frame(read.csv("Suiamp1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Suiame1_TMHMM<- data.frame(read.csv("Suiame1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Suifus1_TMHMM<- data.frame(read.csv("Suifus1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))

#non-Suillus set
Rhivul1_TMHMM<- data.frame(read.csv("Rhivul1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Rhitru1_TMHMM<- data.frame(read.csv("Rhitru1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Amamu1_TMHMM<- data.frame(read.csv("Amamu1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Hebcy2_TMHMM<- data.frame(read.csv("Hebcy2_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Lacbi2_TMHMM<- data.frame(read.csv("Lacbi2_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Paxin1_TMHMM<- data.frame(read.csv("Paxin1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Pilcr1_TMHMM<- data.frame(read.csv("Pilcr1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Pismi1_TMHMM<- data.frame(read.csv("Pismi1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Sclci1_TMHMM<- data.frame(read.csv("Sclci1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Theter1_TMHMM<- data.frame(read.csv("Theter1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Thega1_TMHMM<- data.frame(read.csv("Thega1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Rhivi1_TMHMM<- data.frame(read.csv("Rhivi1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Rhives1_TMHMM<- data.frame(read.csv("Rhives1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Rhisa1_TMHMM<- data.frame(read.csv("Rhisa1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Ruscom1_TMHMM<- data.frame(read.csv("Ruscom1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Rusbre1_TMHMM<- data.frame(read.csv("Rusbre1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Pisti1_TMHMM<- data.frame(read.csv("Pisti1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Lacam2_TMHMM<- data.frame(read.csv("Lacam2_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Hydru2_TMHMM<- data.frame(read.csv("Hydru2_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Gaumor1_TMHMM<- data.frame(read.csv("Gaumor1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Gyrli1_TMHMM<- data.frame(read.csv("Gyrli1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Cananz1_TMHMM<- data.frame(read.csv("Cananz1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Sclci1_TMHMM<- data.frame(read.csv("Sclci1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Hyssto1_TMHMM<- data.frame(read.csv("Hyssto1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))


#isolate relevent lines 
#get only lines starting with a hash
Suivar1_TMHMM_sm<- data.frame(Suivar1_TMHMM[ grep("# jgi", Suivar1_TMHMM$header),])
Suitom1_TMHMM_sm<- data.frame(Suitom1_TMHMM[ grep("# jgi", Suitom1_TMHMM$header),])
Suisub1_TMHMM_sm<- data.frame(Suisub1_TMHMM[ grep("# jgi", Suisub1_TMHMM$header),])
Suisu1_TMHMM_sm<- data.frame(Suisu1_TMHMM[ grep("# jgi", Suisu1_TMHMM$header),])
Suipla1_TMHMM_sm<- data.frame(Suipla1_TMHMM[ grep("# jgi", Suipla1_TMHMM$header),])
Suipic1_TMHMM_sm<- data.frame(Suipic1_TMHMM[ grep("# jgi", Suipic1_TMHMM$header),])
Suipal1_TMHMM_sm<- data.frame(Suipal1_TMHMM[ grep("# jgi", Suipal1_TMHMM$header),])
Suiocc1_TMHMM_sm<- data.frame(Suiocc1_TMHMM[ grep("# jgi", Suiocc1_TMHMM$header),])
Suilu4_TMHMM_sm<- data.frame(Suilu4_TMHMM[ grep("# jgi", Suilu4_TMHMM$header),])
Suilak1_TMHMM_sm<- data.frame(Suilak1_TMHMM[ grep("# jgi", Suilak1_TMHMM$header),])
Suihi1_TMHMM_sm<- data.frame(Suihi1_TMHMM[ grep("# jgi", Suihi1_TMHMM$header),])
Suigr1_TMHMM_sm<- data.frame(Suigr1_TMHMM[ grep("# jgi", Suigr1_TMHMM$header),])
Suidec1_TMHMM_sm<- data.frame(Suidec1_TMHMM[ grep("# jgi", Suidec1_TMHMM$header),])
Suicot1_TMHMM_sm<- data.frame(Suicot1_TMHMM[ grep("# jgi", Suicot1_TMHMM$header),])
Suicli1_TMHMM_sm<- data.frame(Suicli1_TMHMM[ grep("# jgi", Suicli1_TMHMM$header),])
Suibr2_TMHMM_sm<- data.frame(Suibr2_TMHMM[ grep("# jgi", Suibr2_TMHMM$header),])
Suibov1_TMHMM_sm<- data.frame(Suibov1_TMHMM[ grep("# jgi", Suibov1_TMHMM$header),])
Suiamp1_TMHMM_sm<- data.frame(Suiamp1_TMHMM[ grep("# jgi", Suiamp1_TMHMM$header),])
Suiame1_TMHMM_sm<- data.frame(Suiame1_TMHMM[ grep("# jgi", Suiame1_TMHMM$header),])
Suifus1_TMHMM_sm<- data.frame(Suifus1_TMHMM[ grep("# jgi", Suifus1_TMHMM$header),])
#non-Suillus set
Rhivul1_TMHMM_sm<- data.frame(Rhivul1_TMHMM[ grep("# jgi", Rhivul1_TMHMM$header),])
Rhitru1_TMHMM_sm<- data.frame(Rhitru1_TMHMM[ grep("# jgi", Rhitru1_TMHMM$header),])
Amamu1_TMHMM_sm<- data.frame(Amamu1_TMHMM[ grep("# jgi", Amamu1_TMHMM$header),])
Hebcy2_TMHMM_sm<- data.frame(Hebcy2_TMHMM[ grep("# jgi", Hebcy2_TMHMM$header),])
Lacbi2_TMHMM_sm<- data.frame(Lacbi2_TMHMM[ grep("# jgi", Lacbi2_TMHMM$header),])
Paxin1_TMHMM_sm<- data.frame(Paxin1_TMHMM[ grep("# jgi", Paxin1_TMHMM$header),])
Pilcr1_TMHMM_sm<- data.frame(Pilcr1_TMHMM[ grep("# jgi", Pilcr1_TMHMM$header),])
Pismi1_TMHMM_sm<- data.frame(Pismi1_TMHMM[ grep("# jgi", Pismi1_TMHMM$header),])
Sclci1_TMHMM_sm<- data.frame(Sclci1_TMHMM[ grep("# jgi", Sclci1_TMHMM$header),])
Theter1_TMHMM_sm<- data.frame(Theter1_TMHMM[ grep("# jgi", Theter1_TMHMM$header),])
Thega1_TMHMM_sm<- data.frame(Thega1_TMHMM[ grep("# jgi", Thega1_TMHMM$header),])
Rhivi1_TMHMM_sm<- data.frame(Rhivi1_TMHMM[ grep("# jgi", Rhivi1_TMHMM$header),])
Rhives1_TMHMM_sm<- data.frame(Rhives1_TMHMM[ grep("# jgi", Rhives1_TMHMM$header),])
Rhisa1_TMHMM_sm<- data.frame(Rhisa1_TMHMM[ grep("# jgi", Rhisa1_TMHMM$header),])
Ruscom1_TMHMM_sm<- data.frame(Ruscom1_TMHMM[ grep("# jgi", Ruscom1_TMHMM$header),])
Rusbre1_TMHMM_sm<- data.frame(Rusbre1_TMHMM[ grep("# jgi", Rusbre1_TMHMM$header),])
Pisti1_TMHMM_sm<- data.frame(Pisti1_TMHMM[ grep("# jgi", Pisti1_TMHMM$header),])
Lacam2_TMHMM_sm<- data.frame(Lacam2_TMHMM[ grep("# jgi", Lacam2_TMHMM$header),])
Hydru2_TMHMM_sm<- data.frame(Hydru2_TMHMM[ grep("# jgi", Hydru2_TMHMM$header),])
Gaumor1_TMHMM_sm<- data.frame(Gaumor1_TMHMM[ grep("# jgi", Gaumor1_TMHMM$header),])
Gyrli1_TMHMM_sm<- data.frame(Gyrli1_TMHMM[ grep("# jgi", Gyrli1_TMHMM$header),])
Cananz1_TMHMM_sm<- data.frame(Cananz1_TMHMM[ grep("# jgi", Cananz1_TMHMM$header),])
Hyssto1_TMHMM_sm<- data.frame(Hyssto1_TMHMM[ grep("# jgi", Hyssto1_TMHMM$header),])


#combine the size and predicted TMD# by columns 
Suivar1_TMHMM_sm.by.col<- data.frame(cbind(Suivar1_TMHMM_sm[c(TRUE, FALSE), ],
                                           Suivar1_TMHMM_sm[c(FALSE, TRUE), ]))
Suitom1_TMHMM_sm.by.col<- data.frame(cbind(Suitom1_TMHMM_sm[c(TRUE, FALSE), ],
                                           Suitom1_TMHMM_sm[c(FALSE, TRUE), ]))
Suisub1_TMHMM_sm.by.col<- data.frame(cbind(Suisub1_TMHMM_sm[c(TRUE, FALSE), ],
                                           Suisub1_TMHMM_sm[c(FALSE, TRUE), ]))
Suisu1_TMHMM_sm.by.col<- data.frame(cbind(Suisu1_TMHMM_sm[c(TRUE, FALSE), ],
                                          Suisu1_TMHMM_sm[c(FALSE, TRUE), ]))
Suipla1_TMHMM_sm.by.col<- data.frame(cbind(Suipla1_TMHMM_sm[c(TRUE, FALSE), ],
                                           Suipla1_TMHMM_sm[c(FALSE, TRUE), ]))
Suipic1_TMHMM_sm.by.col<- data.frame(cbind(Suipic1_TMHMM_sm[c(TRUE, FALSE), ],
                                           Suipic1_TMHMM_sm[c(FALSE, TRUE), ]))
Suipal1_TMHMM_sm.by.col<- data.frame(cbind(Suipal1_TMHMM_sm[c(TRUE, FALSE), ],
                                           Suipal1_TMHMM_sm[c(FALSE, TRUE), ]))
Suiocc1_TMHMM_sm.by.col<- data.frame(cbind(Suiocc1_TMHMM_sm[c(TRUE, FALSE), ],
                                           Suiocc1_TMHMM_sm[c(FALSE, TRUE), ]))
Suilu4_TMHMM_sm.by.col<- data.frame(cbind(Suilu4_TMHMM_sm[c(TRUE, FALSE), ],
                                          Suilu4_TMHMM_sm[c(FALSE, TRUE), ]))
Suilak1_TMHMM_sm.by.col<- data.frame(cbind(Suilak1_TMHMM_sm[c(TRUE, FALSE), ],
                                           Suilak1_TMHMM_sm[c(FALSE, TRUE), ]))
Suihi1_TMHMM_sm.by.col<- data.frame(cbind(Suihi1_TMHMM_sm[c(TRUE, FALSE), ],
                                          Suihi1_TMHMM_sm[c(FALSE, TRUE), ]))
Suigr1_TMHMM_sm.by.col<- data.frame(cbind(Suigr1_TMHMM_sm[c(TRUE, FALSE), ],
                                          Suigr1_TMHMM_sm[c(FALSE, TRUE), ]))
Suidec1_TMHMM_sm.by.col<- data.frame(cbind(Suidec1_TMHMM_sm[c(TRUE, FALSE), ],
                                           Suidec1_TMHMM_sm[c(FALSE, TRUE), ]))
Suicot1_TMHMM_sm.by.col<- data.frame(cbind(Suicot1_TMHMM_sm[c(TRUE, FALSE), ],
                                           Suicot1_TMHMM_sm[c(FALSE, TRUE), ]))
Suicli1_TMHMM_sm.by.col<- data.frame(cbind(Suicli1_TMHMM_sm[c(TRUE, FALSE), ],
                                           Suicli1_TMHMM_sm[c(FALSE, TRUE), ]))
Suibr2_TMHMM_sm.by.col<- data.frame(cbind(Suibr2_TMHMM_sm[c(TRUE, FALSE), ],
                                          Suibr2_TMHMM_sm[c(FALSE, TRUE), ]))
Suibov1_TMHMM_sm.by.col<- data.frame(cbind(Suibov1_TMHMM_sm[c(TRUE, FALSE), ],
                                           Suibov1_TMHMM_sm[c(FALSE, TRUE), ]))
Suiamp1_TMHMM_sm.by.col<- data.frame(cbind(Suiamp1_TMHMM_sm[c(TRUE, FALSE), ],
                                           Suiamp1_TMHMM_sm[c(FALSE, TRUE), ]))
Suiame1_TMHMM_sm.by.col<- data.frame(cbind(Suiame1_TMHMM_sm[c(TRUE, FALSE), ],
                                           Suiame1_TMHMM_sm[c(FALSE, TRUE), ]))
Suifus1_TMHMM_sm.by.col<- data.frame(cbind(Suifus1_TMHMM_sm[c(TRUE, FALSE), ],
                                           Suifus1_TMHMM_sm[c(FALSE, TRUE), ]))


#non-Suillus set
Rhivul1_TMHMM_sm.by.col<- data.frame(cbind(Rhivul1_TMHMM_sm[c(TRUE, FALSE), ],
                                           Rhivul1_TMHMM_sm[c(FALSE, TRUE), ]))
Rhitru1_TMHMM_sm.by.col<- data.frame(cbind(Rhitru1_TMHMM_sm[c(TRUE, FALSE), ],
                                           Rhitru1_TMHMM_sm[c(FALSE, TRUE), ]))
Amamu1_TMHMM_sm.by.col<- data.frame(cbind(Amamu1_TMHMM_sm[c(TRUE, FALSE), ],
                                          Amamu1_TMHMM_sm[c(FALSE, TRUE), ]))
Hebcy2_TMHMM_sm.by.col<- data.frame(cbind(Hebcy2_TMHMM_sm[c(TRUE, FALSE), ],
                                          Hebcy2_TMHMM_sm[c(FALSE, TRUE), ]))
Lacbi2_TMHMM_sm.by.col<- data.frame(cbind(Lacbi2_TMHMM_sm[c(TRUE, FALSE), ],
                                          Lacbi2_TMHMM_sm[c(FALSE, TRUE), ]))
Paxin1_TMHMM_sm.by.col<- data.frame(cbind(Paxin1_TMHMM_sm[c(TRUE, FALSE), ],
                                          Paxin1_TMHMM_sm[c(FALSE, TRUE), ]))
Pilcr1_TMHMM_sm.by.col<- data.frame(cbind(Pilcr1_TMHMM_sm[c(TRUE, FALSE), ],
                                          Pilcr1_TMHMM_sm[c(FALSE, TRUE), ]))
Pismi1_TMHMM_sm.by.col<- data.frame(cbind(Pismi1_TMHMM_sm[c(TRUE, FALSE), ],
                                          Pismi1_TMHMM_sm[c(FALSE, TRUE), ]))
Sclci1_TMHMM_sm.by.col<- data.frame(cbind(Sclci1_TMHMM_sm[c(TRUE, FALSE), ],
                                          Sclci1_TMHMM_sm[c(FALSE, TRUE), ]))
Theter1_TMHMM_sm.by.col<- data.frame(cbind(Theter1_TMHMM_sm[c(TRUE, FALSE), ],
                                           Theter1_TMHMM_sm[c(FALSE, TRUE), ]))
Thega1_TMHMM_sm.by.col<- data.frame(cbind(Thega1_TMHMM_sm[c(TRUE, FALSE), ],
                                          Thega1_TMHMM_sm[c(FALSE, TRUE), ]))
Rhivi1_TMHMM_sm.by.col<- data.frame(cbind(Rhivi1_TMHMM_sm[c(TRUE, FALSE), ],
                                           Rhivi1_TMHMM_sm[c(FALSE, TRUE), ]))
Rhives1_TMHMM_sm.by.col<- data.frame(cbind(Rhives1_TMHMM_sm[c(TRUE, FALSE), ],
                                           Rhives1_TMHMM_sm[c(FALSE, TRUE), ]))
Rhisa1_TMHMM_sm.by.col<- data.frame(cbind(Rhisa1_TMHMM_sm[c(TRUE, FALSE), ],
                                          Rhisa1_TMHMM_sm[c(FALSE, TRUE), ]))
Ruscom1_TMHMM_sm.by.col<- data.frame(cbind(Ruscom1_TMHMM_sm[c(TRUE, FALSE), ],
                                           Ruscom1_TMHMM_sm[c(FALSE, TRUE), ]))
Rusbre1_TMHMM_sm.by.col<- data.frame(cbind(Rusbre1_TMHMM_sm[c(TRUE, FALSE), ],
                                           Rusbre1_TMHMM_sm[c(FALSE, TRUE), ]))
Pisti1_TMHMM_sm.by.col<- data.frame(cbind(Pisti1_TMHMM_sm[c(TRUE, FALSE), ],
                                          Pisti1_TMHMM_sm[c(FALSE, TRUE), ]))
Lacam2_TMHMM_sm.by.col<- data.frame(cbind(Lacam2_TMHMM_sm[c(TRUE, FALSE), ],
                                          Lacam2_TMHMM_sm[c(FALSE, TRUE), ]))
Hydru2_TMHMM_sm.by.col<- data.frame(cbind(Hydru2_TMHMM_sm[c(TRUE, FALSE), ],
                                          Hydru2_TMHMM_sm[c(FALSE, TRUE), ]))
Gaumor1_TMHMM_sm.by.col<- data.frame(cbind(Gaumor1_TMHMM_sm[c(TRUE, FALSE), ],
                                           Gaumor1_TMHMM_sm[c(FALSE, TRUE), ]))
Gyrli1_TMHMM_sm.by.col<- data.frame(cbind(Gyrli1_TMHMM_sm[c(TRUE, FALSE), ],
                                          Gyrli1_TMHMM_sm[c(FALSE, TRUE), ]))
Cananz1_TMHMM_sm.by.col<- data.frame(cbind(Cananz1_TMHMM_sm[c(TRUE, FALSE), ],
                                           Cananz1_TMHMM_sm[c(FALSE, TRUE), ]))
Hyssto1_TMHMM_sm.by.col<- data.frame(cbind(Hyssto1_TMHMM_sm[c(TRUE, FALSE), ],
                                           Hyssto1_TMHMM_sm[c(FALSE, TRUE), ]))


#only proteins with no TMD
Suivar1_TMHMM_no_TMD<- data.frame(Suivar1_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Suivar1_TMHMM_sm.by.col$X2), ])
Suitom1_TMHMM_no_TMD<- data.frame(Suitom1_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Suitom1_TMHMM_sm.by.col$X2),])
Suisub1_TMHMM_no_TMD<- data.frame(Suisub1_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Suisub1_TMHMM_sm.by.col$X2),])
Suisu1_TMHMM_no_TMD<- data.frame(Suisu1_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Suisu1_TMHMM_sm.by.col$X2),])
Suipla1_TMHMM_no_TMD<- data.frame(Suipla1_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Suipla1_TMHMM_sm.by.col$X2),])
Suipic1_TMHMM_no_TMD<- data.frame(Suipic1_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Suipic1_TMHMM_sm.by.col$X2),])
Suipal1_TMHMM_no_TMD<- data.frame(Suipal1_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Suipal1_TMHMM_sm.by.col$X2),])
Suiocc1_TMHMM_no_TMD<- data.frame(Suiocc1_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Suiocc1_TMHMM_sm.by.col$X2),])
Suilu4_TMHMM_no_TMD<- data.frame(Suilu4_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Suilu4_TMHMM_sm.by.col$X2),])
Suilak1_TMHMM_no_TMD<- data.frame(Suilak1_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Suilak1_TMHMM_sm.by.col$X2),])
Suihi1_TMHMM_no_TMD<- data.frame(Suihi1_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Suihi1_TMHMM_sm.by.col$X2),])
Suigr1_TMHMM_no_TMD<- data.frame(Suigr1_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Suigr1_TMHMM_sm.by.col$X2),])
Suidec1_TMHMM_no_TMD<- data.frame(Suidec1_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Suidec1_TMHMM_sm.by.col$X2),])
Suicot1_TMHMM_no_TMD<- data.frame(Suicot1_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Suicot1_TMHMM_sm.by.col$X2),])
Suicli1_TMHMM_no_TMD<- data.frame(Suicli1_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Suicli1_TMHMM_sm.by.col$X2),])
Suibr2_TMHMM_no_TMD<- data.frame(Suibr2_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Suibr2_TMHMM_sm.by.col$X2),])
Suibov1_TMHMM_no_TMD<- data.frame(Suibov1_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Suibov1_TMHMM_sm.by.col$X2),])
Suiamp1_TMHMM_no_TMD<- data.frame(Suiamp1_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Suiamp1_TMHMM_sm.by.col$X2),])
Suiame1_TMHMM_no_TMD<- data.frame(Suiame1_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Suiame1_TMHMM_sm.by.col$X2),])
Suifus1_TMHMM_no_TMD<- data.frame(Suifus1_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Suifus1_TMHMM_sm.by.col$X2),])

#non-Suillus set
Rhivul1_TMHMM_no_TMD<- data.frame(Rhivul1_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Rhivul1_TMHMM_sm.by.col$X2),])
Rhitru1_TMHMM_no_TMD<- data.frame(Rhitru1_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Rhitru1_TMHMM_sm.by.col$X2),])
Amamu1_TMHMM_no_TMD<- data.frame(Amamu1_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Amamu1_TMHMM_sm.by.col$X2),])
Hebcy2_TMHMM_no_TMD<- data.frame(Hebcy2_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Hebcy2_TMHMM_sm.by.col$X2),])
Lacbi2_TMHMM_no_TMD<- data.frame(Lacbi2_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Lacbi2_TMHMM_sm.by.col$X2),])
Paxin1_TMHMM_no_TMD<- data.frame(Paxin1_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Paxin1_TMHMM_sm.by.col$X2),])
Pilcr1_TMHMM_no_TMD<- data.frame(Pilcr1_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Pilcr1_TMHMM_sm.by.col$X2),])
Pismi1_TMHMM_no_TMD<- data.frame(Pismi1_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Pismi1_TMHMM_sm.by.col$X2),])
Sclci1_TMHMM_no_TMD<- data.frame(Sclci1_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Sclci1_TMHMM_sm.by.col$X2),])
Theter1_TMHMM_no_TMD<- data.frame(Theter1_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Theter1_TMHMM_sm.by.col$X2),])
Thega1_TMHMM_no_TMD<- data.frame(Thega1_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Thega1_TMHMM_sm.by.col$X2),])
Rhivi1_TMHMM_no_TMD<- data.frame(Rhivi1_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Rhivi1_TMHMM_sm.by.col$X2),])
Rhives1_TMHMM_no_TMD<- data.frame(Rhives1_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Rhives1_TMHMM_sm.by.col$X2),])
Rhisa1_TMHMM_no_TMD<- data.frame(Rhisa1_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Rhisa1_TMHMM_sm.by.col$X2),])
Ruscom1_TMHMM_no_TMD<- data.frame(Ruscom1_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Ruscom1_TMHMM_sm.by.col$X2),])
Rusbre1_TMHMM_no_TMD<- data.frame(Rusbre1_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Rusbre1_TMHMM_sm.by.col$X2),])
Pisti1_TMHMM_no_TMD<- data.frame(Pisti1_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Pisti1_TMHMM_sm.by.col$X2),])
Lacam2_TMHMM_no_TMD<- data.frame(Lacam2_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Lacam2_TMHMM_sm.by.col$X2),])
Hydru2_TMHMM_no_TMD<- data.frame(Hydru2_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Hydru2_TMHMM_sm.by.col$X2),])
Gaumor1_TMHMM_no_TMD<- data.frame(Gaumor1_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Gaumor1_TMHMM_sm.by.col$X2),])
Gyrli1_TMHMM_no_TMD<- data.frame(Gyrli1_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Gyrli1_TMHMM_sm.by.col$X2),])
Cananz1_TMHMM_no_TMD<- data.frame(Cananz1_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Cananz1_TMHMM_sm.by.col$X2),])
Hyssto1_TMHMM_no_TMD<- data.frame(Hyssto1_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Hyssto1_TMHMM_sm.by.col$X2),])


#function to isolate the aa length
isolate_aa <- function(in_df) { 
  a<- data.frame(str_split(in_df$X1, "Length:"))
  b<- a[2,]
  c<- data.table::transpose(b)
  out_df<- cbind(in_df, c)
  return= out_df
}



#run the function
Suivar1_df<- isolate_aa(Suivar1_TMHMM_no_TMD)
Suitom1_df<- isolate_aa(Suitom1_TMHMM_no_TMD)
Suisub1_df<- isolate_aa(Suisub1_TMHMM_no_TMD)
Suisu1_df<- isolate_aa(Suisu1_TMHMM_no_TMD)
Suipla1_df<- isolate_aa(Suipla1_TMHMM_no_TMD)
Suipic1_df<- isolate_aa(Suipic1_TMHMM_no_TMD)
Suipal1_df<- isolate_aa(Suipal1_TMHMM_no_TMD)
Suiocc1_df<- isolate_aa(Suiocc1_TMHMM_no_TMD)
Suilu4_df<- isolate_aa(Suilu4_TMHMM_no_TMD)
Suilak1_df<- isolate_aa(Suilak1_TMHMM_no_TMD)
Suihi1_df<- isolate_aa(Suihi1_TMHMM_no_TMD)
Suigr1_df<- isolate_aa(Suigr1_TMHMM_no_TMD)
Suidec1_df<- isolate_aa(Suidec1_TMHMM_no_TMD)
Suicot1_df<- isolate_aa(Suicot1_TMHMM_no_TMD)
Suicli1_df<- isolate_aa(Suicli1_TMHMM_no_TMD)
Suibr2_df<- isolate_aa(Suibr2_TMHMM_no_TMD)
Suibov1_df<- isolate_aa(Suibov1_TMHMM_no_TMD)
Suiamp1_df<- isolate_aa(Suiamp1_TMHMM_no_TMD)
Suiame1_df<- isolate_aa(Suiame1_TMHMM_no_TMD)
Suifus1_df<- isolate_aa(Suifus1_TMHMM_no_TMD)

#non-Suillus set
Rhivul1_df<- isolate_aa(Rhivul1_TMHMM_no_TMD)
Rhitru1_df<- isolate_aa(Rhitru1_TMHMM_no_TMD)
Amamu1_df<- isolate_aa(Amamu1_TMHMM_no_TMD)
Hebcy2_df<- isolate_aa(Hebcy2_TMHMM_no_TMD)
Lacbi2_df<- isolate_aa(Lacbi2_TMHMM_no_TMD)
Paxin1_df<- isolate_aa(Paxin1_TMHMM_no_TMD)
Pilcr1_df<- isolate_aa(Pilcr1_TMHMM_no_TMD)
Pismi1_df<- isolate_aa(Pismi1_TMHMM_no_TMD)
Sclci1_df<- isolate_aa(Sclci1_TMHMM_no_TMD)
Theter1_df<- isolate_aa(Theter1_TMHMM_no_TMD)
Thega1_df<- isolate_aa(Thega1_TMHMM_no_TMD)
Rhivi1_df<- isolate_aa(Rhivi1_TMHMM_no_TMD)
Rhives1_df<- isolate_aa(Rhives1_TMHMM_no_TMD)
Rhisa1_df<- isolate_aa(Rhisa1_TMHMM_no_TMD)
Ruscom1_df<- isolate_aa(Ruscom1_TMHMM_no_TMD)
Rusbre1_df<- isolate_aa(Rusbre1_TMHMM_no_TMD)
Pisti1_df<- isolate_aa(Pisti1_TMHMM_no_TMD)
Lacam2_df<- isolate_aa(Lacam2_TMHMM_no_TMD)
Hydru2_df<- isolate_aa(Hydru2_TMHMM_no_TMD)
Gaumor1_df<- isolate_aa(Gaumor1_TMHMM_no_TMD)
Gyrli1_df<- isolate_aa(Gyrli1_TMHMM_no_TMD)
Cananz1_df<- isolate_aa(Cananz1_TMHMM_no_TMD)
Hyssto1_df<- isolate_aa(Hyssto1_TMHMM_no_TMD)


#gotta convert the aa number into a a numeric type or R thinks it's a char. string.
Suivar1_df.num<-transform(Suivar1_df, V1 = as.numeric(V1))
Suitom1_df.num<-transform(Suitom1_df, V1 = as.numeric(V1))
Suisub1_df.num<-transform(Suisub1_df, V1 = as.numeric(V1))
Suisu1_df.num<-transform(Suisu1_df, V1 = as.numeric(V1))
Suipla1_df.num<-transform(Suipla1_df, V1 = as.numeric(V1))
Suipic1_df.num<-transform(Suipic1_df, V1 = as.numeric(V1))
Suipal1_df.num<-transform(Suipal1_df, V1 = as.numeric(V1))
Suiocc1_df.num<-transform(Suiocc1_df, V1 = as.numeric(V1))
Suilu4_df.num<-transform(Suilu4_df, V1 = as.numeric(V1))
Suilak1_df.num<-transform(Suilak1_df, V1 = as.numeric(V1))
Suihi1_df.num<-transform(Suihi1_df, V1 = as.numeric(V1))
Suigr1_df.num<-transform(Suigr1_df, V1 = as.numeric(V1))
Suidec1_df.num<-transform(Suidec1_df, V1 = as.numeric(V1))
Suicot1_df.num<-transform(Suicot1_df, V1 = as.numeric(V1))
Suicli1_df.num<-transform(Suicli1_df, V1 = as.numeric(V1))
Suibr2_df.num<-transform(Suibr2_df, V1 = as.numeric(V1))
Suibov1_df.num<-transform(Suibov1_df, V1 = as.numeric(V1))
Suiamp1_df.num<-transform(Suiamp1_df, V1 = as.numeric(V1))
Suiame1_df.num<-transform(Suiame1_df, V1 = as.numeric(V1))
Suifus1_df.num<-transform(Suifus1_df, V1 = as.numeric(V1))
#non-Suillus set
Rhivul1_df.num<-transform(Rhivul1_df, V1 = as.numeric(V1))
Rhitru1_df.num<-transform(Rhitru1_df, V1 = as.numeric(V1))
Amamu1_df.num<-transform(Amamu1_df, V1 = as.numeric(V1))
Hebcy2_df.num<-transform(Hebcy2_df, V1 = as.numeric(V1))
Lacbi2_df.num<-transform(Lacbi2_df, V1 = as.numeric(V1))
Paxin1_df.num<-transform(Paxin1_df, V1 = as.numeric(V1))
Pilcr1_df.num<-transform(Pilcr1_df, V1 = as.numeric(V1))
Pismi1_df.num<-transform(Pismi1_df, V1 = as.numeric(V1))
Sclci1_df.num<-transform(Sclci1_df, V1 = as.numeric(V1))
Theter1_df.num<-transform(Theter1_df, V1 = as.numeric(V1))
Thega1_df.num<-transform(Thega1_df, V1 = as.numeric(V1))
Rhivi1_df.num<-transform(Rhivi1_df, V1 = as.numeric(V1))
Rhives1_df.num<-transform(Rhives1_df, V1 = as.numeric(V1))
Rhisa1_df.num<-transform(Rhisa1_df, V1 = as.numeric(V1))
Ruscom1_df.num<-transform(Ruscom1_df, V1 = as.numeric(V1))
Rusbre1_df.num<-transform(Rusbre1_df, V1 = as.numeric(V1))
Pisti1_df.num<-transform(Pisti1_df, V1 = as.numeric(V1))
Lacam2_df.num<-transform(Lacam2_df, V1 = as.numeric(V1))
Hydru2_df.num<-transform(Hydru2_df, V1 = as.numeric(V1))
Gaumor1_df.num<-transform(Gaumor1_df, V1 = as.numeric(V1))
Gyrli1_df.num<-transform(Gyrli1_df, V1 = as.numeric(V1))
Cananz1_df.num<-transform(Cananz1_df, V1 = as.numeric(V1))
Hyssto1_df.num<-transform(Hyssto1_df, V1 = as.numeric(V1))

#only proteins < 300 aa in length
Suivar1_TMHMM_no_TMD_300<- Suivar1_df.num[Suivar1_df.num$V1 < 300, ]
Suitom1_TMHMM_no_TMD_300<- Suitom1_df.num[Suitom1_df.num$V1 < 300, ]
Suisub1_TMHMM_no_TMD_300<- Suisub1_df.num[Suisub1_df.num$V1 < 300, ]
Suisu1_TMHMM_no_TMD_300<- Suisu1_df.num[Suisu1_df.num$V1 < 300, ]
Suipla1_TMHMM_no_TMD_300<- Suipla1_df.num[Suipla1_df.num$V1 < 300, ]
Suipic1_TMHMM_no_TMD_300<- Suipic1_df.num[Suipic1_df.num$V1 < 300, ]
Suipal1_TMHMM_no_TMD_300<- Suipal1_df.num[Suipal1_df.num$V1 < 300, ]
Suiocc1_TMHMM_no_TMD_300<- Suiocc1_df.num[Suiocc1_df.num$V1 < 300, ]
Suilu4_TMHMM_no_TMD_300<- Suilu4_df.num[Suilu4_df.num$V1 < 300, ]
Suilak1_TMHMM_no_TMD_300<- Suilak1_df.num[Suilak1_df.num$V1 < 300, ]
Suihi1_TMHMM_no_TMD_300<- Suihi1_df.num[Suihi1_df.num$V1 < 300, ]
Suigr1_TMHMM_no_TMD_300<- Suigr1_df.num[Suigr1_df.num$V1 < 300, ]
Suidec1_TMHMM_no_TMD_300<- Suidec1_df.num[Suidec1_df.num$V1 < 300, ]
Suicot1_TMHMM_no_TMD_300<- Suicot1_df.num[Suicot1_df.num$V1 < 300, ]
Suicli1_TMHMM_no_TMD_300<- Suicli1_df.num[Suicli1_df.num$V1 < 300, ]
Suibr2_TMHMM_no_TMD_300<- Suibr2_df.num[Suibr2_df.num$V1 < 300, ]
Suibov1_TMHMM_no_TMD_300<- Suibov1_df.num[Suibov1_df.num$V1 < 300, ]
Suiamp1_TMHMM_no_TMD_300<- Suiamp1_df.num[Suiamp1_df.num$V1 < 300, ]
Suiame1_TMHMM_no_TMD_300<- Suiame1_df.num[Suiame1_df.num$V1 < 300, ]
Suifus1_TMHMM_no_TMD_300<- Suifus1_df.num[Suifus1_df.num$V1 < 300, ]

#non-Suillus set
Rhivul1_TMHMM_no_TMD_300<- Rhivul1_df.num[Rhivul1_df.num$V1 < 300, ]
Rhitru1_TMHMM_no_TMD_300<- Rhitru1_df.num[Rhitru1_df.num$V1 < 300, ]
Amamu1_TMHMM_no_TMD_300<- Amamu1_df.num[Amamu1_df.num$V1 < 300, ]
Hebcy2_TMHMM_no_TMD_300<- Hebcy2_df.num[Hebcy2_df.num$V1 < 300, ]
Lacbi2_TMHMM_no_TMD_300<- Lacbi2_df.num[Lacbi2_df.num$V1 < 300, ]
Paxin1_TMHMM_no_TMD_300<- Paxin1_df.num[Paxin1_df.num$V1 < 300, ]
Pilcr1_TMHMM_no_TMD_300<- Pilcr1_df.num[Pilcr1_df.num$V1 < 300, ]
Pismi1_TMHMM_no_TMD_300<- Pismi1_df.num[Pismi1_df.num$V1 < 300, ]
Sclci1_TMHMM_no_TMD_300<- Sclci1_df.num[Sclci1_df.num$V1 < 300, ]
Theter1_TMHMM_no_TMD_300<- Theter1_df.num[Theter1_df.num$V1 < 300, ]
Thega1_TMHMM_no_TMD_300<- Thega1_df.num[Thega1_df.num$V1 < 300, ]
Rhivi1_TMHMM_no_TMD_300<- Rhivi1_df.num[Rhivi1_df.num$V1 < 300, ]
Rhives1_TMHMM_no_TMD_300<- Rhives1_df.num[Rhives1_df.num$V1 < 300, ]
Rhisa1_TMHMM_no_TMD_300<- Rhisa1_df.num[Rhisa1_df.num$V1 < 300, ]
Ruscom1_TMHMM_no_TMD_300<- Ruscom1_df.num[Ruscom1_df.num$V1 < 300, ]
Rusbre1_TMHMM_no_TMD_300<- Rusbre1_df.num[Rusbre1_df.num$V1 < 300, ]
Pisti1_TMHMM_no_TMD_300<- Pisti1_df.num[Pisti1_df.num$V1 < 300, ]
Lacam2_TMHMM_no_TMD_300<- Lacam2_df.num[Lacam2_df.num$V1 < 300, ]
Hydru2_TMHMM_no_TMD_300<- Hydru2_df.num[Hydru2_df.num$V1 < 300, ]
Gaumor1_TMHMM_no_TMD_300<- Gaumor1_df.num[Gaumor1_df.num$V1 < 300, ]
Gyrli1_TMHMM_no_TMD_300<- Gyrli1_df.num[Gyrli1_df.num$V1 < 300, ]
Cananz1_TMHMM_no_TMD_300<- Cananz1_df.num[Cananz1_df.num$V1 < 300, ]
Hyssto1_TMHMM_no_TMD_300<- Hyssto1_df.num[Hyssto1_df.num$V1 < 300, ]


#get #SSP's per genome
#get summary numbers for each of these and attach them to the previous output file. 
Suivar1<- nrow(Suivar1_TMHMM_no_TMD_300)
Suitom1<- nrow(Suitom1_TMHMM_no_TMD_300)
Suisub1<- nrow(Suisub1_TMHMM_no_TMD_300)
Suisu1<- nrow(Suisu1_TMHMM_no_TMD_300)
Suipla1<- nrow(Suipla1_TMHMM_no_TMD_300)
Suipic1<- nrow(Suipic1_TMHMM_no_TMD_300)
Suipal1<- nrow(Suipal1_TMHMM_no_TMD_300)
Suiocc1<- nrow(Suiocc1_TMHMM_no_TMD_300)
Suilu4<- nrow(Suilu4_TMHMM_no_TMD_300)
Suilak1<- nrow(Suilak1_TMHMM_no_TMD_300)
Suihi1<- nrow(Suihi1_TMHMM_no_TMD_300)
Suigr1<- nrow(Suigr1_TMHMM_no_TMD_300)
Suidec1<- nrow(Suidec1_TMHMM_no_TMD_300)
Suicot1<- nrow(Suicot1_TMHMM_no_TMD_300)
Suicli1<- nrow(Suicli1_TMHMM_no_TMD_300)
Suibr2<- nrow(Suibr2_TMHMM_no_TMD_300)
Suibov1<- nrow(Suibov1_TMHMM_no_TMD_300)
Suiamp1<- nrow(Suiamp1_TMHMM_no_TMD_300)
Suiame1<- nrow(Suiame1_TMHMM_no_TMD_300)
Suifus1<- nrow(Suifus1_TMHMM_no_TMD_300)

#non-Suillus set
Rhivul1<- nrow(Rhivul1_TMHMM_no_TMD_300)
Rhitru1<- nrow(Rhitru1_TMHMM_no_TMD_300)
Amamu1<- nrow(Amamu1_TMHMM_no_TMD_300)
Hebcy2<- nrow(Hebcy2_TMHMM_no_TMD_300)
Lacbi2<- nrow(Lacbi2_TMHMM_no_TMD_300)
Paxin1<- nrow(Paxin1_TMHMM_no_TMD_300)
Pilcr1<- nrow(Pilcr1_TMHMM_no_TMD_300)
Pismi1<- nrow(Pismi1_TMHMM_no_TMD_300)
Sclci1<- nrow(Sclci1_TMHMM_no_TMD_300)
Theter1<- nrow(Theter1_TMHMM_no_TMD_300)
Thega1<- nrow(Thega1_TMHMM_no_TMD_300)
Rhivi1<- nrow(Rhivi1_TMHMM_no_TMD_300)
Rhives1<- nrow(Rhives1_TMHMM_no_TMD_300)
Rhisa1<- nrow(Rhisa1_TMHMM_no_TMD_300)
Ruscom1<- nrow(Ruscom1_TMHMM_no_TMD_300)
Rusbre1<- nrow(Rusbre1_TMHMM_no_TMD_300)
Pisti1<- nrow(Pisti1_TMHMM_no_TMD_300)
Lacam2<- nrow(Lacam2_TMHMM_no_TMD_300)
Hydru2<- nrow(Hydru2_TMHMM_no_TMD_300)
Gaumor1<- nrow(Gaumor1_TMHMM_no_TMD_300)
Gyrli1<- nrow(Gyrli1_TMHMM_no_TMD_300)
Cananz1<- nrow(Cananz1_TMHMM_no_TMD_300)
Hyssto1<- nrow(Hyssto1_TMHMM_no_TMD_300)

#add the ttotals 
totals.2<- data.frame(cbind(Suivar1, 
                            Suitom1, 
                            Suisub1, 
                            Suisu1, 
                            Suipla1, 
                            Suipic1, 
                            Suipal1, 
                            Suiocc1, 
                            Suilu4, 
                            Suilak1, 
                            Suihi1, 
                            Suidec1, 
                            Suicot1, 
                            Suicli1, 
                            Suibr2, 
                            Suibov1, 
                            Suiamp1, 
                            Suiame1,
                            Suifus1,
                            Rhivul1,
                            Rhitru1,
                            Amamu1,
                            Hebcy2,
                            Lacbi2,
                            Paxin1,
                            Pilcr1,
                            Pismi1,
                            Sclci1, 
                            Theter1,
                            Thega1,
                            Rhivi1,
                            Rhives1,
                            Rhisa1,
                            Ruscom1,
                            Rusbre1,
                            Pisti1,
                            Lacam2,
                            Hydru2,
                            Gaumor1,
                            Gyrli1,
                            Cananz1,
                            Hyssto1
                            ),
                      row.names = "#SSPs_signalP,TMHMM,lt_300aa")
totals.2[1,]

#write out totals to a file
#results<- rbind(totals, totals.2[1,])
write.csv(totals.2, quote = FALSE, file = "SSP_totals.csv")


#now parse out the SSP's from the original aa sequences. 
#read in the whole proteome fastas
Suivar1_in<- seqinr::read.fasta(file = "Suivar1_wo_stops.fasta", 
                                seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suitom1_in<- seqinr::read.fasta(file = "Suitom1_wo_stops.fasta", 
                                seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suisub1_in<- seqinr::read.fasta(file = "Suisub1_wo_stops.fasta", 
                                seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suisu1_in<- seqinr::read.fasta(file = "Suisu1_wo_stops.fasta", 
                               seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suipla1_in<- seqinr::read.fasta(file = "Suipla1_wo_stops.fasta", 
                                seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suipic1_in<- seqinr::read.fasta(file = "Suipic1_wo_stops.fasta", 
                                seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suipal1_in<- seqinr::read.fasta(file = "Suipal1_wo_stops.fasta", 
                                seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suiocc1_in<- seqinr::read.fasta(file = "Suiocc1_wo_stops.fasta", 
                                seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suilu4_in<- seqinr::read.fasta(file = "Suilu4_wo_stops.fasta", 
                               seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suilak1_in<- seqinr::read.fasta(file = "Suilak1_wo_stops.fasta", 
                                seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suihi1_in<- seqinr::read.fasta(file = "Suihi1_wo_stops.fasta", 
                               seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suigr1_in<- seqinr::read.fasta(file = "Suigr1_wo_stops.fasta", 
                               seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suidec1_in<- seqinr::read.fasta(file = "Suidec1_wo_stops.fasta", 
                                seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suicot1_in<- seqinr::read.fasta(file = "Suicot1_wo_stops.fasta", 
                                seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suicli1_in<- seqinr::read.fasta(file = "Suicli1_wo_stops.fasta", 
                                seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suibr2_in<- seqinr::read.fasta(file = "Suibr2_wo_stops.fasta", 
                               seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suibov1_in<- seqinr::read.fasta(file = "Suibov1_wo_stops.fasta", 
                                seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suiamp1_in<- seqinr::read.fasta(file = "Suiamp1_wo_stops.fasta", 
                                seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suiame1_in<- seqinr::read.fasta(file = "Suiame1_wo_stops.fasta", 
                                seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suifus1_in<- seqinr::read.fasta(file = "Suifus1_wo_stops.fasta", 
                                seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
#non-Suillus set
Rhivul1_in<- seqinr::read.fasta(file = "Rhivul1_wo_stops.fasta", 
                                seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Rhitru1_in<- seqinr::read.fasta(file = "Rhitru1_wo_stops.fasta", 
                                seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Amamu1_in<- seqinr::read.fasta(file = "Amamu1_wo_stops.fasta", 
                               seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Hebcy2_in<- seqinr::read.fasta(file = "Hebcy2_wo_stops.fasta", 
                               seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Lacbi2_in<- seqinr::read.fasta(file = "Lacbi2_wo_stops.fasta", 
                               seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Paxin1_in<- seqinr::read.fasta(file = "Paxin1_wo_stops.fasta", 
                               seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Pilcr1_in<- seqinr::read.fasta(file = "Pilcr1_wo_stops.fasta", 
                               seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Pismi1_in<- seqinr::read.fasta(file = "Pismi1_wo_stops.fasta", 
                               seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Sclci1_in<- seqinr::read.fasta(file = "Sclci1_wo_stops.fasta", 
                               seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Theter1_in<- seqinr::read.fasta(file = "Theter1_wo_stops.fasta", 
                               seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Thega1_in<- seqinr::read.fasta(file = "Thega1_wo_stops.fasta", 
                               seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Rhivi1_in<- seqinr::read.fasta(file = "Rhivi1_wo_stops.fasta", 
                               seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Rhives1_in<- seqinr::read.fasta(file = "Rhives1_wo_stops.fasta", 
                               seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Rhisa1_in<- seqinr::read.fasta(file = "Rhisa1_wo_stops.fasta", 
                               seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Ruscom1_in<- seqinr::read.fasta(file = "Ruscom1_wo_stops.fasta", 
                               seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Rusbre1_in<- seqinr::read.fasta(file = "Rusbre1_wo_stops.fasta", 
                               seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Pisti1_in<- seqinr::read.fasta(file = "Pisti1_wo_stops.fasta", 
                               seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Lacam2_in<- seqinr::read.fasta(file = "Lacam2_wo_stops.fasta", 
                               seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Hydru2_in<- seqinr::read.fasta(file = "Hydru2_wo_stops.fasta", 
                               seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Gaumor1_in<- seqinr::read.fasta(file = "Gaumor1_wo_stops.fasta", 
                               seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Gyrli1_in<- seqinr::read.fasta(file = "Gyrli1_wo_stops.fasta", 
                               seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Cananz1_in<- seqinr::read.fasta(file = "Cananz1_wo_stops.fasta", 
                               seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Hyssto1_in<- seqinr::read.fasta(file = "Hyssto1_wo_stops.fasta", 
                               seqtype = "AA",as.string = TRUE, set.attributes = FALSE)


#function to isolate the gene ID
isolate_ID <- function(in_df) { 
  a<- data.frame(strsplit(as.character(in_df[,1]), split=" "))
  b<- data.frame(a[1,])
  c<- data.table::transpose(b)
  return= c
}



#run function to isolate gene IDs  
Suivar1_subset <- data.frame(substr(Suivar1_TMHMM_no_TMD_300$X1, 3,nchar(Suivar1_TMHMM_no_TMD_300$X1)))
Suivar1.2_subset<- isolate_ID(Suivar1_subset) 
Suivar1.3_subset<- data.frame(lapply(Suivar1.2_subset, gsub, pattern=' ', replacement=''))

Suitom1_subset <- data.frame(substr(Suitom1_TMHMM_no_TMD_300$X1, 3,nchar(Suitom1_TMHMM_no_TMD_300$X1)))
Suitom1.2_subset<- isolate_ID(Suitom1_subset) 
Suitom1.3_subset<- data.frame(lapply(Suitom1.2_subset, gsub, pattern=' ', replacement=''))

Suisub1_subset <- data.frame(substr(Suisub1_TMHMM_no_TMD_300$X1, 3,nchar(Suisub1_TMHMM_no_TMD_300$X1)))
Suisub1.2_subset<- isolate_ID(Suisub1_subset) 
Suisub1.3_subset<- data.frame(lapply(Suisub1.2_subset, gsub, pattern=' ', replacement=''))

Suisu1_subset <- data.frame(substr(Suisu1_TMHMM_no_TMD_300$X1, 3,nchar(Suisu1_TMHMM_no_TMD_300$X1)))
Suisu1.2_subset<- isolate_ID(Suisu1_subset) 
Suisu1.3_subset<- data.frame(lapply(Suisu1.2_subset, gsub, pattern=' ', replacement=''))

Suipla1_subset <- data.frame(substr(Suipla1_TMHMM_no_TMD_300$X1, 3,nchar(Suipla1_TMHMM_no_TMD_300$X1)))
Suipla1.2_subset<- isolate_ID(Suipla1_subset) 
Suipla1.3_subset<- data.frame(lapply(Suipla1.2_subset, gsub, pattern=' ', replacement=''))

Suipic1_subset <- data.frame(substr(Suipic1_TMHMM_no_TMD_300$X1, 3,nchar(Suipic1_TMHMM_no_TMD_300$X1)))
Suipic1.2_subset<- isolate_ID(Suipic1_subset) 
Suipic1.3_subset<- data.frame(lapply(Suipic1.2_subset, gsub, pattern=' ', replacement=''))

Suipal1_subset <- data.frame(substr(Suipal1_TMHMM_no_TMD_300$X1, 3,nchar(Suipal1_TMHMM_no_TMD_300$X1)))
Suipal1.2_subset<- isolate_ID(Suipal1_subset) 
Suipal1.3_subset<- data.frame(lapply(Suipal1.2_subset, gsub, pattern=' ', replacement=''))

Suiocc1_subset <- data.frame(substr(Suiocc1_TMHMM_no_TMD_300$X1, 3,nchar(Suiocc1_TMHMM_no_TMD_300$X1)))
Suiocc1.2_subset<- isolate_ID(Suiocc1_subset) 
Suiocc1.3_subset<- data.frame(lapply(Suiocc1.2_subset, gsub, pattern=' ', replacement=''))

Suilu4_subset <- data.frame(substr(Suilu4_TMHMM_no_TMD_300$X1, 3,nchar(Suilu4_TMHMM_no_TMD_300$X1)))
Suilu4.2_subset<- isolate_ID(Suilu4_subset) 
Suilu4.3_subset<- data.frame(lapply(Suilu4.2_subset, gsub, pattern=' ', replacement=''))

Suilak1_subset <- data.frame(substr(Suilak1_TMHMM_no_TMD_300$X1, 3,nchar(Suilak1_TMHMM_no_TMD_300$X1)))
Suilak1.2_subset<- isolate_ID(Suilak1_subset) 
Suilak1.3_subset<- data.frame(lapply(Suilak1.2_subset, gsub, pattern=' ', replacement=''))

Suihi1_subset <- data.frame(substr(Suihi1_TMHMM_no_TMD_300$X1, 3,nchar(Suihi1_TMHMM_no_TMD_300$X1)))
Suihi1.2_subset<- isolate_ID(Suihi1_subset) 
Suihi1.3_subset<- data.frame(lapply(Suihi1.2_subset, gsub, pattern=' ', replacement=''))

Suigr1_subset <- data.frame(substr(Suigr1_TMHMM_no_TMD_300$X1, 3,nchar(Suigr1_TMHMM_no_TMD_300$X1)))
Suigr1.2_subset<- isolate_ID(Suigr1_subset) 
Suigr1.3_subset<- data.frame(lapply(Suigr1.2_subset, gsub, pattern=' ', replacement=''))

Suidec1_subset <- data.frame(substr(Suidec1_TMHMM_no_TMD_300$X1, 3,nchar(Suidec1_TMHMM_no_TMD_300$X1)))
Suidec1.2_subset<- isolate_ID(Suidec1_subset) 
Suidec1.3_subset<- data.frame(lapply(Suidec1.2_subset, gsub, pattern=' ', replacement=''))

Suicot1_subset <- data.frame(substr(Suicot1_TMHMM_no_TMD_300$X1, 3,nchar(Suicot1_TMHMM_no_TMD_300$X1)))
Suicot1.2_subset<- isolate_ID(Suicot1_subset) 
Suicot1.3_subset<- data.frame(lapply(Suicot1.2_subset, gsub, pattern=' ', replacement=''))

Suicli1_subset <- data.frame(substr(Suicli1_TMHMM_no_TMD_300$X1, 3,nchar(Suicli1_TMHMM_no_TMD_300$X1)))
Suicli1.2_subset<- isolate_ID(Suicli1_subset) 
Suicli1.3_subset<- data.frame(lapply(Suicli1.2_subset, gsub, pattern=' ', replacement=''))

Suibr2_subset <- data.frame(substr(Suibr2_TMHMM_no_TMD_300$X1, 3,nchar(Suibr2_TMHMM_no_TMD_300$X1)))
Suibr2.2_subset<- isolate_ID(Suibr2_subset) 
Suibr2.3_subset<- data.frame(lapply(Suibr2.2_subset, gsub, pattern=' ', replacement=''))

Suibov1_subset <- data.frame(substr(Suibov1_TMHMM_no_TMD_300$X1, 3,nchar(Suibov1_TMHMM_no_TMD_300$X1)))
Suibov1.2_subset<- isolate_ID(Suibov1_subset) 
Suibov1.3_subset<- data.frame(lapply(Suibov1.2_subset, gsub, pattern=' ', replacement=''))

Suiamp1_subset <- data.frame(substr(Suiamp1_TMHMM_no_TMD_300$X1, 3,nchar(Suiamp1_TMHMM_no_TMD_300$X1)))
Suiamp1.2_subset<- isolate_ID(Suiamp1_subset) 
Suiamp1.3_subset<- data.frame(lapply(Suiamp1.2_subset, gsub, pattern=' ', replacement=''))

Suiame1_subset <- data.frame(substr(Suiame1_TMHMM_no_TMD_300$X1, 3,nchar(Suiame1_TMHMM_no_TMD_300$X1)))
Suiame1.2_subset<- isolate_ID(Suiame1_subset) 
Suiame1.3_subset<- data.frame(lapply(Suiame1.2_subset, gsub, pattern=' ', replacement=''))

Suifus1_subset <- data.frame(substr(Suifus1_TMHMM_no_TMD_300$X1, 3,nchar(Suifus1_TMHMM_no_TMD_300$X1)))
Suifus1.2_subset<- isolate_ID(Suifus1_subset) 
Suifus1.3_subset<- data.frame(lapply(Suifus1.2_subset, gsub, pattern=' ', replacement=''))


#non-Suillus set
Rhivul1_subset <- data.frame(substr(Rhivul1_TMHMM_no_TMD_300$X1, 3,nchar(Rhivul1_TMHMM_no_TMD_300$X1)))
Rhivul1.2_subset<- isolate_ID(Rhivul1_subset) 
Rhivul1.3_subset<- data.frame(lapply(Rhivul1.2_subset, gsub, pattern=' ', replacement=''))

Rhitru1_subset <- data.frame(substr(Rhitru1_TMHMM_no_TMD_300$X1, 3,nchar(Rhitru1_TMHMM_no_TMD_300$X1)))
Rhitru1.2_subset<- isolate_ID(Rhitru1_subset) 
Rhitru1.3_subset<- data.frame(lapply(Rhitru1.2_subset, gsub, pattern=' ', replacement=''))

Amamu1_subset <- data.frame(substr(Amamu1_TMHMM_no_TMD_300$X1, 3,nchar(Amamu1_TMHMM_no_TMD_300$X1)))
Amamu1.2_subset<- isolate_ID(Amamu1_subset) 
Amamu1.3_subset<- data.frame(lapply(Amamu1.2_subset, gsub, pattern=' ', replacement=''))

Hebcy2_subset <- data.frame(substr(Hebcy2_TMHMM_no_TMD_300$X1, 3,nchar(Hebcy2_TMHMM_no_TMD_300$X1)))
Hebcy2.2_subset<- isolate_ID(Hebcy2_subset) 
Hebcy2.3_subset<- data.frame(lapply(Hebcy2.2_subset, gsub, pattern=' ', replacement=''))

Lacbi2_subset <- data.frame(substr(Lacbi2_TMHMM_no_TMD_300$X1, 3,nchar(Lacbi2_TMHMM_no_TMD_300$X1)))
Lacbi2.2_subset<- isolate_ID(Lacbi2_subset) 
Lacbi2.3_subset<- data.frame(lapply(Lacbi2.2_subset, gsub, pattern=' ', replacement=''))

Paxin1_subset <- data.frame(substr(Paxin1_TMHMM_no_TMD_300$X1, 3,nchar(Paxin1_TMHMM_no_TMD_300$X1)))
Paxin1.2_subset<- isolate_ID(Paxin1_subset) 
Paxin1.3_subset<- data.frame(lapply(Paxin1.2_subset, gsub, pattern=' ', replacement=''))

Pilcr1_subset <- data.frame(substr(Pilcr1_TMHMM_no_TMD_300$X1, 3,nchar(Pilcr1_TMHMM_no_TMD_300$X1)))
Pilcr1.2_subset<- isolate_ID(Pilcr1_subset) 
Pilcr1.3_subset<- data.frame(lapply(Pilcr1.2_subset, gsub, pattern=' ', replacement=''))

Pismi1_subset <- data.frame(substr(Pismi1_TMHMM_no_TMD_300$X1, 3,nchar(Pismi1_TMHMM_no_TMD_300$X1)))
Pismi1.2_subset<- isolate_ID(Pismi1_subset) 
Pismi1.3_subset<- data.frame(lapply(Pismi1.2_subset, gsub, pattern=' ', replacement=''))

Sclci1_subset <- data.frame(substr(Sclci1_TMHMM_no_TMD_300$X1, 3,nchar(Sclci1_TMHMM_no_TMD_300$X1)))
Sclci1.2_subset<- isolate_ID(Sclci1_subset) 
Sclci1.3_subset<- data.frame(lapply(Sclci1.2_subset, gsub, pattern=' ', replacement=''))

Theter1_subset <- data.frame(substr(Theter1_TMHMM_no_TMD_300$X1, 3,nchar(Theter1_TMHMM_no_TMD_300$X1)))
Theter1.2_subset<- isolate_ID(Theter1_subset) 
Theter1.3_subset<- data.frame(lapply(Theter1.2_subset, gsub, pattern=' ', replacement=''))

Thega1_subset <- data.frame(substr(Thega1_TMHMM_no_TMD_300$X1, 3,nchar(Thega1_TMHMM_no_TMD_300$X1)))
Thega1.2_subset<- isolate_ID(Thega1_subset) 
Thega1.3_subset<- data.frame(lapply(Thega1.2_subset, gsub, pattern=' ', replacement=''))

Rhivi1_subset <- data.frame(substr(Rhivi1_TMHMM_no_TMD_300$X1, 3,nchar(Rhivi1_TMHMM_no_TMD_300$X1)))
Rhivi1.2_subset<- isolate_ID(Rhivi1_subset) 
Rhivi1.3_subset<- data.frame(lapply(Rhivi1.2_subset, gsub, pattern=' ', replacement=''))

Rhives1_subset <- data.frame(substr(Rhives1_TMHMM_no_TMD_300$X1, 3,nchar(Rhives1_TMHMM_no_TMD_300$X1)))
Rhives1.2_subset<- isolate_ID(Rhives1_subset) 
Rhives1.3_subset<- data.frame(lapply(Rhives1.2_subset, gsub, pattern=' ', replacement=''))

Rhisa1_subset <- data.frame(substr(Rhisa1_TMHMM_no_TMD_300$X1, 3,nchar(Rhisa1_TMHMM_no_TMD_300$X1)))
Rhisa1.2_subset<- isolate_ID(Rhisa1_subset) 
Rhisa1.3_subset<- data.frame(lapply(Rhisa1.2_subset, gsub, pattern=' ', replacement=''))

Ruscom1_subset <- data.frame(substr(Ruscom1_TMHMM_no_TMD_300$X1, 3,nchar(Ruscom1_TMHMM_no_TMD_300$X1)))
Ruscom1.2_subset<- isolate_ID(Ruscom1_subset) 
Ruscom1.3_subset<- data.frame(lapply(Ruscom1.2_subset, gsub, pattern=' ', replacement=''))

Rusbre1_subset <- data.frame(substr(Rusbre1_TMHMM_no_TMD_300$X1, 3,nchar(Rusbre1_TMHMM_no_TMD_300$X1)))
Rusbre1.2_subset<- isolate_ID(Rusbre1_subset) 
Rusbre1.3_subset<- data.frame(lapply(Rusbre1.2_subset, gsub, pattern=' ', replacement=''))

Pisti1_subset <- data.frame(substr(Pisti1_TMHMM_no_TMD_300$X1, 3,nchar(Pisti1_TMHMM_no_TMD_300$X1)))
Pisti1.2_subset<- isolate_ID(Pisti1_subset) 
Pisti1.3_subset<- data.frame(lapply(Pisti1.2_subset, gsub, pattern=' ', replacement=''))

Lacam2_subset <- data.frame(substr(Lacam2_TMHMM_no_TMD_300$X1, 3,nchar(Lacam2_TMHMM_no_TMD_300$X1)))
Lacam2.2_subset<- isolate_ID(Lacam2_subset) 
Lacam2.3_subset<- data.frame(lapply(Lacam2.2_subset, gsub, pattern=' ', replacement=''))

Hydru2_subset <- data.frame(substr(Hydru2_TMHMM_no_TMD_300$X1, 3,nchar(Hydru2_TMHMM_no_TMD_300$X1)))
Hydru2.2_subset<- isolate_ID(Hydru2_subset) 
Hydru2.3_subset<- data.frame(lapply(Hydru2.2_subset, gsub, pattern=' ', replacement=''))

Gaumor1_subset <- data.frame(substr(Gaumor1_TMHMM_no_TMD_300$X1, 3,nchar(Gaumor1_TMHMM_no_TMD_300$X1)))
Gaumor1.2_subset<- isolate_ID(Gaumor1_subset) 
Gaumor1.3_subset<- data.frame(lapply(Gaumor1.2_subset, gsub, pattern=' ', replacement=''))

Gyrli1_subset <- data.frame(substr(Gyrli1_TMHMM_no_TMD_300$X1, 3,nchar(Gyrli1_TMHMM_no_TMD_300$X1)))
Gyrli1.2_subset<- isolate_ID(Gyrli1_subset) 
Gyrli1.3_subset<- data.frame(lapply(Gyrli1.2_subset, gsub, pattern=' ', replacement=''))

Cananz1_subset <- data.frame(substr(Cananz1_TMHMM_no_TMD_300$X1, 3,nchar(Cananz1_TMHMM_no_TMD_300$X1)))
Cananz1.2_subset<- isolate_ID(Cananz1_subset) 
Cananz1.3_subset<- data.frame(lapply(Cananz1.2_subset, gsub, pattern=' ', replacement=''))

Hyssto1_subset <- data.frame(substr(Hyssto1_TMHMM_no_TMD_300$X1, 3,nchar(Hyssto1_TMHMM_no_TMD_300$X1)))
Hyssto1.2_subset<- isolate_ID(Hyssto1_subset) 
Hyssto1.3_subset<- data.frame(lapply(Hyssto1.2_subset, gsub, pattern=' ', replacement=''))


#YOU ARE HERE:
#change all "|" to "_" in the input files
names(Suivar1_in)<-lapply(names(Suivar1_in), gsub, pattern="\\|", replacement='_')
names(Suitom1_in)<-lapply(names(Suitom1_in), gsub, pattern="\\|", replacement='_')
names(Suisub1_in)<-lapply(names(Suisub1_in), gsub, pattern="\\|", replacement='_')
names(Suisu1_in)<-lapply(names(Suisu1_in), gsub, pattern="\\|", replacement='_')
names(Suipla1_in)<-lapply(names(Suipla1_in), gsub, pattern="\\|", replacement='_')
names(Suipic1_in)<-lapply(names(Suipic1_in), gsub, pattern="\\|", replacement='_')
names(Suipal1_in)<-lapply(names(Suipal1_in), gsub, pattern="\\|", replacement='_')
names(Suiocc1_in)<-lapply(names(Suiocc1_in), gsub, pattern="\\|", replacement='_')
names(Suilu4_in)<-lapply(names(Suilu4_in), gsub, pattern="\\|", replacement='_')
names(Suilak1_in)<-lapply(names(Suilak1_in), gsub, pattern="\\|", replacement='_')
names(Suihi1_in)<-lapply(names(Suihi1_in), gsub, pattern="\\|", replacement='_')
names(Suigr1_in)<-lapply(names(Suigr1_in), gsub, pattern="\\|", replacement='_')
names(Suidec1_in)<-lapply(names(Suidec1_in), gsub, pattern="\\|", replacement='_')
names(Suicot1_in)<-lapply(names(Suicot1_in), gsub, pattern="\\|", replacement='_')
names(Suicli1_in)<-lapply(names(Suicli1_in), gsub, pattern="\\|", replacement='_')
names(Suibr2_in)<-lapply(names(Suibr2_in), gsub, pattern="\\|", replacement='_')
names(Suibov1_in)<-lapply(names(Suibov1_in), gsub, pattern="\\|", replacement='_')
names(Suiamp1_in)<-lapply(names(Suiamp1_in), gsub, pattern="\\|", replacement='_')
names(Suiame1_in)<-lapply(names(Suiame1_in), gsub, pattern="\\|", replacement='_')
names(Suifus1_in)<-lapply(names(Suifus1_in), gsub, pattern="\\|", replacement='_')
#non-suillus set
names(Rhivul1_in)<-lapply(names(Rhivul1_in), gsub, pattern="\\|", replacement='_')
names(Rhitru1_in)<-lapply(names(Rhitru1_in), gsub, pattern="\\|", replacement='_')
names(Amamu1_in)<-lapply(names(Amamu1_in), gsub, pattern="\\|", replacement='_')
names(Hebcy2_in)<-lapply(names(Hebcy2_in), gsub, pattern="\\|", replacement='_')
names(Lacbi2_in)<-lapply(names(Lacbi2_in), gsub, pattern="\\|", replacement='_')
names(Paxin1_in)<-lapply(names(Paxin1_in), gsub, pattern="\\|", replacement='_')
names(Pilcr1_in)<-lapply(names(Pilcr1_in), gsub, pattern="\\|", replacement='_')
names(Pismi1_in)<-lapply(names(Pismi1_in), gsub, pattern="\\|", replacement='_')
names(Sclci1_in)<-lapply(names(Sclci1_in), gsub, pattern="\\|", replacement='_')
names(Theter1_in)<-lapply(names(Theter1_in), gsub, pattern="\\|", replacement='_')
names(Thega1_in)<-lapply(names(Thega1_in), gsub, pattern="\\|", replacement='_')
names(Rhivi1_in)<-lapply(names(Rhivi1_in), gsub, pattern="\\|", replacement='_')
names(Rhives1_in)<-lapply(names(Rhives1_in), gsub, pattern="\\|", replacement='_')
names(Rhisa1_in)<-lapply(names(Rhisa1_in), gsub, pattern="\\|", replacement='_')
names(Ruscom1_in)<-lapply(names(Ruscom1_in), gsub, pattern="\\|", replacement='_')
names(Rusbre1_in)<-lapply(names(Rusbre1_in), gsub, pattern="\\|", replacement='_')
names(Pisti1_in)<-lapply(names(Pisti1_in), gsub, pattern="\\|", replacement='_')
names(Lacam2_in)<-lapply(names(Lacam2_in), gsub, pattern="\\|", replacement='_')
names(Hydru2_in)<-lapply(names(Hydru2_in), gsub, pattern="\\|", replacement='_')
names(Gaumor1_in)<-lapply(names(Gaumor1_in), gsub, pattern="\\|", replacement='_')
names(Gyrli1_in)<-lapply(names(Gyrli1_in), gsub, pattern="\\|", replacement='_')
names(Cananz1_in)<-lapply(names(Cananz1_in), gsub, pattern="\\|", replacement='_')
names(Hyssto1_in)<-lapply(names(Hyssto1_in), gsub, pattern="\\|", replacement='_')

#get fastas of only positive hits for EffectorP analysis 
Suivar1_SSP_fastas<- Suivar1_in[c(which(names(Suivar1_in) %in% Suivar1.3_subset$V1))]
Suitom1_SSP_fastas<- Suitom1_in[c(which(names(Suitom1_in) %in% Suitom1.3_subset$V1))]
Suisub1_SSP_fastas<- Suisub1_in[c(which(names(Suisub1_in) %in% Suisub1.3_subset$V1))]
Suisu1_SSP_fastas<- Suisu1_in[c(which(names(Suisu1_in) %in% Suisu1.3_subset$V1))]
Suipla1_SSP_fastas<- Suipla1_in[c(which(names(Suipla1_in) %in% Suipla1.3_subset$V1))]
Suipic1_SSP_fastas<- Suipic1_in[c(which(names(Suipic1_in) %in% Suipic1.3_subset$V1))]
Suipal1_SSP_fastas<- Suipal1_in[c(which(names(Suipal1_in) %in% Suipal1.3_subset$V1))]
Suiocc1_SSP_fastas<- Suiocc1_in[c(which(names(Suiocc1_in) %in% Suiocc1.3_subset$V1))]
Suilu4_SSP_fastas<- Suilu4_in[c(which(names(Suilu4_in) %in% Suilu4.3_subset$V1))]
Suilak1_SSP_fastas<- Suilak1_in[c(which(names(Suilak1_in) %in% Suilak1.3_subset$V1))]
Suihi1_SSP_fastas<- Suihi1_in[c(which(names(Suihi1_in) %in% Suihi1.3_subset$V1))]
Suigr1_SSP_fastas<- Suigr1_in[c(which(names(Suigr1_in) %in% Suigr1.3_subset$V1))]
Suidec1_SSP_fastas<- Suidec1_in[c(which(names(Suidec1_in) %in% Suidec1.3_subset$V1))]
Suicot1_SSP_fastas<- Suicot1_in[c(which(names(Suicot1_in) %in% Suicot1.3_subset$V1))]
Suicli1_SSP_fastas<- Suicli1_in[c(which(names(Suicli1_in) %in% Suicli1.3_subset$V1))]
Suibr2_SSP_fastas<- Suibr2_in[c(which(names(Suibr2_in) %in% Suibr2.3_subset$V1))]
Suibov1_SSP_fastas<- Suibov1_in[c(which(names(Suibov1_in) %in% Suibov1.3_subset$V1))]
Suiamp1_SSP_fastas<- Suiamp1_in[c(which(names(Suiamp1_in) %in% Suiamp1.3_subset$V1))]
Suiame1_SSP_fastas<- Suiame1_in[c(which(names(Suiame1_in) %in% Suiame1.3_subset$V1))]
Suifus1_SSP_fastas<- Suifus1_in[c(which(names(Suifus1_in) %in% Suifus1.3_subset$V1))]


#non-Suillus set
Rhivul1_SSP_fastas<- Rhivul1_in[c(which(names(Rhivul1_in) %in% Rhivul1.3_subset$V1))]
Rhitru1_SSP_fastas<- Rhitru1_in[c(which(names(Rhitru1_in) %in% Rhitru1.3_subset$V1))]
Amamu1_SSP_fastas<- Amamu1_in[c(which(names(Amamu1_in) %in% Amamu1.3_subset$V1))]
Hebcy2_SSP_fastas<- Hebcy2_in[c(which(names(Hebcy2_in) %in% Hebcy2.3_subset$V1))]
Lacbi2_SSP_fastas<- Lacbi2_in[c(which(names(Lacbi2_in) %in% Lacbi2.3_subset$V1))]
Paxin1_SSP_fastas<- Paxin1_in[c(which(names(Paxin1_in) %in% Paxin1.3_subset$V1))]
Pilcr1_SSP_fastas<- Pilcr1_in[c(which(names(Pilcr1_in) %in% Pilcr1.3_subset$V1))]
Pismi1_SSP_fastas<- Pismi1_in[c(which(names(Pismi1_in) %in% Pismi1.3_subset$V1))]
Sclci1_SSP_fastas<- Sclci1_in[c(which(names(Sclci1_in) %in% Sclci1.3_subset$V1))]
Theter1_SSP_fastas<- Theter1_in[c(which(names(Theter1_in) %in% Theter1.3_subset$V1))]
Thega1_SSP_fastas<- Thega1_in[c(which(names(Thega1_in) %in% Thega1.3_subset$V1))]
Rhivi1_SSP_fastas<- Rhivi1_in[c(which(names(Rhivi1_in) %in% Rhivi1.3_subset$V1))]
Rhives1_SSP_fastas<- Rhives1_in[c(which(names(Rhives1_in) %in% Rhives1.3_subset$V1))]
Rhisa1_SSP_fastas<- Rhisa1_in[c(which(names(Rhisa1_in) %in% Rhisa1.3_subset$V1))]
Ruscom1_SSP_fastas<- Ruscom1_in[c(which(names(Ruscom1_in) %in% Ruscom1.3_subset$V1))]
Rusbre1_SSP_fastas<- Rusbre1_in[c(which(names(Rusbre1_in) %in% Rusbre1.3_subset$V1))]
Pisti1_SSP_fastas<- Pisti1_in[c(which(names(Pisti1_in) %in% Pisti1.3_subset$V1))]
Lacam2_SSP_fastas<- Lacam2_in[c(which(names(Lacam2_in) %in% Lacam2.3_subset$V1))]
Hydru2_SSP_fastas<- Hydru2_in[c(which(names(Hydru2_in) %in% Hydru2.3_subset$V1))]
Gaumor1_SSP_fastas<- Gaumor1_in[c(which(names(Gaumor1_in) %in% Gaumor1.3_subset$V1))]
Gyrli1_SSP_fastas<- Gyrli1_in[c(which(names(Gyrli1_in) %in% Gyrli1.3_subset$V1))]
Cananz1_SSP_fastas<- Cananz1_in[c(which(names(Cananz1_in) %in% Cananz1.3_subset$V1))]
Hyssto1_SSP_fastas<- Hyssto1_in[c(which(names(Hyssto1_in) %in% Hyssto1.3_subset$V1))]


#make sure nothing screwy went on and the fasta are the correct length
#n SSP's 
totals.2[1,]
# length of a few files
length(Suivar1_SSP_fastas)
length(Suitom1_SSP_fastas)
length(Suiamp1_SSP_fastas)
length(Suisu1_SSP_fastas)
#numbers match

#print fastas to files to run EffectorP / Orthofinder
write.fasta(Suivar1_SSP_fastas, names = names(Suivar1_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suivar1_SSPs_for_EffectorP.fasta")
write.fasta(Suitom1_SSP_fastas, names = names(Suitom1_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suitom1_SSPs_for_EffectorP.fasta")
write.fasta(Suisub1_SSP_fastas, names = names(Suisub1_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suisub1_SSPs_for_EffectorP.fasta")
write.fasta(Suisu1_SSP_fastas, names = names(Suisu1_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suisu1_SSPs_for_EffectorP.fasta")
write.fasta(Suipla1_SSP_fastas, names = names(Suipla1_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suipla1_SSPs_for_EffectorP.fasta")
write.fasta(Suipic1_SSP_fastas, names = names(Suipic1_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suipic1_SSPs_for_EffectorP.fasta")
write.fasta(Suipal1_SSP_fastas, names = names(Suipal1_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suipal1_SSPs_for_EffectorP.fasta")
write.fasta(Suiocc1_SSP_fastas, names = names(Suiocc1_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suiocc1_SSPs_for_EffectorP.fasta")
write.fasta(Suilu4_SSP_fastas, names = names(Suilu4_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suilu4_SSPs_for_EffectorP.fasta")
write.fasta(Suilak1_SSP_fastas, names = names(Suilak1_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suilak1_SSPs_for_EffectorP.fasta")
write.fasta(Suihi1_SSP_fastas, names = names(Suihi1_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suihi1_SSPs_for_EffectorP.fasta")
write.fasta(Suigr1_SSP_fastas, names = names(Suigr1_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suigr1_SSPs_for_EffectorP.fasta")
write.fasta(Suidec1_SSP_fastas, names = names(Suidec1_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suidec1_SSPs_for_EffectorP.fasta")
write.fasta(Suicot1_SSP_fastas, names = names(Suicot1_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suicot1_SSPs_for_EffectorP.fasta")
write.fasta(Suicli1_SSP_fastas, names = names(Suicli1_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suicli1_SSPs_for_EffectorP.fasta")
write.fasta(Suibr2_SSP_fastas, names = names(Suibr2_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suibr2_SSPs_for_EffectorP.fasta")
write.fasta(Suibov1_SSP_fastas, names = names(Suibov1_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suibov1_SSPs_for_EffectorP.fasta")
write.fasta(Suiamp1_SSP_fastas, names = names(Suiamp1_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suiamp1_SSPs_for_EffectorP.fasta")
write.fasta(Suiame1_SSP_fastas, names = names(Suiame1_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suiame1_SSPs_for_EffectorP.fasta")
write.fasta(Suifus1_SSP_fastas, names = names(Suifus1_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suifus1_SSPs_for_EffectorP.fasta")

#non-Suillus set
write.fasta(Rhivul1_SSP_fastas, names = names(Rhivul1_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Rhivul1_SSPs_for_EffectorP.fasta")
write.fasta(Rhitru1_SSP_fastas, names = names(Rhitru1_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Rhitru1_SSPs_for_EffectorP.fasta")
write.fasta(Amamu1_SSP_fastas, names = names(Amamu1_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Amamu1_SSPs_for_EffectorP.fasta")
write.fasta(Hebcy2_SSP_fastas, names = names(Hebcy2_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Hebcy2_SSPs_for_EffectorP.fasta")
write.fasta(Lacbi2_SSP_fastas, names = names(Lacbi2_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Lacbi2_SSPs_for_EffectorP.fasta")
write.fasta(Paxin1_SSP_fastas, names = names(Paxin1_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Paxin1_SSPs_for_EffectorP.fasta")
write.fasta(Pilcr1_SSP_fastas, names = names(Pilcr1_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Pilcr1_SSPs_for_EffectorP.fasta")
write.fasta(Pismi1_SSP_fastas, names = names(Pismi1_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Pismi1_SSPs_for_EffectorP.fasta")
write.fasta(Sclci1_SSP_fastas, names = names(Sclci1_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Sclci1_SSPs_for_EffectorP.fasta")
write.fasta(Theter1_SSP_fastas, names = names(Theter1_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Theter1_SSPs_for_EffectorP.fasta")
write.fasta(Thega1_SSP_fastas, names = names(Thega1_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Thega1_SSPs_for_EffectorP.fasta")
write.fasta(Rhivi1_SSP_fastas, names = names(Rhivi1_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Rhivi1_SSPs_for_EffectorP.fasta")
write.fasta(Rhives1_SSP_fastas, names = names(Rhives1_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Rhives1_SSPs_for_EffectorP.fasta")
write.fasta(Rhisa1_SSP_fastas, names = names(Rhisa1_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Rhisa1_SSPs_for_EffectorP.fasta")
write.fasta(Ruscom1_SSP_fastas, names = names(Ruscom1_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Ruscom1_SSPs_for_EffectorP.fasta")
write.fasta(Rusbre1_SSP_fastas, names = names(Rusbre1_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Rusbre1_SSPs_for_EffectorP.fasta")
write.fasta(Pisti1_SSP_fastas, names = names(Pisti1_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Pisti1_SSPs_for_EffectorP.fasta")
write.fasta(Lacam2_SSP_fastas, names = names(Lacam2_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Lacam2_SSPs_for_EffectorP.fasta")
write.fasta(Hydru2_SSP_fastas, names = names(Hydru2_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Hydru2_SSPs_for_EffectorP.fasta")
write.fasta(Gaumor1_SSP_fastas, names = names(Gaumor1_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Gaumor1_SSPs_for_EffectorP.fasta")
write.fasta(Gyrli1_SSP_fastas, names = names(Gyrli1_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Gyrli1_SSPs_for_EffectorP.fasta")
write.fasta(Cananz1_SSP_fastas, names = names(Cananz1_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Cananz1_SSPs_for_EffectorP.fasta")
write.fasta(Hyssto1_SSP_fastas, names = names(Hyssto1_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Hyssto1_SSPs_for_EffectorP.fasta")


########
#now run EffectorP version2 on the web interface
#######


#read in EffectorP output files
Suivar1_effectors<- seqinr::read.fasta(file = "Suivar1_EffectorCandidates.fasta", 
                                       seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suitom1_effectors<- seqinr::read.fasta(file = "Suitom1_EffectorCandidates.fasta", 
                                       seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suisub1_effectors<- seqinr::read.fasta(file = "Suisub1_EffectorCandidates.fasta", 
                                       seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suisu1_effectors<- seqinr::read.fasta(file = "Suisu1_EffectorCandidates.fasta", 
                                      seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suipla1_effectors<- seqinr::read.fasta(file = "Suipla1_EffectorCandidates.fasta", 
                                       seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suipic1_effectors<- seqinr::read.fasta(file = "Suipic1_EffectorCandidates.fasta", 
                                       seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suipal1_effectors<- seqinr::read.fasta(file = "Suipal1_EffectorCandidates.fasta", 
                                       seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suiocc1_effectors<- seqinr::read.fasta(file = "Suiocc1_EffectorCandidates.fasta", 
                                       seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suilu4_effectors<- seqinr::read.fasta(file = "Suilu4_EffectorCandidates.fasta", 
                                      seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suilak1_effectors<- seqinr::read.fasta(file = "Suilak1_EffectorCandidates.fasta", 
                                       seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suihi1_effectors<- seqinr::read.fasta(file = "Suihi1_EffectorCandidates.fasta", 
                                      seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suigr1_effectors<- seqinr::read.fasta(file = "Suigr1_EffectorCandidates.fasta", 
                                      seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suidec1_effectors<- seqinr::read.fasta(file = "Suidec1_EffectorCandidates.fasta", 
                                       seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suicot1_effectors<- seqinr::read.fasta(file = "Suicot1_EffectorCandidates.fasta", 
                                       seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suicli1_effectors<- seqinr::read.fasta(file = "Suicli1_EffectorCandidates.fasta", 
                                       seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suibr2_effectors<- seqinr::read.fasta(file = "Suibr2_EffectorCandidates.fasta", 
                                      seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suibov1_effectors<- seqinr::read.fasta(file = "Suibov1_EffectorCandidates.fasta", 
                                       seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suiamp1_effectors<- seqinr::read.fasta(file = "Suiamp1_EffectorCandidates.fasta", 
                                       seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suiame1_effectors<- seqinr::read.fasta(file = "Suiame1_EffectorCandidates.fasta", 
                                       seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suifus1_effectors<- seqinr::read.fasta(file = "Suifus1_EffectorCandidates.fasta", 
                                       seqtype = "AA",as.string = TRUE, set.attributes = FALSE)


#non-Suillus set
Rhivul1_effectors<- seqinr::read.fasta(file = "Rhivul1_EffectorCandidates.fasta", 
                                       seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Rhitru1_effectors<- seqinr::read.fasta(file = "Rhitru1_EffectorCandidates.fasta", 
                                       seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Amamu1_effectors<- seqinr::read.fasta(file = "Amamu1_EffectorCandidates.fasta", 
                                      seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Hebcy2_effectors<- seqinr::read.fasta(file = "Hebcy2_EffectorCandidates.fasta", 
                                      seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Lacbi2_effectors<- seqinr::read.fasta(file = "Lacbi2_EffectorCandidates.fasta", 
                                      seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Paxin1_effectors<- seqinr::read.fasta(file = "Paxin1_EffectorCandidates.fasta", 
                                      seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Pilcr1_effectors<- seqinr::read.fasta(file = "Pilcr1_EffectorCandidates.fasta", 
                                      seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Pismi1_effectors<- seqinr::read.fasta(file = "Pismi1_EffectorCandidates.fasta", 
                                      seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Sclci1_effectors<- seqinr::read.fasta(file = "Sclci1_EffectorCandidates.fasta", 
                                      seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Theter1_effectors<- seqinr::read.fasta(file = "Theter1_EffectorCandidates.fasta", 
                                      seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Thega1_effectors<- seqinr::read.fasta(file = "Thega1_EffectorCandidates.fasta", 
                                       seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Rhivi1_effectors<- seqinr::read.fasta(file = "Rhivi1_EffectorCandidates.fasta", 
                                      seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Rhives1_effectors<- seqinr::read.fasta(file = "Rhives1_EffectorCandidates.fasta", 
                                      seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Rhisa1_effectors<- seqinr::read.fasta(file = "Rhisa1_EffectorCandidates.fasta", 
                                       seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Ruscom1_effectors<- seqinr::read.fasta(file = "Ruscom1_EffectorCandidates.fasta", 
                                      seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Rusbre1_effectors<- seqinr::read.fasta(file = "Rusbre1_EffectorCandidates.fasta", 
                                       seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Pisti1_effectors<- seqinr::read.fasta(file = "Pisti1_EffectorCandidates.fasta", 
                                       seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Lacam2_effectors<- seqinr::read.fasta(file = "Lacam2_EffectorCandidates.fasta", 
                                      seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Hydru2_effectors<- seqinr::read.fasta(file = "Hydru2_EffectorCandidates.fasta", 
                                      seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Gaumor1_effectors<- seqinr::read.fasta(file = "Gaumor1_EffectorCandidates.fasta", 
                                      seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Gyrli1_effectors<- seqinr::read.fasta(file = "Gyrli1_EffectorCandidates.fasta", 
                                       seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Cananz1_effectors<- seqinr::read.fasta(file = "Cananz1_EffectorCandidates.fasta", 
                                      seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Hyssto1_effectors<- seqinr::read.fasta(file = "Hyssto1_EffectorCandidates.fasta", 
                                       seqtype = "AA",as.string = TRUE, set.attributes = FALSE)



#get totals
#total number of putative effectors
Suivar1<- length(Suivar1_effectors)
Suitom1<- length(Suitom1_effectors)
Suisub1<- length(Suisub1_effectors)
Suisu1<- length(Suisu1_effectors)
Suipla1<- length(Suipla1_effectors)
Suipic1<- length(Suipic1_effectors)
Suipal1<- length(Suipal1_effectors)
Suiocc1<- length(Suiocc1_effectors)
Suilu4<- length(Suilu4_effectors)
Suilak1<- length(Suilak1_effectors)
Suihi1<- length(Suihi1_effectors)
Suigr1<- length(Suigr1_effectors)
Suidec1<- length(Suidec1_effectors)
Suicot1<- length(Suicot1_effectors)
Suicli1<- length(Suicli1_effectors)
Suibr2<- length(Suibr2_effectors)
Suibov1<- length(Suibov1_effectors)
Suiamp1<- length(Suiamp1_effectors)
Suiame1<- length(Suiame1_effectors)
Suifus1<- length(Suifus1_effectors)
#non-Suillus set
Rhivul1<- length(Rhivul1_effectors)
Rhitru1<- length(Rhitru1_effectors)
Amamu1<- length(Amamu1_effectors)
Hebcy2<- length(Hebcy2_effectors)
Lacbi2<- length(Lacbi2_effectors)
Paxin1<- length(Paxin1_effectors)
Pilcr1<- length(Pilcr1_effectors)
Pismi1<- length(Pismi1_effectors)
Sclci1<- length(Sclci1_effectors)
Theter1<- length(Theter1_effectors)
Thega1<- length(Thega1_effectors)
Rhivi1<- length(Rhivi1_effectors)
Rhives1<- length(Rhives1_effectors)
Rhisa1<- length(Rhisa1_effectors)
Ruscom1<- length(Ruscom1_effectors)
Rusbre1<- length(Rusbre1_effectors)
Pisti1<- length(Pisti1_effectors)
Lacam2<- length(Lacam2_effectors)
Hydru2<- length(Hydru2_effectors)
Gaumor1<- length(Gaumor1_effectors)
Gyrli1<- length(Gyrli1_effectors)
Cananz1<- length(Cananz1_effectors)
Hyssto1<- length(Hyssto1_effectors)




#touare here 
total_effectors<- data.frame(cbind(Suivar1, 
                                   Suitom1, 
                                   Suisub1, 
                                   Suisu1, 
                                   Suipla1, 
                                   Suipic1, 
                                   Suipal1, 
                                   Suiocc1, 
                                   Suilu4, 
                                   Suilak1, 
                                   Suihi1, 
                                   Suidec1, 
                                   Suicot1, 
                                   Suicli1, 
                                   Suibr2, 
                                   Suibov1, 
                                   Suiamp1, 
                                   Suiame1,
                                   Suifus1,
                                   Rhivul1,
                                   Rhitru1,
                                   Amamu1,
                                   Hebcy2,
                                   Lacbi2,
                                   Paxin1,
                                   Pilcr1,
                                   Pismi1,
                                   Sclci1, 
                                   Theter1,
                                   Thega1,
                                   Rhivi1,
                                   Rhives1,
                                   Rhisa1,
                                   Ruscom1,
                                   Rusbre1,
                                   Pisti1,
                                   Lacam2,
                                   Hydru2,
                                   Gaumor1,
                                   Gyrli1,
                                   Cananz1,
                                   Hyssto1),
                             row.names = "n_putative_effectors_from_EffectorP")



##SSP's per genome
Suivar1<- nrow(Suivar1_TMHMM_no_TMD_300)
Suitom1<- nrow(Suitom1_TMHMM_no_TMD_300)
Suisub1<- nrow(Suisub1_TMHMM_no_TMD_300)
Suisu1<- nrow(Suisu1_TMHMM_no_TMD_300)
Suipla1<- nrow(Suipla1_TMHMM_no_TMD_300)
Suipic1<- nrow(Suipic1_TMHMM_no_TMD_300)
Suipal1<- nrow(Suipal1_TMHMM_no_TMD_300)
Suiocc1<- nrow(Suiocc1_TMHMM_no_TMD_300)
Suilu4<- nrow(Suilu4_TMHMM_no_TMD_300)
Suilak1<- nrow(Suilak1_TMHMM_no_TMD_300)
Suihi1<- nrow(Suihi1_TMHMM_no_TMD_300)
Suigr1<- nrow(Suigr1_TMHMM_no_TMD_300)
Suidec1<- nrow(Suidec1_TMHMM_no_TMD_300)
Suicot1<- nrow(Suicot1_TMHMM_no_TMD_300)
Suicli1<- nrow(Suicli1_TMHMM_no_TMD_300)
Suibr2<- nrow(Suibr2_TMHMM_no_TMD_300)
Suibov1<- nrow(Suibov1_TMHMM_no_TMD_300)
Suiamp1<- nrow(Suiamp1_TMHMM_no_TMD_300)
Suiame1<- nrow(Suiame1_TMHMM_no_TMD_300)
Suifus1<- nrow(Suifus1_TMHMM_no_TMD_300)
#non-Suillus set
Rhivul1<- nrow(Rhivul1_TMHMM_no_TMD_300)
Rhitru1<- nrow(Rhitru1_TMHMM_no_TMD_300)
Amamu1<- nrow(Amamu1_TMHMM_no_TMD_300)
Hebcy2<- nrow(Hebcy2_TMHMM_no_TMD_300)
Lacbi2<- nrow(Lacbi2_TMHMM_no_TMD_300)
Paxin1<- nrow(Paxin1_TMHMM_no_TMD_300)
Pilcr1<- nrow(Pilcr1_TMHMM_no_TMD_300)
Pismi1<- nrow(Pismi1_TMHMM_no_TMD_300)
Sclci1<- nrow(Sclci1_TMHMM_no_TMD_300)
Theter1<- nrow(Theter1_TMHMM_no_TMD_300)
Thega1<- nrow(Thega1_TMHMM_no_TMD_300)
Rhivi1<- nrow(Rhivi1_TMHMM_no_TMD_300)
Rhives1<- nrow(Rhives1_TMHMM_no_TMD_300)
Rhisa1<- nrow(Rhisa1_TMHMM_no_TMD_300)
Ruscom1<- nrow(Ruscom1_TMHMM_no_TMD_300)
Rusbre1<- nrow(Rusbre1_TMHMM_no_TMD_300)
Pisti1<- nrow(Pisti1_TMHMM_no_TMD_300)
Lacam2<- nrow(Lacam2_TMHMM_no_TMD_300)
Hydru2<- nrow(Hydru2_TMHMM_no_TMD_300)
Gaumor1<- nrow(Gaumor1_TMHMM_no_TMD_300)
Gyrli1<- nrow(Gyrli1_TMHMM_no_TMD_300)
Cananz1<- nrow(Cananz1_TMHMM_no_TMD_300)
Hyssto1<- nrow(Hyssto1_TMHMM_no_TMD_300)

total_SSPs<- data.frame(cbind(Suivar1, 
                              Suitom1, 
                              Suisub1, 
                              Suisu1, 
                              Suipla1, 
                              Suipic1, 
                              Suipal1, 
                              Suiocc1, 
                              Suilu4, 
                              Suilak1, 
                              Suihi1, 
                              Suidec1, 
                              Suicot1, 
                              Suicli1, 
                              Suibr2, 
                              Suibov1, 
                              Suiamp1, 
                              Suiame1,
                              Suifus1,
                              Rhivul1,
                              Rhitru1,
                              Amamu1,
                              Hebcy2,
                              Lacbi2,
                              Paxin1,
                              Pilcr1,
                              Pismi1,
                              Sclci1, 
                              Theter1,
                              Thega1,
                              Rhivi1,
                              Rhives1,
                              Rhisa1,
                              Ruscom1,
                              Rusbre1,
                              Pisti1,
                              Lacam2,
                              Hydru2,
                              Gaumor1,
                              Gyrli1,
                              Cananz1,
                              Hyssto1),
                        row.names = "#SSPs_signalP,TMHMM,lt_300aa")


#Proteins per genome
Suivar1<- length(Suivar1_in)
Suitom1<- length(Suitom1_in)
Suisub1<- length(Suisub1_in)
Suisu1<- length(Suisu1_in)
Suipla1<- length(Suipla1_in)
Suipic1<- length(Suipic1_in)
Suipal1<- length(Suipal1_in)
Suiocc1<- length(Suiocc1_in)
Suilu4<- length(Suilu4_in)
Suilak1<- length(Suilak1_in)
Suihi1<- length(Suihi1_in)
Suigr1<- length(Suigr1_in)
Suidec1<- length(Suidec1_in)
Suicot1<- length(Suicot1_in)
Suicli1<- length(Suicli1_in)
Suibr2<- length(Suibr2_in)
Suibov1<- length(Suibov1_in)
Suiamp1<- length(Suiamp1_in)
Suiame1<- length(Suiame1_in)
Suifus1<- length(Suifus1_in)
#non-Suillus set
Rhivul1<- length(Rhivul1_in)
Rhitru1<- length(Rhitru1_in)
Amamu1<- length(Amamu1_in)
Hebcy2<- length(Hebcy2_in)
Lacbi2<- length(Lacbi2_in)
Paxin1<- length(Paxin1_in)
Pilcr1<- length(Pilcr1_in)
Pismi1<- length(Pismi1_in)
Sclci1<- length(Sclci1_in)
Theter1<- length(Theter1_in)
Thega1<- length(Thega1_in)
Rhivi1<- length(Rhivi1_in)
Rhives1<- length(Rhives1_in)
Rhisa1<- length(Rhisa1_in)
Ruscom1<- length(Ruscom1_in)
Rusbre1<- length(Rusbre1_in)
Pisti1<- length(Pisti1_in)
Lacam2<- length(Lacam2_in)
Hydru2<- length(Hydru2_in)
Gaumor1<- length(Gaumor1_in)
Gyrli1<- length(Gyrli1_in)
Cananz1<- length(Cananz1_in)
Hyssto1<- length(Hyssto1_in)

total_proteins<-data.frame(cbind(Suivar1, 
                                 Suitom1, 
                                 Suisub1, 
                                 Suisu1, 
                                 Suipla1, 
                                 Suipic1, 
                                 Suipal1, 
                                 Suiocc1, 
                                 Suilu4, 
                                 Suilak1, 
                                 Suihi1, 
                                 Suidec1, 
                                 Suicot1, 
                                 Suicli1, 
                                 Suibr2, 
                                 Suibov1, 
                                 Suiamp1, 
                                 Suiame1,
                                 Suifus1,
                                 Rhivul1,
                                 Rhitru1,
                                 Amamu1,
                                 Hebcy2,
                                 Lacbi2,
                                 Paxin1,
                                 Pilcr1,
                                 Pismi1,
                                 Sclci1, 
                                 Theter1,
                                 Thega1,
                                 Rhivi1,
                                 Rhives1,
                                 Rhisa1,
                                 Ruscom1,
                                 Rusbre1,
                                 Pisti1,
                                 Lacam2,
                                 Hydru2,
                                 Gaumor1,
                                 Gyrli1,
                                 Cananz1,
                                 Hyssto1),
                           row.names = "n_putative_proteins_from_gene_cat")  

#bind results and transform
results<- rbind(total_proteins[1,], total_SSPs[1,], total_effectors[1,])
results.1<- t(results)

#get % SSP's out of total genes
percent_SSPs_out_of_total_genes<- (results.1[,2] / results.1[,1])*100
percent_SSPs_out_of_total_genes<- signif(percent_SSPs_out_of_total_genes, digits = 3)
#get % of Effectors out of total SSP's 
percent_effectors_out_of_total_genes<- (results.1[,3] / results.1[,1])*100
percent_effectors_out_of_total_genes<- signif(percent_effectors_out_of_total_genes, digits = 2)
#get % of Effectors out of total genes
percent_effectors_out_of_SSPs<- (results.1[,3] / results.1[,2])*100
percent_effectors_out_of_SSPs<-signif(percent_effectors_out_of_SSPs, digits = 4)

#slap it together
results_with_percentages<- cbind(results.1, percent_SSPs_out_of_total_genes, percent_effectors_out_of_total_genes, percent_effectors_out_of_SSPs)

#write it up
write.csv(results, quote = FALSE, file = "SSP_and_Effector_totals.csv")

#read in genome size file and run stats for table
genome_size_df<- read.csv("genome_size_data.csv")


#split the two data categories for genome size and get stats for each
Suillus_genome_size<-data.frame(genome_size_df[ grep("Sui", genome_size_df[,1],), ])
Suillus_genome_size_mean<- mean(Suillus_genome_size$genome_size)
std <- function(x) sd(x)/sqrt(length(x))
Suillus_genome_size_SD<- sd(Suillus_genome_size$genome_size)
Suillus_genome_size_SE<- std(Suillus_genome_size$genome_size)  

Other_genome_size<-data.frame(genome_size_df[ grep("Sui", genome_size_df[,1], invert = TRUE), ])
Other_genome_size_mean<- mean(Other_genome_size$genome_size)

Other_genome_size_SD<- sd(Other_genome_size$genome_size)
Other_genome_size_SE<- std(Other_genome_size$genome_size)




#you are here:::

#test normality 
shapiro.test(Suillus_genome_size$genome_size)
#not normal
shapiro.test(Other_genome_size$genome_size)
#not normal


#test for equal variance 
var.test(Suillus_genome_size$genome_size, Other_genome_size$genome_size)
#variance is not sif. different, but can't trust becasue not nomal, be conservatice and use unequal var in test

#t-test 
t.test(Suillus_genome_size$genome_size, Other_genome_size$genome_size, var.equal = FALSE)



#split and get stats for n proteins 
results_with_percentages_for_split<- cbind(results_with_percentages, row.names(results_with_percentages))
Suillus_ssps<-data.frame(results_with_percentages_for_split[ grep("Sui", results_with_percentages_for_split[,7],), ])
Suillus_ssps<- Suillus_ssps[,1:6]
Other_ssps<- data.frame(results_with_percentages_for_split[ grep("Sui", results_with_percentages_for_split[,7], invert = TRUE), ])
Other_ssps<- Other_ssps[,1:6]

#need to convert these to numeric form aparently . . . 
Suillus_ssps[] <- lapply(Suillus_ssps, function(x) {
  if(is.character(x)) as.numeric(as.character(x)) else x
})
sapply(Suillus_ssps, class)

#and the other one too...
Other_ssps[] <- lapply(Other_ssps, function(x) {
  if(is.character(x)) as.numeric(as.character(x)) else x
})
sapply(Other_ssps, class)

Suillus_n_prot_mean<- mean(Suillus_ssps$n_putative_proteins_from_gene_cat)
Suillus_n_prot_SD<- sd(Suillus_ssps$n_putative_proteins_from_gene_cat)
Suillus_n_prot_SE<- std(Suillus_ssps$n_putative_proteins_from_gene_cat)

Other_n_prot_mean<- mean(Other_ssps$n_putative_proteins_from_gene_cat)
Other_n_prot_SD<- sd(Other_ssps$n_putative_proteins_from_gene_cat)
Other_n_prot_SE<- std(Other_ssps$n_putative_proteins_from_gene_cat)


#test normality 
shapiro.test(Suillus_ssps$n_putative_proteins_from_gene_ca)
#normal
shapiro.test(Other_ssps$n_putative_proteins_from_gene_ca)
#normal

#test for equal variance 
var.test(Suillus_ssps$n_putative_proteins_from_gene_cat, Other_ssps$n_putative_proteins_from_gene_cat)
#variance is not sif. different 

#t-test 
t.test(Suillus_ssps$n_putative_proteins_from_gene_cat, Other_ssps$n_putative_proteins_from_gene_cat, var.equal = TRUE)


#get stats for SSP's 
Suillus_ssp_mean<- mean(Suillus_ssps$X.SSPs_signalP.TMHMM.lt_300aa)
Suillus_ssp_SD<- sd(Suillus_ssps$X.SSPs_signalP.TMHMM.lt_300aa)
Suillus_ssp_SE<- std(Suillus_ssps$X.SSPs_signalP.TMHMM.lt_300aa)

Other_ssp_mean<- mean(Other_ssps$X.SSPs_signalP.TMHMM.lt_300aa)
Other_ssp_SD<- sd(Other_ssps$X.SSPs_signalP.TMHMM.lt_300aa)
Other_ssp_SE<- std(Other_ssps$X.SSPs_signalP.TMHMM.lt_300aa)

#test normality 
shapiro.test(Suillus_ssps$X.SSPs_signalP.TMHMM.lt_300aa)
#normal-ish
shapiro.test(Other_ssps$X.SSPs_signalP.TMHMM.lt_300aa)
#normal

#test for equal variance 
var.test(Suillus_ssps$X.SSPs_signalP.TMHMM.lt_300aa, Other_ssps$X.SSPs_signalP.TMHMM.lt_300aa)
#variance is not sif. different 

#t-test 
t.test(Suillus_ssps$X.SSPs_signalP.TMHMM.lt_300aa, Other_ssps$X.SSPs_signalP.TMHMM.lt_300aa, var.equal = TRUE)



#get stats for effectors
Suillus_effectors_mean<- mean(Suillus_ssps$n_putative_effectors_from_EffectorP)
Suillus_effectors_SD<- sd(Suillus_ssps$n_putative_effectors_from_EffectorP)
Suillus_effectors_SE<- std(Suillus_ssps$n_putative_effectors_from_EffectorP)

Other_effectors_mean<- mean(Other_ssps$n_putative_effectors_from_EffectorP)
Other_effectors_SD<- sd(Other_ssps$n_putative_effectors_from_EffectorP)
Other_effectors_SE<- std(Other_ssps$n_putative_effectors_from_EffectorP)

#test normality 
shapiro.test(Suillus_ssps$n_putative_effectors_from_EffectorP)
#not normal
shapiro.test(Other_ssps$n_putative_effectors_from_EffectorP)
#not normal

#test for equal variance 
var.test(Suillus_ssps$n_putative_effectors_from_EffectorP, Other_ssps$n_putative_effectors_from_EffectorP)
#variance is not sif. different, but we can't trust it becaseu it's not normal, use conservatie var = FALSE in test

#t-test 
t.test(Suillus_ssps$n_putative_effectors_from_EffectorP, Other_ssps$n_putative_effectors_from_EffectorP, var.equal = FALSE)
#not different


#parse the SSP's from the SSSP's
#to do this, upload the outout fasta files of SSP's and run them in orthofinder. 
ortho_finder_SSP_results<- read.table("Statistics_PerSpecies_all_ECM.csv", sep = ",", row.names = 1, header=TRUE)

#the perninant info: 
SSP_vs_SSSP1<- ortho_finder_SSP_results[1:3,]
SSP_vs_SSSP<- t(SSP_vs_SSSP1)


#get percent SSSP's of SSP's
percent_SSPs_out_of_SSSPs<- (100*(SSP_vs_SSSP[,3]) / SSP_vs_SSSP[,1])

percent_SSPs_out_of_SSSPs<- data.frame(signif(percent_SSPs_out_of_SSSPs, digits = 3))




#######START HERE WITH GRAPHS######## 
######### % SSP's out of total proteins 
med1= mean(SSPs_out_of_prot_df$`%` [SSPs_out_of_prot_df$group == "a_Suillus"])
med2= mean(SSPs_out_of_prot_df$`%` [SSPs_out_of_prot_df$group == "b_Other"])
par(mar = c(6.5, 8.5, 3, 3.5), mgp = c(6, 2.5, 0))
stripchart(SSPs_out_of_prot_df$`%` ~ SSPs_out_of_prot_df$group,
           vertical = TRUE,
           method = "jitter", jitter = 0.2, 
           pch = 16, 
           col = c("#F0502B",  "#F79552"),
           bg = rep(c("#F0502B", "#F79552")),
           cex.axis = 0.7,
           ylim=c(0,4), 
           ylab = "% SSPs out of all proteins", 
           axes = FALSE, 
           cex = 1.3)
box()
axis(2)
mtext(text = c("Suillus", "Other ECM"),side=1,at=c(1,2),line = 1, font = 3)
segments(x0 = .7, y0 =  med1, x1 = 1.3, y1=med1, lwd = 2, col = "black" )
segments(x0 = 1.7, y0 =  med2, x1 = 2.3, y1=med2, lwd = 2, col = "black" )
#add boxplot around the data? 
#boxplot(`%` ~ group, data = SSPs_out_of_SSSPs_df, add=TRUE, range=0, whisklty = 0, staplelty = 0)




#######SSSP's as a percentage of SSP's########
percent_SSSPs_out_of_SSPs_Suillus<- (100*SSP_vs_SSSP_raw_Suillus[,3] / SSP_vs_SSSP_raw_Suillus[,1])
percent_SSSPs_out_of_SSPs_Suillus<- data.frame(signif(percent_SSSPs_out_of_SSPs_Suillus, digits = 3))

percent_SSSPs_out_of_SSPs_Other<- (100*SSP_vs_SSSP_raw_Other[,3] / SSP_vs_SSSP_raw_Other[,1])
percent_SSSPs_out_of_SSPs_Other<- data.frame(signif(percent_SSSPs_out_of_SSPs_Other, digits = 3))

Suillus2<- cbind(percent_SSSPs_out_of_SSPs_Suillus, rep("a_Suillus", length(percent_SSSPs_out_of_SSPs_Suillus)))
colnames(Suillus2)<- c("%", "group")

Other_ECM2<- cbind(percent_SSSPs_out_of_SSPs_Other, rep("b_Other", length(percent_SSSPs_out_of_SSPs_Other)))
colnames(Other_ECM2)<- c("%", "group")

SSPs_out_of_SSSPs_df<- rbind(Suillus2, Other_ECM2)

med3= mean(SSPs_out_of_SSSPs_df$`%` [SSPs_out_of_SSSPs_df$group == "a_Suillus"])
med4= mean(SSPs_out_of_SSSPs_df$`%` [SSPs_out_of_SSSPs_df$group == "b_Other"])
par(mar = c(6.5, 8.5, 3, 3.5), mgp = c(6, 2.5, 0))
stripchart(SSPs_out_of_SSSPs_df$`%` ~ SSPs_out_of_SSSPs_df$group,
           vertical = TRUE,
           method = "jitter", jitter = 0.2, 
           pch = 19, 
           col = c("#F0502B",  "#F79552"),
           bg = rep(c("#F0502B", "#F79552")),
           cex.axis = 0.7,
           ylim=c(0,70), 
           ylab = "% SSSPs out of SSPs", 
           axes = FALSE, 
           cex = 1.3)
box()
axis(2)
mtext(text = c("Suillus", "Other ECM"),side=1,at=c(1,2),line = 1, font = 3)
segments(x0 = .7, y0 =  med3, x1 = 1.3, y1=med3, lwd = 2, col = "black" )
segments(x0 = 1.7, y0 =  med4, x1 = 2.3, y1=med4, lwd = 2, col = "black" )



#######Effectors as a percentage of SSP's########

percent_effector_Suillus<- as.numeric(Suillus_ssps$percent_effectors_out_of_SSPs)
percent_effector_Other<- as.numeric(Other_ssps$percent_effectors_out_of_SSPs)

Suillus3<- cbind(percent_effector_Suillus, rep("a_Suillus", length(percent_effector_Suillus)))
colnames(Suillus3)<- c("%", "group")

Other_ECM3<- cbind(percent_effector_Other, rep("b_Other", length(percent_effector_Other)))
colnames(Other_ECM3)<- c("%", "group")

effectors_out_of_SSSPs_df<- data.frame(rbind(Suillus3, Other_ECM3))


med3= mean(Suillus_ssps$percent_effectors_out_of_SSPs)
med4= mean(Other_ssps$percent_effectors_out_of_SSPs)
par(mar = c(6.5, 8.5, 3, 3.5), mgp = c(6, 2.5, 0))
stripchart(effectors_out_of_SSSPs_df$X. ~ effectors_out_of_SSSPs_df$group,
           vertical = TRUE,
           method = "jitter", jitter = 0.2, 
           pch = 19, 
           col = c("#F0502B",  "#F79552"),
           bg = rep(c("#F0502B", "#F79552")),
           cex.axis = 0.7,
           ylim=c(0,50), 
           ylab = "% Effectors out of SSPs", 
           axes = FALSE, 
           cex = 1.3)
box()
axis(2)
mtext(text = c("Suillus", "Other ECM"),side=1,at=c(1,2),line = 1, font = 3)
segments(x0 = .7, y0 =  med3, x1 = 1.3, y1=med3, lwd = 2, col = "black" )
segments(x0 = 1.7, y0 =  med4, x1 = 2.3, y1=med4, lwd = 2, col = "black" )




##########same as above but with raw numbers rather than percentages############

###### form Raw data (no %): SSP's #######
Suillus_ssps[,2] #total Suillus SSPs
Other_ssps[,2] #total Other ECM SSPs

SSPs_Suillus<- as.numeric(Suillus_ssps[,2])
SSPs_Suillus_df<- data.frame(cbind(Suillus_ssps[,2], rep("a_Suillus", length(Suillus_ssps[,2]))))

SSPs_Other<- as.numeric(Other_ssps[,2])
SSPs_Other_df<- data.frame(cbind(Other_ssps[,2], rep("b_Other", length(Other_ssps[,2]))))

SSP_df<- data.frame(rbind(SSPs_Suillus_df, SSPs_Other_df))

med3= mean(SSPs_Suillus)
med4= mean(SSPs_Other)

par(mar = c(6.5, 8.5, 3, 3.5), mgp = c(6, 2.5, 0))
stripchart(SSP_df$X1 ~ SSP_df$X2,
           vertical = TRUE,
           method = "jitter", jitter = 0.2, 
           pch = 16, 
           col = c("#F0502B",  "#F79552"),
           bg = rep(c("#F0502B", "#F79552")),
           cex.axis = 0.7,
           ylim=c(0,550), 
           ylab = "SSPs", 
           axes = FALSE, 
           cex = 1.3)
box()
axis(2)
mtext(text = c("Suillus", "Other ECM"),side=1,at=c(1,2),line = 1, font = 3)
segments(x0 = .7, y0 =  med3, x1 = 1.3, y1=med3, lwd = 2, col = "black" )
segments(x0 = 1.7, y0 =  med4, x1 = 2.3, y1=med4, lwd = 2, col = "black" )



###### form Raw data (no %): SSSP's #######
SSSPs_Suillus<- as.numeric(SSP_vs_SSSP_raw_Suillus[,3])
SSSPs_Suillus_df<- data.frame(cbind(SSP_vs_SSSP_raw_Suillus[,3], rep("a_Suillus", length(SSP_vs_SSSP_raw_Suillus[,3]))))

SSSPs_Other<- as.numeric(SSP_vs_SSSP_raw_Other[,3])
SSSPs_Other_df<- data.frame(cbind(SSP_vs_SSSP_raw_Other[,3], rep("b_Other", length(SSP_vs_SSSP_raw_Other[,3]))))

SSSP_df<- data.frame(rbind(SSSPs_Suillus_df, SSSPs_Other_df))

med3= mean(SSSPs_Suillus)
med4= mean(SSSPs_Other)

par(mar = c(6.5, 8.5, 3, 3.5), mgp = c(6, 2.5, 0))
stripchart(SSSP_df$X1 ~ SSSP_df$X2,
           vertical = TRUE,
           method = "jitter", jitter = 0.2, 
           pch = 16, 
           col = c("#F0502B",  "#F79552"),
           bg = rep(c("#F0502B", "#F79552")),
           cex.axis = 0.7,
           ylim=c(0,360), 
           ylab = "SSSPs", 
           axes = FALSE, 
           cex = 1.3)
box()
axis(2)
mtext(text = c("Suillus", "Other ECM"),side=1,at=c(1,2),line = 1, font = 3)
segments(x0 = .7, y0 =  med3, x1 = 1.3, y1=med3, lwd = 2, col = "black" )
segments(x0 = 1.7, y0 =  med4, x1 = 2.3, y1=med4, lwd = 2, col = "black" )

###### form Raw data (no %): Effectors #######
effector_Suillus<- as.numeric(Suillus_ssps$n_putative_effectors_from_EffectorP)
effector_Suillus_df<- data.frame(cbind(Suillus_ssps$n_putative_effectors_from_EffectorP, rep("a_Suillus", length(Suillus_ssps$n_putative_effectors_from_EffectorP))))

effector_Other<- as.numeric(Other_ssps$n_putative_effectors_from_EffectorP)
effector_Other_df<- data.frame(cbind(Other_ssps$n_putative_effectors_from_EffectorP, rep("b_Other", length(Other_ssps$n_putative_effectors_from_EffectorP))))

effector_df<- data.frame(rbind(effector_Suillus_df, effector_Other_df))

med3<- mean(Suillus_ssps$n_putative_effectors_from_EffectorP)
med4<- mean(Other_ssps$n_putative_effectors_from_EffectorP)

par(mar = c(6.5, 8.5, 3, 3.5), mgp = c(6, 2.5, 0))
stripchart(effector_df$X1 ~ effector_df$X2,
           vertical = TRUE,
           method = "jitter", jitter = 0.2, 
           pch = 16, 
           col = c("#F0502B",  "#F79552"),
           bg = rep(c("#F0502B", "#F79552")),
           cex.axis = 0.7,
           ylim=c(0,280), 
           ylab = "Effectors", 
           axes = FALSE, 
           cex = 1.3)
box()
axis(2)
mtext(text = c("Suillus", "Other ECM"),side=1,at=c(1,2),line = 1, font = 3)
segments(x0 = .7, y0 =  med3, x1 = 1.3, y1=med3, lwd = 2, col = "black" )
segments(x0 = 1.7, y0 =  med4, x1 = 2.3, y1=med4, lwd = 2, col = "black" )