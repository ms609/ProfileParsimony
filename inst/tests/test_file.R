path = "C:/work/fun/implied_weight/goloboff/clM_bin33exp/"; repl=4;
require(devtools)
install(tsPath <- "C:/Research/R/TreeSearch")
install(profPath <- "C:/work/information/ProfileParsimony")
require(phangorn)
require(TreeSearch)
require(ProfileParsimony)
PP_SUBOPTIMAL_VALUES <- c(1e-08, round(0.73^((19:0) - 10), 5))
ROOT <- paste0(path, '../../') ## DELETE once we can call library('SlowQuartet') in tree_statistics
source(paste0(ROOT, 'R/treegen_functions.R'))
source(paste0(ROOT, 'R/tree_statistics.R'))
data(referenceTree)

  fileRoot <- paste0(path, 'PP/R', repl, '/', repl)
  matrixFile <- paste0(fileRoot, '_PP.nex')
  if (!file.exists(matrixFile)) stop("\n File not found: ", matrixFile)
  nexusData <- read.nexus.data(matrixFile)
  tokens <- unique(unlist(nexusData))
  if (sum(as.character(0:9) %in% tokens) != 2) stop("Currently, Profile Parsimony requires exactly two tokens per TS")
  dataset <- PhyDat(read.nexus.data(matrixFile), 0:9)
  
  ready2 <- PrepareDataProfile(dataset, precision=12000)
  ready3 <- PrepareDataProfile(dataset, precision=4e+05)
  ready4 <- PrepareDataProfile(dataset, precision=8e+05)
  ready5 <- PrepareDataProfile(dataset, precision=1.6e+06)
  
info2 <- attr(ready2, 'info.amounts') # Quick, but far from brilliant I must confess!
info3 <- attr(ready3, 'info.amounts')
info4 <- attr(ready4, 'info.amounts')
info5 <- attr(ready5, 'info.amounts')


diff32 <- as.double(info3 - info2)
diff42 <- as.double(info4 - info2)
diff43 <- as.double(info4 - info3)
diff54 <- as.double(info5 - info4)
nonzero <- info4 > 0.00001

hist (diff32)
hist (diff43)
hist (thisDiff <- diff54); quantile(thisDiff, probs=c(0, 5, 10, 50, 90, 95, 100)/100)
hist (diff42)
hist(100*(diff32 / info4)[nonzero])
hist(100*(diff42 / info4)[nonzero])
hist(100*(diff43 / info4)[nonzero])




hist(info1 - info2)
hist(info3 - info2)
hist(info3 - info1)
hist(info4 - info2)

  sillyData <- lapply(1:22, function (i) c( rep(0, i - 1), rep(1, 22 - i), rep(1, 22 - i), rep(0, i - 1)))#, sample(2, 20, replace=TRUE)-1))
  names(sillyData) <- as.character(1:22)
  dataset <- PhyDat(sillyData, 0:1)
  readyData <- PrepareDataProfile(dataset, 12000)
  
  rTree <- randomTree <- RandomTree(dataset, '1')
  FitchScore(rTree, dataset, TipsAreNames)
  FitchScore(rTree, readyData, TipsAreColumns)
  FitchScore(referenceTree, dataset, TipsAreNames)
  FitchScore(referenceTree, readyData, TipsAreColumns)
  
  ProfileScore(rTree, readyData)
  ProfileScore(referenceTree, readyData)


  quickTS <- TreeSearch(rTree, readyData, TreeScorer = ProfileScore, Rearrange=RootedNNI, maxIter=1000, maxHits=40) # then without the Do
  quickTS <- TreeSearch(rTree, dataset, TreeScorer = FitchScore  , Rearrange=RootedNNI, maxIter=1000, maxHits=40) # then without the Do
  
  quickFitch <- Ratchet(rTree, dataset, TreeScorer = FitchScore, suboptimal=max(PP_SUBOPTIMAL_VALUES),
                        ratchHits=6, searchHits=25, searchIter=500, ratchIter=500,
                        rearrangements='TBR')
                   
  quick <- Ratchet(rTree, readyData, TreeScorer = ProfileScore, returnAll = FALSE, rooted=TRUE,
                   ratchHits=5, searchHits=30, searchIter=100, ratchIter=50,
                   rearrangements='TBR')
    
    
    
    
   C_Fitch(characters = as.integer(data[whichChar<-1, tipLabel]), 
   nChar = length(whichChar), parent, child, nEdge = length(parent), weight = charWeights, maxNode = max(parent), nTip = length(tipLabel))
    