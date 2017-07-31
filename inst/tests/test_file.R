
path = "C:/work/fun/implied_weight/goloboff/clM_bin33exp/"; repl=4;
require('ProfileParsimony')
require('TreeSearch')
PP_SUBOPTIMAL_VALUES <- c(1e-08, round(0.73^((19:0) - 10), 5))
ROOT <- paste0(path, '../../') ## DELETE once we can call library('SlowQuartet') in tree_statistics
source(paste0(ROOT, 'R/treegen_functions.R'))
source(paste0(ROOT, 'R/tree_statistics.R'))
source(paste0(ROOT, '../../Information/ProfileParsimony/R/info_extra_step.R'))
source(paste0(ROOT, '../../Information/ProfileParsimony/R/ProfileScore.R'))
source(paste0(ROOT, '../../Information/ProfileParsimony/R/data_manipulation.R'))
data(referenceTree)
source(paste0(ROOT, '../../Information/ProfileParsimony/R/info_extra_step.R'))
source(paste0(ROOT, '../../Information/ProfileParsimony/R/ProfileScore.R'))
source(paste0(ROOT, '../../Information/ProfileParsimony/R/data_manipulation.R'))


  sillyData <- lapply(1:22, function (i) c( rep(0, i - 1), rep(1, 22 - i), rep(1, 22 - i), rep(0, i - 1)))#, sample(2, 20, replace=TRUE)-1))
  names(sillyData) <- as.character(1:22)
  dataset <- PhyDat(sillyData, 0:1)
  readyData <- PrepareDataProfile(dataset, 12000)
  
  rTree <- RandomTree(dataset)
  FitchScore(rTree, dataset, TipsAreNames)
  FitchScore(rTree, readyData, TipsAreColumns)
  FitchScore(referenceTree, dataset, TipsAreNames)
  FitchScore(referenceTree, readyData, TipsAreColumns)
  
  ProfileScore(RandomTree(dataset), readyData)
  ProfileScore(referenceTree, readyData)
  
  quickFitch <- Ratchet(RandomTree(dataset), dataset, TreeScorer = FitchScore,
                        keepAll=TRUE, outgroup='1', suboptimal=max(PP_SUBOPTIMAL_VALUES),
                        ratchHits=6, searchHits=25, searchIter=500, ratchIter=500,
                        rearrangements='TBR')
                   
  quick <- Ratchet(RandomTree(dataset), readyData, TreeScorer = ProfileScore,
                   keepAll=TRUE, outgroup='1', suboptimal=max(PP_SUBOPTIMAL_VALUES),
                   ratchHits=9, searchHits=90, searchIter=5000, ratchIter=500,
                   rearrangements='TBR')
    