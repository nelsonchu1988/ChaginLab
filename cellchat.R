library(CellChat)
library(ggplot2)                  
library(patchwork)
library(igraph)
library(ComplexHeatmap)
options(stringsAsFactors = FALSE)

CellC <- function(cellchat){
  CellChatDB <- CellChatDB.mouse
  cellchat@DB <- CellChatDB.use
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  future::plan("multiprocess", workers = 4)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat)
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  selectK(cellchat, pattern = "outgoing")
  nPatterns = 5
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
  ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
  ht1
  cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
  netAnalysis_river(cellchat, pattern = "outgoing")
  dev.off()
  selectK(cellchat, pattern = "incoming")
  nPatterns = 5
  cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
  netAnalysis_river(cellchat, pattern = "incoming")
  ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
  ht2
  return (cellchat)
}

CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

CellChatDB <- CellChatDB.mouse # set CellChatDB <- CellChatDB.human if working on the human dataset

library(NMF)
library(ggalluvial)

data.input <- GetAssayData(tmpGH, assay = "RNA", slot = "data") # normalized data matrix

levels(tmp)
new.cluster.ids <- c("epSC/5", "epSC/1", "CC/4", "CC/6", "CC/2", "PMC/3",
                     "PMC/7", "HC/8", "PHC/0")
names(new.cluster.ids) <- levels(tmp)
tmp <- RenameIdents(tmp, new.cluster.ids)

labels <- Idents(tmpGH)

meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels

cellchat.GH <- createCellChat(object = data.input, meta = meta, group.by = "group")
#> Create a CellChat object from a data matrix


cellchat.GH <- CellC(cellchat.GH)
