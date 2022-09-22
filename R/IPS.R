#' @title IPS
#' @description
#' A function is used to compute the importance pair score (IPS) of the network returned by the TarNet() function.
#'
#' @param TarnetPairs The result returned by TarNet() function.
#' @param IIS The result returned by IIS() function.
#' @return IPS_res A data frame including the importance pair score (IPS) of the target network.
#'
#' @export IPS
#' @examples
#' IIScore <- IIS(microApath = micro.eg, metaApath = metabo.eg,
#'                conf = confounder.eg, groupInfo = groupInfo.eg)
#' MMfunc_res <- MMfunc(IIS = IIScore)
#' TarNet_res <- TarNet(IIS = IIScore, mm2path = MMfunc_res)
#' IPScore <- IPS(TarnetPairs = TarNet_res,IIS = IIScore)

IPS<-function(TarnetPairs,IIS){
  IISdf1 = IIS[["Node_Score"]][,1:2]
  colnames(IISdf1) <- c("to","to_IIS")
  IISdf2 = IIS[["Node_Score"]][,1:2]
  colnames(IISdf2) <- c("from","from_IIS")
  IPS_res <- merge(TarnetPairs,IISdf1,all.x = T)
  IPS_res <- merge(IPS_res,IISdf2,all.x = T)
  IPS_res$r_score <- abs(IPS_res$r)
  IPS_res$r_score <-  seq(nrow(IPS_res),1,-1)
  IPS_res$r_score <- (IPS_res$r_score - 1)/(nrow(IPS_res) - 1) * 2 + 1
  IPS_res$ips <- IPS_res$to_IIS + IPS_res$from_IIS + IPS_res$r_score
  IPS_res <- IPS_res[order(IPS_res$ips,decreasing = T),]
  return(IPS_res)
}
