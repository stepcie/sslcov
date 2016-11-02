#'Cumulative indicator computation ?
#'
#'@return a vector of length(yy)
#'
#'@keywords internal
#'
#'
sum.I=function(yy,FUN,Yi,Vi=NULL){

#   if(FUN == "<" | FUN == "<="){
#     yy=-yy
#     Yi=-Yi
#   }
#   if(substring(FUN,2,2)=="="){
#     pos=rank(c(yy,Yi-1e-8))[1:length(yy)] - rank(yy)
#   }else{
    pos <- rank(c(yy,Yi+1e-8))[1:length(yy)] - rank(yy)
  # }

   if(!is.null(Vi)){
     Vi <- Vi[order(Yi)]
     res <- (c(0,cumsum(Vi))[pos+1])
   }else{
    res <- pos
   }

  return(res)
}
