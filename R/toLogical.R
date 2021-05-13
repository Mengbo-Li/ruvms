toLogical <- function(ctl, n) {

   tmp <- rep(FALSE, n)
   tmp[ctl] <- TRUE
   return(tmp)

}
