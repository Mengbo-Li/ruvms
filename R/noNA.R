noNA <- function(Yc) {

   if (sum(is.na(Yc)) == 0) {    return(TRUE)    } else {    return(FALSE)    }

}
