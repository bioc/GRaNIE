#' @param communities \code{NULL} or numeric vector or character vector. Default \code{NULL}. 
#' If set to \code{NULL}, all community enrichments that have been calculated before are plotted. 
#' If a numeric vector is specified (when \code{selection = "byRank"}), the rank of the communities is specified.
#' For example, \code{communities = c(1,4)} then denotes the first and fourth largest community.
#' If a character vector is specified (when \code{selection = "byLabel"}), the name of the communities is specified instead.
#' For example, \code{communities = c("1","4")} then denotes the communities with the names "1" and "4", which may or may not be the largest and fourth largest communities among all.
