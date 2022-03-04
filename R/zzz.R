#' @importFrom utils packageVersion
.onAttach <- function(libname, pkgname) {
    
    message = paste0("\n",
                     
                     "Welcome to the ", utils::packageName(), " package and thank you for using our software. This is ", utils::packageName(), " version ", utils::packageVersion(utils::packageName()),".\n",
                     
                     "All project-related information, guidance, help, documentation and ways to contact us can be found in the R help or here:\nhttps://grp-zaugg.embl-community.io/GRaNIE\n"
                     )
    
    packageStartupMessage(message)
    
    # Turn off scientific notation
    options(scipen = 999) 

}
