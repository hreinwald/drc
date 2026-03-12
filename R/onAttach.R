#' @title Package attach hook
#' @keywords internal
.onAttach <- function(libname, pkgname)
{
    packageStartupMessage("\n===============================================")
    packageStartupMessage("  'drc' has been loaded")
    packageStartupMessage("  Analysis of Dose-Response Data")
    packageStartupMessage("  Version 3.3.0")
    packageStartupMessage("===============================================\n")

    packageStartupMessage("Developers:")
    packageStartupMessage("  - Christian Ritz (ritz@bioassay.dk)")
    packageStartupMessage("  - Jens C. Streibig (streibig@bioassay.dk)")
    packageStartupMessage("  - Hannes Reinwald (hannes.reinwald@bayer.com)\n")

    packageStartupMessage("Please cite 'drc' if used for a publication:\n")
    packageStartupMessage("  Ritz, C., Jensen, S. M., Gerhard, D., Streibig, J. C. (2019)")
    packageStartupMessage("  Dose-Response Analysis Using R. CRC Press\n")
    packageStartupMessage("Additional references:")
    packageStartupMessage("  - Ritz, C., et al. (2015). Dose-Response Analysis Using R.")
    packageStartupMessage("    PLOS ONE, 10(12), e0146021.")
    packageStartupMessage("  - Ritz, C. and Streibig, J. C. (2005). Bioassay Analysis using R.")
    packageStartupMessage("    Journal of Statistical Software, 12(5), 1-22.\n")

    packageStartupMessage("For citation formats, type: citation('drc')")
    packageStartupMessage("For R citation, type: citation()\n")

    packageStartupMessage("Bug reports and issues:")
    packageStartupMessage("  https://github.com/hreinwald/drc/issues/\n")
    packageStartupMessage("===============================================\n")
}
