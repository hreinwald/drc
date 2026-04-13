#' @title Package attach hook
#' @keywords internal
.onAttach <- function(libname, pkgname)
{
    packageStartupMessage(paste0(
        "\n===============================================\n",
        "  'drc' has been loaded\n",
        "  Analysis of Dose-Response Data\n",
        "  Version 3.3.1\n",
        "===============================================\n\n",
        "Developers:\n",
        "  - Christian Ritz (ritz@bioassay.dk)\n",
        "  - Jens C. Streibig (streibig@bioassay.dk)\n",
        "  - Hannes Reinwald (hannes.reinwald@bayer.com)\n\n",
        "Please cite 'drc' if used for a publication:\n\n",
        "  Ritz, C., Jensen, S. M., Gerhard, D., Streibig, J. C. (2019)\n",
        "  Dose-Response Analysis Using R. CRC Press\n\n",
        "Additional references:\n",
        "  - Ritz, C., et al. (2015). Dose-Response Analysis Using R.\n",
        "    PLOS ONE, 10(12), e0146021.\n",
        "  - Ritz, C. and Streibig, J. C. (2005). Bioassay Analysis using R.\n",
        "    Journal of Statistical Software, 12(5), 1-22.\n\n",
        "For citation formats, type: citation('drc')\n",
        "For R citation, type: citation()\n\n",
        "Bug reports and issues:\n",
        "  https://github.com/hreinwald/drc/issues/\n",
        "===============================================\n"
    ))
}