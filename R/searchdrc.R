"searchdrc" <- function(object, which, range, len = 50)
{
   sv <- object$start
 
   parNames <- object$parNames[[2]]
   whichInd <- regexpr(paste("^", which, ":", sep = ""), parNames)
   whichInd <- ((1:length(parNames))[whichInd>0])[1]
   
   if (length(whichInd)<1) {stop(paste("No such parameter ", which, sep = ""))}

   found <- FALSE
   oldWarn <- getOption("warn")
   on.exit(options(warn = oldWarn), add = TRUE)
   for (i in seq(range[1], range[2], length.out = len))
   {
       sv[whichInd] <- i 

       options(warn = -1)
       modelFit <- try(update(object, start = sv, control = drmc(noMessage = TRUE)), silent=TRUE)
       options(warn = oldWarn)
       if (!inherits(modelFit, "try-error")) {found <- TRUE; break}
   }
   
   if (found) 
   {
       return(modelFit)
   } else {
       warning("Convergence failed.", call. = FALSE)
   }
}
