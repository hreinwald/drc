#' @title Convert parameter vectors to matrices
#' @keywords internal
"drmConvertParm" <- 
function(startVec, startMat, factor1, colList)
{
    startMat2 <- startMat
    if (length(unique(factor1)) == 1) {return(startVec)}
    
    mmat <- model.matrix(~factor(factor1) - 1)
        
    pm <- list()
    for (i in 1:length(colList))
    {
        clElt <- colList[[i]]
        ncclElt <- dim(clElt)[2]   
            
        indVec <- !is.na(startMat2[, i, drop = FALSE])
        indVal <- min(c(sum(indVec), dim(clElt)[2]))
        
        indVec2 <- (1:ncclElt)[indVec]
        if (length(indVec2) > ncclElt) {indVec2 <- 1:ncclElt}

        pm[[i]] <- (ginv(t(clElt)%*%clElt)%*%t(clElt))[1:indVal, , drop = FALSE]%*%mmat[, indVec]%*%startMat2[indVec, i, drop = FALSE]
    }  
    tempVec <- unlist(pm)
    tempVec <- tempVec[!is.na(tempVec)]
    
    ## Checking whether the intercept column has been removed
    indVec3 <- ( abs(tempVec)<1e-10 )
    if (any(indVec3)) 
    {
        tempVec <- startVec
    }

    return(tempVec)
}
