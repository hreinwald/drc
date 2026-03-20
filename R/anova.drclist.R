#' @title ANOVA for list of drc objects
#' @keywords internal
"anova.drclist" <- function(object, ...,  details = TRUE, test = NULL)
{
    objects <- list(object, ...)
    if (length(objects) > 2) {stop("Only two models can be compared")}

    ## Testing two models against each other
    obj1 <- objects[[1]]
    obj2 <- objects[[2]]
    rowNames <- c("1st model", "2nd model")
 
    if ( !(obj1$"type"==obj2$"type") ) {stop("The two models are based on different types on data")}
    if (obj1$"type" == "binomial" && (is.null(test)) ) {test <- "Chisq"}
    if (obj1$"type" == "continuous"  && (is.null(test)) ) {test <- "F"}

    if (!(test == "F"))  # chis-square based test
    {
        loglik <- c(logLik(obj1), logLik(obj2))  # c(sumObj1[[4]][1], sumObj2[[4]][1])
        dfModel <- c(attr(logLik(obj1), "df"),  attr(logLik(obj2), "df"))  # c(sumObj1[[4]][2], sumObj2[[4]][2])
        testStat <- (2*abs(loglik[1] - loglik[2]))
        dfDiff <- c(NA, abs(dfModel[1] - dfModel[2]))

        pVal <- c(NA, 1 - pchisq(testStat, dfDiff[2]))
        testStat <- c(NA, testStat)

        headName <- "ANOVA-like table\n"
        colNames <- c("ModelDf", "Loglik", "Df", "LR value", "p value")
        
    } else {  # F-test       
        df1 <- df.residual(obj1)
        df2 <- df.residual(obj2)
        if (df2 > df1) 
        {
            objTemp <- obj1
            obj1 <- obj2
            obj2 <- objTemp
            df1 <- df.residual(obj1)
            df2 <- df.residual(obj2)
            
            rowNames <- rowNames[c(2, 1)]
        }
        loglik <- c(obj1$"summary"[4], obj2$"summary"[4])
        dfModel <- c(df1, df2)

        dfDiff <- c(NA, df1 - df2)

        testStat <- ((loglik[1] - loglik[2]) / dfDiff[2]) / (loglik[2] / df2)
        # Handle edge cases for the F-test p-value
        if (is.nan(testStat)) {
            # Undefined test statistic (e.g. 0/0), test is meaningless
            pVal <- c(NA, NA)
        } else if (testStat < 0) {
            # Negative F (including -Inf): the simpler model fits better,
            # so there is no evidence for the more complex model
            pVal <- c(NA, 1)
        } else if (is.infinite(testStat)) {
            # +Inf: test is undefined (e.g. models have equal df)
            pVal <- c(NA, NA)
        } else {
            pVal <- c(NA, 1 - pf(testStat, dfDiff[2], df2))
        }
        testStat <- c(NA, testStat)

        headName <- "ANOVA table\n"
        colNames <- c("ModelDf", "RSS", "Df", "F value", "p value")
    }
    if (details)
    {
        ## Specifying the models
        cat("\n")

        collapse1 <- obj1[[8]]$collapse
        if (is.null(collapse1)) {collapse1 <- obj1[[8]]$pmodels}
        if (!is.null(obj1$"pmodelsText")) {collapse1 <- obj1$"pmodelsText"}
        if (is.null(collapse1)) 
        {
            if (is.null(obj1[[8]]$curve))
            {
                collapse1 <- "1 (for all parameters)"
            } else {
                collapse1 <- paste(deparse(obj1[[8]]$curve), "(for all parameters)")
            }
        } else {
            if (!is.character(collapse1)) 
            {
                collapse1 <- paste(deparse(collapse1), collapse = "")
            }
            collapse1 <- gsub("    ", "", collapse1, fixed = TRUE)  # removing extra spaces
        }

        pos <- regexpr("data.frame(", collapse1, fixed = TRUE)
        if (pos > 0)     
        {
            collapse1 <- substring(collapse1, 12, nchar(collapse1)-1)
        }
        pos <- regexpr("list(", collapse1, fixed=TRUE)
        if (pos > 0) 
        {
             collapse1 <- substring(collapse1, 6, nchar(collapse1)-1)
        }      

        collapse2 <- obj2[[8]]$collapse
        if (is.null(collapse2)) {collapse2 <- obj2[[8]]$pmodels}
        if (!is.null(obj2$"pmodelsText")) {collapse2 <- obj2$"pmodelsText"}
        if (is.null(collapse2)) 
        {
            if (is.null(obj2[[8]]$curve))
            {
                collapse2 <- "1 (for all parameters)"
            } else {
                collapse2 <- paste(deparse(obj2[[8]]$curve), "(for all parameters)")
            }
        } else {
            if (!is.character(collapse2)) 
            {
                collapse2 <- paste(deparse(collapse2), collapse = "")
            }
            collapse2 <- gsub("    ", "", collapse2, fixed = TRUE)  # removing extra spaces
        }
    
        pos <- regexpr("data.frame(", collapse2, fixed = TRUE)    
        if (pos > 0)
        {
            collapse2 <- substring(collapse2, 12, nchar(collapse2) - 1)
        }

        pos <- regexpr("list(", collapse2, fixed = TRUE)
        if (pos > 0)
        {
            collapse2 <- substring(collapse2, 6, nchar(collapse2) - 1)
        }               
    
    
        ## Omitting collapse line if content is the same in both lines
        if (identical(collapse1, collapse2)) {colLine <- FALSE} else {colLine <- TRUE}
    
        if ( (!is.null(obj1[[8]]$pmodels)) || (!is.null(obj2[[8]]$pmodels)) )
        {
            fctStart <- " fct:     "
            colStart <- " pmodels: "
        } else {
            fctStart <- " fct:      "
            colStart <- " pmodels: "
            # changed July 4 2011
        }
        cat("1st model\n")
        fctInfo <- ifelse(is.null(obj1$"text"), deparse(obj1[[8]]$fct), obj1$"text")
        cat(paste(fctStart, fctInfo, "\n", sep = ""))        
        if (colLine) {cat(paste(colStart, collapse1, "\n", sep = ""))}
        cat("2nd model\n")
        fctInfo <- ifelse(is.null(obj2$"text"), deparse(obj2[[8]]$fct), obj2$"text")
        cat(paste(fctStart, fctInfo, "\n", sep = ""))        
        if (colLine) {cat(paste(colStart, collapse2, "\n", sep = ""))}
        cat("\n")
        }

    dataFra <- data.frame(dfModel, loglik, dfDiff, testStat, pVal)
    dimnames(dataFra) <- list(rowNames, colNames)
    structure(dataFra, heading = headName, class = c("anova", "data.frame"))
}
