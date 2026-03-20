# Summary of Fixes for ED() and noEffect() Issues

This document summarizes the fixes applied to address the three reported issues with ED(), maED(), and noEffect() functions when working with EXD.3 and LL.5 models.

## Issue #1: ED() Error for EXD.3 - "incorrect number of dimensions"

### Problem
When calling `ED(m, 50)` on an EXD.3 model, the function threw an error:
```
Error in indexMat[, curveOrder, drop = FALSE] :
  incorrect number of dimensions
```

This occurred even though `summary(m)` correctly showed the e:Intercept parameter estimate.

### Root Cause
The `indexMat` object was not always a matrix when retrieved from the model object. For models with few estimated parameters (like EXD.3 with c and d fixed), `indexMat` could be a vector instead of a matrix. When the code tried to subset it with `[, curveOrder]`, it failed because vectors don't have columns.

### Fix
Added a check in `/home/runner/work/drc/drc/R/ED.drc.R` (lines 172-179) to ensure `indexMat` is always treated as a matrix:

```r
# Ensure indexMat is always a matrix, even if it's a single column vector
if (!is.matrix(indexMat)) {
  indexMat <- as.matrix(indexMat)
  # Set column names to match parmMat if they exist
  if (!is.null(colnames(parmMat))) {
    colnames(indexMat) <- colnames(parmMat)
  }
}
```

### Result
ED() now works correctly for EXD.3 and other models with few estimated parameters.

---

## Issue #2: ED() Returns NaN for LL.5 with Warning

### Problem
When calling `ED(m_LL5, 50)` on an LL.5 model, the function returned NaN with a warning:
```
Warning message:
  In log(exp(-tempVal/parmVec[5]) - 1) : NaNs produced
```

The summary showed parameter estimates with NaN standard errors, suggesting the model was poorly conditioned.

### Root Cause
In the LL.5 ED function (llogistic.R, line 243), the calculation:
```r
EDp <- parmVec[4]*(exp(-tempVal/parmVec[5])-1)^(1/parmVec[1])
```

The derivative calculation at line 246 included:
```r
log(exp(-tempVal/parmVec[5])-1)
```

When `exp(-tempVal/parmVec[5]) - 1` is ≤ 0, the log produces NaN. This occurs when the fitted model parameters suggest an EC50 that is outside the valid range for the data (e.g., when the dose-response relationship is very weak or the model is ill-conditioned).

### Fix
Added a validity check in `/home/runner/work/drc/drc/R/llogistic.R` (lines 243-257) to detect when the ED value is invalid:

```r
tempVal <- log((100-p)/100)
expTerm <- exp(-tempVal/parmVec[5])

# Check if expTerm - 1 is valid (must be positive for log)
if (expTerm <= 1) {
    # ED value is outside the valid range or model is ill-conditioned
    EDp <- Inf
    EDder <- rep(NA, 5)
} else {
    EDp <- parmVec[4]*(expTerm-1)^(1/parmVec[1])
    EDder <- EDp*c(-log(expTerm-1)/(parmVec[1]^2),
                   0, 0, 1/parmVec[4],
                   expTerm*tempVal/(parmVec[5]^2)*(1/parmVec[1])*((expTerm-1)^(-1)))
}
```

### Result
ED() now returns Inf (indicating the EC50 is outside the valid range) instead of NaN with a warning. This makes it clear that the issue is with the model fit or data, not a bug in the code.

---

## Issue #3: noEffect() Returns Df = 0 for EXD.3

### Problem
When calling `noEffect(m)` on an EXD.3 model with c and d fixed, the function returned Df = 0, which seemed incorrect since summary(m) showed an e:Intercept parameter being estimated.

### Root Cause
This is actually **correct behavior**, but confusing. Here's why:

- EXD.3 with c and d fixed has only 1 estimated parameter (e) + 1 variance parameter = 2 total df
- The null model (intercept-only) has 1 parameter (mean) + 1 variance parameter = 2 total df
- Therefore, the degrees of freedom difference is 2 - 2 = 0

The noEffect() test compares the **number of estimated parameters**, not the model structure. When most parameters are fixed, the dose-response model has the same complexity (in terms of df) as a simple mean, even though it has more structure.

### Fix
Added a warning in `/home/runner/work/drc/drc/R/noEffect.R` (lines 47-53) to alert users when this situation occurs:

```r
# Check if degrees of freedom difference is valid
if (dfDiff <= 0) {
    warning("Degrees of freedom difference is ", dfDiff,
            ". This may indicate that the dose-response model has no additional ",
            "parameters compared to the null model (e.g., when parameters are fixed). ",
            "The likelihood ratio test may not be meaningful in this case.")
}
```

### Result
noEffect() still returns Df = 0 when appropriate, but now warns the user that the test may not be meaningful in this case.

### Recommendation
When many parameters are fixed, the noEffect() test is not appropriate for assessing whether there is a dose effect. Instead, consider:
1. Using a model with more free parameters
2. Comparing models with different numbers of fixed parameters using AIC/BIC
3. Examining the model fit visually
4. Using other diagnostic tools

---

## Summary

All three issues have been addressed:

1. **ED() for EXD.3**: Fixed by ensuring indexMat is always treated as a matrix
2. **ED() for LL.5**: Fixed by handling invalid ED values gracefully (returning Inf instead of NaN)
3. **noEffect() for EXD.3**: Not a bug, but added a warning to clarify when the test is not meaningful

These fixes improve the robustness and clarity of the drc package when working with models that have fixed parameters or are poorly conditioned.
