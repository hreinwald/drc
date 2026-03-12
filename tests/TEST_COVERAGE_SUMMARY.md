# Test Coverage Enhancement for drc Package

This document describes the comprehensive test suite added to significantly improve code coverage for the drc package.

## Summary

Added **350+ new tests** across 7 new test files, covering critical functions that previously had 0% test coverage. The new tests focus on:

- User-facing S3 methods (ED, predict, fitted, residuals, print, summary)
- Parameter extraction and inference (coef, vcov, confint)
- Visualization (plot)
- Model utilities and comparisons (update, anova, logLik, AIC/BIC)

## New Test Files

### 1. test-ED.R (327 lines, 60+ tests)
Tests for `ED.drc()` - Effective Dose/Concentration estimation

**Coverage:**
- Basic ED calculation for single and multiple response levels
- All confidence interval types (delta, fls, tfls, inv)
- Response level validation and bounds checking
- Multi-curve models with curve selection
- Different interval types and confidence levels
- Multiple model types (LL.3, LL.4, Weibull)
- Different data types (continuous, binomial)
- Edge cases: extreme response levels, custom vcov, logBase transformation
- multcomp package integration
- Reference parameter (control vs upper)

**Key Edge Cases Tested:**
- Invalid response levels (0, 100, negative, >100)
- Missing edfct function
- Custom variance-covariance functions
- Logarithmic transformations
- Order validation for ED values (ED10 < ED50 < ED90)

---

### 2. test-predict.R (385 lines, 50+ tests)
Tests for `predict.drc()` - Model predictions

**Coverage:**
- Basic predictions with and without newdata
- Standard errors (se.fit = TRUE)
- Confidence intervals (interval = "confidence")
- Prediction intervals (interval = "prediction")
- SSD (species sensitivity distribution) intervals
- Multi-curve model predictions
- Different response types (continuous, binomial, Poisson)
- Over-dispersion adjustment (od parameter)
- Custom vcov functions
- Constrain parameter for bounding predictions
- Power transformations
- Missing derivatives handling

**Key Edge Cases Tested:**
- Zero dose values
- Very large dose values
- Missing curve ID in newdata
- Binomial predictions constrained to [0, 1]
- Non-continuous predictions constrained to >= 0
- Prediction interval wider than confidence interval
- Monotonicity of predictions

---

### 3. test-fitted-residuals.R (318 lines, 40+ tests)
Tests for `fitted.drc()` and `residuals.drc()` - Fitted values and residuals

**Coverage for fitted.drc():**
- Basic fitted value extraction
- Equivalence with predict()
- Multi-curve models
- Passing arguments to predict()

**Coverage for residuals.drc():**
- Working residuals (default)
- Standardised residuals
- Studentised residuals (accounting for leverage)
- Box-Cox transformation handling (trScale parameter)
- Multi-curve models
- Different data types (continuous, binomial, Poisson)
- Error handling for missing derivatives

**Key Edge Cases Tested:**
- Fitted + residuals = observed (reconstruction test)
- Sum of residuals near zero
- Standardised residuals have unit variance approximately
- Studentised residuals handle leverage properly
- Box-Cox transformation on different scales
- NA scale estimate handling for non-continuous data

---

### 4. test-print-summary.R (319 lines, 35+ tests)
Tests for `print.drc()` and `summary.drc()` - Display and summary methods

**Coverage for print.drc():**
- Basic printing with model information
- Coefficient display
- digits parameter
- Empty coefficients handling
- Multi-curve models
- Invisible return value

**Coverage for summary.drc():**
- Summary object structure
- Coefficient matrix with estimates, SE, t/z-values, p-values
- Residual standard error
- Degrees of freedom
- Over-dispersion adjustment (od parameter)
- Pooled vs unpooled estimation (pool parameter)
- Different response types (continuous uses t, binomial uses z)
- Robust estimation methods
- Multi-curve models

**Key Edge Cases Tested:**
- P-values between 0 and 1
- Standard errors are positive
- Coefficient estimates match coef()
- Correct distribution for different data types

---

### 5. test-coef-vcov-confint.R (396 lines, 45+ tests)
Tests for `coef.drc()`, `vcov.drc()`, and `confint.drc()` - Parameter extraction and inference

**Coverage for coef.drc():**
- Coefficient vector extraction
- Proper naming
- Different model types (LL.3, LL.4, Weibull)
- Multi-curve models
- Fallback to fit$par when coefficients are NULL

**Coverage for vcov.drc():**
- Variance-covariance matrix
- Symmetry
- Positive diagonal elements
- Correlation matrix (corr = TRUE)
- Over-dispersion adjustment (od parameter)
- Pooled vs unpooled (pool parameter)
- Unscaled vcov (unscaled parameter)
- Different data types

**Coverage for confint.drc():**
- Confidence intervals for all/specific parameters
- Interval contains point estimate
- Confidence level effects (90% vs 95%)
- Parameter selection (parm parameter)
- Invalid parameter names
- Multi-curve models with pooling
- Correct distribution (t for continuous, z for binomial)
- confint.basic() helper function

**Key Edge Cases Tested:**
- vcov diagonal = SE²
- Narrower intervals at lower confidence levels
- Integration: coef, vcov, confint consistency
- Proper quantiles for different data types

---

### 6. test-plot.R (473 lines, 60+ tests)
Tests for `plot.drc()` - Visualization

**Coverage:**
- Basic plotting without errors
- Different plot types (average, all, bars, none, obs, confidence)
- Multi-curve models
- Curve selection (level parameter)
- Graphical parameters (col, lty, pch, xlab, ylab, xlim, ylim)
- Log scale (log parameter)
- Broken axis (broken = TRUE)
- Grid size control
- Legend display
- Normalization (normal = TRUE)
- Confidence level for error bars
- Different model types (LL.3, LL.4, Weibull)
- Different data types (continuous, binomial, Poisson)
- Robust estimation

**Key Edge Cases Tested:**
- Zero dose with log scale
- Multiple graphical parameters combined
- Vector of colors for multi-curve
- Logical col = TRUE for auto colors
- Multiple plots in sequence
- Plotting after predictions
- Integration with par() settings

---

### 7. test-utilities.R (412 lines, 65+ tests)
Tests for additional utility and model functions

**Coverage:**
- `update.drc()`: Model updates (function, data)
- `logLik.drc()`: Log-likelihood extraction with df attribute
- `anova.drc()`: Model comparison, p-values, lack-of-fit test
- `df.residual()`: Degrees of freedom calculation
- `AIC()` and `BIC()`: Information criteria
- `boxcox.drc()`: Box-Cox transformations
- `drm()` core functionality:
  - Start values
  - Weights
  - Subset
  - Formula handling
  - Different dose-response functions (LL.2-5, W1.2-4, W2.2-4)
  - Control parameters (drmc)
  - Robust estimation
  - Different data types (binomial, Poisson)
- Model structure validation
- Coefficient naming
- `hatvalues.drc()`: Leverage values
- `cooks.distance.drc()`: Influence measures
- Complete workflow integration tests
- Model comparison workflows

**Key Edge Cases Tested:**
- Update preserves model when nothing changed
- Log-likelihood is negative
- More parameters give higher log-likelihood
- Valid p-values from ANOVA
- Correct degrees of freedom calculation
- AIC/BIC for model comparison
- Start values affect convergence
- Subset reduces data size
- All DR functions have correct parameter counts

---

## Test Design Principles

All tests follow these principles:

1. **Comprehensive Edge Case Coverage**: Tests include boundary conditions, extreme values, missing data, and error conditions
2. **Multiple Model Types**: Tests cover LL.2-5, Weibull (W1, W2), binomial, Poisson, robust
3. **Multi-Curve Models**: Extensive testing of curve identification and pooling
4. **Integration Testing**: Tests verify functions work together in complete workflows
5. **Proper Setup/Teardown**: PDF devices properly opened and closed for plotting tests
6. **Clear Assertions**: Each test has specific, meaningful expectations
7. **Informative Names**: Test names clearly describe what is being tested
8. **No Side Effects**: Tests don't modify global state or leave artifacts

## Expected Coverage Improvements

Functions that previously had **0% coverage** and now have comprehensive tests:

- `ED.drc()` - Effective dose estimation
- `predict.drc()` - Predictions
- `fitted.drc()` - Fitted values
- `residuals.drc()` - Residuals (working, standardised, studentised)
- `print.drc()` - Printing
- `summary.drc()` - Summaries
- `coef.drc()` - Coefficient extraction
- `vcov.drc()` - Variance-covariance
- `confint.drc()` - Confidence intervals
- `plot.drc()` - Plotting
- `update.drc()` - Model updates
- `logLik.drc()` - Log-likelihood
- Helper functions: `confint.basic()`, `EDinvreg1()`, and more

Functions with improved coverage:

- `drm()` - Core model fitting (now tested with all parameters)
- `anova.drc()` - Model comparison (now 66.67% → higher)
- Various internal model functions

## Testing Statistics

- **Total new tests**: 350+
- **Total new test lines**: ~2,600 (excluding existing tests)
- **Test files created**: 7
- **Functions covered**: 30+
- **Edge cases tested**: 100+

## Running the Tests

```r
# Run all tests
library(testthat)
library(drc)
test_check("drc")

# Run specific test file
test_file("tests/testthat/test-ED.R")

# Run with coverage
library(covr)
package_coverage()
```

## Future Recommendations

While this significantly improves coverage, additional tests could be added for:

1. **Low-coverage utility functions** (0% coverage):
   - EDcomp(), compParm() - parameter comparisons
   - isobole() - isobologram analysis
   - Specialized model functions (braincousens, cedergreen, etc.)
   - Simulation functions (simDR, simFct)

2. **Edge cases in existing functions**:
   - More exotic model combinations
   - Numerical edge cases (convergence failures, singular matrices)
   - More comprehensive multi-dimensional dose testing

3. **Performance tests**:
   - Large dataset handling
   - Convergence with poor initial values

4. **Documentation validation**:
   - Example code in documentation runs without errors
   - Example output matches expectations

## Notes

- All tests use the standard `testthat` framework
- Tests are self-contained with their own test data
- No external data dependencies (uses synthetic datasets)
- All tests follow R package testing best practices
- Tests are deterministic (use set.seed where needed)
