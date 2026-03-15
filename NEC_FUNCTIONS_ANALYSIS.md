# Analysis of NEC Functions in the drc Package

## Executive Summary

The `drc` R package contains 4 NEC (No Effect Concentration) functions: `NEC`, `NEC.2`, `NEC.3`, and `NEC.4`. After thorough analysis, **all functions are necessary and serve distinct purposes**. There is no redundancy.

## Function Overview

### Base Implementation: NEC (Not Exported)

**Location:** `R/nec.R` (lines 38-91)

**Purpose:** The core implementation function that provides the flexible NEC dose-response model.

**Model Equation:**
```
f(x) = c + (d-c) * exp(-b(x-e) * I(x-e))
```
where `I(x-e)` is an indicator function (0 when x≤e, 1 when x>e).

**Parameters:**
- `b`: Slope/rate parameter
- `c`: Lower limit (control response)
- `d`: Upper limit (maximum response)
- `e`: NEC threshold (no effect concentration)

**Key Features:**
- Accepts a `fixed` argument to specify which parameters should be fixed
- Uses log-logistic self-starter function for initialization
- Returns a model list with nonlinear function, self starter, and parameter names
- **Not exported** in NAMESPACE - serves as an internal implementation engine

---

### NEC.2 (2-parameter variant)

**Location:** `R/nec.R` (lines 109-122)

**Purpose:** Convenience wrapper for highly constrained scenarios where both lower and upper limits are known.

**Implementation:**
```r
NEC(fixed = c(fixed[1], 0, upper, fixed[2]),
    names = c(names[1], "c", "d", names[2]), ...)
```

**Free Parameters:** 2
- `b`: Slope parameter (free)
- `e`: NEC threshold (free)

**Fixed Parameters:**
- `c`: Fixed at 0
- `d`: Fixed at user-specified value (default 1)

**Use Cases:**
- Response bounded on a known scale (e.g., 0-1 for proportions, 0-100 for percentages)
- Both bounds are well-defined from experimental design
- Focus estimation on slope and threshold only
- Reduces model complexity and improves parameter identifiability

---

### NEC.3 (3-parameter variant)

**Location:** `R/nec.R` (lines 138-151)

**Purpose:** Most common variant - assumes zero baseline response with variable maximum.

**Implementation:**
```r
NEC(fixed = c(fixed[1], 0, fixed[2:3]),
    names = c(names[1], "c", names[2:3]), ...)
```

**Free Parameters:** 3
- `b`: Slope parameter (free)
- `d`: Upper limit (free)
- `e`: NEC threshold (free)

**Fixed Parameters:**
- `c`: Fixed at 0

**Use Cases:**
- Standard toxicological/biological scenarios
- Baseline response is zero (no treatment/exposure)
- Maximum response varies by treatment
- Balances flexibility with model stability
- Reduces overfitting compared to NEC.4

---

### NEC.4 (4-parameter variant)

**Location:** `R/nec.R` (lines 167-177)

**Purpose:** Full flexibility - all parameters estimated from data.

**Implementation:**
```r
NEC(fixed = fixed, names = names, ...)
```

**Free Parameters:** 4
- `b`: Slope parameter (free)
- `c`: Lower limit (free)
- `d`: Upper limit (free)
- `e`: NEC threshold (free)

**Use Cases:**
- No biological constraints on parameters
- Both baseline and maximum responses vary
- Model selection and comparison workflows
- Maximum flexibility when data supports it
- Cases where control/baseline response is non-zero and unknown

---

## Comparison Matrix

| Aspect | NEC (base) | NEC.2 | NEC.3 | NEC.4 |
|--------|-----------|-------|-------|-------|
| **Exported** | No | Yes | Yes | Yes |
| **Free Parameters** | Configurable | 2 (b, e) | 3 (b, d, e) | 4 (b, c, d, e) |
| **Fixed c (lower)** | Configurable | 0 | 0 | Free |
| **Fixed d (upper)** | Configurable | User-defined | Free | Free |
| **Model Complexity** | Depends | Lowest | Medium | Highest |
| **When to Use** | Internal only | Known bounds | Zero baseline | Full flexibility |
| **Identifiability** | Depends | Excellent | Good | May be challenging |

---

## Design Pattern Analysis

This follows the **standard drc package design pattern** used consistently across all model families:

### Examples of Similar Patterns in drc:

1. **Log-logistic models:** `llogistic`, `LL.2`, `LL.3`, `LL.3u`, `LL.4`, `LL.5`
2. **Weibull type 1:** `weibull1`, `W1.2`, `W1.3`, `W1.3u`, `W1.4`
3. **Weibull type 2:** `weibull2`, `W2.2`, `W2.3`, `W2.3u`, `W2.4`
4. **Gompertz:** `gompertz`, `G.2`, `G.3`, `G.3u`, `G.4`
5. **Log-normal:** `lnormal`, `LN.2`, `LN.3`, `LN.3u`, `LN.4`
6. **Brain-Cousens:** `braincousens`, `BC.4`, `BC.5`
7. **Cedergreen-Ritz-Streibig:** `cedergreen`, `CRS.4a`, `CRS.4b`, `CRS.4c`, `CRS.5`, `CRS.5a`, `CRS.5b`, `CRS.5c`, `CRS.6`

### Pattern Structure:

1. **Base function** (e.g., `llogistic`, `NEC`)
   - Provides core implementation with full parameter flexibility
   - Often not exported (used internally)
   - Accepts `fixed` argument for parameter constraints

2. **Numbered variants** (e.g., `LL.2`, `LL.3`, `LL.4`, `LL.5`)
   - Convenience wrappers with common parameter combinations
   - Exported for user convenience
   - Number indicates count of free parameters
   - Each serves specific biological/experimental scenarios

### Benefits of This Design:

- **User convenience**: Common cases are easy to specify
- **Parameter identifiability**: Constraining parameters when appropriate improves estimation
- **Model selection**: Easy to compare nested models
- **Biological meaning**: Parameter constraints reflect experimental knowledge
- **Backwards compatibility**: Adding variants doesn't break existing code
- **Documentation clarity**: Each variant can have specific use-case documentation

---

## Redundancy Assessment

### Question: Are any NEC functions redundant?

**Answer: NO - All functions are necessary.**

### Reasoning:

1. **NEC (base function)**
   - **Cannot be removed**: It contains the actual mathematical implementation
   - All other functions are wrappers that call `NEC` with specific constraints
   - Removing it would break NEC.2, NEC.3, and NEC.4

2. **NEC.2**
   - **Unique purpose**: Only variant with both upper and lower limits fixed
   - **Distinct use case**: Bounded response scales (proportions, percentages)
   - **Cannot be replicated**: NEC.3 fixes only lower limit, NEC.4 fixes nothing
   - **Statistical benefit**: Reduces parameters from 4 to 2, greatly improving identifiability

3. **NEC.3**
   - **Most common scenario**: Standard toxicology with zero baseline
   - **Optimal balance**: More flexible than NEC.2, more stable than NEC.4
   - **Common convention**: Matches typical experimental designs where control = 0
   - **Unique constraint**: Only variant fixing lower limit while freeing upper limit

4. **NEC.4**
   - **Essential for flexibility**: Only way to estimate all 4 parameters
   - **Model selection**: Needed for comparing against constrained models
   - **Non-zero baselines**: Only option when control response is unknown and non-zero
   - **Diagnostic tool**: Helps determine if constraints are appropriate

### If Functions Were Combined:

Users would need to manually specify constraints:
```r
# Current (user-friendly):
drm(y ~ x, data = mydata, fct = NEC.3())

# If combined (cumbersome):
drm(y ~ x, data = mydata, fct = NEC(fixed = c(NA, 0, NA, NA)))
```

This would:
- Reduce usability
- Increase errors (wrong constraint specifications)
- Eliminate helpful documentation for common cases
- Break backwards compatibility
- Deviate from established drc package conventions

---

## Verification of Usage

All three variants are properly:
- **Exported** in NAMESPACE (line 38)
- **Documented** with individual .Rd files in `man/`
- **Cross-referenced** in each other's documentation
- **Example provided** in NEC.Rd showing NEC.4 usage

The base `NEC` function is:
- **Not exported** (correct - it's an internal implementation)
- **Referenced** in documentation via `\code{\link{NEC}}`
- **Called by** all three numbered variants

---

## Conclusion

**All 4 NEC functions should be retained.**

The design represents:
1. **Sound software architecture**: Internal implementation separated from user interface
2. **Statistical best practice**: Providing appropriate model complexity for different scenarios
3. **User experience optimization**: Common cases are simple, complex cases are possible
4. **Package consistency**: Matches the established pattern used for all other model families

**Recommendation: No changes needed.** The current implementation is well-designed, follows package conventions, and serves distinct user needs.

---

## References

**Source Code:**
- Implementation: `R/nec.R`
- Exports: `NAMESPACE` (line 38)
- Documentation: `man/NEC.Rd`, `man/NEC.2.Rd`, `man/NEC.3.Rd`, `man/NEC.4.Rd`

**Scientific Reference:**
Pires, A. M., Branco, J. A., Picado, A., Mendonca, E. (2002)
Models for the estimation of a 'no effect concentration',
*Environmetrics*, **13**, 15-27.
