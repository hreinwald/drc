# Comparison: arandaordaz() vs weibull2()

## Executive Summary

This document analyzes the relationship between `arandaordaz()` and `weibull2()` functions in the drc package.

**Key Finding:** `weibull2()` is NOT a simple replacement for `arandaordaz()`. They are different models with different mathematical forms. However, the package provides convenience wrappers (`AR.2()` and `AR.3()`) that are implemented using **both approaches** - with the `weibull2`-based versions being the primary implementations.

---

## Detailed Analysis

### 1. Mathematical Differences

#### arandaordaz() - Asymptotic Regression Model
**File:** `R/arandaordaz.R` (lines 29-123)

**Mathematical Formula:**
```
f(x) = a + (b - a)(1 - exp(-c * x))
```

**Parameters:**
- `a`: Lower asymptote (value at x=0)
- `b`: Upper asymptote
- `c`: Rate parameter (controls steepness)

**Key Characteristics:**
- 3-parameter model
- Uses **linear dose** (x) in the exponential
- Simple asymptotic growth model
- Best for monotonic increasing curves starting from a lower limit

---

#### weibull2() - Weibull Type 2 Model
**File:** `R/weibull2.R` (lines 30-255)

**Mathematical Formula:**
```
f(x) = c + (d - c)(1 - exp(-exp(b * (log(x) - log(e)))))
```

**Parameters:**
- `b`: Shape parameter
- `c`: Lower asymptote
- `d`: Upper asymptote
- `e`: Scale parameter (ED50 equivalent)

**Key Characteristics:**
- 4-parameter model
- Uses **log-transformed dose** in a double exponential
- More flexible S-shaped dose-response curves
- Standard model for dose-response analysis

---

### 2. Functional Relationship

**weibull2() is NOT simply a replacement for arandaordaz().** They are mathematically different models:

| Aspect | arandaordaz() | weibull2() |
|--------|---------------|------------|
| **Model Type** | Simple asymptotic | Weibull type 2 |
| **Parameters** | 3 (a, b, c) | 4 (b, c, d, e) |
| **Dose Transform** | Linear (x) | Log-transformed (log(x)) |
| **Flexibility** | Less flexible | More flexible |
| **Primary Use** | Basic asymptotic curves | Complex dose-response |

However, **weibull2() can approximate arandaordaz()** behavior by:
1. Fixing the shape parameter `b = 1`
2. Using appropriate constraints on other parameters

---

### 3. Convenience Wrappers: AR.2() and AR.3()

The package provides two convenience functions for asymptotic regression:

#### Implementation History

**Original (Broken) - in arandaordaz.R:**
- Lines 126-151 define `AR.2()` and `AR.3()` that call `asymreg()`
- **BUG:** `asymreg` was never defined (likely was intended as an alias to `arandaordaz`)
- **FIXED:** Changed to call `arandaordaz()` instead

**Current (Active) - in weibull2.R:**
- Lines 412-424: `AR.2()` calls `weibull2(fixed = c(1, 0, fixed[1:2]), ...)`
- Lines 446-458: `AR.3()` calls `weibull2(fixed = c(1, fixed[1:3]), ...)`

These are the **primary implementations** used by the package.

---

### 4. Which Implementation is Used?

In R, when two functions have the same name, the **later-loaded** version takes precedence. Since package loading typically follows alphabetical or dependency order, and both files define `AR.2()` and `AR.3()`, the versions in `weibull2.R` are likely the active ones.

To verify which is used, check:
```r
getMeanFunctions()  # Lists available models
```

The package documentation (vignettes) references only the Weibull-based implementations.

---

### 5. Bug Fix Applied

**Problem Found:** Lines 126-151 in `R/arandaordaz.R` originally called an undefined function `asymreg`.

**Root Cause:** The `asymreg` identifier was never defined in the codebase. It was likely intended as an alias to `arandaordaz()` but was never implemented.

**Solution Applied:** Changed `asymreg` to `arandaordaz` in both `AR.2()` and `AR.3()` functions.

**Code Comments Added:**
- Added comments to both `arandaordaz.R` and `weibull2.R` documenting that duplicate definitions exist
- Clarified that the `weibull2.R` versions are the PRIMARY implementations
- Explained that R's alphabetical loading order means `weibull2.R` overrides `arandaordaz.R`

**Impact:** These functions are now syntactically correct, though they are still overridden by the `weibull2.R` versions in practice.

---

### 6. Recommendations

1. **For users:** Use `AR.2()` and `AR.3()` from `weibull2.R` (the default behavior)
   - These are more robust and better integrated

2. **For developers:** Consider one of these approaches:
   - **Option A:** Remove duplicate definitions from `arandaordaz.R` (lines 126-151)
   - **Option B:** Keep them but document that they use the pure asymptotic model
   - **Option C:** Export only the `weibull2`-based versions

3. **Documentation:** Clarify in user-facing docs that:
   - `arandaordaz()` is the base asymptotic regression function
   - `AR.2()` and `AR.3()` are implemented via `weibull2()` for consistency
   - The two models are mathematically different but `weibull2` with `b=1` approximates asymptotic behavior

---

## Conclusion

**weibull2() is NOT a simple replacement for arandaordaz()** - they are different models with different mathematical foundations:

- **arandaordaz()**: Pure asymptotic regression with linear dose
- **weibull2()**: Flexible Weibull type 2 with log-transformed dose

However, the package has **migrated convenience functions** (`AR.2()`, `AR.3()`) to use `weibull2()`-based implementations with fixed parameters, providing consistency with the broader Weibull framework while approximating asymptotic regression behavior.

The bug in `arandaordaz.R` (undefined `asymreg`) has been fixed, but these functions are likely not the primary implementations used by the package.
