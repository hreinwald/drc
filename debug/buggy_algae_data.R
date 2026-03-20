library(dplyr)
library(data.table)
devtools::load_all()

# Import example data
dt = fread("./tests/buggy_algae_data_light.tsv")

# FCT.ls for running mselect()
lowerl = 1e-9
upperl = dt %>% filter(conc == 0) %>% .[["yield"]] %>% mean(.) # 649565.1
FCT.ls <- list(
  LL.4(fixed = c(NA, lowerl, upperl, NA)),
  LL.5(fixed = c(NA, lowerl, upperl, NA, NA)),
  LN.4(fixed = c(NA, lowerl, upperl, NA)),
  W1.4(fixed = c(NA, lowerl, upperl, NA)),
  W2.4(fixed = c(NA, lowerl, upperl, NA)),
  #EXD.3(fixed = c(lowerl, upperl, NA)),
  # Hormesis specific models - c flexibel
  BC.5(fixed = c(NA, lowerl, upperl, NA, NA)),
  CRS.5(fixed = c(NA, lowerl, upperl, NA, NA), alpha_type = "a"), #return_error = T),
  CRS.5(fixed = c(NA, lowerl, upperl, NA, NA), alpha_type = "b"), #return_error = T),
  CRS.5(fixed = c(NA, lowerl, upperl, NA, NA), alpha_type = "c") #return_error = T) 
)

# Fit initial model with EXD.3
m = drm(
  formula = yield ~ conc,
  weights = NULL,
  data = dt,
  fct = EXD.3(fixed = c(lowerl, upperl, NA)),
  #fct = LL.4(fixed = c(NA, lowerl, upperl, NA)),
  type = "continuous"
)

# Inspect drm results
plot(m, type = "all")
m$fct$name
noEffect(m) # <- why is Df 0 here??? is that a bug?
summary(m)
modelFit(m)
AIC(m)

# Find best fitting model
mselect(m, fctList = FCT.ls, nested = TRUE)

# According to mselect EXD.3 is the best fitting model following the IC score. 
# Hence it would make sense to retrieve the EC50 for this model.
summary(m) # This gave us an e intercept at 16193 -> Hence an EC50 value > maxconc would be expected.
try( ED(m, 50) ) # This returns an error which I did not expect ... 

# This fails completely now! 
try( maED(m, fctList = FCT.ls, respLev = 50) )
maED_robust(m, fctList = FCT.ls, respLev = 50) 

# Let's see what is happening with LL.5 function
m_LL5 = update(m, fct = LL.5(fixed = c(NA, lowerl, upperl, NA, NA)))
ED(m_LL5, 50)

# Ok I understand a view issues now! 

# 1. With EXD.3 no EC50 value can be computed -> is there a reason or might this be a bug? 
# > ED(m, 50)
# > Error in indexMat[, curveOrder, drop = FALSE] : 
# >  incorrect number of dimensions

# However this makes litle sense as summary() returns the e:Intercept ... 
# > summary(m)
# > Estimate Std. Error t-value p-value
# > e:(Intercept)    16193      13838  1.1702  0.2455
# Hence, to me the error in the ED(m, 50) call seems like a bug ... 

# 2. With LL.5 no EC50 value can be computed but ED() function returns NaN with a warning
# -> is there a reason or might this be a bug?
# > m_LL5 = update(m, fct = LL.5(fixed = c(NA, lowerl, upperl, NA, NA)))
# > ED(m_LL5, 50)
# > Warning message:
# >  In log(exp(-tempVal/parmVec[5]) - 1) : NaNs produced

# For summary(m_LL5) we get the following output:
# Estimate Std. Error t-value p-value
# b:(Intercept)  0.131147        NaN     NaN     NaN
# e:(Intercept)  0.261586        NaN     NaN     NaN
# f:(Intercept) -0.036225        NaN     NaN     NaN
# Hencem the e intercept is estimated here as well. 

# 3. I don't understand why noEffect(m) returns Df  with 0 ??? Is that a bug?

# Why is it impossible to compute an EC50 via ED() for EXD.3 and LL.5 for this particular dataset?
# Is there a bug? 