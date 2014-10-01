
f1 <- lFormula(formula = Reaction ~ Days + (Days|Subject), data = sleepstudy)

d1 <- do.call(mkLmerDevfun, f1)

p1 <- optimizeLmer(d1)

m1 <- mkMerMod(rho = environment(d1),
  opt = p1,
  reTrms = f1$reTrms,
  fr = f1$fr)
  
###

f2 <- lFormulaLwts(formula = Reaction ~ Days + (Days|Subject), data = sleepstudy)

d2 <- do.call(mkLmerDevfunLwts, f2)
  
