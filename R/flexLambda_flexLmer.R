#' @export
flexLmer <- function(formula, data, specials=c("cs", "d", "ar1d"), verbose=0L){
	#split off reGenerator terms:
	frmlterms <- terms(formula, specials=specials)
	termnames <- attr(frmlterms, "variables")
	where.specials <- unlist(attr(frmlterms, "specials"))-1 #offset by one since response is counted
	lmerformula <- formula(paste(formula[[2]], "~", 
								 paste(attr(frmlterms, "term.labels")[-where.specials], collapse="+")))
	reGenerators <- as.formula(paste("~",
							paste(attr(frmlterms, "term.labels")[where.specials], collapse="+")))
	# see example(modular)
	lmod <- lFormula(lmerformula, data, reGenerators=reGenerators)
	devfun <- do.call(mkLmerDevfun, lmod)
	opt <- optimizeLmer(devfun, verbose=verbose)
	#FIXME: this should not be necessary...
	environment(devfun)$pp$theta <- opt$par
	list(model=mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr),
		 opt=opt, devfun=devfun)
}
