##' @rdname modular
##' @param fr A model frame containing the variables needed to create an
##'   \code{\link{lmerResp}} or \code{\link{glmResp}} instance
##' @param X fixed-effects design matrix
##' @param reTrms information on random effects structure (see \code{\link{mkReTrms}})
##' @param REML (logical) fit restricted maximum likelihood model?
##' @param start starting values
##' @param verbose print output?
##' @return \bold{mkLmerDevfunLwts, mkGlmerDevfun}: A function to calculate deviance
##' (or restricted deviance) as a function of the theta (random-effect) parameters
##' (for GlmerDevfun, of beta (fixed-effect) parameters as well).  These deviance
##' functions have an environment containing objects required for their evaluation.
##' CAUTION: The output object of \code{mk(Gl|L)merDevfun} is an \code{\link{environment}}
##' containing reference class objects (see \code{\link{ReferenceClasses}}, \code{\link{merPredD-class}},
##' \code{\link{lmResp-class}}), which behave in ways that may surprise many users. For example, if the
##' output of \code{mk(Gl|L)merDevfun} is naively copied, then modifications to the original will
##' also appear in the copy (and vice versa). To avoid this behavior one must make a deep copy
##' (see \code{\link{ReferenceClasses}} for details).
##' \cr
##' \cr
##' @export
mkLmerDevfunLwts <- function(fr, X, reTrms, REML = TRUE, start = NULL, verbose=0, control=lmerControl(), ...) {

    ## FIXME: make sure verbose gets handled properly
    #if (missing(fr)) {
    ## reconstitute frame
    #}
    ## pull necessary arguments for making the model frame out of ...
    p <- ncol(X) # maybe also do rank check on X here??
    rho <- new.env(parent=parent.env(environment()))
    rho$pp <- do.call(merPredDLwts$new, c(reTrms[c("Zt","theta","Lambdat","Lind", "Lwts")],
                                      n=nrow(X), list(X=X)))
    REMLpass <- if(REML) p else 0L
    if(missing(fr)) rho$resp <- mkRespMod(REML = REMLpass, ...)
    else rho$resp <- mkRespMod(fr, REML = REMLpass)
    ## note: REML does double duty as rank of X and a flag for using
    ## REML maybe this should be mentioned in the help file for
    ## mkRespMod??  currently that help file says REML is logical.  a
    ## consequence of this double duty is that it is impossible to fit
    ## a model with no fixed effects using REML.
    devfun <- mkdevfun(rho, 0L, verbose, control)
    theta <- getStart(start,reTrms$lower,rho$pp)
    if (length(rho$resp$y) > 0)  ## only if non-trivial y
        devfun(rho$pp$theta) # one evaluation to ensure all values are set
    rho$lower <- reTrms$lower # SCW:  in order to be more consistent with mkLmerDevfunLwts
    return(devfun) # this should pass the rho environment implicitly
}

##' @rdname modular
##' @param control a list giving (for \code{[g]lFormulaLwts}) all options (see \code{\link{lmerControl}} for running the model;
##' (for \code{mkLmerDevfun,mkGlmerDevfun}) options for inner optimization step;
##' (for \code{optimizeLmer} and \code{optimize[Glmer}) control parameters for nonlinear optimizer (typically inherited from the \dots argument to \code{lmerControl})
##' @return \bold{lFormulaLwts, glFormulaLwts}: A list containing components,
##' \item{fr}{model frame}
##' \item{X}{fixed-effect design matrix}
##' \item{reTrms}{list containing information on random effects structure: result of \code{\link{mkReTrms}}}
##' \item{REML}{(lFormulaLwts only): logical flag: use restricted maximum likelihood? (Copy of argument.)}
##' @importFrom Matrix rankMatrix
##' @export
lFormulaLwts <- function(formula, data=NULL, REML = TRUE,
                     subset, weights, na.action, offset, contrasts = NULL,
                     control=lmerControl(), ...)
{
    control <- control$checkControl ## this is all we really need
    mf <- mc <- match.call()

    ignoreArgs <- c("start","verbose","devFunOnly","control")
    l... <- list(...)
    l... <- l...[!names(l...) %in% ignoreArgs]
    do.call("checkArgs",c(list("lmer"),l...))
    if (!is.null(list(...)[["family"]])) {
        ## lmer(...,family=...); warning issued within checkArgs
        mc[[1]] <- quote(lme4::glFormulaLwts)
        if (missing(control)) mc[["control"]] <- glmerControl()
        return(eval(mc, parent.frame()))
    }

    cstr <- "check.formula.LHS"
    checkCtrlLevels(cstr,control[[cstr]])
    denv <- checkFormulaData(formula,data,checkLHS=(control$check.formula.LHS=="stop"))
    #mc$formula <- formula <- as.formula(formula,env=denv) ## substitute evaluated call
    formula <- as.formula(formula,env=denv)
    ## as.formula ONLY sets environment if not already explicitly set ...
    ## ?? environment(formula) <- denv
    # get rid of || terms so update() works as expected
    RHSForm(formula) <- expandDoubleVerts(RHSForm(formula))
    mc$formula <- formula

    m <- match(c("data", "subset", "weights", "na.action", "offset"),
               names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    fr.form <- subbars(formula) # substitute "|" by "+"
    environment(fr.form) <- environment(formula)
    ## (DRY! copied from glFormulaLwts)
    ## model.frame.default looks for these objects in the environment
    ## of the *formula* (see 'extras', which is anything passed in ...),
    ## so they have to be put there ...
    for (i in c("weights", "offset")) {
        if (!eval(bquote(missing(x=.(i)))))
            assign(i,get(i,parent.frame()),environment(fr.form))
    }
    mf$formula <- fr.form
    fr <- eval(mf, parent.frame())
    ## store full, original formula & offset
    attr(fr,"formula") <- formula
    attr(fr,"offset") <- mf$offset
    n <- nrow(fr)
    ## random effects and terms modules
    reTrms <- mkReTrms(findbars(RHSForm(formula)), fr)
    checkNlevels(reTrms$flist, n=n, control)
    checkZdims(reTrms$Ztlist, n=n, control, allow.n=FALSE)
    if (any(is.na(reTrms$Zt))) {
        stop("NA in Z (random-effects model matrix): ",
             "please use ",
             shQuote("na.action='na.omit'"),
             " or ",
             shQuote("na.action='na.exclude'"))
    }
    checkZrank(reTrms$Zt, n=n, control, nonSmall = 1e6)

    ## fixed-effects model matrix X - remove random effect parts from formula:
    fixedform <- formula
    RHSForm(fixedform) <- if(is.null(nb <- nobars(RHSForm(fixedform)))) 1 else nb
    mf$formula <- fixedform
    ## re-evaluate model frame to extract predvars component
    fixedfr <- eval(mf, parent.frame())
    attr(attr(fr,"terms"),"predvars.fixed") <-
        attr(attr(fixedfr,"terms"),"predvars")
    X <- model.matrix(fixedform, fr, contrasts)#, sparse = FALSE, row.names = FALSE) ## sparseX not yet
    ## backward compatibility (keep no longer than ~2015):
    if(is.null(rankX.chk <- control[["check.rankX"]]))
        rankX.chk <- eval(formals(lmerControl)[["check.rankX"]])[[1]]
    X <- chkRank.drop.cols(X, kind=rankX.chk, tol = 1e-7)
    if(is.null(scaleX.chk <- control[["check.scaleX"]]))
        scaleX.chk <- eval(formals(lmerControl)[["check.scaleX"]])[[1]]
    X <- checkScaleX(X, kind=scaleX.chk)
    
    reTrms$Lwts <- reTrms$Lind
    list(fr = fr, X = X, reTrms = reTrms, REML = REML, formula = formula)
}
