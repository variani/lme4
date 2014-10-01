
### FIXME
### shouldn't we have "merPred"  with two *sub* classes "merPredDLwts" and "merPredS"
### for the dense and sparse X cases ?

##' Generator object for the \code{\linkS4class{merPredDLwts}} class
##'
##' The generator object for the \code{\linkS4class{merPredDLwts}} reference class.
##' Such an object is primarily used through its \code{new} method.
##'
##' @param ... List of arguments (see Note).
##' @section Methods:
##' \describe{
##'   \item{new(X, Zt, Lambdat, Lind, theta, n):}{Create a new \code{\linkS4class{merPredDLwts}} object}
##' }
##' @note Arguments to the \code{new} methods must be named arguments:
##' \itemize{
##' \item{X}{ dense model matrix for the fixed-effects parameters, to be stored
##'     in the \code{X} field.}
##' \item{Zt}{ transpose of the sparse model matrix for the random effects.  It
##'     is stored in the \code{Zt} field.}
##' \item{Lambdat}{ transpose of the sparse lower triangular relative variance
##'     factor (stored in the \code{Lambdat} field).}
##' \item{Lind}{ integer vector of the same length as the \code{"x"} slot in the
##'     \code{Lambdat} field.  Its elements should be in the range 1 to the length
##'     of the \code{theta} field.}
##' \item{theta}{ numeric vector of variance component parameters (stored in the
##'     \code{theta} field).}
##' \item{n}{ sample size, usually \code{nrow(X)}.}
##' }
##' @seealso \code{\linkS4class{merPredDLwts}}
##' @keywords classes
##' @export
merPredDLwts <-
    setRefClass("merPredDLwts", # Predictor class for mixed-effects models with dense X
                fields =
                list(Lambdat = "dgCMatrix",   # depends: theta, Lind and Lwts
                     LamtUt  = "dgCMatrix",   # depends: Lambdat and Ut
                     Lind    = "integer",     # depends: nothing
                     Lwts    = "numeric",     # depends: nothing                     
                     Ptr     = "externalptr", # depends: 
                     RZX     = "matrix",      # depends: lots
                     Ut      = "dgCMatrix",   # depends: Zt and weights
                     Utr     = "numeric",     # depends: lots
                     V       = "matrix",      # depends: 
                     VtV     = "matrix",
                     Vtr     = "numeric",
                     X       = "matrix",
                     Xwts    = "numeric",
                     Zt      = "dgCMatrix",
                     beta0   = "numeric",
                     delb    = "numeric",
                     delu    = "numeric",
                     theta   = "numeric",
                     u0      = "numeric"),
                methods =
                list(
                     initialize = function(X, Zt, Lambdat, Lind, Lwts, theta, n, ...) {
                         if (!nargs()) return
                         X <<- as(X, "matrix")
                         Zt <<- as(Zt, "dgCMatrix")
                         Lambdat <<- as(Lambdat, "dgCMatrix")
                         Lind <<- as.integer(Lind)
                         Lwts <<- as.numeric(Lwts)
                         theta <<- as.numeric(theta)
                         N <- nrow(X)
                         p <- ncol(X)
                         q <- nrow(Zt)
                         stopifnot(length(theta) > 0L,
                                   length(Lind) > 0L,
                                   length(Lwts) > 0L,
                                   length(Lind) == length(Lwts),
                                   all(sort(unique(Lind)) == seq_along(theta)))
                         RZX <<- array(0, c(q, p))
                         Utr <<- numeric(q)
                         V <<- array(0, c(n, p))
                         VtV <<- array(0, c(p, p))
                         Vtr <<- numeric(p)
                         b0 <- list(...)$beta0
                         beta0 <<- if (is.null(b0)) numeric(p) else b0
                         delb <<- numeric(p)
                         delu <<- numeric(q)
                         uu <- list(...)$u0
                         u0 <<- if (is.null(uu)) numeric(q) else uu
                         Ut <<- if (n == N) Zt + 0 else
                             Zt %*% sparseMatrix(i=seq_len(N), j=as.integer(gl(n, 1, N)), x=rep.int(1,N))
                         ## The following is a kludge to overcome problems when Zt is square
                         ## by making LamtUt rectangular
                         LtUt <- Lambdat %*% Ut
                         ## if (nrow(LtUt) == ncol(LtUt))
                         ##     LtUt <- cbind2(LtUt,
                         ##                    sparseMatrix(i=integer(0),
                         ##                                 j=integer(0),
                         ##                                 x=numeric(0),
                         ##                                 dims=c(nrow(LtUt),1)))
                         LamtUt <<- LtUt
                         Xw <- list(...)$Xwts
                         Xwts <<- if (is.null(Xw)) rep.int(1, N) else as.numeric(Xw)
                         initializePtr()
                     },
                     CcNumer      = function() {
                         'returns the numerator of the orthogonality convergence criterion'
                         .Call(merPredDLwtsCcNumer, ptr())
                     },
                     L            = function() {
                         'returns the current value of the sparse Cholesky factor'
                         .Call(merPredDLwtsL, ptr())
                     },
                     P            = function() {
                         'returns the permutation vector for the sparse Cholesky factor'
                         .Call(merPredDLwtsPvec, ptr())
                     },
                     RX           = function() {
                         'returns the dense downdated Cholesky factor for the fixed-effects parameters'
                         .Call(merPredDLwtsRX, ptr())
                     },
                     RXi          = function() {
                         'returns the inverse of the dense downdated Cholesky factor for the fixed-effects parameters'
                         .Call(merPredDLwtsRXi, ptr())
                     },
                     RXdiag       = function() {
                         'returns the diagonal of the dense downdated Cholesky factor'
                         .Call(merPredDLwtsRXdiag, ptr())
                     },
                     b            = function(fac) {
                         'random effects on original scale for step factor fac'
                         .Call(merPredDLwtsb, ptr(), as.numeric(fac))
                     },
                     beta         = function(fac) {
                         'fixed-effects coefficients for step factor fac'
                         .Call(merPredDLwtsbeta, ptr(), as.numeric(fac))
                     },
                     copy         = function(shallow = FALSE) {
                         def <- .refClassDef
                         selfEnv <- as.environment(.self)
                         vEnv    <- new.env(parent=emptyenv())
                         
                         for (field in setdiff(names(def@fieldClasses), "Ptr")) {
                             if (shallow)
                                 assign(field, get(field, envir = selfEnv), envir = vEnv)
                             else {
                                 current <- get(field, envir = selfEnv)
                                 if (is(current, "envRefClass"))
                                     current <- current$copy(FALSE)
                                 ## hack (https://stat.ethz.ch/pipermail/r-devel/2014-March/068448.html)
                                 ## ... to ensure real copying
                                 ## forceCopy() does **NOT** work here, but +0 does
                                 ## we can get away with this because all fields other than Ptr are numeric
                                 assign(field, current+0, envir = vEnv)
                             }
                         }
                         do.call(merPredDLwts$new, c(as.list(vEnv), n=nrow(vEnv$V), Class=def))
                     },
                     ldL2         = function() {
                         'twice the log determinant of the sparse Cholesky factor'
                         .Call(merPredDLwtsldL2, ptr())
                     },
                     ldRX2        = function() {
                         'twice the log determinant of the downdated dense Cholesky factor'
                         .Call(merPredDLwtsldRX2, ptr())
                     },
                     unsc         = function() {
                         'the unscaled variance-covariance matrix of the fixed-effects parameters'
                         .Call(merPredDLwtsunsc, ptr())
                     },
                     linPred      = function(fac) {
                         'evaluate the linear predictor for step factor fac'
                         .Call(merPredDLwtslinPred, ptr(), as.numeric(fac))
                     },
                     installPars  = function(fac) {
                         'update u0 and beta0 to the values for step factor fac'
                         .Call(merPredDLwtsinstallPars, ptr(), as.numeric(fac))
                     },
                    initializePtr = function() {
                        Ptr <<- .Call(merPredDLwtsCreate, as(X, "matrix"), Lambdat,
                                      LamtUt, Lind, Lwts, RZX, Ut, Utr, V, VtV, Vtr,
                                      Xwts, Zt, beta0, delb, delu, theta, u0)
                        .Call(merPredDLwtssetTheta, Ptr, theta)
                        .Call(merPredDLwtsupdateXwts, Ptr, Xwts)
                        .Call(merPredDLwtsupdateDecomp, Ptr, NULL)
                     },
                     ptr          = function() {
                         'returns the external pointer, regenerating if necessary'
                         if (length(theta)) {
                             if (.Call(isNullExtPtr, Ptr)) initializePtr()
                         }
                         Ptr
                     },
                     setBeta0     = function(beta0) {
                         'install a new value of beta'
                         .Call(merPredDLwtssetBeta0, ptr(), as.numeric(beta0))
                     },
                     setTheta     = function(theta) {
                         'install a new value of theta'
                         .Call(merPredDLwtssetTheta, ptr(), as.numeric(theta))
                     },
                     setZt        = function(ZtNonZero) {
                         'install new values in Zt'
                         .Call(merPredDLwtssetZt, ptr(), as.numeric(ZtNonZero))
                     },
                     solve        = function() {
                         'solve for the coefficient increments delu and delb'
                         .Call(merPredDLwtssolve, ptr())
                     },
                     solveU       = function() {
                         'solve for the coefficient increment delu only (beta is fixed)'
                         .Call(merPredDLwtssolveU, ptr())
                     },
                     setDelu      = function(val) {
                         'set the coefficient increment delu'
                         .Call(merPredDLwtssetDelu , ptr(), as.numeric(val))
                     },
                     setDelb      = function(val) {
                         'set the coefficient increment delb'
                         .Call(merPredDLwtssetDelb , ptr(), as.numeric(val))
                     },
                     sqrL         = function(fac) {
                         'squared length of u0 + fac * delu'
                         .Call(merPredDLwtssqrL, ptr(), as.numeric(fac))
                     },
                     u            = function(fac) {
                         'orthogonal random effects for step factor fac'
                         .Call(merPredDLwtsu, ptr(), as.numeric(fac))
                     },
                     updateDecomp = function(XPenalty = NULL) {
                         'update L, RZX and RX from Ut, Vt and VtV'
                         invisible(.Call(merPredDLwtsupdateDecomp, ptr(), XPenalty))
                     },
                     updateL = function() {
                         'update LamtUt and L'
                         .Call(merPredDLwtsupdateL, ptr())
                     },
                     updateLamtUt = function() {
                         'update LamtUt and L'
                         .Call(merPredDLwtsupdateLamtUt, ptr())
                     },
                     updateRes    = function(wtres) {
                         'update Vtr and Utr using the vector of weighted residuals'
                         .Call(merPredDLwtsupdateRes, ptr(), as.numeric(wtres))
                     },
                     updateXwts   = function(wts) {
                         'update Ut and V from Zt and X using X weights'
                         .Call(merPredDLwtsupdateXwts, ptr(), wts)
                     }
                     )
                )
merPredDLwts$lock("Lambdat", "LamtUt", "Lind", "Lwts", "RZX", "Ut", "Utr", "V", "VtV", "Vtr",
              "X", "Xwts", "Zt", "beta0", "delb", "delu", "theta", "u0")

##' Class \code{"merPredDLwts"} - a dense predictor reference class
##'
##' A reference class for a mixed-effects model predictor module with a dense
##' model matrix for the fixed-effects parameters.  The reference class is
##' associated with a C++ class of the same name.  As is customary, the
##' generator object, \code{\link{merPredDLwts}}, for the class has the same name as
##' the class.
##' @name merPredDLwts-class
##' @note Objects from this reference class correspond to objects in a C++
##'     class.  Methods are invoked on the C++ class object using the external
##'     pointer in the \code{Ptr} field.  When saving such an object the external
##'     pointer is converted to a null pointer, which is why there are redundant
##'     fields containing enough information as R objects to be able to regenerate
##'     the C++ object.  The convention is that a field whose name begins with an
##'     upper-case letter is an R object and the corresponding field, whose name
##'     begins with the lower-case letter is a method.  References to the
##'     external pointer should be through the method, not directly through the
##'     \code{Ptr} field.
##' @section Extends: All reference classes extend and inherit methods from
##'     \code{"\linkS4class{envRefClass}"}.
##' @seealso \code{\link{lmer}}, \code{\link{glmer}}, \code{\link{nlmer}},
##'     \code{\link{merPredDLwts}}, \code{\linkS4class{merMod}}.
##' @keywords classes
##' @examples
##'
##' showClass("merPredDLwts")
##' str(slot(lmer(Yield ~ 1|Batch, Dyestuff), "pp"))
##'
NULL
