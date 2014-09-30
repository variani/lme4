
#include "predModuleExt.h"

namespace lme4 {
    using    Rcpp::as;

    using     std::invalid_argument;
    using     std::runtime_error;

    using   Eigen::ArrayXd;

    typedef Eigen::Map<MatrixXd>  MMat;
    typedef Eigen::Map<VectorXd>  MVec;
    typedef Eigen::Map<VectorXi>  MiVec;

    merPredDLwts::merPredDLwts(SEXP X, SEXP Lambdat, SEXP LamtUt, SEXP Lind,
                       SEXP Lwts,
                       SEXP RZX, SEXP Ut, SEXP Utr, SEXP V, SEXP VtV,
                       SEXP Vtr, SEXP Xwts, SEXP Zt, SEXP beta0,
                       SEXP delb, SEXP delu, SEXP theta, SEXP u0)
        : merPredD(X, Lambdat, LamtUt, Lind, RZX, Ut, Utr, V, VtV,
                      Vtr, Xwts, Zt, beta0, delb, delu, theta, u0),
          d_Lwts(as<MVec>(Lwts))
    {
    }
}
