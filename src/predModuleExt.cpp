
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
    
    void merPredDLwts::setTheta(const VectorXd& theta) {

        if (theta.size() != d_theta.size()) {
	          Rcpp::Rcout << "(" << theta.size() << "!=" <<
        		d_theta.size() << ")" << std::endl;
      	    // char errstr[100];
      	    // sprintf(errstr,"theta size mismatch (%d != %d)",
      	    //	    theta.size(),d_theta.size());
      	    throw invalid_argument("theta size mismatch");
      	}
      	// update theta
      	std::copy(theta.data(), theta.data() + theta.size(), 
		        d_theta.data());
      	// update Lambdat
      	int    *lipt = d_Lind.data();
      	double *LamX = d_Lambdat.valuePtr(), *thpt = d_theta.data();
        double *lwpt = d_Lwts.data();
      	for (int i = 0; i < d_Lind.size(); ++i) {
      	    LamX[i] = lwpt[i] * thpt[lipt[i] - 1];
      	}
    }
}
