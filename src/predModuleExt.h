// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// This file is part of lme4.

#ifndef LME4_PREDMODULEEXT_H
#define LME4_PREDMODULEEXT_H

#include <RcppEigen.h>
#include "lme4CholmodDecomposition.h"
#include "predModule.h"

namespace lme4 {

    using Eigen::ArrayXd;
    using Eigen::LLT;
    using Eigen::MatrixXd;
    using Eigen::VectorXd;
    using Eigen::VectorXi;

    class merPredDLwts : public merPredD {

    protected:
	MVec         d_Lwts;

    public:	
	merPredDLwts(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, 
                 SEXP, 
		 SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

	void             setTheta(const VectorXd&);
    };
}

#endif // LME4_PREDMODULEEXT_H
