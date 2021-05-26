/*
 * GLM.hpp
 *
 *  Created on: May 26, 2014
 *      Author: fabio
 */

#ifndef GLM_HPP_
#define GLM_HPP_

#include "LinkFunc.h"

#include <gsl/gsl_blas.h>
#include <cmath>

class GLM : public LinkFunc {
public:
	GLM() {}
	~GLM() {}



	/**
	 * mv equals y, except that when it is 0, it becomes 0.01
	 */
	void init_mv(gsl_vector *y, gsl_vector * mv) {
		size_t n = y->size;
		for(size_t i = 0; i < n; ++i) {
			if(gsl_vector_get(y, i) == 0)
				gsl_vector_set(mv, i, 0.01);
			else
				gsl_vector_set(mv, i, gsl_vector_get(y, i));
		}
	}

	/**
	 * Apply link function to determine the relationship between the
	 * linear predictor (y ==> mv) and the mean of the distribution function
	 */
	void compute_z(gsl_vector * y, gsl_vector * mv, gsl_vector * offset, gsl_vector * z) {
		size_t n = y->size;
		double mv_i, y_i, val;
		for(size_t i = 0; i < n; ++i) {
			mv_i = gsl_vector_get(mv, i);
			y_i = gsl_vector_get(y, i);
			val = log(mv_i) + (1.0/mv_i)*(y_i-mv_i) - gsl_vector_get(offset, i);
			gsl_vector_set(z, i, val);
		}
	}

	// needed for p-value
	void compute_weights(gsl_vector * mv, gsl_vector * w) {
		size_t n = mv->size;
		for(size_t i = 0; i < n; ++i)
			gsl_vector_set(w, i, gsl_vector_get(mv, i));
	}

	// needed to compute dispersion
	void compute_mv(gsl_vector * bv, gsl_matrix * Xv, gsl_vector * offset, gsl_vector * mv) {
		size_t n = Xv->size1, p = Xv->size2;

		gsl_matrix * B = gsl_matrix_calloc(p, 1);
		for(size_t i = 0; i < p; ++i)
			gsl_matrix_set(B, i, 0, gsl_vector_get(bv, i));

		gsl_matrix * fit = gsl_matrix_calloc(n, 1);
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, Xv, B, 0, fit);

		for(size_t i = 0; i < n; ++i)
			gsl_vector_set(mv, i, exp(gsl_matrix_get(fit, i, 0)
					+ gsl_vector_get(offset, i)));

		gsl_matrix_free(B);
		gsl_matrix_free(fit);
	}


	// weighted total sum of squares (WTSS)
	double compute_dispersion(gsl_vector * y, gsl_matrix * Xv,
			gsl_vector * bv, gsl_vector * offset,
			gsl_vector * mv, double rank, bool quasi_lik) {
		double psi;
		if(! quasi_lik){
			psi = 1.0;
		}
		else{
			compute_mv(bv, Xv, offset, mv);
			double wtss = 0.0;
			for(size_t i = 0; i < y->size; ++i)
				wtss += pow(gsl_vector_get(y,i) - gsl_vector_get(mv,i), 2) /
				gsl_vector_get(mv,i);
			psi = (1/(y->size-rank)) * wtss;
		}
		return psi;
	}

	inline void setQuasi(bool q) { quasi = q; }
	inline bool checkQuasi() { return quasi; }

private:
	bool quasi;
};

#endif /* GLM_HPP_ */
