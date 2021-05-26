/*
 * FIt.hpp
 *
 *  Created on: May 26, 2014
 *      Author: fabio
 */

#ifndef FIT_HPP_
#define FIT_HPP_

#include "GLM.hpp"
//#include "../Profile.h"
#include <string.h>
#include <vector>
#include <iostream>

#include <gsl/gsl_multifit.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_linalg.h>

class Fit {

private:
	void compute_variance(gsl_vector *w) {
		if(VB != 0)
			gsl_matrix_free(VB);

		VB = gsl_matrix_calloc(p, p);
		gsl_matrix * W = gsl_matrix_calloc(n, n);
		//std::cout << "Compute variance alloc W.\n";
		for(size_t i = 0; i < n; ++i)
			gsl_matrix_set(W, i, i, gsl_vector_get(w, i));

		gsl_matrix * t1 = gsl_matrix_calloc(p, n);
		gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, X, W, 0, t1);

		gsl_matrix * t2 = gsl_matrix_calloc(p, p);
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, t1, X, 0, t2);

		// invert t2
		int ss;
		gsl_permutation * pp = gsl_permutation_alloc(p);
		gsl_linalg_LU_decomp(t2, pp, &ss);
		gsl_linalg_LU_invert(t2, pp, VB);

		gsl_matrix_scale(VB, psi); // cf. quasi-likelihood

		gsl_matrix_free(W);
		gsl_matrix_free(t1);
		gsl_matrix_free(t2);
		gsl_permutation_free(pp);

	}

public:


	Fit() {
		link = new GLM();
		link->setQuasi(true);
		bv = 0;
		VB = 0;
	}

	Fit(const char * link_type) {
		if(strcmp(link_type,"log-link") == 0)
			link = new GLM();
		link->setQuasi(true);
		bv = 0;
		VB = 0;
	}

	~Fit() {
		delete link;

		if(free_data){
			gsl_vector_free(y);
			gsl_matrix_free(X);
			gsl_vector_free(offset);
		}

		if(bv !=0)
			gsl_vector_free(bv);

		if(VB != 0)
			gsl_matrix_free(VB);
	}

	void load_data(const std::vector<double>& yv, const std::vector<std::vector<double> >& Xv,
			const std::vector<double>& offv) {
		free_data = true;

		n = yv.size();
		p = 1 + Xv.size();

		y = gsl_vector_calloc(n);
		X = gsl_matrix_calloc(n, p);
		offset = gsl_vector_calloc(n);

		for(size_t i = 0; i < n; ++i){
			gsl_vector_set(y, i, yv[i]);
			gsl_matrix_set(X, i, 0, 1.0); // intercept
			for(size_t j = 1; j < p; ++j)
				gsl_matrix_set(X, i, j, Xv[j-1][i]);
		}
		if(! offv.empty())
			for(size_t i = 0; i < n; ++i)
				gsl_vector_set(offset, i, offv[i]);
	}

	void set_data(gsl_vector * yv, gsl_matrix * Xv, gsl_vector * offv) {
		free_data = false;

		n = yv->size;
		p = Xv->size2;

		y = yv;
		X = Xv;
		if(offv != NULL)
			offset = offv;
		else{
			offset = gsl_vector_calloc(n);
		}
	}

	void set_data(const std::vector<int>& yv, const std::vector< std::vector<double> >& Xv, const std::vector<double>& offv) {
		free_data = true;

		n = yv.size();
		p = Xv.size() + 1;

		gsl_vector *t_y = gsl_vector_alloc(n);
		gsl_matrix *t_X = gsl_matrix_alloc(n,p);

		for(size_t i = 0; i < n; ++i){
			gsl_vector_set(t_y, i, yv[i]);
			gsl_matrix_set(t_X, i, 0, 1.0); // intercept
			for(size_t j = 1; j < p; ++j)
				gsl_matrix_set(t_X, i, j, Xv[j-1][i]);
		}

		y = t_y;
		X = t_X;

		if(!offv.empty()) {
			offset = gsl_vector_calloc(n);
			for(size_t i=0; i<n; ++i)
				gsl_vector_set(offset, i, offv[i]);
		} else {
			offset = gsl_vector_calloc(n);
		}

#ifdef VERBOSE
		std::cout << "gsl: Data Structures constructed:\n\tn = " << n << ", p = " << p << std::endl;
#endif
	}

	void set_data(int* yv, size_t s_yv, double** Xv, size_t s_Xv, double* offv) {
		free_data = true;

		n = s_yv;
		p = s_Xv + 1;

		gsl_matrix *t_X = gsl_matrix_alloc(n,p);
		gsl_vector *t_y = gsl_vector_alloc(n); // (gsl_vector*) ::malloc(n*sizeof(gsl_vector));

		for(size_t i = 0; i < n; ++i) {
			gsl_vector_set(t_y, i, yv[i]);
			gsl_matrix_set(t_X, i, 0, 1.0); // intercept
			for(size_t j = 1; j < p; ++j)
				gsl_matrix_set(t_X, i, j, Xv[j-1][i]);
		}

		y = t_y;
		X = t_X;

		if(offv != NULL) {
			offset = gsl_vector_alloc(n);
			for(size_t i=0; i<n; ++i)
				gsl_vector_set(offset, i, offv[i]);
		} else {
			offset = gsl_vector_calloc(n);
		}

#ifdef VERBOSE
		std::cout << "gsl: Data Structures constructed:\n\tn = " << n << ", p = " << p << std::endl;
#endif
	}

	void fit_model() {
		gsl_vector *mv = gsl_vector_alloc(n);
		link->init_mv(y, mv);

		double old_chisq = -1.0, chisq;
		gsl_vector *z = gsl_vector_alloc(n);
		gsl_vector *w = gsl_vector_alloc(n);
		bv = gsl_vector_alloc(p);
		gsl_matrix *cov = gsl_matrix_alloc(p, p);

		link->compute_z(y, mv, offset, z); // here apply link-function
		link->compute_weights(mv, w);

		gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc (n, p);

		while(true) {

			// weighted least square fitting
			gsl_multifit_wlinear_svd(X, w, z, GSL_DBL_EPSILON, &rank, bv, cov, &chisq, work);

			if(fabs(chisq - old_chisq) < 1e-6) { // check convergence
				psi = link->compute_dispersion(y, X, bv, offset, mv, rank, link->checkQuasi());

#ifdef PRINTFIT
				compute_variance(w);
				std::cout << "variance computed.\n";
#endif
				break;
			}
			old_chisq = chisq;
			link->compute_mv(bv, X, offset, mv);
		}

		gsl_vector_free(mv);
		gsl_vector_free(z);
		gsl_vector_free(w);
		gsl_matrix_free(cov);
		gsl_multifit_linear_free(work);
	}


	void get_coef(std::vector<double> &cf) {
		cf.reserve(p);
		for(size_t i=0; i < p; ++i)
			cf[i] = gsl_vector_get(bv, i);
	}

	void get_stderr(std::vector<double> &sev) {
		sev.reserve(p);
		for(size_t i = 0; i < p; ++i)
			sev[i] = sqrt(gsl_matrix_get(VB, i, i));
	}

	size_t get_rank_X() 	const { return rank; };
	double get_dispersion() const { return psi; };

	inline void setQuasi(bool q) { link->setQuasi(q); }
	inline bool checkQuasi() { return link->checkQuasi(); }

private:
	GLM * link;
	gsl_vector * y;      // response vector - data count
	gsl_matrix * X;      // design matrix
	gsl_vector * offset; // offset vector
	bool free_data;      // depends on load_data() or set_data()

	size_t n;            // sample size;
	size_t p;            // number of parameters (including intercept: +1) [intercept, len_vec, gcc]
	size_t rank;         // of X (useful for p-values)

	gsl_vector * bv;     // vector of estimated effect sizes
	gsl_matrix * VB;     // covariance matrix of estimated effect sizes
	double psi;          // dispersion
};

#endif /* FIT_HPP_ */
