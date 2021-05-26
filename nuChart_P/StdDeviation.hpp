/*
 * StdDeviation.hpp
 *
 *  Created on: May 7, 2015
 *      Author: fabio
 */

#ifndef STDDEVIATION_HPP_
#define STDDEVIATION_HPP_

#include <cmath>
#include <cstdlib>

template <class T>
class StdDeviation {
private:
	double mean, var, stdev;
	size_t dim;
	T *vals;

	void meanVar() {
		double sum=0;
		for(unsigned i=0; i<dim; ++i)
			sum += vals[i];
		mean = sum/dim;

		sum=0;
		for(unsigned i=0; i<dim; ++i)
			sum += (vals[i]-mean)*(vals[i]-mean);
		var = sum/dim;

		stdev = sqrt(var);
	}

public:
	StdDeviation(size_t size, T *vals_) : mean(0.0), var(0.0), stdev(0.0), dim(size) {
		vals = new T[dim]();
		memcpy(vals, vals_, dim);
		meanVar();
	}
	StdDeviation(std::vector<T>& vals_) : mean(0.0), var(0.0), stdev(0.0), dim(vals_.size()) {
		vals = new T[dim]();
		for(unsigned i=0; i<dim; ++i)
			vals[i] = vals_[i];
		meanVar();
	}

	inline double getMeanValue()	const { return mean; }
	inline double getVariance() 	const { return var;  }
	inline double getStdDev()		const { return stdev;}
	inline size_t getDim()			const { return dim;  }
	inline T getValAt(size_t p)		const { return p < dim ? vals[p] : -1; }
};


template <class C>
class Covariance {
private:
	StdDeviation<C> *x, *y;
	double cov, corr;

	inline void xyCov() {
		double xMean = x->getMeanValue(), yMean = y->getMeanValue(), tot=0;
		size_t dim = x->getDim();

		assert(x->getDim() == y->getDim());

		for(unsigned i=0; i<dim; ++i)
			tot += (x->getValAt(i)-xMean)*(y->getValAt(i)-yMean);
		cov = tot/dim;

		corr = cov / (x->getStdDev()*y->getStdDev());
	}

public:
	Covariance(std::vector<C> s1, std::vector<C> s2) : cov(0.0), corr(0.0) {
		x = new StdDeviation<C>(s1);
		y = new StdDeviation<C>(s2);
		xyCov();
	}
	Covariance(size_t sz1, C *s1, size_t sz2, C *s2) : cov(0.0), corr(0.0) {
		x = new StdDeviation<C>(sz1, s1);
		y = new StdDeviation<C>(sz2, s2);
		xyCov();
	}

	inline double getCovariance() 	const { return cov;  }
	inline double getCorrelation() 	const { return corr; }
};

#endif /* STDDEVIATION_HPP_ */
