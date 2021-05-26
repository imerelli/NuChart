/*
 * CMap.hpp
 *
 *  Created on: Jul 19, 2014
 *      Author: fabio
 */

#ifndef CMAP_HPP_
#define CMAP_HPP_

#include "common.hpp"

/*
 * Class representing a CMap matrix
 * Very simple matrix class
 */
template <class T>
class CMap {
public:

	CMap() : rows_(0), cols_(0), utr_sz(0), ltr_sz(0), max(0), data_(NULL) { }

	/**
	 * Constructor that creates an empty matrix of size r*c
	 * and initializes its values to zero
	 */
	CMap(unsigned r, unsigned c) : rows_ (r), cols_ (c), utr_sz(0), ltr_sz(0), max(0) {
		if (r < 2 || c < 2) {
			std::cerr << "Matrix cannot have dimensions of size less than 2." << std::endl;
			exit(EXIT_FAILURE);
		}

		//std::cerr << "CALLED CMAP CTOR2" << std::endl;

		data_ = new T[r*c];

		for(unsigned i=0; i<(r*c); ++i)
				data_[i] = 0;
	}

	// Copy constructor
	CMap(const CMap<T>& cm) :
		rows_(cm.rows_), cols_(cm.cols_), max(cm.max) {
		if(data_) delete[] data_;
		data_ = new T[cm.rows_*cm.cols_];
		memcpy(data_, cm.data_, cm.rows_*cm.cols_);

		//std::cerr << "CALLED CMAP CCTOR" << std::endl;
	}

	// Move constructor
	CMap(CMap<T>&& cm) :
		rows_( cm.rows_ ), cols_( cm.cols_ ), max(cm.max), data_( NULL ) {
		data_ = cm.data_;
		cm.rows_ = 0;
		cm.cols_ = 0;
		cm.data_ = NULL;

		//std::cerr << "CALLED CMAP MTOR" << std::endl;
	}

	// Copy assignment operator
	CMap<T>& operator=(const CMap<T>& cm) {
		if(this != &cm) {
			if(data_) delete[] data_;
			rows_ = cm.rows_;
			cols_ = cm.cols_;
			max = cm.max;
			data_ = new T[cm.rows_*cm.cols_];
			memcpy(data_, cm.data_, cm.rows_*cm.cols_);

			//std::cerr << "CALLED CMAP CA_OPR" << std::endl;
		}
		return *this;
	}

	// move assignment operator
	CMap<T>& operator=(CMap<T>&& cm) {
		if(this != &cm) {
			if(data_) delete[] data_;
			rows_ = cm.rows_;
			cols_ = cm.cols_;
			max = cm.max;
			data_ = cm.data_;

			cm.rows_ = 0;
			cm.cols_ = 0;
			//cm.data_= NULL;
			//std::cerr << "CALLED CMAP MV_OPR" << std::endl;
		}
		return *this;
	}

	// Destructor
	~CMap() {
		//std::cerr << "CALLED CMAP DTOR" << std::endl;
		if(data_)
			delete[] data_;
	}

	/**
	 * Allocate memory for the matrix and fill it with 'val' (default is zero)
	 */
	void buildMatrix (unsigned r, unsigned c, T val=0) {
		if (r < 2 || c < 2) {
			std::cerr << "Matrix cannot have dimensions of size less than 2." << std::endl;
			exit(EXIT_FAILURE);
		}

		if(data_ == NULL)
			data_ = new T[r*c];

		for(unsigned i=0; i<(r*c); ++i)
				data_[i] = val;
	}


	/**
	 * Write the upperTriangular matrix in vector 'matx' passed by reference
	 *
	 * if diagonal = 0, diagonal is not considered.
	 */
	void upperTriangular(std::vector<T>& matx, int diagonal=0) {
		utr_sz=0;

		if(data_ == NULL) {
			std::cerr << "No Matrix data." << std::endl;
			exit(EXIT_FAILURE);
		}

		if(diagonal == 0) {
			for(unsigned s=1; s<cols_; ++s)
				utr_sz += s;
			matx.reserve( utr_sz );

			for(unsigned i=0; i<rows_; ++i)
				for(unsigned j=1+i; j<cols_; ++j)
					matx.push_back(data_[cols_*i+j]);
		} else {
			for(unsigned s=1; s<=cols_; ++s)
				utr_sz += s;
			matx.reserve( utr_sz );

			for(unsigned i=0; i<rows_; ++i)
				for(unsigned j=0+i; j<cols_; ++j)
					matx.push_back(data_[cols_*i+j]);
		}
	}

	/**
	 * Returns the upper triangular matrix as a pointer to array of double
	 *
	 * if diagonal = 0, diagonal is not considered.
	 */
	T *upperTriangular(int diagonal=0) {
		T *matx = NULL;
		utr_sz=0;

		if(data_ == NULL) {
			std::cerr << "No Matrix data." << std::endl;
			exit(EXIT_FAILURE);
		}

		if(!diagonal) {
			for(unsigned s=1; s<cols_; ++s)
				utr_sz += s;
			matx = new T[utr_sz];

			unsigned f=0;
			for(unsigned i=0; i<rows_; ++i) {
				for(unsigned j=1+i; j<cols_; ++j) {
					matx[f] = (data_[cols_*i+j]);
					++f;
				}
			}
		} else {
			for(unsigned s=1; s<=cols_; ++s)
				utr_sz += s;
			matx = new T[utr_sz];

			unsigned f=0;
			for(unsigned i=0; i<rows_; ++i) {
				for(unsigned j=0+i; j<cols_; ++j) {
					matx[f] = (data_[cols_*i+j]);
					++f;
				}
			}
		}

		return matx;
	}

	/// return the size of the upper triangular matrix
	unsigned int upperTriangularSize(int diagonal=0) const {
		return utr_sz;
	}

	/// return the max value stored in matrix
	T getMaxValue() {
		for(unsigned i=0; i<(rows_*cols_); ++i)
			if(data_[i] > max)
				max = data_[i];
		return max;
	}

	/**
	 * Calculate the mean value of all values in the matrix
	 */
	double meanVal() {
		double sum=0;
		for(unsigned i=0; i<rows_; ++i)
			for(unsigned j=0; j<cols_; ++j)
				sum += data_[cols_*i+j];
		return (sum/(cols_*rows_));
	}

	/**
	 * Calculate the variance of all values in the matrix
	 */
	double variance(double mn=-1) {
		double mean=0;
		if(mn == -1)
			mean = meanVal();
		else mean = mn;

		double temp=0;
		for(unsigned i=0; i<rows_; ++i)
			for(unsigned j=0; j<cols_; ++j)
				temp += (data_[cols_*i+j] - mean) * (data_[cols_*i+j] - mean);
		return ( temp/((cols_*rows_)-1) );
	}

	/**
	 * Calculate the standard deviation of all values in the matrix
	 */
	double standardDeviation(double mn=-1) {
		if(mn == -1)
			return sqrt(variance());
		else
			return sqrt(variance(mn));
	}

	// same as operator(), useful when CMap<T> object is a pointer to object
	inline T& atPoint(unsigned r, unsigned c) { return data_[cols_*r + c]; }

	/// increment value at 'cmap(r,c)' by 'inc' (default: + 1)
	inline void incPoint(unsigned r, unsigned c, int inc=1) { data_[cols_*r + c] += inc; }

	/// return matrix dimension, as a product of rows x columns
	inline unsigned dim() 			const { return rows_*cols_; }
	/// return the allocated space for the matrix
	inline unsigned size()			const { return (rows_*cols_)*sizeof(T); }

	inline unsigned numRows() 		const { return rows_; }
	inline unsigned numColumns() 	const { return cols_; }
	inline T* asArray() 			const { return data_; }


	inline std::vector<double> asVector() {
		std::vector<double> v;
		v.reserve(rows_*cols_);
		for(unsigned i=0; i<dim(); ++i)
			v.push_back(data_[i]);
		return v;
	}

	void printContacts() {
		for(unsigned i=0; i < rows_; ++i) {
			for(unsigned j=0; j < cols_; ++j) {
				if(data_[cols_*i + j] != 0) {
					std::cout << "R" << i+1 << " : ";
					std::cout << data_[cols_*i + j] << " : C" << j+1 << " ";
				}
			}
			std::cout << std::endl;
		}
	}

	T& operator() (unsigned r, unsigned c) {
		if (r >= rows_ || c >= cols_) {
			std::cerr << "Matrix subscripts out of bounds:\n"
					<< "N rows: " << rows_ << ", r: " << r
					<< "; N cols: " << cols_ << ", c: " << c << std::endl;
			exit(EXIT_FAILURE);
		}
		return data_[cols_*r + c];
	}

	T  operator() (unsigned r, unsigned c) const {
		if (r >= rows_ || c >= cols_) {
			std::cerr << "Matrix subscripts out of bounds:\n"
					<< "N rows: " << rows_ << ", r: " << r
					<< "; N cols: " << cols_ << ", c: " << c << std::endl;
			exit(EXIT_FAILURE);
		}
		return data_[cols_*r + c];
	}

	/**
	 * Friendly print a CMap
	 */
	friend std::ostream& operator<<(std::ostream &out, const CMap& m) {
		for(unsigned i=0; i < m.rows_; ++i) {
			out << "R" << i+1 << ":\t";
			for(unsigned j=0; j < m.cols_; ++j) {
					out << m(i,j) << " ";
			}
			out << std::endl;
		}
		out << std::endl;

		return out;
	}

private:
	unsigned int rows_, cols_;  // zero-based rows and columns number
	unsigned int utr_sz, ltr_sz; // upper- and lower- triangular matrix sizes
	T max;
	T* data_;
};


#endif /* CMAP_HPP_ */
