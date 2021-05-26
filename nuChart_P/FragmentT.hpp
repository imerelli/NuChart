/*
 * FragmentT.hpp
 *
 *  Created on: Apr 16, 2014
 *      Author: fabio
 */

#ifndef FRAGMENTT_HPP_
#define FRAGMENTT_HPP_

#include "common.hpp"

class FragmentT {
	friend class Gene;
	friend class Edge;

public:
	FragmentT() :
			start_p(0), stop_p(0)
	{ } // frag_num(0), , hs(0)

	FragmentT(uint_64 st, uint_64 sp) :
		start_p(st), stop_p(sp)	{ } // frag_num(frnum), , hs(h)

	// FragmentT(const FragmentT& fg) = default;

	// // move constructor
	// FragmentT(FragmentT&& fg) :
	// 		start_p(std::move(fg.start_p)), stop_p(std::move(fg.stop_p))
	// { } // frag_num(std::move(fg.frag_num)), , hs(std::move(fg.hs)

	// FragmentT& operator=(const FragmentT& fg) = default;

	// FragmentT& operator=(FragmentT&& fg) {
	// 	if (&fg != this) {
	// 		start_p = std::move(fg.start_p);
	// 		stop_p = std::move(fg.stop_p);
	// 	}
	// 	return *this;
	// }

	// getters
	inline uint_64 getStart() const {
		return start_p;
	}
	inline uint_64 getStop() const {
		return stop_p;
	}

private:
	uint_64 start_p;
	uint_64 stop_p;
};

#endif /* FragmentT_HPP_ */
