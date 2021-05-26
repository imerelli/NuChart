/*
 * samdatat.hpp
 *
 *  Created on: Apr 14, 2014
 *      Author: fabio
 */

#ifndef SAMDATAT_HPP_
#define SAMDATAT_HPP_

#include "common.hpp"

//using namespace ff;

class SamDataT {
	friend class EdgeT;

private:
	size_t hs1, hs2;
	long id;
	uint_64 start1, start2; // cols 4, 8

public:

	SamDataT() : hs1(0), hs2(0), id(-1), start1(0), start2(0) { }
	SamDataT(size_t _hs1, size_t _hs2, long _id, uint_64 st1, uint_64 st2) :
		hs1(_hs1), hs2(_hs2), id(_id), start1(st1), start2(st2) { }
	// SamDataT(const SamDataT &s) { *this = s; }
	// SamDataT(SamDataT&& o) : hs1(std::move(o.hs1)), hs2(std::move(o.hs2)),
	// 		id(o.id), start1(std::move(o.start1)), start2(std::move(o.start2)) { }

	// SamDataT& operator=(const SamDataT& s) = default;
	// SamDataT& operator=(SamDataT&& s) = default;

	// getters
	inline uint_64 getStart1() 	const { return start1; }
	inline uint_64 getStart2() 	const { return start2; }
	inline size_t getHS1() 		const { return hs1; }
	inline size_t getHS2() 		const { return hs2; }
	inline long getId()			const { return id; }

	// order sam data by chr1 - should consider start1 as well
	// inline bool operator<(const SamDataT &sd) const {
	// 	return (hs1 < sd.hs1) || (hs1 == sd.hs1 && start1 < sd.getStart1());
	// }


	// needed to write an ordered SAM file, so that the parser does not
	// have to be modified. may be useful again
	friend std::ostream& operator<<(std::ostream &out, const SamDataT &sd) {
	 	out << 	//sd.chr1 	<< "\t" <<
	 			sd.start1 	<< "\t" <<
	 			//sd.chr2 	<< "\t" <<
	 			sd.start2 	<< "\t" ;
	 			//.seq1;

	 	return out;
	}
};


#endif /* SAMDATAT_HPP_ */
