/*
 * samdata.hpp
 *
 *  Created on: Apr 14, 2014
 *      Author: fabio
 */

#ifndef SAMDATA_HPP_
#define SAMDATA_HPP_

#include "common.hpp"
#include "SamDataT.hpp"

//using namespace ff;

class SamData {
	friend class Edge;

private:
	std::string chr1, chr2, seq1, entry; // cols 1, 3, 7, 10, rest...
	size_t hs1, hs2;
	long id;
	uint_64 start1, start2; // cols 4, 8

public:

	SamData() : hs1(0), hs2(0), id(-1), start1(0), start2(0) { }
	SamData(const SamData &s) { *this = s; }
	SamData(SamData&& o) : chr1(std::move(o.chr1)), chr2(std::move(o.chr2)),
			seq1(std::move(o.seq1)), entry(std::move(o.entry)), hs1(std::move(o.hs1)),
			hs2(std::move(o.hs2)), id(o.getId()), start1(std::move(o.start1)),
			start2(std::move(o.start2)) { }

	SamData& operator=(const SamData& s) = default;
	SamData& operator=(SamData&& s) = default;

	~SamData() { }

	// getters
	inline uint_64 getStart1() 		const { return start1; }
	inline uint_64 getStart2() 		const { return start2; }
	inline std::string getSeq() 	const { return seq1; }
	inline std::string getChr1()	const { return chr1; }
	inline const char* getCChr1()	const { return chr1.c_str(); }
	inline std::string getChr2()	const { return chr2; }
	inline const char* getCChr2()	const { return chr2.c_str(); }
	inline std::string getEntry() 	const { return entry; }
	inline size_t getHS1() 			const { return hs1; }
	inline size_t getHS2()	 		const { return hs2; }
	inline long getId()				const {return id; }

	// setters
	void setEntry(std::string en) 	{ entry  = en; }
	void setChr1(std::string chr1_) { chr1 = chr1_; }
	void setStart1(uint_64 st) 		{ start1 = st; }
	void setChr2(std::string chr2_) { chr2 = chr2_; }
	void setStart2(uint_64 st) 		{ start2 = st; }
	void setHS1(size_t h)			{ hs1 = h; }
	void setHS2(size_t h)			{ hs2 = h; }
	void setId(long l)				{ id = l; }
	void setSeq(std::string s)		{ seq1 = s; }

	SamDataT* toSamDataT() {
		SamDataT *st = new SamDataT(hs1, hs2, id, start1, start2);
		return st;
	}

	// order sam data by chr1 - should consider start1 as well
	inline bool operator<(const SamData &sd) const {
		return (hs1 < sd.hs1) || (hs1 == sd.hs1 && start1 < sd.getStart1());
	}

//	friend std::ostream& operator<<(std::ostream &out, const SamData &sd) {
//		out << sd.rowname 	<< " | " <<
//				sd.chr1 	<< " | " <<
//				sd.start1 	<< " | " <<
//				sd.chr2 	<< " | " <<
//				sd.start2 	<< " | " <<
//				sd.seq1 	<< " | ";
//
//		return out;
//	}

	// needed to write an ordered SAM file, so that the parser does not
	// have to be modified. may be useful again
	friend std::ostream& operator<<(std::ostream &out, const SamData &sd) {
		out << sd.entry 	<< "\t" <<
				sd.chr1 	<< "\t" <<
				sd.start1 	<< "\t" <<
				sd.chr2 	<< "\t" <<
				sd.start2 	<< "\t" <<
				sd.seq1;

		return out;
	}
};


#endif /* SAMDATA_HPP_ */
