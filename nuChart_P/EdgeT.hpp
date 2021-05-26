/*
 * EdgeT.hpp
 *
 *  Created on: Jul 18, 2014
 *      Author: fabio
 */

#ifndef EDGET_HPP_
#define EDGET_HPP_

#include "common.hpp"
#include "Edg.hpp"
#include "GeneT.hpp"
#include "SamDataT.hpp"

struct miniEdge {
	miniEdge(size_t _HS1, size_t _HS2, uint_64 _st1, uint_64 _st2, bool intr, double w) :
		HS1(_HS1), HS2(_HS2), st1(_st1), st2(_st2),
		inter(intr), score(w), prob(0.0) { }

	size_t HS1, HS2;
	uint_64 st1, st2;
	bool inter;
	double score, prob;
};

/*
 * Class representing an Edge of the graph
 */
class EdgeT : public Edg {

public:
	EdgeT() :
			Edg(), conn1(NULL), id1(-1), id2(-1), hs1(0), st1(0), end1(0), hs2(0), st2(0), end2(0),
			hashed(0), /*mirror(0),*/ inter(false) {	}
	EdgeT(GeneT *geneA, GeneT *geneB, SamDataT *connA, /*SamDataT *connB,*/ uint_64 stA,
			uint_64 endA, uint_64 stB, uint_64 endB, bool intr) {
		Edg::v1 = new Vertex(*geneA);
		Edg::v2 = new Vertex(*geneB);
		conn1 = new SamDataT(*connA);
		setStart1 (stA);
		setEnd1 (endA);
		setStart2 (stB);
		setEnd2(endB);
		if(intr)
			setIntergenic();

		hashed = pairZ(geneA->getId(), geneB->getId());
		//mirror = hashThemAll(geneB, geneA);
		id1 = geneA->getId();
		id2 = geneB->getId();
	}

	EdgeT(const EdgeT& e) : Edg(e), conn1(new SamDataT(*e.conn1)), id1(e.getVertex1()->getId()), id2(e.getVertex2()->getId()),
			hs1(e.getHS1()), st1(e.getStart1()),
			end1(e.getEnd1()), hs2(e.getHS2()), st2(e.getStart2()), end2(e.getEnd2()),
			hashed(e.getHashedValue()), /*mirror(e.getMirrorValue()),*/ inter(e.hasIntergenic()) { }
	EdgeT(EdgeT&& e) : Edg(std::move(e)), id1(e.getVertex1()->getId()), id2(e.getVertex2()->getId()), hs1(e.getHS1()), st1(e.getStart1()),
			end1(e.getEnd1()), hs2(e.getHS2()), st2(e.getStart2()), end2(e.getEnd2()),
			hashed(e.getHashedValue()), /*mirror(e.getMirrorValue()),*/ inter(e.hasIntergenic()) {
		conn1 = e.conn1;
		e.conn1 = NULL;
	}

	EdgeT& operator=(const EdgeT& e) {
		if(this != &e) {
			Edg::operator=(e);
			if(conn1) delete conn1;
			conn1 = new SamDataT(*e.conn1);
			hs1 = e.getHS1();
			st1 = e.getStart1();
			end1 = e.getEnd1();
			hs2 = e.getHS2();
			st2 = e.getStart2();
			end2 = e.getEnd2();
			hashed = e.getHashedValue();
			inter = e.hasIntergenic();
			//mirror = e.getMirrorValue();
			id1 = e.getVertex1()->getId();
			id2 = e.getVertex2()->getId();
		}
		return *this;
	}
	EdgeT& operator=(EdgeT&& e) {
		if(this != &e) {
			Edg::operator=(std::move(e));
			if(conn1) delete conn1;
			conn1 = e.conn1;
			hs1 = e.getHS1();
			st1 = e.getStart1();
			end1 = e.getEnd1();
			hs2 = e.getHS2();
			st2 = e.getStart2();
			end2 = e.getEnd2();
			hashed = e.getHashedValue();
			inter = e.hasIntergenic();
			//mirror = e.getMirrorValue();
		}
		return *this;
	}

	~EdgeT() {
		//if(conn1) delete conn1;
	}

	// setters
	inline void setHS1(uint_64 h) 		{ hs1 = h; }
	inline void setStart1(uint_64 v) 	{ st1 = v; }
	inline void setEnd1(uint_64 v) 		{ end1 = v; }
	inline void setHS2(uint_64 h) 		{ hs2 = h; }
	inline void setStart2(uint_64 v) 	{ st2 = v; }
	inline void setEnd2(uint_64 v) 		{ end2 = v; }
	inline void setHashValue(uint_64 v) { hashed = v; }
	inline void setIntergenic()			{
		inter = true;
		setWeight(-1.1);
	}

	//getters
	inline uint_64 getStart1() 			const { return st1; }
	inline uint_64 getEnd1() 			const { return end1; }
	inline uint_64 getStart2() 			const { return st2; }
	inline uint_64 getEnd2() 			const { return end2; }
	inline uint_64 getHS1() 			const { return hs1; }
	inline uint_64 getHS2() 			const {	return hs2; }
	inline uint_64 getHashedValue() 	const {	return hashed; }
	inline bool hasIntergenic()			const { return inter; }
	//inline uint_64 getMirrorValue() 	const { return mirror; }

	inline miniEdge toMiniEdge() {
		miniEdge me(hs1, hs2, st1, st2, inter, getWeight());
		return me;
	}

	bool operator==(const EdgeT& e) const {
		if( (id1 == e.id1 && id2 == e.id2) ||
				(id1 == e.id2 && id2 == e.id1) )
			return  true;
		else return false;
	}
	bool operator<(const EdgeT& e) const { return hashed < e.getHashedValue(); }

	int fill() {
		setHS1(conn1->getHS1());
		setHS2(conn1->getHS2());

#ifdef TEST3
		if( getVertex2()->getGraphRoot() != getVertex1()->getGraphRoot() ) {
			//std::cout << "TOUCH" << std::endl;
			return -1; // nodes belong to different roots, i.e. different graphs
		}
#endif

		return 0;
	}

	SamDataT *conn1;
	long id1, id2;

private:

	uint_64 hashThemAll(Vertex *g1, Vertex *g2) {
		uint_64 res = (uint_64) (g1->getId()) + 2 * (uint_64) (g2->getId());
		return res;
	}

	uint_64 hs1;
	uint_64 st1;
	uint_64 end1;
	uint_64 hs2;
	uint_64 st2;
	uint_64 end2;
	uint_64 hashed;	// hashed value for identifying EdgeT , mirror
	bool inter;
};

namespace {
bool sortMiniEdge(miniEdge e1, miniEdge e2) { return e1.prob > e2.prob; }
bool sortEdgeT(EdgeT *e1, EdgeT *e2) 		{ return e1->getProb() > e2->getProb(); }
}

#endif /* EdgeT_HPP_ */
