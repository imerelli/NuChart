/*
 * GeneT.hpp
 *
 *  Created on: Apr 16, 2014
 *      Author: fabio
 */

#ifndef GeneT_HPP_
#define GeneT_HPP_

#include "common.hpp"
#include "Vertex.hpp"
#include "SamDataT.hpp"

/**
 * This class describes a Gene.
 */
class GeneT : public Vertex {
	friend class Finder;

public:

	/**
	 * Default constructor
	 */
	GeneT() :
		Vertex(), start_p(0), stop_p(0), hs(0), conn(NULL), cis(0), trans(0) {
	}
	// GeneT(uint_64 st, uint_64 sp, int lv, SamDataT *c) :
	// 	Vertex(), start_p(st), stop_p(sp), conn(c) {
	// 	Vertex::level = lv;
	// }
  GeneT(uint_64 st, uint_64 sp, long eid, int lv, size_t h, uint_64 est, uint_64 esp,
	long pos, int cis_, int trans_) :
    Vertex(),start_p(st), stop_p(sp), exst(est), exsp(esp), hs(h), conn(NULL), cis(cis_), trans(trans_) {
    setId(pos);
    Vertex::level = lv;
  }
	// GeneT(const GeneT& g) : Vertex(g) {
	// 		start_p = g.start_p;
	// 		stop_p = g.stop_p;
	// 		exsp = g.exsp;
	// 		exst = g.exst;
	// 		//entrez_ID = g.entrez_ID;
	// 		hs = g.hs;
	// 		conn = g.getConnection();
	// 		cis = g.cis;
	// 		trans = g.trans;
	// }
	// GeneT(GeneT&& g) :
	// 	Vertex(std::move(g)), start_p(std::move(g.start_p)), stop_p(std::move(g.stop_p)),
	// 	exst(std::move(g.exst)), exsp(std::move(g.exsp)),
	// 	hs(std::move(g.hs)), conn(NULL) {
	// 	conn = g.getConnection();
	// 	g.setConnection(NULL);
	// 	cis = g.cis;
	// 	trans = g.trans;
	// }

	// GeneT& operator=(const GeneT& g) {
	// 	if (this != &g) {
	// 		Vertex::operator =(g);
	// 		start_p = g.start_p;
	// 		stop_p = g.stop_p;
	// 		exsp = g.exsp;
	// 		exst = g.exst;
	// 		hs = g.hs;
	// 		conn = g.getConnection();
	// 		cis = g.cis;
	// 		trans = g.trans;
	// 	}
	// 	return *this;
	// }
	// GeneT& operator=(GeneT&& g) {
	// 	if (this != &g) {
	// 		Vertex::operator =(std::move(g));
	// 		start_p = std::move(g.start_p);
	// 		stop_p = std::move(g.stop_p);
	// 		exsp = std::move(g.exsp);
	// 		exst = std::move(g.exst);
	// 		hs = std::move(g.hs);
	// 		conn = g.getConnection();
	// 		g.setConnection(NULL);
	// 		cis = g.cis;
	// 		trans = g.trans;
	// 	}
	// 	return *this;
	// }

	// getters
	//inline long getID() 				const { return entrez_ID; }
	inline uint_64 getStart() 			const { return start_p; }
	inline uint_64 getStop() 			const { return stop_p; }
	inline long getPosition() 			const { return Vertex::id; }
	inline size_t getHS() 				const { return hs; }
	inline uint_64 getExtendedStart() 	const { return exst; }
	inline uint_64 getExtendedStop() 	const { return exsp; }
	inline SamDataT* getConnection()	const { return conn; }
	inline long getNumCis() 			const { return cis; }
	inline long getNumTrans() 			const { return trans; }

	// setters
	//void setID(long id)					{ entrez_ID = id; }
	void setStart(uint_64 _start) 		{ start_p = _start;	}
	void setStop(uint_64 _stop) 		{ stop_p = _stop;	}
	void setPosition(long p) 			{ Vertex::id = p; }
	void setHS(size_t h)				{ hs = h; }
	void setExtendedStart(uint_64 st)	{ exst = st; }
	void setExtendedStop(uint_64 sp)	{ exsp = sp; }
	void setConnection(SamDataT* c)		{ conn = c; }
	inline void setCis()				{ ++cis; }
	inline void setTrans()				{ ++trans; }

	// bool operator==(const GeneT& g) const { return getPosition() == g.getPosition(); }
	// bool operator<(const GeneT& g) const {
	// 	return getPosition() < g.getPosition();
	// }


private:
	uint_64 start_p;
	uint_64 stop_p;
	uint_64 exst;
	uint_64 exsp;
	size_t hs;		// hashed chromosome
	SamDataT *conn;
	int cis, trans;
};

#endif /* GeneT_HPP_ */
