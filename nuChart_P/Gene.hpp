/*
 * Gene.hpp
 *
 *  Created on: Apr 16, 2014
 *      Author: fabio
 */

#ifndef GENE_HPP_
#define GENE_HPP_

#include "common.hpp"
#include "Vertex.hpp"
#include "GeneT.hpp"

struct Node {
  Node(std::string lbl, long id, long dg, long intrs, size_t hs_, int cis_, int trans_, int l_,
			bool rt, bool it, double gxp) :
    label(lbl), n_id(id), deg(dg), intrs(intrs), cis(cis_), trans(trans_), lvl(l_), hs(hs_),
		isRoot(rt), isInter(it), gExp(gxp) { }
	Node(const Node& t) = default;
    Node(Node&& t) = default;

    Node& operator=(const Node& other) = default;
    Node& operator=(Node&& other) = default;

    bool operator==(const Node& g) const { return n_id == g.n_id; }
    bool operator<(const Node& g) const { return n_id < g.n_id; }

	std::string label;
	long n_id, deg, intrs;
	int cis, trans, lvl;
	size_t hs;
	bool isRoot, isInter;
	double gExp;
};

/**
 * This class describes a Gene.
 */
class Gene : public Vertex {

public:

	/**
	 * Default constructor
	 */
	Gene() :
		Vertex(), start_p(0), stop_p(0), entrez_ID(-1), hs(0), intr(0), cis(-1), trans(-1) {

	}
	Gene(std::string sym, std::string chr, uint_64 st_, uint_64 sp_, size_t h, int lv, long ps) :
		symbol(sym), chrom(chr), start_p(st_), stop_p(sp_), entrez_ID(-1), hs(h), intr(0), cis(-1), trans(-1) {
		setLevel(lv);
		setId(ps);
	}
	Gene(const Gene& g) : Vertex(g),symbol(g.symbol), chrom(g.chrom),
			start_p(g.start_p), stop_p(g.stop_p), entrez_ID(g.entrez_ID), hs(g.hs), intr(g.intr),
			exp_t(g.exp_t),	gene_contacts(g.gene_contacts) { }
	Gene(Gene&& g): Vertex(std::move(g)), symbol(std::move(g.symbol)), chrom(std::move(g.chrom)),
			start_p(std::move(g.start_p)), stop_p(std::move(g.stop_p)),
			entrez_ID(std::move(g.entrez_ID)), hs(std::move(g.hs)), intr(g.intr), cis(g.cis), trans(g.trans),
			exp_t(std::move(g.exp_t)), gene_contacts(g.gene_contacts){ };

	Gene& operator=(const Gene& g) {
		if(this != &g) {
			Vertex::operator =(g);
			symbol = g.getSymbol();
			chrom = g.chrom;
			start_p = g.start_p;
			stop_p = g.stop_p;
			entrez_ID = g.entrez_ID;
			id = g.getPosition();
			level = g.getLevel();
			exp_t = g.exp_t;
			hs = g.hs;
			intr = g.intr;
			cis = g.cis;
			trans = g.trans;
			gene_contacts = g.gene_contacts;
		}
		return *this;
	}
	Gene& operator=(Gene&& g){
		if(this != &g){
			Vertex::operator =(std::move(g));
			symbol = std::move(g.getSymbol());
			chrom = std::move(g.chrom);
			start_p = std::move(g.start_p);
			stop_p = std::move(g.stop_p);
			entrez_ID = std::move(g.entrez_ID);
			id = std::move(g.getPosition());
			level = std::move(g.getLevel());
			exp_t = std::move(g.exp_t);
			hs = std::move(g.hs);
			intr = g.intr;
			cis = g.cis;
			trans = g.trans;
			gene_contacts = std::move(g.gene_contacts);
		}
		return *this;
	}


	~Gene() { }

	// getters
	inline long getEID()				const { return entrez_ID; }
	inline uint_64 getStart() 			const { return start_p; }
	inline uint_64 getStop() 			const { return stop_p; }
	inline std::string getChr() 		const { return chrom; }
	inline uint_64 getExtendedStart() 	const { return gene_contacts.start; }
	inline uint_64 getExtendedStop() 	const { return gene_contacts.stop; }
	inline std::string getSymbol() 		const { return symbol; }
	inline size_t getHS() 				const { return hs; }
	inline long getIntr() 				const { return intr; }
	inline long getNumCis() 			const { return cis; }
	inline long getNumTrans() 			const { return trans; }
	inline long getPosition()			const { return getId(); }
	inline double getGeneExp()			const { return exp_t.empty() ? 0.0 : exp_t[0].getLogFC(); }

	// setters
	inline void setChromosomeName(std::string chr) 	{ chrom = chr; }
	inline void setStartCoordinate(uint_64 _start) 	{ start_p = _start; }
	inline void setStopCoordinate(uint_64 _stop)   	{ stop_p  = _stop;  }
	inline void setEID(long id)						{ entrez_ID = id; }
	inline void setHS(size_t h)						{ hs = h; }
	inline void setIntr()							{ ++intr; }
	inline void setCis(int c)						{ cis = c; }
	inline void setTrans(int t)						{ trans = t; }
	inline void setSymbol(std::string _s)  			{ symbol = _s; }
	inline void setPosition(long p)					{ setId(p); }
	inline void setExtendedStart(uint_64 st)		{ gene_contacts.start = st; }
	inline void setExtendedStop(uint_64 sp)			{ gene_contacts.stop = sp; }
	inline void pushExpression(Expression e)		{ exp_t.push_back(e); }
	inline void pushBindings(BindingSite b)			{ bindings.emplace_back(b); }

	inline void sortExpressions() {
		if(!exp_t.empty())
			std::sort(exp_t.begin(), exp_t.end());
	}

	GeneT* toGeneT() {
		GeneT *gt = new GeneT(start_p, stop_p, entrez_ID, getLevel(), hs,
				gene_contacts.start, gene_contacts.stop, getPosition(), cis, trans);
		return gt;
	}

	bool operator==(const Gene& g) const { return entrez_ID == g.getEID(); }
	bool operator<(const Gene& g) const {
		if( (compareChar(chrom.c_str(), g.getChr().c_str()) < 0) ||
//				(hs == g.getHS() && entrez_ID < g.getEID()) ||
//				(hs == g.getHS() && entrez_ID == g.getEID() && start_p < g.getStart()) ||
//				(hs == g.getHS() && entrez_ID == g.getEID() && start_p == g.getStart() && stop_p < g.getStop())
				(hs == g.getHS() && start_p < g.getStart()) ||
				(hs == g.getHS() && start_p == g.getStart() && stop_p < g.getStop())
		)
			return true;
		else
			return false;
	}

	/*
	 * Friendly print a Gene entry
	 */
	friend std::ostream& operator<<(std::ostream &out, const Gene& gn) {
		// TODO: consider headers?
		out << " Chrom: "  << gn.chrom 	<< " | " <<
				//"Strand: " << gn.strand 	<< " | " <<
				"Start: "  << gn.start_p 	<< " | " <<
				"Stop: "   << gn.stop_p 	<< " | " <<
				"Symbol: " << gn.getSymbol()<< " | " <<
				"ID: "     << gn.entrez_ID 	<< " | ";

		return out;
	}

private:

	struct GeneFrags {
		uint_64 start, stop;
		GeneFrags() : start(0), stop(0) { }
		GeneFrags(uint_64 st, uint_64 sp) : start(st), stop(sp) { }
		~GeneFrags() {}
	};

	std::string symbol;	// hgnc symbol
	std::string chrom;
	//std::string strand;
	uint_64 start_p;
	uint_64 stop_p;
	long entrez_ID;
	size_t hs;		// hashed chromosome
	long intr;
	int cis, trans;

	std::vector<Expression> exp_t;  // at most 10?
	std::vector<BindingSite> bindings;
	GeneFrags gene_contacts;
};





#endif /* GENE_HPP_ */
