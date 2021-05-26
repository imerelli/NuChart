/*
 * Edge.hpp
 *
 *  Created on: Jul 18, 2014
 *      Author: fabio
 */

#ifndef EDGE_HPP_
#define EDGE_HPP_

#include "common.hpp"
#include "Vertex.hpp"
#include "Edg.hpp"
#include "EdgeT.hpp"

/*
 * Class representing an Edge of the graph
 */
class Edge : public Edg {

public:
	Edge() : Edg(), hashed(0)/*, mirror(0)*/ { }
	Edge(Vertex *v1, Vertex *v2, std::string sym1, std::string sym2, std::string s1, std::string s2, std::string c1, std::string c2,
			size_t h1, size_t h2, uint_64 t1, uint_64 t2, uint_64 e1, uint_64 e2, uint_64 hsh, /*uint_64 mr,*/ double w, double p, bool ntr) {
		symbol1 = sym1;
		symbol2 = sym2;
		seq1 = s1;
		seq2 = s2;
		chr1 = c1;
		chr2 = c2;
		hs1 = h1;
		hs2 = h2;
		st1 = t1;
		st2 = t2;
		end1 = e1;
		end2 = e2;
		hashed = hsh;
		//mirror = mr;
		inter = ntr;
		Edg::prob = p;
		Edg::weight = w;
		Edg::v1 = v1;
		Edg::v2 = v2;

//		if(ntr) {
//			std::stringstream igstring;
//			igstring << chr1 << "_" << st2
//					<< "_" << end2;
//			inter_symb = igstring.str();
//		} else inter_symb = "NONE";
	}

	Edge(const Edge& e) : Edg(e) { *this = e; }
	Edge(Edge&& e) : Edg(std::move(e)), inter(e.hasIntergenic()), symbol1(std::move(e.getGeneSymbol1())), symbol2(std::move(e.getGeneSymbol2())),
			seq1(std::move(e.getSeq1())), seq2(std::move(e.getSeq2())), chr1(std::move(e.getChromosome1())),
			hs1(std::move(e.getHS1())), st1(std::move(e.getStart1())), end1(std::move(e.getEnd1())), chr2(std::move(e.getChromosome2())),
			hs2(std::move(e.getHS2())), st2(std::move(e.getStart2())), end2(std::move(e.getEnd2())),
			hashed(e.getHashedValue())/*, mirror(e.getHashedValue2())*/	{ }

	Edge& operator=(const Edge& e) { // shallow-copy is fine here!
		if(this != &e) {
			Edg::operator =(e);
			inter = e.hasIntergenic();
			symbol1 = e.getGeneSymbol1();
			symbol2 = e.getGeneSymbol2();
			seq1 = e.getSeq1();
			seq2 = e.getSeq2();
			chr1 = e.getChromosome1();
			hs1 = e.getHS1();
			st1 = e.getStart1();
			end1 = e.getEnd1();
			chr2 = e.getChromosome2();
			hs2 = e.getHS2();
			st2 = e.getStart2();
			end2 = e.getEnd2();
			hashed = e.getHashedValue();
			//mirror = e.getHashedValue2();
		}
		return *this;
	}

	Edge& operator=(Edge&& e) {
		if(this != &e) {
			Edg::operator=(std::move(e));
			inter = e.hasIntergenic();
			symbol1 = std::move(e.getGeneSymbol1());
			symbol2 = std::move(e.getGeneSymbol2());
			seq1 = std::move(e.getSeq1());
			seq2 = std::move(e.getSeq2());
			chr1 = std::move(e.getChromosome1());
			hs1 = e.getHS1();
			st1 = e.getStart1();
			end1 = e.getEnd1();
			chr2 = std::move(e.getChromosome2());
			hs2 = e.getHS2();
			st2 = e.getStart2();
			end2 = e.getEnd2();
			hashed = e.getHashedValue();
			//mirror = e.getHashedValue2();
		}
		return *this;
	}

	~Edge() {}

	// setters
	inline void setGeneSymbol1(std::string s) 	{ symbol1 = s; }
	inline void setGeneSymbol2(std::string s) 	{ symbol2 = s; }
	inline void setSeq1(std::string seq) 		{ seq1 = seq; }
	inline void setSeq2(std::string seq) 		{ seq2 = seq; }
	inline void setChromosome1(std::string f) 	{ chr1 = f; }
	inline void setHS1(size_t h) 				{ hs1 = h; }
	inline void setStart1(uint_64 v)			{ st1 = v; }
	inline void setEnd1(uint_64 v) 				{ end1 = v; }
	inline void setChromosome2(std::string f) 	{ chr2 = f; }
	inline void setHS2(size_t h)	 			{ hs2 = h; }
	inline void setStart2(uint_64 v) 			{ st2 = v; }
	inline void setEnd2(uint_64 v) 				{ end2 = v; }
	inline void setHashValue(uint_64 v) 		{ hashed = v; }
	inline void setIntergenic()					{ inter = true; }
	//inline void setHashValue2(uint_64 v) 		{ mirror = v; }

	//getters
	inline std::string getGeneSymbol1() const { return symbol1; }
	inline std::string getGeneSymbol2() const { return symbol2; }
	inline std::string getSeq1() 		const { return seq1; }
	inline std::string getSeq2() 		const { return seq2; }
	inline std::string getChromosome1() const { return chr1; }
	inline uint_64 getStart1() 			const { return st1; }
	inline uint_64 getEnd1() 			const { return end1; }
	inline std::string getChromosome2() const { return chr2; }
	inline uint_64 getStart2() 			const { return st2; }
	inline uint_64 getEnd2() 			const { return end2; }
	inline uint_64 getHashedValue() 	const { return hashed; }
	inline size_t getHS1() 				const { return hs1; }
	inline size_t getHS2() 				const { return hs2; }
	inline bool hasIntergenic()			const { return inter; }
	//inline uint_64 getHashedValue2() 	const { return mirror; }


//	bool operator==(const Edge& e) const { return hashed == e.getHashedValue(); }
//	bool operator<(const Edge& e) const { return hashed < e.getHashedValue(); }

	/*
	 * Friendly print an Edge.
	 */
	friend std::ostream& operator<<(std::ostream &out, const Edge &ed) {
		// TODO: consider headers?
		out << " Gene1: "   << ed.getGeneSymbol1() 	<< " | " <<
				"Gene2: " 	<< ed.getGeneSymbol2() 	<< " | " <<
				"Chr1: "  	<< ed.getChromosome1() 	<< " | " <<
				"Seq1: "   	<< ed.getSeq1() 	   	<< " | " <<
				"Start1: " 	<< ed.getStart1() 	   	<< " | " <<
				"End1: "    << ed.getEnd1()	  	   	<< " | " <<
				"Chr2: "  	<< ed.getChromosome2() 	<< " | " <<
				"Seq2: "   	<< ed.getSeq2() 	   	<< " | " <<
				"Start2: " 	<< ed.getStart2() 	   	<< " | " <<
				"End2: "    << ed.getEnd2() 	   	<< " | " <<
				//"InterGenic: " << ed.inter_symb		<< " | " <<
				"Score: " 	<< ed.getWeight()		<< " | " <<
				"Prob: "	<< ed.getProb()			<< " |";

		return out;
	}

private:

	bool inter;
	std::string symbol1;
	std::string symbol2;
	std::string seq1;
	std::string seq2;
	std::string chr1;
	size_t hs1;			// hashed SAM chromosome1
	uint_64 st1;
	uint_64 end1;
	std::string chr2;
	size_t hs2;			// hashed SAM chromosome2
	uint_64 st2;
	uint_64 end2;
	uint_64 hashed;		// hashed value for identifying the edge
	//std::string inter_symb;

	//uint_64 mirror;
};


#endif /* EDGE_HPP_ */
