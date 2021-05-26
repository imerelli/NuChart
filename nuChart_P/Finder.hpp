/*
 * Finder.hpp
 *
 *  Created on: Jun 20, 2014
 *      Author: fabio
 */

#ifndef FINDER_HPP_
#define FINDER_HPP_

#include <ff/parallel_for.hpp>

#include <unordered_map>

#include "Parsers.hpp"
#include "PrintOut.hpp"
#include "Graph.hpp"
#include "Edge.hpp"
#include "EdgeT.hpp"
#include "CMap.hpp"
#include "glm/Fit.hpp"

#include "Timings.hpp"
#include "FileIO.hpp"
#include "Profile.h"

using namespace ff;

typedef std::shared_ptr<CMap<int>> shPtrCM;
typedef std::unordered_map<uint_64, shPtrCM> umCM;

struct nc_task_t {

	nc_task_t() { }
    nc_task_t(GeneT *g, SamDataT *c):
        w_gene(g), conn(c) {}
    nc_task_t(const nc_task_t& t) = default;
    nc_task_t(nc_task_t&& t) = default;

    nc_task_t& operator=(const nc_task_t& other) = default;
    nc_task_t& operator=(nc_task_t&& other) = default;

	GeneT *w_gene;
	SamDataT *conn;
};


class Finder {
private:

	// search and store fragments that compose a gene
	void findLimits	(Gene &git, uint_64 p, uint_64 q, uint_64 *start, uint_64 *stop);

	void genesNumbers() {
		size_t qt=0, i=0, chr;

		gene_qt.reserve(24);  // chr 1 - 22 + x, y
		gene_pos.reserve(24); // chr 1 - 22 + x, y

		while(i < genes.size()) {
			qt = 0;
			chr = genes[i].getHS();
			gene_pos.emplace( chr, i );
			while( i < genes.size() && genes[i].getHS() == chr ) {
				++i;
				++qt;
			}
			gene_qt.emplace( chr, qt );
		}
	}

	void fragmentsNumbers() {
		size_t qt=0, i=0, chr=0;

		frag_qt.reserve(24);  // chr 1 - 22 + x, y
		frag_pos.reserve(24); // chr 1 - 22 + x, y

		while(i < frags.size()) {
			qt = 0;
			chr = frags[i].getHS();
			frag_pos.emplace( chr, i );
			while( i < frags.size() && frags[i].getHS() == chr ) {
				++i;
				++qt;
			}
			frag_qt.emplace( chr, qt );
		}
	}

	void featuresNumbers() {
		size_t qt=0, i=0, chr;

		feat_pos.reserve(24); // chr 1 - 22 + x, y
		feat_qt.reserve(24);  // chr 1 - 22 + x, y

		while(i < uniques.size()) {
			qt = 0;
			chr = uniques[i]->hs;
			feat_pos.emplace( chr, i );
			while( i < uniques.size() && uniques[i]->hs == chr ) {
				++i;
				++qt;
			}
			feat_qt.emplace( chr, qt );
		}
	}

	void samNumbers() {
		size_t qt=0, i=0, chr;

		sam_qt.reserve(24);  // chr 1 - 22 + x, y
		sam_pos.reserve(24); // chr 1 - 22 + x, y

		while(i < sam->size()) {
			qt = 0;
			chr = sam->at(i).getHS1();
			sam_pos.emplace( chr, i );
			while( i < sam->size() && sam->at(i).getHS1() == chr ) {
				++i;
				++qt;
			}
			sam_qt.emplace( chr, qt );
		}
	}

	void bedNumbers() {
		size_t qt=0, i=0, chr;

		bed_qt.reserve(24);  // chr 1 - 22 + x, y
		bed_pos.reserve(24); // chr 1 - 22 + x, y

		while(i < beds.size()) {
			qt = 0;
			chr = beds[i]->getHS();
			bed_pos.emplace( chr, i );
			while( i < beds.size() && beds[i]->getHS() == chr ) {
				++i;
				++qt;
			}
			bed_qt.emplace( chr, qt );
		}
	}

	void locNumbers() {
		size_t qt=0, i=0, chr;

		loc_qt.reserve(24);  // chr 1 - 22 + x, y
		loc_pos.reserve(24); // chr 1 - 22 + x, y

		while(i < locations.size()) {
			qt = 0;
			chr = locations[i].hs;
			loc_pos.emplace( chr, i );
			while( i < locations.size() && locations[i].hs == chr ) {
				++i;
				++qt;
			}
			loc_qt.emplace( chr, qt );
		}
	}

	// perform only intra-chromosome search (chr1 == chr2) - needed for testing
	void useOnlyChrom(std::string chr) {
		size_t f_p, f_q, stop, i, hs;
		std::vector<SamData> ns;

		std::cout << "Using only fragments on Chromosome " << chr << std::endl;

		hs = DJBHash(chr);
		i=0;
		f_p = sam_pos[hs];
		f_q = sam_qt[hs];
		stop = f_p + f_q;
		ns.reserve(stop);
		while (f_p < stop /*|| i < sam.size()*/) {
			if( sam->at(f_p).getHS2() == hs) {
				ns.push_back( static_cast<SamData&&>( sam->at(f_p) ));
			}
			++f_p;
			++i;
		}

		sam->swap(ns);
		sam->shrink_to_fit();

		sam_qt.clear(); sam_pos.clear();
		samNumbers();
	}

	void buildContactMap(size_t chr1, size_t chr2, long bin=1000000);
	void findGeneExpression(Gene& g);
	void findBindings(Gene& g);

	// Average Degree Rewiring - Extract a representative instance of the graph
	void ADR(Graph *g, int p);
	void kCores(Graph *meg);
	void ketaCores(Graph *meg, double eta=0.4);

    void intersections(std::string inf);
    std::string genesCisTrans(Graph *g, std::vector<Node>& vertices, std::string start, std::string dts,
    		int sc_limit, int wrks, std::string out_fold="csv/");

public:
	Finder() {}

	// copy constructor
	Finder(const Finder& f) {
		*this = f;
	}

	// destructor
	~Finder() {
		if(!uniques.empty())
			free_pointed_to(uniques);

		sam->clear();
		if(pnp) {
#ifdef USE_NUMA
			if(numa_available >= 0)
				numa_free(pnp, (sizeof(SamData)+24));
			else
				::free(pnp);
#else
			::free(pnp);
#endif
		}
		//malloc_trim(0);
	}

	// parsers
	void parseFiles(std::string genefile, std::string fragfile, std::string samfile, std::string exprfile);
	void detectIntersections(std::string inf);

	// Create a contact map between two chromosomes.
	void buildCMaps();
	void calculateEdgeScore(miniEdge& e, /*std::vector<SamDataT*>& smt,*/
			double *gcc, double *map, uint_64 *len, float limit=0, long bin=1000000);

	// Get a Gene entry
	inline Gene getSingleGeneByEntrezId(int id);
	inline Gene getSingleGeneBySymbol(const std::string sym);
	inline Gene* getPtrSingleGeneByEntrezId(int id);
	inline Gene* getPtrSingleGeneBySymbol(const std::string sym);
	std::vector<Gene> getGenesByCoordinates(std::string chr_, uint_64 gst, uint_64 gend);

	// routine for discovering neighbor connected genes
	void parFindConnections(std::vector<std::string> gvec, std::string dts,
			int sc_limit, int w, std::string clstr, int chunk, float e_lm, int plot,
			int colls, std::string locs="",	bool active=false, int only_chr=0);
	void csStdDev(std::string dirPath);

private:
	std::vector<Gene>      		genes;
	std::vector<Fragment>  		frags;
	std::vector<SamData>   		*sam;
	std::vector<Features*> 		uniques;
	std::vector<Expression*> 	expressions;
	std::vector<BindingSite*> 	beds;
	std::vector<Intr> 			locations;

	void *pnp;	// pointer for SamData vector's placement new

	// these are needed to speed up searches over parsed data
	std::unordered_map<size_t, size_t> frag_qt;
	std::unordered_map<size_t, size_t> frag_pos;
	std::unordered_map<size_t, size_t> gene_qt;
	std::unordered_map<size_t, size_t> gene_pos;
	std::unordered_map<size_t, size_t> sam_pos;
	std::unordered_map<size_t, size_t> sam_qt;
	std::unordered_map<size_t, size_t> feat_pos;
	std::unordered_map<size_t, size_t> feat_qt;
	std::unordered_map<size_t, size_t> bed_pos;
	std::unordered_map<size_t, size_t> bed_qt;
	std::unordered_map<size_t, size_t> loc_pos;
	std::unordered_map<size_t, size_t> loc_qt;
	umCM cmaps; //std::unordered_map<uint_64, CMap<int> > cmaps;
};

// anonymous namespace for static compare functions
namespace {

bool maxFunc(Features *a, Features *b) { return a->stop < b->stop; }
bool sortEdges(Edge e1, Edge e2 ) { return e1.getHashedValue() < e2.getHashedValue(); }

}

#endif /* FINDER_HPP_ */
