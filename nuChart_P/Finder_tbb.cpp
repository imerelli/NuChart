/*
 * Finder.cpp
 *
 *  Created on: Jun 20, 2014
 *      Author: fabio
 *
 */

#include "Finder.hpp"

#ifdef USE_OPENMP
#include <omp.h>
#endif

#ifdef USE_TBB
#include <tbb/tbb.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/task_scheduler_init.h>
#endif

#ifdef USE_NUMA
#include <numa.h>
#endif


// parse all needed file
void Finder::parseFiles(std::string genefile, std::string fragfile, std::string samfile) {
	std::cout << "\nPARSING INPUT FILES..." << std::endl;

#ifdef USE_NUMA
	if(numa_available() >= 0) {
		pnp = numa_alloc_interleaved(sizeof(SamData)+24); // std::vector uses 24 bytes for the vector object itself
		sam = new (pnp) std::vector<SamData>();
	} else {
		pnp = ::malloc(sizeof(SamData)+24);
		sam = new (pnp) std::vector<SamData>();
	}
#else
	pnp = ::malloc(sizeof(SamData)+24);
	sam = new (pnp) std::vector<SamData>();
#endif

	Parser ps(genefile, fragfile, samfile);
	ps.parseFiles(genes, frags, *sam); //, uniques, expressions, beds);

	fragmentsNumbers();
	genesNumbers();
	samNumbers();

#ifdef WEIGHIT
	ps.parseFeaturesFile(uniques);
	featuresNumbers();
#endif

	// bedNumbers();
}

void Finder::intersections(std::string infile) {
  Parser ps;
  ps.parseIntersecFile(locations, infile);
  locNumbers();
}

/* Build a contact map for the two chromosomes and store pointer in an unordered map
 * key is obtained pairing the two chromosomes' hashes with Cantor pairing function */
void Finder::buildContactMap(size_t chr1, size_t chr2, long bin) {
	int tot_size=0, r=0, c=0, count=0;
	size_t f_p=0, f_q=0, stop=0;
	size_t chr_list[2];
	uint_64 m_key=0;

	chr_list[0] = chr1;
	chr_list[1] = chr2;
	tot_size = feat_qt[chr1] + feat_qt[chr2];

	CMap<int> *cmap = new CMap<int>(tot_size, tot_size);

	for(unsigned i=0; i<2; ++i) {
		f_p = sam_pos[chr_list[i]];
		f_q = sam_qt[chr_list[i]];
		stop = f_p + f_q;
		while(f_p < stop) {
			if(	sam->at(f_p).getHS2() == chr_list[i] ) {
				c = ceil( sam->at(f_p).getStart1()/bin );
				r = ceil( sam->at(f_p).getStart2()/bin );

				cmap->incPoint(c,r);
				cmap->incPoint(r,c);

				++count;
			}
			++f_p;
		}
	}

	m_key = pairZ(chr1,chr2);
	shPtrCM sh_CM(cmap);
	cmaps[m_key] = sh_CM;
}


void Finder::buildCMaps() {
	if(uniques.empty()) {
		std::cerr << "ERROR - Features file is empty!" << std::endl;
		exit(EXIT_FAILURE);
	}

	size_t i=0, chr=0, sz_uniques;
	std::vector<size_t> chrs;
	chrs.reserve(24);
	cmaps.reserve(576);
	sz_uniques = uniques.size();

	while(i < sz_uniques) {
		chr = uniques[i]->hs;
		chrs.push_back( uniques[i]->hs );
		while( i < sz_uniques && uniques[i]->hs == chr ) ++i;
	}
	size_t sz_chrs = chrs.size();

	for(unsigned j=0; j<sz_chrs; ++j)
		for(unsigned k=0; k<sz_chrs; ++k)
			buildContactMap(chrs[j], chrs[k]);
}


void Finder::detectIntersections(std::string infile) {
	size_t l_p=0, l_q=0, stop=0;
	std::vector<Gene*> hotspots;
	hotspots.reserve(genes.size());

	intersections(infile);

	for(unsigned j=0; j<genes.size(); ++j) {
		l_p = loc_pos[genes[j].getHS()];
		l_q = loc_qt[genes[j].getHS()];
		stop = l_p + l_q;
		while (l_p < stop) {
			if( locations[l_p].pos >= genes[j].getStart() &&
					locations[l_p].pos <= genes[j].getStop() )
				genes[j].setIntr();
			++l_p;
		}
		if(genes[j].getIntr() > 0) hotspots.push_back(&genes[j]);
	}

	std::sort(hotspots.begin(), hotspots.end(), [](const Gene* g1, const Gene* g2) {
		return g1->getIntr() > g2->getIntr();
	});

	std::string::size_type found = infile.find("/mld");
	std::string mld = infile.substr(found+1, 5);

	// write output file
	std::stringstream ess;
	ess << mld << "_TopTenHotSpots.dat";
	std::ofstream e_file(ess.str().c_str(),
				std::ios_base::out | std::ios_base::trunc );
	e_file << "Number of Hotspots: " << hotspots.size() << std::endl;
	for(long j=0; j<10; ++j)
		e_file << hotspots[j]->getSymbol() << " := " << hotspots[j]->getIntr() << std::endl;
	e_file.close();

}

void Finder::findGeneExpression(Gene& g) {
	if(!expressions.empty()) {
		for(unsigned u=0; u<expressions.size(); ++u)
			if( expressions[u]->getHS() == DJBHash(g.getSymbol()) ||
					(expressions[u]->getID() == g.getEID()) )
				g.pushExpression(*expressions[u]);
		g.sortExpressions();
	} else {
		std::cerr << "ERROR - Expression file is empty." << std::endl;
		return; //exit(EXIT_FAILURE);
	}
}

void Finder::findBindings(Gene& g) {
	if(!beds.empty()) {
		size_t f_p=0, f_q=0, stop=0;

		f_p = bed_pos[g.getHS()];
		f_q = bed_qt[g.getHS()];
		stop = f_p + f_q;
		while(f_p < stop) {
			if(beds[f_p]->getStart() >= g.getExtendedStart() &&
					beds[f_p]->getStop() <= g.getExtendedStop() )
				g.pushBindings(*beds[f_p]);

			++f_p;
		}
	} else {
		std::cerr << "ERROR - Bed file is empty." << std::endl;
		return; //exit(EXIT_FAILURE);
	}
}

// get the score of an edge between two fragments using the contact map
void Finder::calculateEdgeScore(miniEdge& edge, /*std::vector<SamDataT*>& samt,*/
		double *gcc, double *map, uint_64 *len, float limit, long bin) {

	size_t f_p=0, f_q=0, stop=0, m_f=0, f=0;

	size_t chr_list[2];
	double *rhand[2];
	double *map_u;
	uint_64 m_key = pairZ(edge.HS1, edge.HS2);

	chr_list[0] = (edge.HS1);
	chr_list[1] = (edge.HS2);
	m_f = feat_qt[chr_list[0]] + feat_qt[chr_list[1]];

	CMap<int> *cmap = cmaps[m_key].get();

	for(unsigned i=0; i<2; ++i) {
		f_p = feat_pos[chr_list[i]];
		f_q = feat_qt[chr_list[i]];
		stop = f_p + f_q;
		for(; f_p < stop && f < m_f; ++f_p, ++f) {
			gcc[f] = uniques[f_p]->gcc;
			map[f] = uniques[f_p]->map;
			len[f] = uniques[f_p]->len;
		}
	}

	int *cmap_vec = cmap->upperTriangular();

	// get covariance matrices ------------------------------------
	CMap<double> len_m(m_f, m_f), map_m(m_f, m_f), gcc_m(m_f, m_f);
	for(unsigned i=0; i<m_f; ++i) {
		for(unsigned j=0; j<m_f; ++j) {
			len_m(i,j) = log(len[i] * len[j]);
			map_m(i,j) = log(map[i] * map[j]);
			gcc_m(i,j) = log(gcc[i] * gcc[j]);
		}
	}

	// mean and standard deviation ---------------
	double l_mean = len_m.meanVal();
	double g_mean = gcc_m.meanVal();
	double l_sd = len_m.standardDeviation(l_mean);
	double g_sd = gcc_m.standardDeviation(g_mean);

	// normalize len and gcc using mean and standard deviation
	for(unsigned i=0; i<m_f; ++i) {
		for(unsigned j=0; j<m_f; ++j) {
			len_m(i,j) = (len_m(i,j) - l_mean) / l_sd;
			gcc_m(i,j) = (gcc_m(i,j) - g_mean) / g_sd;
		}
	}

	// get upper triangular part of the matrix with no diagonal
	// so that vectors will be of the same length as the contact map (transformed into vector)

	map_u = map_m.upperTriangular();

	// build X matrix, the linear predictor (a linear combination of unknown parameters)
	rhand[0] = len_m.upperTriangular(); //rhand.push_back( len_u );
	rhand[1] = gcc_m.upperTriangular(); //rhand.push_back(gcc_u );

	size_t cm_size = cmap->upperTriangularSize();
	size_t tr_size = len_m.upperTriangularSize();

	assert(cm_size == tr_size);

	// Fit model with Iterated Weighted Least Square
	Fit fit;
	fit.set_data(cmap_vec, cm_size, rhand, 2, map_u); // fit.set_data(cmap_vec, rhand, map_u);
	fit.fit_model();

	std::vector<double> coef;
	fit.get_coef(coef);

#ifdef PRINTFIT
	std::cout << "\n --- Summary from Iterated Weighted Least Square with Poisson Regression --- \n"
			<< std::endl;
	std::cout << "dispersion = " << fit.get_dispersion() << std::endl;
	std::vector<double> serr;
	fit.get_stderr(serr);

	printf("%10s%12s%12s%15s\n", "", "Estimate", "Std.Error", "p-value");
	for(size_t i = 0; i < 3; ++i) {
		printf("X%-9zu%12.9f%12.8f", i, coef[i], serr[i]);
		if(!fit.checkQuasi())
			printf("%15.6e\n", 2 * gsl_cdf_gaussian_P(-fabs(coef[i]/serr[i]), 1.0));
		else
			printf("%15.6e\n", 2 * gsl_cdf_tdist_P(-fabs(coef[i]/serr[i]), cm_size-fit.get_rank_X()));
	}
#endif

	double score=0, l=0, g=0;
	int r = ceil( edge.st1/bin );
	int c = ceil( edge.st2/bin );

	l = coef[2]*len_m(r,c);
	g = coef[3]*gcc_m(r,c);
	double tmp_res = coef[1] + l + g + map_m(r,c);

	score = ( 1 / exp(tmp_res) );

	edge.score = score;

	//delete[] coef;
	delete[] map_u; delete[] cmap_vec;
	delete[] rhand[0]; delete[] rhand[1];
}


/**
 * Search for Genes according to their entrezID
 *
 * Returns a Gene with the requested entrezID
 */
Gene Finder::getSingleGeneByEntrezId(int id) {
	if(genes.empty()) {
		std::cerr << "ERROR: Genes vector is empty! " << std::endl;
		exit(EXIT_FAILURE);
	}

	Gene g;
	std::vector<Gene> found;
	std::vector<Gene>::iterator it;
	for(it=genes.begin(); it!=genes.end(); ++it){
		if( (*it).getEID() == id )
			found.push_back((*it));
	}

	if(found.size() == 0) {
		std::cerr << "Gene " << id << " not found. Aborting..." << std::endl;
		exit(EXIT_FAILURE);
	}

	found.erase( std::unique(found.begin(), found.end()), found.end() );
	g = found[0];

	return g;
}

/**
 * Search for Genes according to their entrezID
 *
 * Returns a pointer to the Gene with the requested entrezID
 */
Gene* Finder::getPtrSingleGeneByEntrezId(int id) {
	if(genes.empty()) {
		std::cerr << "ERROR: Genes vector is empty!" << std::endl;
		exit(EXIT_FAILURE);
	}

	std::vector<Gene> found;
	std::vector<Gene>::iterator it;
	for(it=genes.begin(); it!=genes.end(); ++it){
		if( (*it).getEID() == id )
			found.push_back((*it));
	}

	if(found.size() == 0) {
		std::cerr << "Gene " << id << " not found. Aborting..." << std::endl;
		exit(EXIT_FAILURE);
	}

	found.erase( std::unique(found.begin(), found.end()), found.end() );
	if(found.size() == 1)
		return &genes[found[0].getPosition()];
	else
		return NULL;

}

/**
 * Search for a Gene according to its hgnc symbol
 * Returns a Gene with the requested hgnc_symbol
 */
Gene Finder::getSingleGeneBySymbol(const std::string sym) {
	if(genes.empty()) {
		std::cerr << "ERROR: Genes vector is empty!" << std::endl;
		exit(EXIT_FAILURE);
	}

	std::vector<Gene> found;
	std::vector<Gene>::iterator it;
	Gene g;

	for(it=genes.begin(); it!=genes.end(); ++it) {
		if(sym.compare((*it).getSymbol()) == 0)
			found.push_back( (*it) );
	}

	if(found.size() == 0) {
		std::cerr << "Gene " << sym << " not found. Aborting..." << std::endl;
		exit(EXIT_FAILURE);
	}

	found.erase( std::unique(found.begin(), found.end()), found.end() );
	g = found[0];

	return g;
}

/**
 * Search for a Gene according to its hgnc symbol
 * Returns a pointer to the Gene with the requested hgnc_symbol
 */
Gene* Finder::getPtrSingleGeneBySymbol(const std::string sym) {
	if(genes.empty()) {
		std::cerr << "ERROR: Genes vector is empty!" << std::endl;
		exit(EXIT_FAILURE);
	}

	std::vector<Gene> found;
	std::vector<Gene>::iterator it;
	Gene g;

	for(it=genes.begin(); it!=genes.end(); ++it) {
		if(sym.compare((*it).getSymbol()) == 0)
			found.push_back( (*it) );
	}

	if(found.size() == 0) {
		std::cerr << "Gene " << sym << " not found. Aborting..." << std::endl;
		exit(EXIT_FAILURE);
	}

	found.erase( std::unique(found.begin(), found.end()), found.end() );
	if(found.size() == 1)
		return &genes[found[0].getPosition()];
	else
		return NULL;
}

/**
 * Search for Genes within a range of coordinates
 * Returns genes within the specified coordinates
 */
std::vector<Gene> Finder::getGenesByCoordinates(std::string chr_, uint_64 gst, uint_64 gend) {
	if(genes.empty()) {
		std::cerr << "ERROR: Genes vector is empty!" << std::endl;
		exit(EXIT_FAILURE);
	}
	size_t f_p=0, f_q=0, stop=0;
	std::vector<Gene> found;

	f_p = gene_pos[DJBHash(chr_)];
	f_q = gene_qt[DJBHash(chr_)];
	stop = f_p + f_q;
	while(f_p < stop) {
		if( genes[f_p].getStart() >= gst && genes[f_p].getStop() <= gend )
			found.push_back(genes[f_p]);
		++f_p;
	}

	if(found.size() == 0) {
		std::cerr << "Genes in gene cluster not found. Aborting..." << std::endl;
		exit(EXIT_FAILURE);
	}

	found.erase( std::unique(found.begin(), found.end()), found.end() );
	return found;
}

/**
 * Look for neighbour genes over paired-ends Hi-C reads
 *
 * From the given initial gene(s), it builds a graph representing the chromosome
 * conformation, showing their disposition in relation with other genes.
 * The graph construction is a BFS-like procedure, which iteratively builds the
 * graph by visiting all nodes (genes) at each graph level.
 *
 * This function uses the PARFOR data-parallel skeleton to iterate over Hi-C reads.
 *
 * Grpah edges contain information related to the chromosome location where two
 * genes have been found. A normalisation step is used to reduce experimental biases.
 */
void Finder::parFindConnections(std::vector<std::string> gvec, std::string dts,
		int sc_limit, int w, std::string clstr, int chunk, float e_lm, int plot,
		int colls, std::string locs, bool activeSched, int only) {
	size_t genes_size = genes.size();
	size_t frags_size = frags.size();
	size_t sam_size = sam->size();
	std::vector<nc_task_t*> tasklist;
	std::vector<bool> v_mask(genes_size, false);
	std::vector<GeneT*> gqueue;
	Graph* bigGraph = new Graph(genes_size);

	int level=0;
	size_t nwrks = w;
	size_t m_s;
	m_s = maxMapValue(sam_qt);   	// max number of reads - per chr

	// PREALLOCATIONS -------------------------------------------------
	std::vector< std::vector<size_t> > 		positions(nwrks);
	std::vector< std::vector<nc_task_t*> >	nutasks(nwrks);
	std::vector< std::vector<EdgeT*> >		localGraphs(nwrks);

	for(unsigned f=0; f<nwrks; ++f) {
		nutasks[f].reserve(m_s/nwrks);
		positions[f].reserve(genes_size/nwrks);
		localGraphs[f].reserve(frags_size/nwrks);
	}
	tasklist.reserve(sam_size);
	gqueue.reserve(genes_size);
	// ----------------------------------------------------------------

	ParallelFor ffpf(nwrks, true /*spin-wait*/, true /*spin-barrier*/);
	ffpf.disableScheduler(activeSched);

	// extend genes to fragments coordinates
	if(frags.empty()) {
		std::cerr << "ERROR: Fragments vector is empty!" << std::endl;
		exit(EXIT_FAILURE);
	}
	ffpf.parallel_for(0, genes_size, 1, chunk, [&] (const unsigned s) {
		uint_64 start=0, stop=0;
		findLimits(genes[s], frag_pos[genes[s].getHS()], frag_qt[genes[s].getHS()], &start, &stop);
		if (start != 0) {
			genes[s].setExtendedStart(start);
			genes[s].setExtendedStop(stop);
		}
	}, nwrks);

	ffpf.threadPause(); // put spinning threads to sleep

	// -------------- DUPLICATE----------------------------------------
	std::vector<SamDataT*> *samt;
	void *pt;

#ifdef USE_NUMA
	if(numa_available >= 0) {
		pt = numa_alloc_interleaved(sizeof(SamDataT*)+24);
		samt = new (pt) std::vector<SamDataT*>();
	} else {
		pt = ::malloc(sizeof(SamDataT*)+24);
		samt = new (pt) std::vector<SamDataT*>();
	}
#else
	pt = ::malloc(sizeof(SamDataT*)+24);
	samt = new (pt) std::vector<SamDataT*>();
#endif

	samt->reserve(sam_size);
	for(unsigned i=0; i<sam_size; ++i)
		samt->push_back(sam->at(i).toSamDataT());

	std::vector<FragmentT*> fragst;
	fragst.reserve(frags_size);
	for(unsigned i=0; i<frags_size; ++i)
		fragst.push_back(frags[i].toFragmentT());

	std::vector<GeneT*> genest;
	genest.reserve(genes_size);
	for(unsigned i=0; i<genes_size; ++i)
		genest.push_back(genes[i].toGeneT());
	// ----------------------------------------------------------------

#ifdef WEIGHIT
	std::cout << " [ContactMaps] Building Contact Maps for the given SAM file..." << std::endl;
	double cmt = tmn::elapsedTime(0);
	buildCMaps();
	cmt = tmn::elapsedTime(1);
	std::cout << " [ContactMaps] Done. " << cmaps.size() << " ContactMaps built (" << cmt << "ms) | ";
	long rsz = ( (48+sizeof(umCM))+(cmaps.size()*(sizeof(shPtrCM)+sizeof(uint_64)+32)) ) >> 10;
	std::cout << rsz << " KB" << std::endl;
#endif

	//std::cout << "\n*Current Mem. usage: " << (getCurrentRSS()/1024) << " KB" << std::endl;

	// prepare root node(s) in order to start building the graph ------
	std::vector<std::string> roots;
	roots.reserve(gvec.size());
	for(unsigned gv=0; gv<gvec.size(); ++gv) {
		Gene *root;
		if(isEntrezID(gvec[gv]))
			root = getPtrSingleGeneByEntrezId( atoi(gvec[gv].c_str()) );
		else
			root = getPtrSingleGeneBySymbol(gvec[gv]);

		if(root->getExtendedStart()) {
			v_mask[root->getPosition()] = true;
			genes[root->getPosition()].setLevel(0);
			genest[root->getPosition()]->setLevel(0);
			genes[root->getPosition()].setGraphRoot(&genes[root->getPosition()]);
			genest[root->getPosition()]->setGraphRoot(genest[root->getPosition()]);
			gqueue.push_back( genest[root->getPosition()] );
			roots.push_back(root->getSymbol());
		} else {
			std::cerr << "WARNING: fragments not found for Gene "
					<< gvec[gv] << ". Skipping gene..." << std::endl;
		}
	}
	// ----------------------------------------------------------------

	if(gqueue.empty()) {
		std::cerr << "ERROR: chromosome fragments not found for given genes. Aborting..." << std::endl;
		exit(EXIT_FAILURE);
	}

	double red=0, red2=0, epf=0, epf2=0, esc=0, esc2=0;

	std::cout << "\n-- START FINDING NEIGHBOURS..." << std::endl;
	while(!gqueue.empty() && level < sc_limit) {
		size_t gqsize = gqueue.size();
		// search for connections using fragments sequences from sam file -------------------------
		std::cout << " Level " << level << " with " << gqsize << " node(s)" << std::endl;
		esc = tmn::elapsedTime(0);

		//		tbb::task_scheduler_init init(nwrks);
		//		tbb::affinity_partitioner ap;
		//		tbb::parallel_for(tbb::blocked_range<long>(0, gqsize, chunk),
		//				[&] (const tbb::blocked_range<long>& r) {
		//					for (long s=r.begin();s!=r.end();++s) {
		//					// create thread-local storage

		ffpf.parallel_for_thid(0, gqsize, 1, chunk, [&] (const unsigned s, const int thid) {
			size_t f_p, f_q, stop;
			f_p = sam_pos[gqueue[s]->getHS()];
			f_q = sam_qt[gqueue[s]->getHS()];
			stop = f_p + f_q;
			while(f_p < stop ) {
				if( (samt->at(f_p)->getStart1() >= gqueue[s]->getExtendedStart()) &&
						(samt->at(f_p)->getStart1() <= gqueue[s]->getExtendedStop()) ) {
					nutasks[thid].push_back( new nc_task_t((gqueue[s]), samt->at(f_p)) );
				}
				++f_p;
			}
		}, nwrks);

		// });

		esc = tmn::elapsedTime(1);
		esc2 += esc;

		gqueue.clear();

		for(unsigned u=0; u<nwrks; ++u) {
			for(unsigned w=0; w<nutasks[u].size(); ++w)
				tasklist.push_back( nutasks[u][w] );
			nutasks[u].clear();
		}

		if(tasklist.empty()) {
			std::cerr << "ERROR: zero connections found to start the computation. Aborting..." << std::endl;
			exit(EXIT_FAILURE);
		}
		// ----------------------------------------------------------------------------------------

		std::cout << " Found " << tasklist.size() << " connections ("
				<< esc << "ms)"<< std::endl; //" in " << qq << " ms." << std::endl;

		size_t conn_size  = tasklist.size();

		epf = tmn::elapsedTime(0);

//		tbb::task_scheduler_init init(nwrks);
//		tbb::affinity_partitioner ap;
//		tbb::parallel_for(tbb::blocked_range<long>(0, conn_size, chunk),
//				[&] (const tbb::blocked_range<long>& r) {
//					for (long i=r.begin();i!=r.end();++i) {
//					// create thread-local storage

		// kernel core from here, for each connection search neighbors in parallel
		ffpf.parallel_for_thid( 0, conn_size, 1, chunk, [&, this] (const unsigned i, const int thid) {

			SamDataT *r1 = tasklist[i]->conn;
			size_t r1HS2 = r1->getHS2();
			uint_64 F1st=0, F1end=0, F2st=0, F2end=0;

			long f_pgene = gene_pos[r1HS2];
			long f_qgene = gene_qt[r1HS2];

#if defined(WEIGHIT) || defined(EXPINTER)
			size_t r1HS1 = r1->getHS1();
			long f_p1 = frag_pos[r1HS1];
			long f_q1 = frag_qt[r1HS1];
			long f_p2 = frag_pos[r1HS2];
			long f_q2 = frag_qt[r1HS2];

			bool found_F1 = false;
			bool found_F2 = false;

			// search F1 coordinates ----------------------------------
			for (; f_q1 > 0 && !found_F1; ++f_p1, --f_q1)
				if( (fragst[f_p1]->getStart() <= r1->getStart1()) &&
						(fragst[f_p1]->getStop() >= r1->getStart1()) ) {
					F1st = fragst[f_p1]->getStart();
					F1end = fragst[f_p1]->getStop();
					found_F1 = true;
				}
			// search F2 coordinates ----------------------------------
			for (; f_q2 > 0 && !found_F2; ++f_p2, --f_q2)
				if( (fragst[f_p2]->getStart() <= r1->getStart2()) &&
						(fragst[f_p2]->getStop() >= r1->getStart2()) ) {
					F2st = fragst[f_p2]->getStart();
					F2end = fragst[f_p2]->getStop();
					found_F2 = true;
				}
#endif

			// Search for genes possibly contained in connection fragments
			size_t intra=0;
			for (; f_qgene > 0; ++f_pgene, --f_qgene) {
				if (genest[f_pgene]->getExtendedStart() <= r1->getStart2()
						&& genest[f_pgene]->getExtendedStop() >= r1->getStart2()) {
					++intra;
					positions[thid].push_back(genest[f_pgene]->getPosition());
					localGraphs[thid].push_back( new EdgeT((tasklist[i]->w_gene), genest[f_pgene], tasklist[i]->conn,
							/*NULL,*/ F1st, F1end, F2st, F2end, false) );
				}
			} //end while f_qgene

#ifdef EXPINTER
			if(intra == 0) { // intergenic case

				// before genes -----------------------------------------------------------------------------------
				long f_pbg = gene_pos[r1HS2];
				long f_qbg = gene_qt[r1HS2];

				long stop = std::max(f_pbg + f_qbg - 1, (long)0);
				bool found_bg = false;
				for(; stop >= f_pbg && !found_bg; --stop) {
					if( (genest[stop]->getExtendedStart() > (F2st-1000000)) &&
							(genest[stop]->getExtendedStart() < F2st) ) {
						//if (/*!v_mask[genest[stop]->getPosition()] &&*/ addAB /*&& genest[stop]->getExtendedStart() > 0*/)
						positions[thid].push_back(genest[stop]->getPosition());

						localGraphs[thid].push_back( new EdgeT((tasklist[i]->w_gene), genest[stop],
								tasklist[i]->conn, /*NULL,*/ F1st, F1end, F2st, F2end, true) );
						found_bg = true; //break; //stop on nearest
					}
				}
				// ------------------------------------------------------------------------------------------------

				// after genes ------------------------------------------------------------------------------------
				long f_pag = gene_pos[r1HS2];
				long f_qag = gene_qt[r1HS2];
				bool found_ag = false;
				for(; f_qag > 0 && !found_ag; ++f_pag, --f_qag) {
					if( (genest[f_pag]->getExtendedStart() > F2end) &&
							(genest[f_pag]->getSExtendedtart() < (F2end+1000000)) ) {
						//if (/*!v_mask[genest[stop]->getPosition()] &&*/ addAB /*&& genest[stop]->getExtendedStart() > 0*/)
						positions[thid].push_back(genest[f_pag]->getPosition());

						localGraphs[thid].push_back( new EdgeT((tasklist[i]->w_gene), genest[f_pag],
								tasklist[i]->conn, /*NULL,*/ F1st, F1end, F2st, F2end, true) );
						found_ag = true; //break; //stop on nearest
					}
				}
				// ------------------------------------------------------------------------------------------------
			} // end-intergenic
#endif

		}, nwrks);

		// }); // now modify reduce

		epf = tmn::elapsedTime(1);
		epf2 += epf;

		ffpf.threadPause();

		red = tmn::elapsedTime(0);
		int t3=0;
		for(unsigned l=0; l<nwrks; ++l) {
			// handle normal vertices
			for(unsigned k=0; k<positions[l].size(); ++k) {
				if(!v_mask[positions[l][k]]) {
					v_mask[positions[l][k]] = true;
					genest[positions[l][k]]->setLevel(level+1);
					gqueue.push_back( genest[positions[l][k]] );
				}
			}

			for (unsigned z=0; z<localGraphs[l].size(); ++z)
				t3 += bigGraph->add_edge(localGraphs[l][z], true);

			positions[l].clear();
			localGraphs[l].clear();
			//edges.clear();
		}

		red = tmn::elapsedTime(1);
		red2 += red;
		std::cout << " Built partial Graph at level " << level << " (" << epf << " + "<< red << "ms)" << std::endl;

		tasklist.clear();

#ifdef TEST3
		if(t3<0) level = sc_limit;
#endif

		level++;
	} // -- while

	double st1 = epf2+red2+esc2;
	std::cout << " -- STOP FINDING NEIGHBOURS: " << st1 << "ms." << std::endl;
	std::cout << "\n-- TIMINGS:\n"
			<< " Connections: " << esc2 << "ms (tot)\n"
			<< " Final Graph in " << epf2 << "ms (+ " << red2 << "ms of reduce phase)\n" << std::endl;

	// PRINT GRAPH INFO --------------------------
	std::cout << " Nodes: " << bigGraph->getNumNodes();
	std::cout << "; Edges: " << bigGraph->edgesSize()
						<<  "\n MaxDegree: " << bigGraph->getMaxDegree()
						<< "; MinDegree: " << bigGraph->getMinDegree()
						<< "; AvgDegree: " << (double) bigGraph->getAverageDegree() << std::endl;
	// -------------------------------------------

	// FREE VECTORS ------------------------------
	for(unsigned l=0; l<nwrks; ++l) {
		freeVector(positions[l]);
		freeVector(localGraphs[l]);
	}
	freeVector(gqueue);
	freeVector(tasklist);
	freeVector(positions);
	freeVector(localGraphs);

	for(unsigned i=0; i<fragst.size(); ++i)
	  free(fragst[i]);
	
	freeVector(frags);
	malloc_trim(0);
	// -------------------------------------------

	//std::cout << "\n*Current Mem. usage: " << (getCurrentRSS()/1024) << " KB" << std::endl;

	std::vector<EdgeT*> t_edges(bigGraph->getEdges());
	std::vector<miniEdge> me;
	size_t sz_edges = t_edges.size();

	me.reserve(sz_edges);
	for(unsigned i=0; i<sz_edges; ++i)
		me.push_back( t_edges[i]->toMiniEdge() );

#ifdef WEIGHIT
	if(uniques.empty()) {
		std::cerr << "\nFeatures file not found! Can't compute edge's weight" << std::endl;
		exit(EXIT_FAILURE);
	} else {
		double st2=0;
		size_t m_f = maxMapValue(feat_qt);		// max qt of features - per chr

		double **gcc = new double*[nwrks];
		double **map = new double*[nwrks];
		uint_64 **len = new uint_64*[nwrks];

		for(unsigned f=0; f<nwrks; ++f) {
			gcc[f] = new double[2*m_f]();
			map[f] = new double[2*m_f]();
			len[f] = new uint_64[2*m_f]();
		}

		std::cout << "\n-- START SCORE COMPUTATION FOR " << sz_edges << " EDGES..." << std::endl;

#if defined(USE_OPENMP)

		st2 = tmn::elapsedTime(0);
		std::cout << " USE_OMP" << std::endl;
		if(chunk == 0) chunk = 1;

#pragma omp parallel for num_threads(nwrks) default(shared) schedule(dynamic,chunk)
		for(unsigned long j=0; j<me.size(); ++j) {
			long thid = omp_get_thread_num();

#ifdef EXPINTER
			if( !me[j].inter ) // do not calculate af/bf edges
#endif
				calculateEdgeScore(me[j], /*samt,*/, gcc[thid], map[thid], len[thid], e_lm);
		}

		st2 = tmn::elapsedTime(1);

#elif defined(USE_TBB)


		std::cout << " USE_TBB" << std::endl;
		tbb::task_scheduler_init init(nwrks);
		//tbb::affinity_partitioner ap;
		st2 = tmn::elapsedTime(0);
		tbb::parallel_for(tbb::blocked_range<long>(0, me.size(), chunk),
				[&] (const tbb::blocked_range<long>& r) {
			double *t_gcc = new double[2*m_f];
			double *t_map = new double[2*m_f];
			uint_64 *t_len = new uint_64[2*m_f];

			for (long j=r.begin();j!=r.end();++j) {

#ifdef EXPINTER
			if( !me[j].inter ) // do not calculate af/bf edges
#endif
				calculateEdgeScore(me[j], /*samt,*/, t_gcc, t_map, t_len, e_lm);
			}

			delete[] t_gcc;
			delete[] t_map;
			delete[] t_len;
		});

		st2 = tmn::elapsedTime(1);

#else

		std::cout << " USE_FF" << std::endl;
		st2 = tmn::elapsedTime(0);
		ffpf.parallel_for_thid(0, sz_edges, 1, chunk, [&] (const unsigned j, const int thid) {

#ifdef EXPINTER
			if( !me[j].inter ) // do not calculate af/bf edges
#endif
				calculateEdgeScore(me[j], /**samt,*/ gcc[thid], map[thid], len[thid], e_lm);
		}, nwrks);

		st2 = tmn::elapsedTime(1);
		ffpf.threadPause();

#endif // OMP-TBB-FF
		std::cout << "-- STOP SCORE CALCULATION: " << st2 << " ms." << std::endl;

		// --- scale score to range (0,1] -----------------------------------------------
		double max_sc=0, min_sc=0;
		for(unsigned i=0; i<sz_edges; ++i) {
			if(me[i].score>max_sc)
				max_sc = me[i].score;
			if(min_sc==0 || me[i].score<min_sc)
				min_sc = me[i].score;
		}
		std::cout << " \nMax Score: " << max_sc << "; Min Score: " << min_sc << std::endl;

		double ead=0.0, psum=0.0;
		for(unsigned i=0; i<sz_edges; ++i) {
			me[i].prob = (me[i].score - min_sc) / (max_sc - min_sc);

#ifdef TEST_EDGES
			me[i].prob = me[i].prob > e_lm ? me[i].prob : -0.01;
#endif

			psum += me[i].prob;
		}

		long probSum = std::round(psum);
		ead = ((double) 2 / (double) bigGraph->getNumNodes()) * psum; // expected average degree

		bigGraph->setEdgesProb(me);

		//for(unsigned i=0; i<5; ++i)
//		ADR(bigGraph, probSum);
		printf(" Cumulative Edges Probability: %ld;\n Expected Average Degree: %f\n", probSum, ead);

		kCores(bigGraph);
		//ketaCores(bigGraph);

		for(unsigned w=0; w<nwrks; ++w) {
			delete[] gcc[w];
			delete[] map[w];
			delete[] len[w];
		}
		delete[] gcc; delete[] map; delete[] len;
	}
		// ------------------------------------------------------------------------------
#endif // weighit


	//std::cout << "\n*Current Mem. Usage: " << (getCurrentRSS() / 1024) << " KB" << std::endl;
	malloc_trim(0);

	// FIX VERTICES AND EDGES --------------------
	//size_t tsz = bigGraph->edgesSize();
	size_t *ids = new size_t[2*(sz_edges+1)]();
	std::vector<Node> t_vertex;   	// list of nodes, with attributes for plotting

	long id1=-1, id2=-1;
	for (unsigned i = 0; i < sz_edges; ++i) {
		id1 = t_edges[i]->getVertex1()->getId();
		id2 = t_edges[i]->getVertex2()->getId();
		ids[2*i]   = id1;
		ids[2*i+1] = id2;

#ifdef TEST_EDGES
		if(t_edges[i]->getProb() > e_lm) {
#endif

			if(t_edges[i]->getHS1() == t_edges[i]->getHS2()) {
				genest[id1]->setCis(); genest[id2]->setCis();
			} else {
				genest[id1]->setTrans(); genest[id2]->setTrans();
			}

#ifdef TEST_EDGES
		}
#endif

	}

	Parser *p = new Parser();
	p->parseExpressionFile(expressions);
	delete p;

	t_vertex.reserve(bigGraph->getNumNodes());
	for(unsigned v=0; v<genes_size; ++v) {
		if(v_mask[v]) {
			//findBindings(genes[v]);
			findGeneExpression(genes[v]);
			genes[v].setCis(genest[v]->getNumCis());
			genes[v].setTrans(genest[v]->getNumTrans());
			if(genes[v].getLevel() == 0) {
				t_vertex.push_back( Node( genes[v].getSymbol(), genes[v].getPosition(),
							  bigGraph->getNodeDegree(genes[v].getPosition()), genes[v].getIntr(), genes[v].getHS(),
							  genes[v].getNumCis(), genes[v].getNumTrans(), 0, true, false, genes[v].getGeneExp()) );
			} else {
				t_vertex.push_back( Node(genes[v].getSymbol(), genes[v].getPosition(),
							 bigGraph->getNodeDegree(genes[v].getPosition()), genes[v].getIntr(), genes[v].getHS(),
							 genes[v].getNumCis(), genes[v].getNumTrans(), genest[v]->getLevel(), false, false, genes[v].getGeneExp()) );
				genes[v].setLevel( genest[v]->getLevel() );
			}
		}
	}
	for(unsigned i=0; i<genes_size; ++i)
		::operator delete(genest[i]);
	genest.clear();
	// -------------------------------------------

	// RESTORE FULL EDGES ---------------------------------------------
	std::cout << "\n-- RESTORING FULL EDGES... " << std::endl;
	//std::vector<Edge*> edgs(sz_edges);
	std::vector<Edge*> *edgs;
	void *ept;

#ifdef USE_NUMA
	if(numa_available >= 0) {
		ept = numa_alloc_interleaved(sizeof(EdgeT*)+24);
		edgs = new (ept) std::vector<Edge*>(sz_edges);
	} else {
		ept = ::malloc(sizeof(Edge*)+24);
		edgs = new (ept) std::vector<Edge*>(sz_edges);
	}
#else
	ept = ::malloc(sizeof(Edge*)+24);
	edgs = new (ept) std::vector<Edge*>(sz_edges);
#endif

	double rest = tmn::elapsedTime(0);
	ffpf.parallel_for(0, sz_edges, 1, chunk, [&] (const unsigned i) {
	    //for(unsigned i=0; i<tsz; ++i) {
		uint_64 st1 = t_edges[i]->conn1->getStart1();
		uint_64 st2 = t_edges[i]->conn1->getStart2();
		size_t r2HS2 = t_edges[i]->conn1->getHS2();
		bool found_r2=false;
		long s2id=-1;
		long f_p = sam_pos[r2HS2];
		long f_q = sam_qt[r2HS2];

		for(; f_q > 0 && !found_r2; ++f_p, --f_q) {
			if( (st1 == samt->at(f_p)->getStart2()) &&
					(st2 == samt->at(f_p)->getStart1()) ) {
				s2id = samt->at(f_p)->getId();
				found_r2=true;

				try {
					edgs->at(i) = new Edge( t_edges[i]->getVertex1(), t_edges[i]->getVertex2(),
							genes[ids[2*i]].getSymbol(), genes[ids[2*i+1]].getSymbol(),
							sam->at(t_edges[i]->conn1->getId()).getSeq(), sam->at(s2id).getSeq(),
							sam->at(t_edges[i]->conn1->getId()).getChr1(), sam->at(t_edges[i]->conn1->getId()).getChr2(),
							t_edges[i]->getHS1(), t_edges[i]->getHS2(), t_edges[i]->getStart1(),
							t_edges[i]->getStart2(), t_edges[i]->getEnd1(), t_edges[i]->getEnd2(),
							t_edges[i]->getHashedValue(), me[i].score, me[i].prob, t_edges[i]->hasIntergenic()
					);
				} catch(const std::out_of_range& e) {
					std::cout << "[ERROR] " << e.what() << std::endl;
				}
			}
		}
		//}
	}, nwrks);
	rest = tmn::elapsedTime(1);
	delete[] ids;
	std::cout << "-- DONE IN " << rest << "ms." << std::endl;
	// ----------------------------------------------------------------

	ffpf.threadPause();
	std::cout << "\n-- WRITING OUTPUTS..." << std::endl;

	// needed for plotting edges and vertices files -------------------
	std::string fst;
	std::stringstream strts;
	if(roots.size() == 1)
		fst = roots[0];
	else if(!clstr.empty()) {
			strts << "Cl_" << clstr;
			fst = strts.str();
		} else {
		strts << roots[0];
		if(roots.size() > 1)
			for(unsigned i=1; i<roots.size(); ++i)
				strts << "--" << roots[i];
		fst = strts.str();
	}
	// ----------------------------------------------------------------

	int d=-1;
#ifdef PRINTLOG
	// PRINT VERTICES AND EDGES TO LOG FILE
	d = mkdir("dats", 0777);
	std::stringstream sst;
	sst << "dats/par_NuC" << sc_limit << "_" << fst << "_"<< w << "_" << timeStamp() << ".dat";
	std::ofstream log_file( sst.str().c_str(),
			std::ios_base::out | std::ios_base::trunc );

	log_file << "Starting with " << gqueue.size() << " genes:\n";
	for(unsigned r=0; r<roots.size(); ++r)
		log_file << roots[r] << "; ";
	log_file << std::endl;

	log_file << "Used " << nwrks << " workers.\n";
	log_file << bigGraph->getNumNodes() << " neighbour genes found in " << st1 << " ms." << std::endl;
	log_file << "Score for " << edgs.size() << " edges calculated in " << st2 << " ms.\n" << std::endl;

	log_file << "## " << bigGraph->getNumNodes()<< " Vertices (Genes) --- " << std::endl;
	for(size_t i=0; i<v_mask.size(); ++i)
		if(v_mask[i])
			log_file << genes[i].getSymbol() << std::endl;

	log_file << "## " << edgs.size() << " Edges --- " << std::endl;
	for(size_t i=0; i<edgs.size(); ++i)
		log_file << *(edgs[i]) << std::endl;

	log_file.close();
#endif

	std::string f_out("csv"), tm(tmn::timeStamp());
	mkdir(f_out.c_str(), 0777);
	f_out.append("/").append(tm);
	d = mkdir(f_out.c_str(), 0777);
	f_out.append("/");

#ifdef WEIGHIT
	//printTimesForCharts(t_edges.size(), st2, fst, dts, sc_limit, e_lm, nwrks);
#endif

#if defined(TEST2)
	d = mkdir("test2", 0777);
	printNodesWithLevel(v_mask, genes, fst, dts, sc_limit, nwrks, e_lm, "test2/");
	printEdgeList(edgs, fst, dts, sc_limit, nwrks, "test2/");
#elif defined(TEST3)
	d = mkdir("test3", 0777);
	printCSV1(edgs, fst, dts, sc_limit, nwrks, e_lm, "test3/");
	printEdgeList(edgs, fst, dts, sc_limit, nwrks, "test3/");
#elif defined(TEST4)
	d = mkdir("test4", 0777);
	printEdgeList(edgs, fst, dts, sc_limit, nwrks, "test4/");
	printEdgeListWeighted(edgs, fst, dts, sc_limit, nwrks, "test4/");
#endif

	std::string plt;
	if(colls) {
		plt = plotGVizIntr(bigGraph, *edgs, t_vertex, fst, dts, locs, sc_limit, nwrks);
		printEdgesStat(*edgs, genes, fst, dts, locs, sc_limit, nwrks, plt);
	}
	if(plot) {
		plt = plotGViz(bigGraph, *edgs, t_vertex, fst, dts, sc_limit, nwrks);
		printEdgesStat(*edgs, genes, fst, dts, "", sc_limit, nwrks, plt);
		//plt = plotGraphML(bigGraph, *edgs, t_vertex, fst, dts, sc_limit, nwrks);
	}

	printDegreeDistribution(bigGraph, bigGraph->getNumNodes(), fst, dts, sc_limit, nwrks, f_out);
	printNodesDegree(bigGraph, t_vertex, bigGraph->getNumNodes(), fst, dts, sc_limit, nwrks, f_out);
	//printNodesWithLevel(v_mask, t_vertex, bigGraph->getNumNodes(), fst, dts, sc_limit, nwrks, f_out);
	printCSV1(*edgs, fst, dts, sc_limit, nwrks, e_lm, f_out);
	std::string csdir = genesCisTrans(bigGraph, t_vertex, fst, dts, sc_limit, nwrks, f_out);

	//std::cout << "\n*Peak Mem. Usage: " << (getPeakRSS()/1024) << " KB" << std::endl;

	//csStdDev(csdir);

	++d; // avoid that annoying warning!

	std::cout << "\n-- CLEANING UP..." << std::endl;

	// free containers --------
	// disambiguation for new/delete
	for(unsigned i=0; i<edgs->size(); ++i)
		::operator delete(edgs->at(i));
	edgs->~vector();

	for(unsigned i=0; i<samt->size(); ++i)
		::operator delete(samt->at(i));
	samt->~vector();
	samt->clear();

#ifdef USE_NUMA
	if(numa_available >= 0) {
		if(pt)
			numa_free(pt, (sizeof(SamDataT*)+24));
		if(ept)
		numa_free(ept, (sizeof(Edge*)+24));
	} else {
		if(pt)
			::free(pt);
		if(ept)
			::free(ept);
	}
#else
	if(pt)
		::free(pt);
	if(ept)
		::free(ept);
#endif

	bigGraph->~Graph();
	::operator delete(bigGraph);
	// -------------------------------------

 /*
  * //normal version ----------------> calls these two functions
  *	MyClass *pMemory = new MyClass();  void *pMemory = operator new(sizeof(MyClass));
  *	                                   MyClass *pMyClass = new (pMemory) MyClass();
  *
  *	//normal version ----------------> calls these two functions
  *	delete pMemory;                    pMyClass->~MyClass();
  *	                                   operator delete(pMemory);
  */
}

// search and store fragments that compose a gene
void Finder::findLimits	(Gene &git, uint_64 p, uint_64 q, uint_64 *start, uint_64 *stop) {
	uint_64 stp = p + q;
	while(p < stp) {
		if( ((frags[p].getStart() <= git.getStart()) &&
				(frags[p].getStop() >= git.getStart()))
		) {
			*start = frags[p].getStart();
			while(p < stp) {
				if( ( (frags[p].getStop() >= git.getStop()) &&
						 (frags[p].getStart() <= git.getStop()))
				) {
					*stop = frags[p].getStop();
					return;
				}
				++p;
			}
		}
		++p;
	}
}

void Finder::csStdDev(std::string dirPath) {
	Parser p;
	std::vector<double> cis_V, trans_V;
	std::vector<CisTrans> stdev_C, stdev_T;
	std::vector<std::string> files;

	getdirFiles(dirPath, files);
	size_t files_sz = files.size();
	if(files_sz < 2) {
		std::cerr << "ERROR: needed at least two files to compute "
				<< "standard deviation among cis and trans" << std::endl;
		return;
	}

	std::vector< std::vector<CisTrans> > csts(files_sz);
	for(unsigned i=0; i<files_sz; ++i)
		p.buildCisTransMatrix(csts[i], files[i]);

	size_t csts_sz = csts[0].size();
	cis_V.reserve(files_sz); trans_V.reserve(files_sz);
	stdev_C.reserve(csts_sz); stdev_T.reserve(csts_sz);

	std::cout << std::fixed;
	std::cout.precision(8);
	for(unsigned j=0; j<csts_sz; ++j) {
		for(unsigned i=0; i<files_sz; ++i) {
			cis_V[i] = (double) csts[i][j].cis;
			trans_V[i] = (double) csts[i][j].trans;
		}
		stdev_C[j].nid = stdev_T[j].nid = csts[0][j].nid;
		stdev_C[j].stdev = vecStdDev(cis_V);
		std::cout << "stDev_C_" << stdev_C[j].nid << " = " << stdev_C[j].stdev << "; ";
		stdev_T[j].stdev = vecStdDev(trans_V);
		std::cout << "stDev_T_" << stdev_T[j].nid << " = " << stdev_T[j].stdev << "\n";
	}

	std::sort(stdev_C.begin(), stdev_C.end(), [&] (const CisTrans& c1, const CisTrans& c2) {
		return c1.stdev < c2.stdev;
	});

	std::sort(stdev_T.begin(), stdev_T.end(), [&] (const CisTrans& c1, const CisTrans& c2) {
		return c1.stdev < c2.stdev;
	});

//	for(unsigned i=0; i<10; ++i) {
//		std::cout << stdev_C[i].nid << ": " << stdev_C[i].stdev << std::endl;
//		std::cout << stdev_T[i].nid << ": " << stdev_T[i].stdev << std::endl;
//	}

}

std::string Finder::genesCisTrans(Graph *g, std::vector<Node>& vertices, std::string start, std::string dts,
		int sc_limit, int wrks, std::string out_folder) {
	std::stringstream ess;
	std::string fold("_csv");
	std::ofstream e_file;
	size_t vx_size = vertices.size();

	start.append(fold);
	out_folder.append(start).append("/cisTrans");
	mkdir(out_folder.c_str(), 0777);

	ess << out_folder << "/nodesCisTrans_" << dts << "_L" << sc_limit << "_" << tmn::timeStamp() << ".csv";
	e_file.open(ess.str().c_str(),
			std::ios_base::out | std::ios_base::trunc );

	e_file << "#NodeId\tDegree\tNumCis\tNumTrans" << std::endl;
	for(unsigned i=0; i<vx_size; ++i)
		e_file << vertices[i].n_id << "\t" << vertices[i].deg << "\t"
		<< vertices[i].cis << "\t" << vertices[i].trans << "\n";
	e_file << std::endl;

	e_file.close();
	out_folder.append("/");
	return fold;
}

void Finder::ADR(Graph *g, int p) {
	Graph* meg = g;
	std::vector<EdgeT*> dgs(meg->getEdges());
	std::vector<long> nids;
	nids.reserve(p);
	size_t r2int=0;
	double d1=0.0, d2=0.0, disA=0, disB=0;

	// sort by probability in non-increasing order
	std::sort(dgs.begin(), dgs.end(), sortEdgeT);

	for(long i=0; i<p; ++i) {
		double r = randomR();
		if(r <= dgs[i]->getProb()) {
			nids.push_back(i);
			meg->add_edge(dgs[i], false, true);
		}
	}

	for(unsigned j=0; j<genes.size(); ++j) {
		if(meg->isNode(j)) {
			EdgeT *e1 = meg->getRandomEdge((long) j);
			if(!e1) continue;

			r2int = randomI(dgs.size()-1);
			while( std::binary_search(nids.begin(), nids.end(), r2int) )
				r2int = randomI(dgs.size()-1);
			EdgeT *e2 = dgs[r2int];
			//				if(!e2) continue;

			disA = ( meg->getNodeDegree(e1->id1) ) - ( g->getNodeUDegree(e1->id1) );
			disB = ( meg->getNodeDegree(e1->id2) ) - ( g->getNodeUDegree(e1->id2) );
			d1 = (std::abs(disA - 1)+std::abs(disB - 1)) - (std::abs(disA)+std::abs(disB));

			disA = ( meg->getNodeDegree(e2->id1) ) - ( g->getNodeUDegree(e2->id1) );
			disB = ( meg->getNodeDegree(e2->id2) ) - ( g->getNodeUDegree(e2->id2) );
			d2 = (std::abs(disA - 1)+std::abs(disB - 1)) - (std::abs(disA)+std::abs(disB));

			if(d1+d2 < 0) {
				meg->remove_edge(e1);
				meg->add_edge(e2, false, true);
			}
		}
	}

	meg->printGraph(std::cout);
	std::cout << "MaxDegree: " << meg->getMaxDegree()
						<< "; MinDegree: " << meg->getMinDegree()
						<< "; AvgDegree: " << (double) meg->getAverageDegree() << std::endl;
}

void Finder::kCores(Graph *meg) {
	Graph *g = meg;

	long *d = new long[g->nodesSpace()]();
	long *c = new long[g->nodesSpace()]();
	std::vector<long> uv;
	std::deque<long>::iterator ituv;
	std::vector< std::deque<long> > bigD(g->getMaxDegree()+1);

	for(long i=0; i<g->nodesSpace(); ++i) {
		d[i] = g->getNodeDegree(i);
		bigD[d[i]].push_back(i);
	}

	long v=0;
	double elp = tmn::elapsedTime(0);
	for(unsigned k=0; k<=g->getMaxDegree(); ++k) {
		while(!bigD[k].empty()) {
			v = bigD[k].back();
			bigD[k].pop_back();
			c[v] = k;

			g->getNodesList(v, uv);
			for(unsigned i=0; i<uv.size(); ++i) {
				if( d[uv[i]] > k && !bigD[d[uv[i]]].empty() ) {
					ituv = std::find(bigD[d[uv[i]]].begin(), bigD[d[uv[i]]].end(), uv[i]);
					if(ituv != bigD[d[uv[i]]].end()) {
						bigD[d[uv[i]]].erase(ituv);
						bigD[d[uv[i]]-1].push_front(uv[i]);
						--d[uv[i]];
					}
				}
			}
			uv.clear();
			//g->remove_node(v);
		}
	}
	elp = tmn::elapsedTime(1);

	for(long i=0; i<g->nodesSpace(); ++i) {
		g->setNodeKcore(i, c[i]);
		meg->setNodeKcore(i, c[i]);
	}

	std::cout << " K-cores computed in " << elp << "ms." << std::endl;

	// print kcores to file
	std::string kc_folder("kcores");
	mkdir(kc_folder.c_str(), 0777);
	std::stringstream ess;
	ess << kc_folder << "/graph_kCore_" << tmn::timeStamp() << ".dat";
	std::ofstream e_file(ess.str().c_str(),
			std::ios_base::out | std::ios_base::app);

	g->printNodesKcore(e_file);
	e_file << std::endl;
	g->printNodesInKcore(e_file);

	e_file.close();
	delete[] d; delete[] c;
}

void Finder::ketaCores(Graph *meg, double eta) {
	Graph *g = meg;

	// eta_d := max k € [0..d_v] | Pr[deg(v) >= k] >= eta
	long *eta_d = new long[g->nodesSpace()]();
	long *c = new long[g->nodesSpace()]();
	std::vector<long> uv;
	std::deque<long>::iterator ituv;
	std::vector< std::deque<long> > bigD(g->getMaxDegree()+1);

	for(long i=0; i<g->nodesSpace(); ++i) {
		eta_d[i] = g->computeEtaDegree(i, eta);
		bigD[eta_d[i]].push_back(i);
	}

	long v=0;
	double elp = tmn::elapsedTime(0);
	for(unsigned k=0; k<=g->getMaxDegree(); ++k) {
		while(!bigD[k].empty()) {
			v = bigD[k].back();
			bigD[k].pop_back();
			c[v] = k;

			g->getNodesList(v, uv);
			for(unsigned i=0; i<uv.size(); ++i) {
				if( eta_d[uv[i]] > k && !bigD[eta_d[uv[i]]].empty() ) {
					ituv = std::find(bigD[eta_d[uv[i]]].begin(), bigD[eta_d[uv[i]]].end(), uv[i]);
					if(ituv != bigD[eta_d[uv[i]]].end()) {
						long etadegU = g->updateEtaDegree(uv[i], v, eta);
						bigD[eta_d[uv[i]]].erase(ituv);
						bigD[etadegU].push_front(uv[i]);
						eta_d[uv[i]] = etadegU;
					}
				}
			}
			uv.clear();
			//g->remove_node(v);
		}
	}
	elp = tmn::elapsedTime(1);

	for(long i=0; i<g->nodesSpace(); ++i) {
		g->setNodeKetacore(i, c[i]);
		meg->setNodeKetacore(i, c[i]);
	}

	std::cout << " Keta-cores computed in " << elp << "ms." << std::endl;

	//g->printGraph(std::cout);

	// print kcores to file
	std::string kc_folder("ketacores");
	mkdir(kc_folder.c_str(), 0777);
	std::stringstream ess;
	ess << kc_folder << "/graph_ketaCore_eta_" << eta << "_" << tmn::timeStamp() << ".dat";
	std::ofstream e_file(ess.str().c_str(),
			std::ios_base::out | std::ios_base::trunc );

	g->printNodesKetacore(e_file);
	e_file << std::endl;
	g->printNodesInKetacore(e_file);

	e_file.close();
	delete[] eta_d; delete[] c;
}



