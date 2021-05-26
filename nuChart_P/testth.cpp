
#include <getopt.h>

#include "common.hpp"
#include "FileIO.hpp"
#include "Finder.hpp"
#include "Fragment.hpp"
#include "SamData.hpp"

void usage() {
	std::cerr << "\n -- USAGE: ./nuchart {-G [gene(s)] | -C [coords]} -S[sam-file] -N[genes-file]\n -F[frag-file] -T[interactions-file] -X[expression-file] -W[num-workers] \n-L[search-level] -g[grain] "
			<< "-E[edge-prob] -P[plot] -H[help]\n\n -- DESCRIPTION:\n"
			<< "\t[help]:\t\tprints this help message and exits.\n"
			<< "\t[genes]:\tentrezID or symbol of the starting gene. There can be more than one gene. (REQUIRED)\n"
		  << "\t[coords]:\tcluster of Genes. use all the genes found within given coordinates. Chromosome must be specified\n"
		  << "\t\t\t(es: -C chr7:27126219,27250333. Some cluster name are available, e.g., -C KRAB)\n"
			<< "\t[sam-file]:\tpath to the file of nucleotide sequence alignments (DEFAULT is \'~/datasets/LiebermanAiden.sam\')\n"
		  << "\t[genes-file]:\tpath to the species gene (DEFAULT is \'./extdata/human/human_genes.txt\')\n"
		  << "\t[frag-file]:\tpath to the chromosme fragments digeste by an enzime (DEFAULT is \'./extdata/human/DIG/HindIII.txt\')\n"
			<< "\t[interactions-file]: file containing virus interactions\n"
			<< "\t[expression-file]: path to the expression file (DEFAULT is NULL) \n"
			<< "\t[num-workers]:\tnumber of worker threads to be used (DEFAULT is 1, i.e. sequential)\n"
			<< "\t[search-level]:\tdeepness of the search - e.g.: 1 = direct neighbours; 2 = neighbours of neighbours, etc. (DEFAULT is 1)\n"
			<< "\t[grain]:\tdecide the grain for the parallel execution (DEFAULT IS 0, which means \'num_tasks/num_workers\')\n."
			<< "\t[edge-prob]:\tedges with probability below this threshold will not be considered (DEFAULT is 0, which means\n"
			<< "\t\t\tcompute score for ALL edges. Allowed probability values range from 0.2 to 1)\n"
			<< "\t[plot]:\t\tplot graphs using GraphViz\n\n"
			<< " -- EXAMPLES:\n"
			<< "\t./nuchart -S ~/datasets/SRR65470.sam -W 16 -L 3 -G 7157 2035 4622\n"
			<< "\tThis will run NuChart over the given SAM file, using 16 worker and iterating until\n"
			<< "\tneighbours of neighbours of neighbours of the starting genes are found.\n"
			<< "\tGraph construction begins from the 3 specified gene IDs.\n\n"
			<< "\t./nuchart -L 2 -W 8 -C chr7:27126219,27250333\n"
			<< "\tThis will run NuChart over the default SAM file, using 8 workers and iterating until\n"
			<< "\tneighbours of neighbours are found. Graph construction begins with the genes found within the given coordinates.\n\n"
			<< "\t./nuchart -G TP53 -L 2 -E 0.6\n"
			<< "\tThis will run NuChart over the default SAM file in sequential mode, iterating until\n"
			<< "\tneighbours of neighbours are found. Graph construction begins with the gene of given symbol.\n"
			<< "\tEdges with probability lower than 0.6 will not be plotted.\n" << std::endl;
	exit(EXIT_FAILURE);
}


int main(int argc, char **argv) {

	if(argc < 3) {
		usage();
	}

	if(compareChar(argv[1], "cistrans") == 0 ) {
		Finder ff;
		ff.csStdDev(argv[2]);
		return 0;
	}

	int c=0, nw=1, level=1, index=0, chunk=0, plot=0, intrs=0;
	float e_lm=0;
	uint_64 st=0, sp=0;
	std::string coords, chr, crd, clstr, mld="";
	std::string::size_type found;
	opterr = 0;

	std::string genes = "./extdata/human/human_genes.txt";
	std::string frags = "./extdata/human/DIG/HindIII.txt";
	std::string sam = "";
	std::string locs = "";
	std::string exprfile = ""; //DAgo

	std::vector<std::string> gvec;
	std::unordered_map<std::string, std::string> clusters;
	clustersMap(clusters);

	while ((c = getopt (argc, argv, "S:N:F:W:L:T:g:G:E:C:X:HP")) != -1) {
		switch(c) {
		case 'S':
			sam = optarg;
			break;

		case 'N':
			genes = optarg;
			break;

		case 'F':
			frags = optarg;
			break;

		case 'X':
			exprfile = optarg;
			break;

		case 'T':
			locs = optarg;
			intrs = 1;
			break;

		case 'L':
			level = atoi(optarg);
			break;

		case 'W':
			nw = atoi(optarg);
			break;

		case 'G':
			gvec.push_back(optarg);
			for (index = optind; index < argc && *argv[index] != '-'; index++)
				gvec.push_back( argv[index] );
			break;

		case 'g':
			chunk = atoi(optarg);
			break;

		case 'E':
			e_lm = atof(optarg);
			break;

		case 'P':
			plot=1;
			break;

		case 'C':
			coords = optarg;
			clstr = optarg;
			if(coords.compare(0, 3, "chr") != 0)
				coords = clusters[optarg];
			found = coords.rfind(":");
			chr = coords.substr(0, found);
			crd = coords.substr(found+1, coords.length());
			eraseTrailingWhiteSpaces(crd);
			found = crd.rfind(",");
			st = atol( crd.substr(0, found).c_str() );
			sp = atol( crd.substr(found+1, crd.length()).c_str() );
			break;

		case 'H':
			usage();
			break;

		case '?':
			std::cerr << "Unknown option!\n" << std::endl;
			usage();
			break;

		default:
			usage(); break;
		}
	}

	// assert correctness
	if(level < 1) {
		std::cerr << "\nERROR: Search Level must be greater than ZERO. Aborting..." << std::endl;
		usage();
	}
	if(nw < 1) {
		std::cerr << "\nERROR: Number of worker threads must be greater than ZERO. Aborting..." << std::endl;
		usage();
	}
	if(e_lm != 0) {
		if(e_lm < 0.2 || e_lm > 1.0) {
			std::cerr << "\nERROR: Probability threshold for edge plotting must fall within [0.2, 1.0]. Aborting..." << std::endl;
			usage();
		}
	}
	if(gvec.empty() && coords.empty()) {
		std::cerr << "\nERROR: Options -G OR -C are required in order to specify at least one Gene ID. Aborting..." << std::endl;
		usage();
	}
	if(intrs && locs.empty()) {
		std::cerr << "\nERROR: -T options requires a collisions file to be specified. Aborting..." << std::endl;
		usage();
	}

	int cpus = getNumCpus();
	std::cout << "\n System configuration:\n"
				<< "   CPUs: " << cpus << " (physical)\n"
				<< "   RAM:  " << (getTotRAM() >> 30) << " GB (total)" << std::endl;

	bool active=false;
	if(cpus > 8) active=true;

	std::cout << "\n Using " << nw << " worker(s) | Search distance: " << level << "\n"
			<< " Grain of parallel execution: " << chunk << "\n"
			<< " SAM file:\t\t\'" << sam << "\'\n"
			<< " Genes file:\t\t\'" << genes << "\'\n"
			<< " Fragments file:\t\'" << frags << "\'\n"
			<< " Intersections:\t\t\'" << locs << "\'\n\n"
			<< " Starting Gene(s): ";

	for(unsigned i=0; i<gvec.size(); ++i)
		std::cout << gvec[i] << " ";

	std::cout << "\n" << std::endl;

	Finder f;
	f.parseFiles(genes, frags, sam, exprfile);
	if(intrs) {
		f.detectIntersections(locs);
		found = locs.find("/mld");
		mld = locs.substr(found+1, 5);
	}

	if(!coords.empty()) {
		std::vector<Gene> gc = f.getGenesByCoordinates(chr, st, sp);
		for(unsigned g=0; g<gc.size(); ++g)
			gvec.push_back( gc[g].getSymbol() );
	}

	// fix SAM filename for plotting utility
	found = sam.find("ord_");
	if(found == std::string::npos) {
		found = sam.rfind("/");
		sam = sam.substr(found+1, sam.length());
	} else sam = sam.substr(found+4, sam.length());

	found = sam.find(".sam");
	sam = sam.substr(0, found);

	f.parFindConnections(gvec, sam, level, nw, clstr, chunk, e_lm, plot, intrs, mld, active);

	exit(EXIT_SUCCESS);
}
