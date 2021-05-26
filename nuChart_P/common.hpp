/*
 * common.hpp
 *
 *  Created on: Apr 14, 2014
 *      Author: fabio
 */

#ifndef COMMON_HPP_
#define COMMON_HPP_

#include <string>
#include <vector>
#include <utility>
#include <cassert>
#include <set>
#include <unordered_map>
#include <map>
#include <queue>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>

#ifdef _MSC_VER
#include <windows.h>
#else
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/sysinfo.h>
#include <errno.h>
#include <malloc.h>
#endif


namespace {

typedef long int uint_64;

// GENE CLUSTERS
// gene clusters with formatted coordinates.
// add here new gene clusters
void clustersMap(std::unordered_map<std::string, std::string>& cl) {
	cl.reserve(9);

	cl.emplace( "IGH", "chr14:105566277,106879844" );
	cl.emplace( "IGK", "chr2:88857361,89330679" );
	cl.emplace( "IGL", "chr22:22026076,22922913" );
	cl.emplace( "HOXA", "chr7:27126219,27250333" );
	cl.emplace( "HOXB", "chr17:46523832,46827571" );
	cl.emplace( "HOXC", "chr12:54326041,54451948" );
	cl.emplace( "HOXD", "chr2:176955950,177067912" );
	cl.emplace( "HLA", "chr6:29760000,33110000" );
	cl.emplace( "KRAB", "chr19:56300000,59128983" );
}

// random real-numbers generator within range [0, 1]
inline double randomR() {
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(0.0, 1.01);	// works on interval [min, max)

	return dis(gen);
}

// random integers generator within range [0, max]
inline int randomI(long max, long min=0) {
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> dis(min, max);		// works on interval [min, max]

	return dis(gen);
}

// Total physical memory (RAM). This is not the
// available free RAM
unsigned long getTotRAM() {

#ifdef _MSC_VER
	MEMORYSTATUSEX memInfo;
	memInfo.dwLength = sizeof(MEMORYSTATUSEX);
	GlobalMemoryStatusEx(&memInfo);
	DWORDLONG totalPhysMem = memInfo.ullTotalPhys;
	return totalPhysMem;
#else
	struct sysinfo memInfo;

	sysinfo(&memInfo);
	unsigned long totalPhysMem = memInfo.totalram;
	totalPhysMem *= memInfo.mem_unit;
	return totalPhysMem;
#endif
}

// LINUX ONLY
/* Get the number of cpus/cores in the machine
 * if ht==0 checks for hyper-threading flag: in case it is
 * enabled it halves the number of cpus, so that only
 * physical cpus are considered. */
inline int getNumCpus(int ht=0) {

#ifdef _MSC_VER

    SYSTEM_INFO sysinfo;
    GetSystemInfo(&sysinfo);
    return sysinfo.dwNumberOfProcessors;

#else

	FILE *fp;
	int num_cpus=-1;
	char cmd[256];

	// count number of cupus listed in '/sys/devices/system'
	cmd[0] = '\0';
	strcpy(cmd, "ls -d /sys/devices/system/cpu/cpu[0-9]* | wc -l");

	fp = popen(cmd,"r");
	if(fscanf(fp,"%d",&num_cpus)>0)
		fclose(fp);

	if(ht==0) { // exclude hyper-threading
		FILE *f;
		int hyper;
		// check for 'ht' flag in cpus list
		char cmd2[] = "grep '^flags\\b.*: .*\\bht\\b' /proc/cpuinfo | wc -l";

		f = popen(cmd2, "r");
		if(fscanf(f,"%d",&hyper)>0)
			fclose(f);
		if(hyper != 0) // ht enabled
			num_cpus = num_cpus/2;
	}
	return num_cpus;

#endif
}

// check if a string is an EntrezID
inline bool isEntrezID(std::string& s) {
	if( s.empty() || ((!isdigit(s[0]))) ) return false ;

	char * p ;
	long ret = strtol(s.c_str(), &p, 10) ;
	assert(ret > 0);

	return (*p == 0) ; // only digits found
}

// check if a string is an EntrezID
inline bool isEntrezID(const char* s) {
	if( strlen(s)==0 || ((!isdigit(s[0]))) ) return false ;

	char * p ;
	long ret = strtol(s, &p, 10) ;
	assert(ret > 0);

	return (*p == 0) ; // only digits found
}

// Hash function for strings
inline size_t DJBHash(const std::string& str)
{
   size_t hash = 5381;

   for(size_t i = 0; i < str.length(); i++) {
      hash = ((hash << 5) + hash) + str[i];
   }

   return hash;
}

// hash integers using Cantor pairing function
template<typename valT>
inline uint_64 pairZ(const valT v1, const valT v2) {
	return (uint_64) ( (((v1+v2)*(v1+v2+1))/2)+v2+1 );
}

struct hash_symbol {
	size_t operator() (const std::string& str) {
		size_t hash = 5381;

		for(size_t i = 0; i < str.length(); i++)
			hash = ((hash << 5) + hash) + str[i];
		return hash;
	}
};


/* some utilities to be used with containers
 * (free vector of pointers, free memory of vector, max value on map).
 * This should be effective in order to avoid (or minimize) memory leaks
 */

template <class C>
inline size_t maxMapValue(C &umap) {
	typename C::iterator it, end;

	it = umap.begin();
	end= umap.end();

	size_t max_value = it->second;
	for( ; it != end; ++it) {
		if((*it).second > max_value)
			max_value = it->second;
	}
	return max_value;
}

template <typename C>
inline void free_pointed_to(C &cntr) {
	typename C::iterator it;
	for( it = cntr.begin(); it != cntr.end(); ++it )
		delete *it;

	cntr.clear();
}

// clear a vector and free its memory
template <typename C>
inline void freeVector(C &vec) {
	vec.clear();
	C().swap(vec);
}

template <typename T>
inline double vecMean(std::vector<T>& vec) {
	double sum = std::accumulate(vec.begin(), vec.end(), 0.0);
	size_t vec_size = vec.size();
	double m =  sum / (double) vec_size;

	return m;
}

template <typename T>
inline double vecStdDev(std::vector<T>& vec) {
	double sum = std::accumulate(vec.begin(), vec.end(), 0.0);
	size_t vec_size = vec.size();
	double m =  sum / (double) vec_size;

	double accum = 0.0;
	std::for_each(vec.begin(), vec.end(), [&](const T d) {
	    accum += (d - m) * (d - m);
	});
	double stdev = sqrt(accum / (double) (vec_size-1));

	return stdev;
}
// ------------------------------------------------------------------


//// -------- FUNCTIONS FOR LEXICOGRAPHICAL ORDERING IN HUMAN-LIKE MANNER ---------

/**
 * Check if the character is a digit
 */
bool isDigit(const char c) { return c>='0' && c<='9'; }


/**
 * compare l and r with strcmp() semantics. This function
 * is designed to read through the l and r strings only one
 * time, for maximum performance. It does not allocate memory
 * for substrings.
 *
 * Code based on the "Alphanum Algorithm" implementation by
 * Dirk Jagdmann. The algorithm is discussed at http://www.DaveKoelle.com
 *
 * @param l NULL-terminated C-style string
 * @param r NULL-terminated C-style string
 * @return negative if l<r, 0 if l equals r, positive if l>r
 */
inline int compareChar(const char *l, const char *r) {
	enum mode_t {
		STRING, NUMBER
	} mode;

	mode = STRING;

	while(*l && *r)	{
		if(mode == STRING) {
			char l_char, r_char;
			while((l_char=*l) && (r_char=*r)) {
				// check if this are digit characters
				const bool l_digit=isDigit(l_char), r_digit=isDigit(r_char);
				// if both characters are digits, we continue in NUMBER mode
				if(l_digit && r_digit) {
					mode=NUMBER;
					break;
				}
				// if only the left character is a digit
				if(l_digit) return -1;
				// if only the right character is a digit
				if(r_digit) return +1;
				// compute the difference of both characters
				const int diff=l_char - r_char;
				// if they differ we have a result
				if(diff != 0) return diff;
				// otherwise process the next characters
				++l;
				++r;
			}
		} else { // mode==NUMBER
			// get the left number
			char *end;
			unsigned long l_int= strtoul(l, &end, 0);
			l=end;

			// get the right number
			unsigned long r_int= strtoul(r, &end, 0);
			r=end;

			while(*l && isDigit(*l)) {
				l_int=l_int*10 + *l-'0';
				++l;
			}

			while(*r && isDigit(*r)) {
				r_int=r_int*10 + *r-'0';
				++r;
			}

			// if the difference is not equal to zero, we have a comparison result
			const long diff=l_int-r_int;
			if(diff != 0)
				return diff;

			// otherwise we process the next substring in STRING mode
			mode=STRING;
		}
	}

	if(*r) return -1;
	if(*l) return +1;
	return 0;
}


/**
 * Functor class to compare two objects. To be used
 * with std::sort algorithm
 */
template<class Ty>
struct charLess : public std::binary_function<Ty, Ty, bool> {
	bool operator()(const Ty& left, const Ty& right) const	{
		return compareChar(left, right) < 0;
	}
};
// --------------------------------------------------------------------------------
} // namespace

// ------------------------------------------------------------------------------------------------------------------

/*
 * Expression data
 */
struct Expression {

	Expression() : eid(0), logfc(0), pval(0), hs(0) { }
	Expression(size_t _hs, int id, double lg, double pv) :
		eid(id), logfc(lg), pval(pv), hs(_hs) { }
	Expression(const Expression& e) { *this = e; }

	~Expression() {}

	//inline std::string getSymbol() 	const 	{ return symbol; }
	inline long getID()				const 	{ return eid; }
	inline double getPvalue() 		const	{ return pval; }
	inline double getLogFC() 		const	{ return logfc; }
	inline size_t getHS() 			const	{ return hs; }

	bool operator<(const Expression& e) const {	return pval < e.getPvalue(); }
	bool operator=(const Expression& e) const {	return eid == e.getID(); }

	friend std::ostream& operator<<(std::ostream &out, const Expression& e) {
		out << std::fixed;
		out.precision(8);
		out << " EntrezID: "  	<< e.eid 	<< " | " <<
				"logFC: "   	<< e.logfc 	<< " | " <<
				"pval: " 		<< e.pval	<< " | ";

		return out;
	}

	long eid;
	double logfc, pval;
	size_t hs;
};


struct Features {
	//std::string chr_number;
	double start, stop, gcc, map;
	uint_64 len;
	size_t hs;

	Features() {}
	Features(const Features& f) = default;
	Features(Features&& f) = default;

	Features& operator=(const Features& f) = default;
	Features& operator=(Features&& f) = default;

	bool operator==(const Features& f) const { return (hs == f.hs && start == f.start && stop == f.stop); }
	bool operator<(const Features& f) const {
		if(hs < f.hs) return true;
		else if(hs == f.hs && start < f.start) return true;
		else return false;
	}

	friend std::ostream& operator<<(std::ostream &out, const Features& f) {
		out << "Start: "  << f.start 		<< " | " <<
				"Stop: "   << f.stop 		<< " | " <<
				"HS: " << f.hs;

		return out;
	}
};

struct BindingSite {
	std::string seq;
	uint_64 st, sp;
	size_t hs;

	BindingSite() : st(0), sp(0), hs(0) { }
	BindingSite(size_t _hs, std::string s, uint_64 start, uint_64 stop) :
		seq(s), st(start), sp(stop), hs(_hs) { }
	~BindingSite() { }

	inline size_t getHS() { return hs; }
	inline uint_64 getStart() {return st; }
	inline uint_64 getStop() { return sp; }
	//inline std::string getChrom() { return chrom; }

	friend std::ostream& operator<<(std::ostream &out, const BindingSite& b) {
		out <<  "Start: "  << b.st 		<< " | " <<
				"Stop: "   << b.sp 		<< " | " <<
				"Seq: "	   << b.seq		<< " |";

		return out;
	}
};

struct CisTrans {
	long nid, deg, cis, trans;
	double stdev;

	CisTrans() : nid(0), deg(0), cis(0), trans(0), stdev(0.0) {}
	CisTrans(long _nid, long _deg, long _cis, long _trans) : nid(_nid), deg(_deg), cis(_cis), trans(_trans), stdev(0.0) {}
};


#endif /* COMMON_HPP_ */
