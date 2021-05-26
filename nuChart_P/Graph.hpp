/*
 * Graph.hpp
 *
 *  Created on: Aug 26, 2014
 *      Author: droccom
 */

#ifndef GRAPH_HPP_
#define GRAPH_HPP_

#include "common.hpp"
#include "EdgeT.hpp"

struct item {
	item(EdgeT *e) :
		edge(e), dest(0), next(NULL) { }
	item(const item& it) : edge(it.edge), dest(it.dest), next(it.next) { }
	item& operator=(const item& it) {
		if(this != &it) {
			edge = it.edge;
			dest = it.dest;
			if(next)
				delete next;
			next = new item(it);
		}
		return *this;
	}

	EdgeT *edge;
	long dest;
	item *next;
};

struct itemList {
	item *head; 	// head of the adj list
	long degree;  	// degree of this node
	long degU;		// degree in uncertain graph
	long kcore;		// kcore of the node
	long ketacore;	// ketacore of the node
};


class Graph {
public:
	Graph(long n) :
		ngenes(n), mindeg(0), maxdeg(0), numnodes(0) {
		adjlists = new itemList[ngenes];
		for (long i = 0; i < ngenes; ++i) {
			adjlists[i].head = NULL;
			adjlists[i].degree = 0;
			adjlists[i].degU = 0;
			adjlists[i].kcore = -1;
			adjlists[i].ketacore = -1;
		}
		edges.reserve(2*n);
	}

	// copy constructor
	Graph(const Graph& g) : edges(g.getEdges()), ngenes(g.nodesSpace()), mindeg(g.getMinDegree()),
			maxdeg(g.getMaxDegree()), numnodes(g.numnodes) {
		adjlists = new itemList[g.ngenes];
		memcpy(adjlists, g.adjlists, g.ngenes);
	}

	// copy assignment operator
	Graph& operator=(const Graph& g) {
		if(this != &g) {
			if(adjlists)
				clear();
			delete[] adjlists;

			ngenes = g.nodesSpace();
			mindeg = g.getMinDegree();
			maxdeg = g.getMaxDegree();
			numnodes = g.getNumNodes();
			edges = g.getEdges();

			adjlists = new itemList[g.ngenes];
			memcpy(adjlists, g.adjlists, g.ngenes);
		}
		return *this;
	}

	~Graph() {
		clear();
		delete[] adjlists;
		for (unsigned long i = 0; i < edges.size(); ++i)
			free(edges[i]);
		edges.clear();
	}

	void clear() {
		for (long i=0; i<ngenes; ++i) {
			item *delptr = adjlists[i].head;
			if(delptr) {
				while (delptr != NULL) {
					item *itm = delptr->next;
					delete delptr;
					delptr = itm;
				}
				adjlists[i].head = NULL;
				adjlists[i].degree = 0;
				adjlists[i].degU = 0;
				adjlists[i].kcore = -1;
				adjlists[i].ketacore = -1;
			}
		}
		updateNodes();
	}

	/* add an edge to the graph, if not already present
	 * (symmentry is checked within EdgeT operator==)
	 *
	 * adr prevents dangerous remove of edges when performing
	 * average degree rewiring */
	int add_edge(EdgeT *e, bool reduce, bool adr=false) {
		//search for duplicated
		item *addptr = NULL;
		int r=-1;
		long src_idx = e->id1;
		long dst_idx = e->id2;

		addptr = adjlists[src_idx].head;
		while (addptr != NULL) {
			if ((*addptr->edge) == *e) {
				if(!adr) delete e;
				return 0;
			}
			addptr = addptr->next;

		}

		addptr = new item(e);
		addptr->dest = dst_idx;
		addptr->next = adjlists[src_idx].head;
		adjlists[src_idx].head = addptr;
		++adjlists[src_idx].degree;
		if(adjlists[src_idx].degree > maxdeg)
			maxdeg = adjlists[src_idx].degree;
		if (mindeg == 0 || adjlists[src_idx].degree < mindeg)
			mindeg = adjlists[src_idx].degree;

		addptr = new item(e);
		addptr->dest = src_idx;
		addptr->next = adjlists[dst_idx].head;
		adjlists[dst_idx].head = addptr;
		++adjlists[dst_idx].degree;
		if(adjlists[dst_idx].degree > maxdeg)
			maxdeg = adjlists[dst_idx].degree;
		if (mindeg == 0 || adjlists[dst_idx].degree < mindeg)
			mindeg = adjlists[dst_idx].degree;

#ifdef TEST3
		if(reduce) {
			if(!e->getVertex2()->getGraphRoot()) {
				e->getVertex2()->setGraphRoot(e->getVertex1()->getGraphRoot());
			}
			r = e->fill();
		}
#else
		if(reduce) {
			e->getVertex2()->setGraphRoot(e->getVertex1()->getGraphRoot());
			r = e->fill();
		}
#endif

		edges.push_back(e);
		updateNodes();
		return r;
	}

	// function that removes an edge
	void remove_edge(EdgeT *e) {
		long e_idx = e->id1;
		long l_idx = e->id2;

		item *delptr = adjlists[e_idx].head;
		while (delptr) {
			if(delptr->dest == l_idx) {
				--adjlists[e_idx].degree;
				adjlists[e_idx].degU -= e->getProb();
				if(adjlists[e_idx].degree == 0) {
					adjlists[e_idx].head = NULL;
					adjlists[e_idx].degU = 0;
					adjlists[e_idx].kcore = -1;
					adjlists[e_idx].ketacore = -1;
					delete delptr;
					break;
				}
				if(delptr->next)
					delptr->next = delptr->next->next;
				break;
			} else delptr = delptr->next;
		}

		delptr = adjlists[l_idx].head;
		while (delptr) {
			if(delptr->dest == e_idx) {
				--adjlists[l_idx].degree;
				adjlists[l_idx].degU -= e->getProb();
				if(adjlists[l_idx].degree == 0) {
					adjlists[l_idx].head = NULL;
					adjlists[l_idx].degU = 0;
					adjlists[l_idx].kcore = -1;
					adjlists[l_idx].ketacore = -1;
					delete delptr;
					break;
				}
				if(delptr->next)
					delptr->next = delptr->next->next;
				break;
			} else delptr = delptr->next;
		}

		std::vector<EdgeT*>::iterator it = std::find(edges.begin(), edges.end(), e);
		if(it != edges.end())
			edges.erase(it);

		updateNodes();
	}

	void remove_node(long nid) {
		if(nid<ngenes) {
			item *delptr = adjlists[nid].head;
			while (delptr) {
				remove_edge(delptr->edge);
				if(!adjlists[nid].degree)
					break;
				else delptr = delptr->next;
			}
			adjlists[nid].head = NULL;
			adjlists[nid].degree = 0;
			adjlists[nid].degU = 0;
			adjlists[nid].kcore = -1;
			adjlists[nid].ketacore = -1;
		}
		updateNodes();
	}

	void setEdgesProb(std::vector<miniEdge>& me) {
		for(unsigned i=0; i<me.size(); ++i) {
			edges[i]->setProb(me[i].prob);
			adjlists[edges[i]->getVertex1()->getId()].degU += me[i].prob;
			adjlists[edges[i]->getVertex2()->getId()].degU += me[i].prob;
		}
	}

	void setNodeKcore(long nid, long kc) {
		adjlists[nid].kcore = kc;
	}

	void setNodeKetacore(long nid, long kc) {
		adjlists[nid].ketacore = kc;
	}

	// compute Pr[deg(v) == i] for keta-core
	double computeX(itemList *node, long i) {
		long deg = node->degree;
		if((i >= deg) || (deg != 0 && i == 0)) return 0;
		if((deg == i) ) return 1.;

		item *itm = node->head;
		std::vector<double> probs;

		while(itm) {
			if(itm->edge != NULL) {
				probs.push_back(itm->edge->getProb());
				itm = itm->next;
			} else itm = itm->next;
		}

		double **x = new double*[deg];
		for(unsigned l=0; l<deg; ++l )
			x[l] = new double[i+1]();

		for(long h=0; h<deg; ++h) {
			for(long j=0; j<=i; ++j) {
				if(h==0 && j==0)
					x[h][j] = 1;
				else if(j > h || (h != 0 && j==0))
					x[h][j] = 0;
				else
					x[h][j] = (probs[h] * x[h-1][j-1]) + ((1-probs[h]) * x[h-1][j]);
			}
		}
		return x[deg-1][i];
	}

	// Compute Pr[deg(v) >= k] for keta-core
	long computeEtaDegree(long nid, double eta) {
		if(adjlists[nid].degree == 0) return 0L;

		double *pr = new double[adjlists[nid].degree+1]();
		pr[0] = 1.;

		for(long k=1; k<=adjlists[nid].degree; ++k) {
				pr[k] = pr[k-1] - computeX(&adjlists[nid], k-1);
				if(pr[k] < eta) {
					delete[] pr;
					return k-1;
				}
		}
		delete[] pr;
		return adjlists[nid].degree;
		// eta_d is max k € [0..d_v] | Pr[deg(v) >= k] >= eta
	}

	/*
	 * Compute Pr[deg(v|~e) = i] for each i = 0, ..., eta_d[v].
	 * prNE[0] means Pr[deg(v|~e)=0] for keta-core
	 *
	 * prNE[0] = 1/(1-p_e)*Pr[deg(v)=0]
	 */
	long updateEtaDegree(long nid1, long nid2, double eta) {
		if(adjlists[nid1].degree == 0) return 0L;

		double pe=0.;
		item *itm = adjlists[nid1].head;

		while(itm) {
			if(itm->edge->id2 == nid2) {
				pe = itm->edge->getProb();
				break;
			} else itm = itm->next;
		}

		double *prNE = new double[adjlists[nid1].degree+1]();
		prNE[0] = (1/(1-pe)); //* computeX(&adjlists[nid1], 0);

		for(long k=1; k<=adjlists[nid1].degree; ++k) {
			prNE[k] = ( computeX(&adjlists[nid1], k) - (pe*prNE[k-1]) ) / (1 - pe);
			if(/*prNE[k] > 0. &&*/ std::abs(prNE[k]) < eta ) {
				delete[] prNE;
				return k-1;
			}
		}
		delete[] prNE;
		return adjlists[nid1].degree-1;
	}

	void getNodesList(long nid, std::vector<long>& nds) {
		if(adjlists[nid].head != NULL) {
			nds.reserve(adjlists[nid].degree);
			item* itm = adjlists[nid].head;
			while(itm) {
				nds.push_back(itm->dest);
				itm = itm->next;
			}
		}
	}

	long* getNodesList(long nid) {
		long *nds = NULL;
		if(adjlists[nid].head != NULL) {
			item* itm = adjlists[nid].head;
			nds = new long[adjlists[nid].degree];
			unsigned i=0;
			while(itm) {
				nds[i] = itm->dest;
				itm = itm->next;
				++i;
			}
		}
		return nds;
	}

	// select a random edge from edges set
	EdgeT *getRandomEdge(long nid) {
		if(adjlists[nid].head != NULL ) {
			int i=0, target=0;
			item* itm = adjlists[nid].head;
			EdgeT *ed = NULL;
			if(adjlists[nid].degree > 1)
				target = randomI(adjlists[nid].degree-1);
			else target = 1;

			while(itm && i<target) {
				ed = itm->edge;
				itm = itm->next;
				++i;
			}
			return ed;
		}
		return NULL;	// should never be reached
	}

	long getNodeIdDegreeMax() const {
		for(unsigned i=0; i<ngenes; ++i)
			if(adjlists[i].degree == maxdeg)
				return i;
		return -1; // never reached
	}

	int* getDegreeDistrib() {
		long max_deg = getMaxDegree(), cnt;
		int *ddc = new int[max_deg+1];

		for(unsigned j=0; j<=max_deg; ++j) {
			cnt=0;
			for(unsigned i=0; i<ngenes; ++i)
				if(adjlists[i].head && adjlists[i].degree == j)
					++cnt;
			ddc[j] = cnt;
		}
		return ddc;
	}

	bool isNode(unsigned id)	const {
		if(id < ngenes)
			return adjlists[id].head != NULL;
		return false;
	} // true if 'id' is a valid node

	long getNodeDegree(long n_id) 		const { return adjlists[n_id].degree; }
	double getNodeUDegree(long n_id) 	const { return adjlists[n_id].degU; }
	long getNumNodes() 					const { return numnodes; }

	double getAverageDegree() 			const { return (2*edges.size())/numnodes; }
	long getMaxDegree() 				const {	return maxdeg; }
	long getMinDegree() 				const { return mindeg; }
	double getGraphDensity()			const { return (2*edges.size()) / (numnodes*(numnodes-1)); }

	inline size_t edgesSize() 				const { return edges.size(); }
	inline long nodesSpace() 				const {return ngenes; }
	inline std::vector<EdgeT *> getEdges() 	const { return edges; }
	inline EdgeT* getEdge(unsigned i) 		const { return edges[i]; }


	// ----------------------- PRETTY PRINT ----------------------------

	void printGraph(std::ostream &os) {
		os << "PRINTING GRAPH.....\n";
		os << getNumNodes() << " nodes, " << edges.size() << " edges\n";
		for (long i = 0; i < ngenes; ++i) {
			if(adjlists[i].head) {
				item* itm = adjlists[i].head;
				os << i;
				while(itm) {
					os << " -> " << itm->dest;
					itm = itm->next;
				}
				os << std::endl;
			}
		}
		os << "end GRAPH\n";
	}

	void printEdges(std::ostream &os) {
		for(unsigned i=0; i<edges.size(); ++i)
			os << edges[i]->id1 << " -- " << edges[i]->id2 << std::endl;
	}

	void printNodesKcore(std::ostream &os) {
		os << "K-CORES (!= 0).....\n";
		for(unsigned i=0; i<ngenes; ++i)
				if(adjlists[i].kcore != 0)
					os << "Node " << i << " => k-core: " << adjlists[i].kcore
					<< "; degree: " << adjlists[i].degree << "\n";
		os << std::endl;
	}

	void printNodesInKcore(std::ostream &os) {
		os << "NODES IN K-CORES.....\n";
		for(unsigned i=1; i<=maxdeg; ++i) {
			os << "[k-core " << i << "]:= ";
			for(unsigned j=0; j<ngenes; ++j)
				if(adjlists[j].kcore == i)
					os << j << "; ";
			os << "\n";
		}
		os << std::endl;
	}

	void printNodesKetacore(std::ostream &os) {
		os << "KETA-CORES (!= 0).....\n";
		for(unsigned i=0; i<ngenes; ++i)
			if(adjlists[i].ketacore != 0)
				os << "Node " << i << " => keta-core " << adjlists[i].ketacore << "\n";
		os << std::endl;
	}

	void printNodesInKetacore(std::ostream &os) {
		os << "NODES IN KETA-CORES.....\n";
		for(unsigned i=1; i<=maxdeg; ++i) {
			os << "[keta-core " << i << "]:= ";
			for(unsigned j=0; j<ngenes; ++j)
				if(adjlists[j].ketacore == i)
					os << j << "; ";
			os << "\n";
		}
		os << std::endl;
	}

private:
	void updateNodes() {
		numnodes=0;
		for(unsigned i=0; i<ngenes; ++i)
			if(adjlists[i].head)
				++numnodes;
	}

	itemList *adjlists;
	std::vector<EdgeT *> edges;
	long ngenes, mindeg, maxdeg, numnodes;
};

#endif /* GRAPH_HPP_ */
