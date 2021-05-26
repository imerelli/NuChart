/*
 * PrintOut.hpp
 *
 *  Created on: Sep 15, 2014
 *      Author: fabio
 */

#ifndef PRINTOUT_HPP_
#define PRINTOUT_HPP_

#include "common.hpp"
#include "Timings.hpp"
#include "FileIO.hpp"
#include "Edge.hpp"
#include "Graph.hpp"

namespace {

void printEdgesList(std::vector<Edge*>& edges, std::string start, std::string dts, int sc_limit,
		int wrks, std::string out_folder="") {
	std::stringstream ess;
	std::string fold("_csv");

	start.append(fold);
	out_folder.append(start);
	mkdir(out_folder.c_str(), 0777);

	ess << out_folder << "/edgesLis_" << dts << "_L" << sc_limit << "_" << tmn::timeStamp() << ".csv";
	std::ofstream e_file(ess.str().c_str(),
			std::ios_base::out | std::ios_base::trunc );

	e_file << "#Gene1\tGene2" << std::endl;
	for(unsigned i=0; i<edges.size(); ++i) {
		if( edges[i]->getWeight() != -2  ) {
			e_file << edges[i]->getGeneSymbol1() << " -- "
					<< edges[i]->getGeneSymbol2() << std::endl;
		}
	}
	e_file.close();
}

void printEdgeListWeighted(std::vector<Edge*>& edges, std::string start, std::string dts, int sc_limit,
		int wrks, std::string out_folder="") {
	std::stringstream ess;
	std::string fold("_csv");

	start.append(fold);
	out_folder.append(start);
	mkdir(out_folder.c_str(), 0777);

	ess << out_folder << "/edgesWeight_" << dts << "_L" << sc_limit << "_" << tmn::timeStamp() << ".csv";
	std::ofstream e_file(ess.str().c_str(),
			std::ios_base::out | std::ios_base::trunc );

	e_file << "#Gene1--Gene2\tWeight\tProbability" << std::endl;
	e_file << std::fixed;
	e_file.precision(6);
	for(unsigned i=0; i<edges.size(); ++i) {
		if( edges[i]->getWeight() != -2  ) {
			e_file << edges[i]->getGeneSymbol1() << "--"
					<< edges[i]->getGeneSymbol2() << "\t"
					<< edges[i]->getWeight() << "\t"
					<< edges[i]->getProb() << std::endl;
		}
	}
	e_file.close();
}

// csv with edges and weigth
void printCSV1(std::vector<Edge*>& edges, std::string start, std::string dts, int sc_limit,
		int wrks, float e_lm, std::string out_folder="csv/") {
	std::stringstream ess;
	std::string fold("_csv"), star;

	star = start;
	start.append(fold);
	out_folder.append(start);
	mkdir(out_folder.c_str(), 0777);

	ess << out_folder << "/edges_" << dts << "_L" << sc_limit << "_th" << e_lm << "_" << tmn::timeStamp() << ".csv";
	std::ofstream e_file(ess.str().c_str(),
			std::ios_base::out | std::ios_base::trunc );

	e_file << "#Root Gene: " << star << "\t Tot. Edges: " << edges.size() << std::endl;
	e_file << "#Edge\tWeight\tProb ( > " << e_lm << " )\tCis/Trans" << std::endl;
	e_file << std::fixed;
	e_file.precision(6);
	for(unsigned i=0; i<edges.size(); ++i) {
		e_file << edges[i]->getGeneSymbol1() << "--" << edges[i]->getGeneSymbol2()
						<< "\t" << edges[i]->getWeight()
						<< "\t" << edges[i]->getProb() << "\t";
		if(edges[i]->getHS1() == edges[i]->getHS2())
			e_file << "CIS" << std::endl;
		else
			e_file << "TRANS" << std::endl;
	}
	e_file.close();
}

// csv with nodes and level
void printNodesWithLevel(std::vector<bool>& vm, std::vector<Node>& genes, size_t tot, std::string start, std::string dts,
		int sc_limit, int wrks, std::string out_folder="") {
	std::stringstream ess;
	std::string fold("_csv");

	start.append(fold);
	out_folder.append(start);
	mkdir(out_folder.c_str(), 0777);

	ess << out_folder << "/nodesLevel_" << dts << "_L" << sc_limit << "_" << tmn::timeStamp() << ".csv";
	std::ofstream e_file(ess.str().c_str(),
			std::ios_base::out | std::ios_base::app );

	e_file << "#Root ( " << start << " )" << std::endl;
	e_file << "#GENE\tLEVEL" << std::endl;
	for(unsigned v=0; v<vm.size(); ++v) {
		if(vm[v]) {
			e_file << genes[v].label << "\t"
					<< genes[v].lvl << std::endl;
		}
	}
	e_file.close();
}

void printTimesForCharts(size_t edgs, double time, std::string start, std::string dts,
		int sc_limit, int wrks, float e_lm, std::string out_folder="") {
	std::stringstream ess;
	std::ofstream e_file;
	std::string fold("_csv");

	start.append(fold);
	out_folder.append(start);
	mkdir(out_folder.c_str(), 0777);

	ess << out_folder << "/edgesTime_" << dts << "_L" << sc_limit << "_th" << e_lm << ".csv";
	if(!exist_f(ess.str().c_str())) {
		e_file.open(ess.str().c_str(),
				std::ios_base::out | std::ios_base::app );
		e_file << "#NUM_WRKRS" << "\t" << "TOT_TIME FOR " << edgs << "EDGES" << std::endl;
		e_file << wrks << "\t" << time << " ms." << std::endl;
	} else {
		e_file.open(ess.str().c_str(),
				std::ios_base::out | std::ios_base::app );
		e_file << wrks << "\t" << time << " ms." << std::endl;
	}
	e_file.close();
}

void printDegreeDistribution(Graph *g, size_t tot, std::string start, std::string dts,
		int sc_limit, int wrks, std::string out_folder="") {
	std::stringstream ess;
	std::ofstream e_file;
	std::string fold("_csv");

	start.append(fold);
	out_folder.append(start);
	mkdir(out_folder.c_str(), 0777);

	long degM = g->getMaxDegree();
	double avg = g->getAverageDegree();
	int *ddc = g->getDegreeDistrib();

	ess << out_folder << "/degreeDistrib_" << dts << "_L" << sc_limit << "_" << tmn::timeStamp() << ".csv";
	e_file.open(ess.str().c_str(),
			std::ios_base::out | std::ios_base::trunc );

	e_file << "#Degree max: " << degM << "\t"
			<< "#Average degree: " << avg << "\n"
			<< "#Tot Nodes: " << tot << "\t"
			<< "#Tot Edges: " << g->edgesSize() << std::endl;

	e_file << "#NumNodes" << "\t" << "Degree" << std::endl;
	for(unsigned i=0; i<=degM; ++i)
		if(ddc[i])
			e_file << ddc[i] << "\t" << i << "\n";
	e_file << std::endl;

	e_file.close();
}

void printNodesDegree(Graph *g, std::vector<Node>& vertices, size_t tot, std::string start, std::string dts,
		int sc_limit, int wrks, std::string out_folder="") {
	std::stringstream ess;
	std::ofstream e_file;
	std::string fold("_csv");

	start.append(fold);
	out_folder.append(start);
	mkdir(out_folder.c_str(), 0777);

	long degMax = g->getMaxDegree();
	double avg = g->getAverageDegree();

	ess << out_folder << "/nodesDegree_" << dts << "_L" << sc_limit << "_" << tmn::timeStamp() << ".csv";
	e_file.open(ess.str().c_str(),
			std::ios_base::out | std::ios_base::trunc );

	e_file << "#Degree max: " << degMax << "\t"
			<< "#Average degree: " << avg << "\n"
			<< "#Tot Nodes: " << tot << "\t"
			<< "#Tot Edges: " << g->edgesSize() << std::endl;

	e_file << "#NodeId" << "\t" << "Name" << "\t" << "Degree" << "\t" << "Expression" << std::endl;
	for(unsigned i=0; i<vertices.size(); ++i)
		e_file << vertices[i].n_id << "\t" << vertices[i].label << "\t" << vertices[i].deg << "\t" << vertices[i].gExp << "\n";
	e_file << std::endl;

	e_file.close();
}

void printEdgesStat(std::vector<Edge*>& edges, std::vector<Gene>& vertices, std::string start,
		std::string dts, std::string locs, int sc_limit, int wrks, std::string out_folder="csv/") {
	std::stringstream edgStat;

	edgStat << out_folder << "/edgesStats_" << dts << "_L" << sc_limit << "_" << locs << "_" << tmn::timeStamp() << ".csv";

	std::ofstream eStat_file(edgStat.str().c_str(),
			std::ios_base::out | std::ios_base::trunc );

	eStat_file << "#Edge\tScore\tProb\tCis/Trans\tCollisions" << std::endl;
	eStat_file << std::fixed;
	eStat_file.precision(8);
	for (unsigned i = 0; i < edges.size(); ++i) {
		long id1 = edges[i]->getVertex1()->getId();
		long id2 = edges[i]->getVertex2()->getId();

		eStat_file << edges[i]->getGeneSymbol1() << "--"
				<< edges[i]->getGeneSymbol2() << "\t"
				<< edges[i]->getWeight() << "\t" << edges[i]->getProb() << "\t";

		if (edges[i]->getVertex1()->getGraphRoot() == edges[i]->getVertex2()->getGraphRoot()) {
			if(edges[i]->getHS1() == edges[i]->getHS2()) {// intra-chromosome
				eStat_file << "CIS";
			} else {
				eStat_file << "TRANS";
			}
		}

		if(locs == "")
			eStat_file << "\tNaN" << std::endl;
		else
			eStat_file << "\t" << vertices[id1].getIntr() << "--"
			<< vertices[id2].getIntr() << std::endl;
	}
	eStat_file.close();
}


int plotGraph(std::string dotFile, std::string ext, std::string l) {
	//std::string layoutExe = getGvLayout(l);
	std::string gVpath = "/usr/bin/", tmp;

	std::string::size_type found1, found2;
	found1 = dotFile.rfind(".");
	tmp = dotFile.substr(0, found1);
	tmp.append("_").append(l).append(".").append(ext);

	// FIX FORMAT
	std::stringstream exec_cmd, full_cmd;
	exec_cmd << l << " -T" << ext << " " << dotFile << " -o " << tmp << " -q";

	if (system(exec_cmd.str().c_str())== -1) {
		gVpath.append(exec_cmd.str());
		if (system(gVpath.c_str())==-1) {
			//			full_cmd.clear(); full_cmd.str(std::string());
			//			full_cmd << "/home/tordini/gviz/bin/" << exec_cmd.str();
			//			std::cout << full_cmd.str() << std::endl;
			//			if (system(full_cmd.str().c_str())==-1) {
			std::cerr << "ERROR: cannot find GraphViz. "
					<< "Locate your Graphviz installation and execute:\n"
					<< gVpath << std::endl;
			return -1;
			//			}
		}
	}

	return 0;
}

/*
 * PRINT VERTICES AND EDGES FOR PLOTTING WITH GRAPHVIZ
 */
std::string plotGViz(Graph *g, std::vector<Edge*>& edges, std::vector<Node>& vertices, std::string start,
		std::string dts, int sc_limit, int wrks) {
	std::stringstream ess;

	std::string out_folder("plots"), baseFolder(tmn::timeStamp()), fold("_gViz");
	mkdir(out_folder.c_str(), 0777);
	out_folder.append("/").append(baseFolder);
	mkdir(out_folder.c_str(), 0777);

	start.append(fold);
	out_folder.append("/").append(start);
	mkdir(out_folder.c_str(), 0777);

	ess << out_folder << "/" << dts << "_L" << sc_limit << "_" << tmn::timeStamp() << ".dot";
	std::ofstream e_file(ess.str().c_str(),
			std::ios_base::out | std::ios_base::trunc );

	e_file << "/*****\n" << ess.str() << "\n*****/\n\n" << std::endl;
	e_file << "graph G {" << std::endl;
	e_file << " graph [splines=line, overlap=compress, outputorder=edgesfirst, sep=3]\n"
			<< " node [shape=ellipse, style=filled, fixedsize=true]\n" << std::endl;


	double avgdeg = g->getAverageDegree();	//long maxdeg = g->getDegreeMax();
	int minPlotDegree = vertices.size() > 500 ? 2 : 1;
	minPlotDegree = minPlotDegree > avgdeg ? avgdeg : minPlotDegree;
	// node colors
	for (unsigned i = 0; i < vertices.size(); ++i) {
		if(vertices[i].deg > minPlotDegree) {
			e_file << " " << vertices[i].n_id << " [";
			if(vertices[i].isRoot) {
				e_file << "width=0.5, height=0.5, penwidth=1.0, fillcolor=\"yellow\", label=\""
						<< vertices[i].label << "\"," << " fontsize=10, fontname=\"times bold\"";
			} else if(vertices[i].deg >= avgdeg && vertices[i].deg < 2*avgdeg) { // average degree nodes
				e_file << "width=0.35, height=0.35, penwidth=0.6, fillcolor=\"navajowhite\", label=";
				if(vertices.size() < 500)
					e_file << "\"" << vertices[i].label << "\"," << " fontsize=5, fontname=\"times\"";
				else
					e_file << "\"\"";
			} else if(vertices[i].deg >= 2*avgdeg) { // higher degree nodes
				e_file << "width=0.7, height=0.7, penwidth=0.8, fillcolor=\"orange1\", label=\"" << vertices[i].label << "\", fontsize=8, fontname=\"times\"";
			} else { // less influential
				e_file << "width=0.12, height=0.12, penwidth=0.4, fillcolor=\"grey95\", label=\"\"";
			}
			e_file << "];" << std::endl;
		}
	}

	for (unsigned i = 0; i < edges.size(); ++i) {
		long dg1 = g->getNodeDegree( edges[i]->getVertex1()->getId() );
		long dg2 = g->getNodeDegree( edges[i]->getVertex2()->getId() );

		if(dg1 > minPlotDegree && dg2 > minPlotDegree) {
			e_file << " " << edges[i]->getVertex1()->getId() << " -- " << edges[i]->getVertex2()->getId();
		} else if(dg1 > 1 && dg2 <= minPlotDegree /*&& vertices.size() < 500*/) {
			e_file << " " << edges[i]->getVertex2()->getId() << " ["
					<< "width=0.12, height=0.12, penwidth=0.4, fillcolor=\"grey95\"" << ", label=\"\"];\n";
			e_file << " " << edges[i]->getVertex1()->getId() << " -- " << edges[i]->getVertex2()->getId();
		} else if(dg1 <= minPlotDegree && dg2 > 1 /*&& vertices.size() < 500*/) {
			e_file << " " << edges[i]->getVertex1()->getId() << " ["
					<< "width=0.12, height=0.12, penwidth=0.4, fillcolor=\"grey95\"" << ", label=\"\"];\n";
			e_file << " " << edges[i]->getVertex1()->getId() << " -- " << edges[i]->getVertex2()->getId();
		} else if(dg1 <= minPlotDegree && dg2 <= minPlotDegree) continue;

		if ( edges[i]->hasIntergenic() ) {
			e_file << "[style=dotted, penwidth=0.6, color=\"chartreuse\"];\n";
			continue;
		}

		double pwidth = edges[i]->getProb();
		if( edges[i]->getHS1() == edges[i]->getHS2() ) // cis/trans
			e_file << " [color=\"blueviolet\"";
		else
			e_file << " [color=\"palegreen\"";

		if (pwidth != 0) {
			//double pwidth = std::fmod(edges[i]->getWeight(), 3.0);
			if(pwidth <= 0.6)
				e_file << ", penwidth=" << 2*pwidth << "];\n";
			else
				e_file << ", style=\"bold\", penwidth=" << 2.2*pwidth << "];\n";
		} else {
			e_file << ", penwidth=0.3];\n";
		}
	} // for-edges

	e_file << "\n";

	e_file << "}" << std::endl;
	e_file.close();

	// ----- PLOT IT -----
	std::vector<std::string> lts;
	lts.push_back("fdp"); lts.push_back("neato");
	for(unsigned i=0; i<lts.size(); ++i) {
		if(vertices.size() < 3500) {
			plotGraph(ess.str(), "svg", lts[i]);
			plotGraph(ess.str(), "png", lts[i]);
		}
		else {
			std::stringstream exec_cmd;
			std::string es_file(ess.str()), tmp;
			std::string::size_type found1, found2;
			found1 = es_file.rfind(".");
			tmp = es_file.substr(0, found1);
			tmp.append("_").append(lts[i]).append(".png");
			exec_cmd << lts[i] << " -Tpng" << " " << ess.str() << " -o " << tmp << " -q" << std::endl;
			tmp.clear();
			tmp = es_file.substr(0, found1);
			tmp.append("_").append(lts[i]).append(".svg");
			exec_cmd << lts[i] << " -Tsvg" << " " << ess.str() << " -o " << tmp << " -q" << std::endl;
			std::cerr << "Very big graph: dot file created, here is the "
					<< "plotting command. see dot manual for tuning parameters.\n\n"
					<< exec_cmd.str() << std::endl;
		}
	}

	return out_folder; //baseFolder;
}

/*
 * PRINT GRAPH IN GRAPHML FORMAT (TO BE USED WITH GEPHI)
 * TODO
 */
std::string plotGraphML(Graph *g, std::vector<Edge*>& edges, std::vector<Node>& vertices, std::string start,
		std::string dts, int sc_limit, int wrks) {
	std::stringstream ess;

	std::string out_folder("plots"), baseFolder(tmn::timeStamp()), fold("_graphML");
	mkdir(out_folder.c_str(), 0777);
	out_folder.append("/").append(baseFolder);
	mkdir(out_folder.c_str(), 0777);

	start.append(fold);
	out_folder.append("/").append(start);
	mkdir(out_folder.c_str(), 0777);

	ess << out_folder << "/" << dts << "_L" << sc_limit << "_" << tmn::timeStamp() << ".graphml";
	std::ofstream e_file(ess.str().c_str(),
			std::ios_base::out | std::ios_base::trunc );

	e_file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
	e_file << "<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\" ";
	e_file << "xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" ";
	e_file << "xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd\">" << std::endl;

	e_file  << "  <key id=\"d0\" for=\"node\" attr.name=\"color\" attr.type=\"string\">\n"
			<< "    <default>silver</default>\n"
			<< "  </key>" << std::endl;

	e_file  << "  <key id=\"d1\" for=\"edge\" attr.name=\"weight\" attr.type=\"double\">\n"
			<< "    <default>1.0</default>\n"
			<< "  </key>" << std::endl;

	e_file  << "  <key id=\"d2\" for=\"edge\" attr.name=\"Rsequence\" attr.type=\"string\">\n"
			<< "    <default>ATGC</default>\n"
			<< "  </key>" << std::endl;

	e_file  << "  <key id=\"d3\" for=\"edge\" attr.name=\"Lsequence\" attr.type=\"string\">\n"
			<< "    <default>CGTA</default>\n"
			<< "  </key>" << std::endl;

	e_file  << "  <key id=\"d4\" for=\"edge\" attr.name=\"color\" attr.type=\"string\">\n"
			<< "    <default>black</default>\n"
			<< "  </key>" << std::endl;

	e_file  << "  <graph id=\"G\" edgedefault=\"undirected\"\n"
			<< "         parse.nodes=\"" << g->getNumNodes() << "\" parse.edges=\"" << edges.size() << "\"\n"
			<< "         parse.nodeids=\"free\" parse.edgeids=\"free\"\n"
			<< "         parse.order=\"free\" >" << std::endl;

	double avgdeg = g->getAverageDegree();	//long maxdeg = g->getDegreeMax();
	int minPlotDegree = vertices.size() > 500 ? 2 : 1;
	minPlotDegree = minPlotDegree > avgdeg ? avgdeg : minPlotDegree;

	// nodes first
	for (unsigned i = 0; i < vertices.size(); ++i) {
		if(vertices[i].deg > minPlotDegree) {
			e_file << "    <node id=\"" << vertices[i].n_id << "\">" << std::endl;
			if(vertices[i].isRoot) {
				e_file  << "      <data key=\"d0\">yellow</data>\n"
						<< "      <data key=\"label\">" << vertices[i].label << "</data>\n";
			} else if(vertices[i].deg >= avgdeg && vertices[i].deg < 2*avgdeg) { // average degree nodes
				e_file  << "      <data key=\"d0\">NavajoWhite</data>\n";
				if(vertices.size() < 500)
					e_file << "      <data key=\"label\">" << vertices[i].label << "</data>\n";
			} else if(vertices[i].deg >= 2*avgdeg) { // higher degree nodes
				e_file  << "      <data key=\"d0\">orange</data>\n"
						<< "      <data key=\"label\">" << vertices[i].label << "</data>\n";
			} else { // less influential
				e_file  << "      <data key=\"d0\">silver</data>\n";
			}
			e_file << "    </node>" << std::endl;
		}
	}

	for (unsigned i = 0; i < edges.size(); ++i) {
		long dg1 = g->getNodeDegree( edges[i]->getVertex1()->getId() );
		long dg2 = g->getNodeDegree( edges[i]->getVertex2()->getId() );
		double pwidth = edges[i]->getProb();

		if(dg1 > minPlotDegree && dg2 > minPlotDegree) {
			e_file  << "    <edge id=\"e"<< i << "\" source=\"" << edges[i]->getVertex1()->getId()
					<< "\" target=\"" << edges[i]->getVertex2()->getId() << "\">" << std::endl;

			if( edges[i]->getHS1() == edges[i]->getHS2() ) // cis/trans
				e_file  << "      <data key=\"d4\">blueviolet</data>\n";
			else
				e_file  << "      <data key=\"d4\">palegreen</data>\n";

			e_file  << "      <data key=\"d1\">"<< pwidth << "</data>\n";
			e_file  << "      <data key=\"d2\">"<< edges[i]->getSeq1() << "</data>\n";
			e_file  << "      <data key=\"d3\">"<< edges[i]->getSeq2() << "</data>" << std::endl;

		} else if(dg1 > 1 && dg2 <= minPlotDegree /*&& vertices.size() < 500*/) {
			e_file  << "    <node id=\"" << edges[i]->getVertex2()->getId() << "\">\n";
			e_file  << "      <data key=\"d0\">silver</data>\n";
			e_file  << "    </node>" << std::endl;

			e_file  << "    <edge id=\"e"<< i << "\" source=\"" << edges[i]->getVertex1()->getId()
					<< "\" target=\"" << edges[i]->getVertex2()->getId() << "\">" << std::endl;

			if( edges[i]->getHS1() == edges[i]->getHS2() ) // cis/trans
				e_file  << "      <data key=\"d4\">blueviolet</data>\n";
			else
				e_file  << "      <data key=\"d4\">palegreen</data>\n";

			e_file  << "      <data key=\"d1\">"<< pwidth << "</data>\n";
			e_file  << "      <data key=\"d2\">"<< edges[i]->getSeq1() << "</data>\n";
			e_file  << "      <data key=\"d3\">"<< edges[i]->getSeq2() << "</data>" << std::endl;

		} else if(dg1 <= minPlotDegree && dg2 > 1 /*&& vertices.size() < 500*/) {
			e_file  << "    <node id=\"" << edges[i]->getVertex1()->getId() << "\">\n";
			e_file  << "      <data key=\"d0\">silver</data>\n";
			e_file  << "    </node>" << std::endl;

			e_file  << "    <edge id=\"e"<< i << "\" source=\"" << edges[i]->getVertex2()->getId()
					<< "\" target=\"" << edges[i]->getVertex1()->getId() << "\">" << std::endl;

			if( edges[i]->getHS1() == edges[i]->getHS2() ) // cis/trans
				e_file  << "      <data key=\"d4\">blueviolet</data>\n";
			else
				e_file  << "      <data key=\"d4\">palegreen</data>\n";

			e_file  << "      <data key=\"d1\">"<< pwidth << "</data>\n";
			e_file  << "      <data key=\"d2\">"<< edges[i]->getSeq1() << "</data>\n";
			e_file  << "      <data key=\"d3\">"<< edges[i]->getSeq2() << "</data>" << std::endl;
		} else if(dg1 <= minPlotDegree && dg2 <= minPlotDegree) continue;

		e_file << "    </edge>" << std::endl;
	} // for-edges

	e_file << "  </graph>" << std::endl;
	e_file << "</graphml>" << std::endl;
	e_file.close();

//	// ----- PLOT IT -----
//	std::vector<std::string> lts;
//	lts.push_back("fdp"); lts.push_back("neato");
//	for(unsigned i=0; i<lts.size(); ++i) {
//		if(vertices.size() < 3500) {
//			plotGraph(ess.str(), "svg", lts[i]);
//			plotGraph(ess.str(), "png", lts[i]);
//		}
//		else {
//			std::stringstream exec_cmd;
//			std::string es_file(ess.str()), tmp;
//			std::string::size_type found1, found2;
//			found1 = es_file.rfind(".");
//			tmp = es_file.substr(0, found1);
//			tmp.append("_").append(lts[i]).append(".png");
//			exec_cmd << lts[i] << " -Tpng" << " " << ess.str() << " -o " << tmp << " -q" << std::endl;
//			tmp.clear();
//			tmp = es_file.substr(0, found1);
//			tmp.append("_").append(lts[i]).append(".svg");
//			exec_cmd << lts[i] << " -Tsvg" << " " << ess.str() << " -o " << tmp << " -q" << std::endl;
//			std::cerr << "Very big graph: dot file created, here is the "
//					<< "plotting command. see dot manual for tuning parameters.\n\n"
//					<< exec_cmd.str() << std::endl;
//		}
//	}

	return baseFolder;
}

/*
 * PRINT GRAPH WITH VIRUS INTERSECTIONS INFORMATION WITH GRAPHVIZ
 */
std::string plotGVizIntr(Graph *g, std::vector<Edge*>& edges, std::vector<Node>& vertices, std::string start,
		std::string dts, std::string mld, int sc_limit, int wrks) {
	std::stringstream ess;

	long maxIntrs = ( std::max_element(vertices.begin(), vertices.end(), [](const Node& n1, const Node& n2) {
		return n1.intrs < n2.intrs;
	}) )->intrs;

	std::vector<long> nds;
	nds.reserve(vertices.size());

	std::string  out_folder("plots"), baseFolder(tmn::timeStamp()), fold("_gVizCollisions");
	mkdir(out_folder.c_str(), 0777);
	out_folder.append("/").append(baseFolder);
	mkdir(out_folder.c_str(), 0777);

	start.append(fold);
	out_folder.append("/").append(start);
	mkdir(out_folder.c_str(), 0777);

	ess << out_folder << "/" << dts << "_L" << sc_limit << "_" << mld << "_" << tmn::timeStamp() << ".dot";
	std::ofstream e_file(ess.str().c_str(),
			std::ios_base::out | std::ios_base::trunc );

	e_file << "/*****\n" << ess.str() << "\n*****/\n\n" << std::endl;
	e_file << "graph G {" << std::endl;
	e_file << " graph [splines=true, overlap=";

	if(vertices.size() < 2500) e_file << "prism,";
	else e_file << "compress,";

	e_file << " outputorder=edgesfirst, sep=3]\n"
			<< " node [shape=ellipse, style=filled, fixedsize=true]\n" << std::endl;

	// node colors
	for (unsigned i = 0; i < vertices.size(); ++i) {
		if(vertices[i].intrs > 1) {
			e_file << " " << vertices[i].n_id << " [";
			if(vertices[i].intrs <= maxIntrs/4) {
				e_file << "width=0.12, height=0.12, penwidth=0.5, fillcolor=\"grey95\", label=\"\"";
			} else if(vertices[i].intrs < maxIntrs/2) {
				e_file << "width=" << 1.5*(double)vertices[i].intrs/100 << ", height=" << 1.5*(double)vertices[i].intrs/100
						<< ", penwidth=0.7, fillcolor=\"orange1\", label=\"" << vertices[i].label << "\"," << " fontsize=5, fontname=\"times\"";
			} else {
				e_file << "width=" << 1.6*(double)vertices[i].intrs/100 << ", height=" << 1.6*(double)vertices[i].intrs/100
						<< ", penwidth=1.0, fillcolor=\"red\", label=\"" << vertices[i].label << "\"," << " fontsize=8, fontname=\"times\"";
			}

			e_file << "];" << std::endl;
			nds.push_back(vertices[i].n_id);
		}
	}

	// edges edgesCol: chartreuse, cornflowerblue, yellow, orangered
	int nedg=0;
	for (unsigned i = 0; i < edges.size(); ++i) {
		long dg1 = edges[i]->getVertex1()->getId();
		long dg2 = edges[i]->getVertex2()->getId();

		if( std::binary_search(nds.begin(), nds.end(), dg1) && std::binary_search(nds.begin(), nds.end(), dg2) ) {
			e_file << " " << dg1 << " -- " << dg2;
			e_file << " [penwidth=0.6,";

			if (edges[i]->getVertex1()->getGraphRoot() == edges[i]->getVertex2()->getGraphRoot()) {
				if(edges[i]->getHS1() == edges[i]->getHS2()) { // intra-chromosome
					e_file << " color=\"blueviolet\"];\n";
				} else {
					e_file << " color=\"palegreen\"];\n";
				}
			}
			++nedg;
		}
	}

	if(nedg==0) {
		std::cout << "WARNING: Virus interactions too low in this genome region.\n"
				<< "Graph not plotted!\n"
				<< "Start from a different gene." << std::endl;
		e_file.clear();
		return baseFolder;
	}

	e_file << "\n}" << std::endl;
	e_file.close();

	// ----- PLOT IT -----
	std::vector<std::string> lts;
	lts.push_back("fdp"); lts.push_back("neato");
	for(unsigned i=0; i<lts.size(); ++i) {
		if(vertices.size() < 3500) {
			plotGraph(ess.str(), "png", lts[i]);
			plotGraph(ess.str(), "svg", lts[i]);
		}
		else {
			std::stringstream exec_cmd;
			std::string es_file(ess.str()), tmp;
			std::string::size_type found1, found2;
			found1 = es_file.rfind(".");
			tmp = es_file.substr(0, found1);
			tmp.append("_").append(lts[i]).append(".png");
			exec_cmd << lts[i] << " -Tpng" << " " << ess.str() << " -o " << tmp << " -q" << std::endl;
			tmp.clear();

			tmp = es_file.substr(0, found1);
			tmp.append("_").append(lts[i]).append(".svg");
			exec_cmd << lts[i] << " -Tsvg" << " " << ess.str() << " -o " << tmp << " -q" << std::endl;
			std::cerr << "Very big graph: dot file created, here is the "
					<< "plotting command. see dot manual for tuning parameters.\n\n"
					<< exec_cmd.str() << std::endl;
		}
	}

	return baseFolder;
}

// print an ordered sam file
std::string printOrderedSam(std::string filename, std::vector<SamData>& sm) {
	std::string::size_type found = filename.rfind("/");
	std::string fl1 = filename.substr(0, found);
	std::string fl2 = filename.substr(found+1, filename.length());

	std::stringstream vss;
	vss << "/ord_" << fl2;

	fl1.append(vss.str());

	std::ofstream v_file( fl1, std::ios_base::out | std::ios_base::trunc );

	for(size_t i=0; i<sm.size(); ++i)
		v_file << sm[i] << std::endl;

	v_file.close();

	return fl1;
}

// print a new genes file with ordered genes
std::string printOrderedGeneFile(std::string filename, std::vector<Gene>& sm) {
	std::ofstream v_file( filename, std::ios_base::out | std::ios_base::trunc );

	for(size_t i=0; i<sm.size(); ++i)
		v_file << sm[i] << std::endl;

	v_file.close();

	return filename;
}


} // namespace



#endif /* PRINTOUT_HPP_ */
