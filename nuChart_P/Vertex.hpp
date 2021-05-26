/*
 * Vertex.hpp
 *
 *  Created on: Aug 8, 2014
 *      Author: fabio
 */

#include "common.hpp"

#ifndef VERTEX_HPP_
#define VERTEX_HPP_

// abstract class for a graph node.
class Vertex {
public:
	Vertex() : id(-1), level(-1), root(NULL) {	}
	Vertex(const Vertex &v) : id(v.getId()), level(v.getLevel()), root(v.root) {
	}
	Vertex(Vertex&& v) : id(v.getId()), level(v.getLevel()), root(NULL) {
		root = v.getGraphRoot();
		v.setGraphRoot(NULL);
	}
	Vertex& operator=(const Vertex& v) {
		if(this != &v) {
//			if(root) delete root;
			id = v.getId();
			level = v.getLevel();
			root = v.root;
		}
		return *this;
	}
	Vertex& operator=(Vertex&& v) {
		if(this != &v) {
//			if(root) delete root;
			id = v.getId();
			level = v.getLevel();
			root = v.getGraphRoot();
			v.setGraphRoot(NULL);
		}
		return *this;
	}

	virtual ~Vertex() {
		//if(root) delete root;
	};

	inline int getLevel() 			const { return level; }
	inline long getId()				const { return id; }
	inline Vertex* getGraphRoot()	const { return root; }

	inline void setLevel(int l) 			{ assert(l >= 0 ); level = l; }
	inline void setId(long pos)				{ id = pos; }
	inline void setGraphRoot(Vertex *r)		{ root = r; }


protected:
	long id;			// node ID - its position in nodes list
	int level;			// distance from the root (0 is root)
	Vertex *root;		// pointer to graph root
};



#endif /* VERTEX_HPP_ */
