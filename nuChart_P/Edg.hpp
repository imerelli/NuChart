/*
 * Edg.hpp
 *
 *  Created on: Oct 9, 2014
 *      Author: fabio
 */

#ifndef EDG_HPP_
#define EDG_HPP_

#include "common.hpp"
#include "Vertex.hpp"

class Edg {
public:
	Edg() : v1(NULL), v2(NULL), weight(0.0), prob(0.0) { }
	Edg(Vertex *_v1, Vertex *_v2, double w, double p) : v1(_v1), v2(_v2), weight(w), prob(p) { }
	Edg(const Edg& e) : v1(new Vertex(*e.v1)), v2(new Vertex(*e.v2)), weight(e.getWeight()), prob(e.getProb()) { }
	Edg(Edg&& e) : v1(NULL), v2(NULL), weight(e.getWeight()), prob(e.getProb()) {
		v1 = e.getVertex1();
		v2 = e.getVertex2();

//		e.setVertex1(NULL);
//		e.setVertex2(NULL);
	}

	Edg& operator=(const Edg& e) { // shallow-copy is fine!
		if(this != &e) {
			if(v1) delete v1;
			if(v2) delete v2;

			v1 = new Vertex(*e.v1);
			v2 = new Vertex(*e.v2);
			weight = e.getWeight();
			prob = e.getProb();
		}
		return *this;
	}
	Edg& operator=(Edg&& e) {
		if(this != &e) {
			if(v1) delete v1;
			if(v2) delete v2;

			v1 = e.getVertex1();
			v2 = e.getVertex2();
			weight = e.getWeight();
			prob = e.getProb();

//			e.setVertex1(NULL);
//			e.setVertex2(NULL);
		}
		return *this;
	}

	virtual ~Edg() { };

	inline Vertex* getVertex1() const { return v1; }
	inline Vertex* getVertex2() const { return v2; }
	inline double getWeight()   const { return weight; }
	inline double getProb()   	const { return prob; }

	void setVertex1(Vertex *v)  { v1 = v; }
	void setVertex2(Vertex *v)  { v2 = v; }
	void setWeight(double w)    { weight = w; }
	void setProb(double p)		{ prob = p; }

protected:
	Vertex *v1;
	Vertex *v2;
	double weight, prob;

};


#endif /* EDG_HPP_ */
