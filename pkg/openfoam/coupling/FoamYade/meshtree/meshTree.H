// Deepak Kunhappan (deepak.Kunhappan@3sr-grenoble.fr, deepak.kn1990@gmail.com)
// Part of YADE - OpenFOAM coupling, Build a  k-d tree from the mesh, and
// search for cell centers within a given distance range and return them.
// Inspired from the C version at https://rosettacode.org/wiki/K-d_tree#C
// this method is way faster than the FOAM's search (mesh.findCell(vector) )
// source file : meshTree.C

#ifndef meshTree_H
#define meshTree_H
#include <algorithm>
#include <vector>

// Check development branch, the foundation version reorganized classes after version 6
#if foamVersion < 1000 && foamVersion > 6
#define newFoundationVersion 1
#else
#define isFoundationVersion 0
#endif

#if (newFoundationVersion)
// The content of "fvCFD.H" has been splitted, pick what's necessary below
#include "volFields.H"
// #include "fvMesh.H"
// #include "fvSchemes.H"
// #include "fvSolution.H"
// #include "surfaceFields.H"
// #include "fvm.H"
#else
#include "fvCFD.H"
#endif

namespace Foam {

class meshpt {
public:
	meshpt(const double* x, const double* y, const double* z, int cid)
	{
		pt.push_back(x);
		pt.push_back(y);
		pt.push_back(z);
		id = cid;
	}
	std::vector<const double*> pt;
	int                        id;
	virtual ~meshpt() {};
};

class kdNode {
public:
	kdNode(meshpt point)
	        : p(point)
	{
		left  = NULL;
		right = NULL;
	}
	meshpt  p;
	kdNode* left;
	kdNode* right;
	virtual ~kdNode() {};
};

class cmpnode {
public:
	cmpnode() {};
	bool operator()(const std::pair<kdNode*, double>& n1, const std::pair<kdNode*, double>& n2) { return n1.second < n2.second; }
};

class cmpvec {
public:
	cmpvec(int axis)
	        : a(axis) {};

	bool operator()(const meshpt& p1, const meshpt& p2) { return *(p1.pt[a]) < *(p2.pt[a]); }
	int  a;
};


class pqueue {
public:
	pqueue(unsigned int bound)
	{
		maxbound = bound;
		container.reserve(bound + 1);
	}
	unsigned int                            maxbound;
	double                                  maxdist;
	std::vector<std::pair<kdNode*, double>> container;
	void                                    push_node(std::pair<kdNode*, double> aNode)
	{
		if (container.size() == maxbound) {
			if (container[maxbound].second > aNode.second && !incontainer(aNode)) {
				container.pop_back();
				container.push_back(aNode);
				std::sort(container.begin(), container.end(), cmpnode());
			}
		}
		if (!incontainer(aNode)) {
			container.push_back(aNode);
			std::sort(container.begin(), container.end(), cmpnode());
		}
	}

	bool incontainer(std::pair<kdNode*, double> aNode)
	{
		bool value = false;
		int  c     = 0;
		for (unsigned int i = 0; i != container.size(); ++i) {
			if (aNode.first->p.id == container[i].first->p.id) c += 1;
		}
		if (c > 0) value = true;
		return value;
	}

	virtual ~pqueue() {};
};


class meshTree {
public:
	meshTree(const volVectorField& mshC)
	        : meshm(mshC) {};
	const volVectorField&                  meshm;
	kdNode*                                root;
	void                                   build_tree();
	kdNode*                                recursive_build_tree(std::vector<meshpt>&, int);
	int                                    nearestCell(const vector&);
	std::vector<int>                       nnearestCellsRange(const vector& v, const double& range, const bool&);
	kdNode*                                nnearest(kdNode*, const meshpt&, kdNode*, double&, int, pqueue&);
	const int                              ndim = 3;
	typedef std::vector<meshpt>::size_type vec_sz;
	void                                   get_median(std::vector<meshpt>&, const int&);
	kdNode*                                recursive_nearest_cell(kdNode*, const meshpt&, kdNode*, double&, int);
	double                                 distance(const meshpt&, const meshpt&);
	int                                    numlevels;
	virtual ~meshTree() {};
	void traversTree();
	void _traversTree(kdNode*);
};
}

#endif
