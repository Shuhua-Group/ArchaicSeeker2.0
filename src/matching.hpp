/*
 * matching.hpp
 *
 *  Created on: Oct 10, 2018
 *      Author: yuankai
 */

#ifndef MATCHING_HPP_
#define MATCHING_HPP_

# include <iostream>
# include <map>
# include <vector>
# include <set>
# include <math.h>

# include <nlopt.h>

class matchModel;

class newickFunc;

class newickDirFunc;

class newickNode
{
public:
	newickNode() : label(""), dis(0), nchild(0), child(NULL), parent(NULL) {};
	newickNode(const std::string& tree);
	~newickNode();

	std::string getLab() const { return label; }
	double getDis() const { return dis; }

	void traversal(newickFunc& func);
	void search(newickDirFunc& func);

	std::string _load(const std::string &subTree, newickNode* upper);

	std::string label;
	double dis;
	int nchild;
	newickNode *child, *parent;
};

class newickFunc
{
public:
	virtual ~newickFunc() {};
	virtual void operator() (newickNode& node) = 0;
};

class newickDirFunc
{
public:
	newickDirFunc() : preNode(NULL), upper(true), isFirst(true), upSearchStat(-1) {};
	virtual ~newickDirFunc() {};
	virtual void operator() (newickNode& node) = 0;

	newickNode* preNode;
	bool upper;

	bool isFirst;
	int upSearchStat; // 1 up search, -1 down search, -1 cur
	std::set<newickNode*> check;
};

class labcheck : public newickFunc
{
public:
	void operator() (newickNode& node);
	std::map<std::string, newickNode*> access;
};

class getNodeLab : public newickFunc
{
public:
	getNodeLab(std::vector<std::string>& _labels, std::vector<std::string>& _leafLab, \
			std::vector<std::string>& _midLab) : labels(_labels), leafLab(_leafLab), midLab(_midLab) {}
	void operator() (newickNode& node);
	std::vector<std::string>& labels, & leafLab, & midLab;
};

#ifdef DEBUG

class labprint : public newickFunc
{
public:
	void operator() (newickNode& node);
};

#endif//DEBUG

class disSearch : public newickDirFunc
{
public:
	disSearch() : newickDirFunc() {};
	void operator() (newickNode& node);
	std::map<newickNode*, double> disMem;
	std::map<std::string, double> labDis;
	std::map<std::string, int> isUpSearch; // 1 up search, -1 down search, -1 cur
	std::map<std::string, std::set<std::string> > nodeCheck;
};

class matFuncData
{
public:
	matFuncData() : nleaf(0), nlab(0), dc(NULL), eta(NULL), diff(NULL), mu(0) {}
	~matFuncData();
	void init(int _nlaef, int _nlab);

	int nleaf, nlab;
	double ***dc; //[leaf]testRecNode [lab]testIntorNode [leaf]
	double **eta; //[lab] [leaf]
	double *diff; //[leaf]
	double mu;
};

class curMatFuncData
{
public:
	curMatFuncData(const matFuncData& data, double length, int testRec, int testIntro) :
		nleaf(data.nleaf), inter(0), curDc(data.dc[testRec][testIntro]),
		curEta(data.eta[testIntro]), curDiff(data.diff), mu(data.mu), len(length) {}
	int nleaf, inter;
	double *curDc, *curEta, *curDiff;
	double mu, len;
};

double matLKFunc(unsigned int n, const double *x, double *grad, void *matData);

double globalMatFunc(unsigned int n, const double *x, double *grad, void *matData);

class modelCalibrateData
{
public:
	modelCalibrateData(const std::vector<std::vector<double> >& _pairDiff, \
			const std::vector<std::vector<std::vector<int> > >& _pairLab, \
			const std::vector<std::vector<int> >& _rootLab, \
			double _mu, long _length) : \
			pairDiff(_pairDiff), pairLab(_pairLab), rootLab(_rootLab), \
			nlab(_pairLab[1].front().size()), \
			nleaf(_pairLab.size()), mu(_mu), length(_length){};
	const std::vector<std::vector<double> > &pairDiff;
	const std::vector<std::vector<std::vector<int> > > &pairLab;
	const std::vector<std::vector<int> > &rootLab;
	int nlab, nleaf;
	double mu, length;
};

double modelCalibrateFunc(unsigned int n, const double *x, double *grad, void *matData);

class modelCalibrateConstraintData
{
public:
	modelCalibrateConstraintData(const std::vector<int>& _lab1, \
			const std::vector<int>& _lab2, double _diff) :\
			lab1(_lab1), lab2(_lab2), diff(_diff) {};
	const std::vector<int> &lab1, &lab2;
	double diff;
};

double modelCalibrateConstraintFunc(unsigned int n, const double *x, double *grad, \
		void *cData);

class matchModel
{
public:
	matchModel(const std::string& m);

	~matchModel();

	friend class newickNode;

	//return the best matched lab
	int segMatchLK(double* ndiff, double length, int recp, double* mat, double* llk) const;

	void modelCalibrate(const std::vector<std::vector<double> >& pairDiff, long len);

#ifdef DEBUG

	void print() { labprint p; model.traversal(p); }

#endif//DEBUG

	newickNode model;
	labcheck access;

	std::map<std::string, disSearch> search;
	std::vector<std::string> labels, leafLabs, midLabs;
	std::vector<double> nodeLen, rootLen, nodeTestUpbound;
	std::vector<std::vector<int> > rootLab; //[nleaf]
	std::vector<std::vector<double> > nodeDis; //[labels] [leafs]
	std::vector<std::vector<double> > pairTreeDis; //[leafs] [leafs]
	std::vector<std::vector<std::vector<int> > > pairLab; //[leafs] [leafs] [labels]
	std::vector<double> calibratedLen;
/*
 * m node and n leafs
 * t stands for the time before introgress node's time
 * For a specific introgress node, the time to the test node is d
 * alpha = 1 / ( L * mu )
 * if is upper search eta = 0, or eta = 2
 * diff = ( d + dis_intro2leaf + eta * t ) * alpha
 * x = t * alpha; y = alpha;
 * f(x,y) = sum(1->n){( (d + c) * y + eta * x - diff ) ^ 2}
 * df/dx = sum(1->n){ 2 * eta * ( (d + c) * y + eta * x - diff ) } = 0
 * df/dy = sum(1->n){ 2 * (d + c) * ( (d + c) * y + eta * x - diff) } = 0
 * A = sum(1->n){ eta * (d + c) } 	[n][m] TestRecieveNode TestIntroNode
 * B = sum(1->n){ eta ^ 2 }			[m] TestIntroNode
 * C = sum(1->n){ eta * diff }		eta[m][n] TestIntroNode PairDifference
 * D = sum(1->n){ (d + c) ^ 2 }		[n][m] TestRecieveNode TestIntroNode
 * E = sum(1->n){ (d + c) * diff }	(d + c)[n][m][n] TestRecieveNode TestIntroNode PairDifference
 * F = sum(1->n){ diff ^ 2 }
 * x = (C*D - A*E) / (B*D - A^2)
 * y = (B*E - A*C) / (B*D - A^2)
 * t = x / y = (C*D - A*E) / (B*E - A*C)
 * sqrt = D * y^2 + B * x^2 + F + 2 * A * xy - 2 * E * y - 2 * C *x
 */

/*
 * L(x,y) = sum(1->n){ diff * log( (d + c + eta * t ) * alpha ) - (d + c + eta * t ) * alpha
 * 		- log(diff!) }
 * x[0] = t, x[1] = alpha
 * dL/dt = sum(1->n){ diff * eta / (d + c + eta * t) - alpha * eta}
 * dL/d(alpha) = sum(1->n){ diff / alpha - (d + c + eta * t)}
 */
	matFuncData matData;
	double** initVal, ** lb, ** ub;//[nlab] [2];

	double mu;
};

#endif /* MATCHING_HPP_ */
