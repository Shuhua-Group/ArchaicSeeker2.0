/*
 * matching.cpp
 *
 *  Created on: Oct 10, 2018
 *      Author: yuankai
 */

# include "matching.hpp"

# include <cassert>
# include <cstdlib>
# include <vector>
# include <math.h>

# include "boost/algorithm/string.hpp"

const int MAXEVAL = 10000;

newickNode::newickNode(const std::string& tree) : label(""), dis(0), \
		nchild(0), child(NULL), parent(NULL)
{
	std::string str = boost::algorithm::trim_copy(tree);
	if(str[str.size() - 1] == ';')
		str = str.substr(0, str.size() - 1);
	str += ":0";
	label = _load(str, NULL);
}

newickNode::~newickNode()
{
	if(nchild)
		delete[] child;
}

//Go through all the node from up to bottom
void newickNode::traversal(newickFunc& func)
{
	func(*this);
	for(int i = 0 ; i < nchild; ++i)
		child[i].traversal(func);
}

void newickNode::search(newickDirFunc& func)
{
	if(func.check.count(this))
		return;
	func.check.insert(this);
	func(*this);
	bool curFirst(func.isFirst);
	if(func.isFirst)
		func.isFirst = false;
	if(parent)
	{
		func.preNode = this;
		func.upper = true;
		if(curFirst)
			func.upSearchStat = 1;
		parent->search(func);
	}
	if(nchild)
		for(int i = 0 ; i < nchild; ++i)
		{
			func.preNode = this;
			func.upper = false;
			if(curFirst)
				func.upSearchStat = -1;
			child[i].search(func);
		}
}

void labcheck::operator ()(newickNode& node)
{
	assert(access.count(node.getLab()) == 0);
	access[node.getLab()] = &node;
}

void getNodeLab::operator ()(newickNode& node)
{
	if(node.nchild == 0)
		leafLab.push_back(node.getLab());
	else if(node.parent != NULL)
		midLab.push_back(node.getLab());
	labels.push_back(node.getLab());
}

#ifdef DEBUG

void labprint::operator ()(newickNode& node)
{
	std::cout << node.getLab() << '\t' << node.getDis() << std::endl;
	if(node.parent != NULL)
		std::cout << "\tParent:\t" << node.parent->getLab() << std::endl;
	if(node.nchild != 0)
	{
		std::cout << "\tChild:";
		for(int i = 0 ; i < node.nchild; ++i)
			std::cout << '\t' << node.child[i].getLab() ;
		std::cout << std::endl;
	}
}

#endif//DEBUG

void disSearch::operator ()(newickNode& node)
{
	if(preNode)
	{
		if(upper)
		{
			disMem[&node] = preNode->getDis() + disMem[preNode];
			nodeCheck[node.getLab()] = nodeCheck[preNode->getLab()];
			nodeCheck[node.getLab()].insert(preNode->getLab());
		}
		else
		{
			disMem[&node] = node.getDis() + disMem[preNode];
			nodeCheck[node.getLab()] = nodeCheck[preNode->getLab()];
			nodeCheck[node.getLab()].insert(node.getLab());
		}
	}
	else
		disMem[&node] = 0;
	labDis[node.getLab()] = disMem[&node];
	isUpSearch[node.getLab()] = upSearchStat;
}

std::string newickNode::_load(const std::string& subTree, newickNode* upper)
{
	parent = upper;
	int p(subTree.size() - 1);
	while(p > 0 && subTree[p] != ':')
		--p;
	assert(p != 0);
	assert(p != (int)subTree.size() - 1);

	dis = atof(subTree.substr(p + 1).c_str());

	if(subTree[0] == '(')
	{
		assert(p >= 2);
		assert(subTree[p - 1] == ')');
		std::vector<std::string> sp;
		std::string t = subTree.substr(1, p - 2);
		int bracket(0);
		uint pre(0);
		for(uint i = 0 ; i < t.size(); ++i)
		{
			if(t[i] == '(')
				++bracket;
			else if(t[i] == ')')
				--bracket;
			else if(t[i] == ',' && bracket == 0)
			{
				sp.push_back(t.substr(pre, i - pre));
				pre = i + 1;
			}
		}
		sp.push_back(t.substr(pre));
		nchild = sp.size();
		child = new newickNode[nchild];
		label = child[0]._load(sp[0], this);
		for(int i = 1 ; i < nchild ; ++i)
			label += "_" + child[i]._load(sp[i], this);
	}
	else
		label = subTree.substr(0, p);
	return label;
}

matFuncData::~matFuncData()
{
	if(nleaf)
	{
		for(int i = 0 ; i < nleaf; ++i)
		{
			for(int j = 0 ; j < nlab; ++j)
				delete[] dc[i][j];
			delete[] dc[i];
		}
		for(int i = 0 ; i < nlab ; ++i)
			delete[] eta[i];
		delete[] dc;
		delete[] eta;
		delete[] diff;
	}
}

void matFuncData::init(int _nleaf, int _nlab)
{
	if(nleaf || nlab)
	{
		for(int i = 0 ; i < nleaf; ++i)
		{
			for(int j = 0 ; j < nlab; ++j)
				delete[] dc[i][j];
			delete[] dc[i];
		}
		for(int i = 0 ; i < nlab ; ++i)
			delete[] eta[i];
		delete[] dc;
		delete[] eta;
		delete[] diff;
	}
	assert(_nleaf);
	assert(_nlab);
	nleaf = _nleaf;
	nlab = _nlab;
	dc = new double** [nleaf];
	for(int i = 0 ; i < nleaf; ++i)
	{
		dc[i] = new double* [nlab];
		for(int j = 0 ; j < nlab; ++j)
		{
			dc[i][j] = new double [nleaf];
			memset(dc[i][j], 0, sizeof(double) * nleaf);
		}
	}
	eta = new double* [nlab];
	for(int i = 0 ; i < nlab; ++i)
	{
		eta[i] = new double[nleaf];
		memset(eta[i], 0, sizeof(double) * nleaf);
	}
	diff = new double [nleaf];
	memset(diff, 0, sizeof(double) * nleaf);
}

double matLKFunc(unsigned int n, const double *x, double *grad, void *matData)
{
	curMatFuncData *data = (curMatFuncData *) matData;
	double denon[data->nleaf];
//	std::cout << "mu: " << data->mu << std::endl;
	double mu(data->mu);
	double len(data->len);
	for(int i = 0 ; i < data->nleaf; ++i)
		denon[i] = data->curDc[i] + data->curEta[i] * x[0];
	if(grad)
	{
		grad[0] = 0;
		for(int i = 0 ; i < data->nleaf; ++i)
		{
			grad[0] += data->curDiff[i] / len / mu * data->curEta[i] / denon[i] \
					- data->curEta[i];
		}
	}
	double ret(0);
	for(int i = 0 ; i < data->nleaf; ++i)
	{
		ret += data->curDiff[i] / len / mu * log(denon[i]) - denon[i];
	}
/*	++data->inter;

	std::cout << "iter: " << data->inter << std::endl;
	for(int i = 0 ; i < data->nleaf; ++i)
		std::cout << '\t' << data->curDiff[i];
	std::cout << std::endl;
	for(int i = 0 ; i < data->nleaf; ++i)
		std::cout << '\t' << data->curEta[i];
	std::cout << std::endl;
	std::cout << grad[0] << '\t' << grad[1] << std::endl;
	std::cout <<  x[0] << '\t' << x[1] << std::endl;
*/
	return ret;
}

double modelCalibrateFunc(unsigned int n, const double *x, double *grad, void *matData)
{
	modelCalibrateData *data = (modelCalibrateData *) matData;
	int nlab(data->nlab), nleaf(data->nleaf);
	double mu(data->mu), length(data->length);
	double ret(0);
	if(grad)
	{
		memset(grad, 0, sizeof(double) * ( nlab ));
		for(int i = 0 ; i < nleaf ; ++i)
		{
			for(int j = 0 ; j < i ; ++j)
			{
				const double& diff(data->pairDiff[i][j] / length / mu);
				double t(0);
				const std::vector<int>& curLab(data->pairLab[i][j]);
				for(int k = 0 ; k < nlab; ++k)
				{
					t += curLab[k] * x[k];
				}
				for(int k = 0 ; k < nlab; ++k)
					grad[k] += curLab[k] * ( diff / t - 1.0 );
				ret += diff * log( t ) - t;
			}
			double rt(0);
			const std::vector<int>& curRootLab(data->rootLab[i]);
			for(int k = 0 ; k < nlab; ++k)
				rt += curRootLab[k] * x[k];
			const double& rdiff(data->pairDiff[i][i] / length / mu);
			for(int k = 0 ; k < nlab; ++k)
				grad[k] += curRootLab[k] * ( rdiff / rt - 1.0 );
			ret += rdiff * log(rt) - rt;
		}
	}
	else
	{
		for(int i = 0 ; i < nleaf ; ++i)
		{
			for(int j = 0 ; j < i ; ++j)
			{
				const double& diff(data->pairDiff[i][j] / length / mu);
				double t(0);
				const std::vector<int>& curLab(data->pairLab[i][j]);
				for(int k = 0 ; k < nlab; ++k)
					t += curLab[k] * x[k];
				ret += diff * log( t ) - t;
			}
			double rt(0);
			const std::vector<int>& curRootLab(data->rootLab[i]);
			for(int k = 0 ; k < nlab; ++k)
				rt += curRootLab[k] * x[k];
			const double& rdiff(data->pairDiff[i][i] / length / mu);
			ret += rdiff * log(rt) - rt;
		}
	}
	return ret;
}

double modelCalibrateConstraintFunc(unsigned int n, const double *x, double *grad, \
		void *cData)
{
	modelCalibrateConstraintData *data = (modelCalibrateConstraintData *) cData;
	int nlab(data->lab1.size());
	if(grad)
	{
		for(int i = 0 ; i < nlab; ++i)
			grad[i] = data->lab1[i] - data->lab2[i];
		grad[nlab] = 0;
	}
	double ret(0);
	for(int i = 0 ; i < nlab ; ++i)
	{
		ret += data->lab1[i] * x[i] - data->lab2[i] * x[i];
	}
	ret -= data->diff;
	return ret;
}

matchModel::matchModel(const std::string& m) : model(m)
{
	model.traversal(access);
	getNodeLab gleaf(labels, leafLabs, midLabs);
	model.traversal(gleaf);
	int nlab = labels.size();
	int nleaf = leafLabs.size();
	nodeTestUpbound.resize(nlab);
	for(int i = 0 ; i < nlab; ++i)
	{
		access.access.at(labels[i])->search(search[labels[i]]);
		nodeTestUpbound[i] = access.access.at(labels[i])->dis;
	}
	assert(nlab > 2);
	nodeLen.resize(nleaf);
	rootLen.resize(nlab);
	const disSearch& rootSearch = search[model.getLab()];
	double maxLen(0);
	for(int i = 0; i < nleaf; ++i)
		nodeLen[i] = rootSearch.labDis.at(leafLabs[i]);
	for(int i = 0; i < nlab ; ++i)
	{
		rootLen[i] = rootSearch.labDis.at(labels[i]);
		if(maxLen < rootLen[i])
			maxLen = rootLen[i];
	}
	nodeTestUpbound[0] = maxLen;
	nodeDis.resize(nlab);
	for(int i = 0; i < nlab; ++i)
	{
		std::vector<double>& curDis = nodeDis[i];
		const disSearch& curSearch = search.at(labels[i]);
		curDis.resize(nleaf);
		for(int j = 0; j < nleaf; ++j)
			curDis[j] = curSearch.labDis.at(leafLabs[j]);
	}
	matData.init(nleaf, nlab);
	initVal = new double* [nlab];
	lb = new double* [nlab];
	ub = new double* [nlab];
	for(int i = 0 ; i < nlab; ++i)
	{
		lb[i] = new double[1];
		lb[i][0] = 1e-128;
		ub[i] = new double[1];
		ub[i][0] = nodeTestUpbound[i];
		initVal[i] = new double[1];
		initVal[i][0] = nodeTestUpbound[i] / 2;
	}
	for(int i = 0 ; i < nlab; ++i)
	{
		const disSearch& curSearch = search[labels[i]];
		for(int j = 0; j < nleaf; ++j)
			if(curSearch.isUpSearch.at(leafLabs[j]) != 1)
				matData.eta[i][j] = 2;
	}
	for(int i = 0 ; i < nleaf; ++i)
		for(int j = 0 ; j < nlab; ++j)
		{
			double d(nodeLen[i] - rootLen[j]);
			for(int k = 0 ; k < nleaf; ++k)
				matData.dc[i][j][k] = d + nodeDis[j][k];
		}
	pairTreeDis.resize(nleaf);
	pairLab.resize(nleaf);
	for(int i = 0; i < nlab; ++i)
	{
		std::vector<double>& curDis = nodeDis[i];
		const disSearch& curSearch = search.at(labels[i]);
		curDis.resize(nleaf);
		for(int j = 0; j < nleaf; ++j)
			curDis[j] = curSearch.labDis.at(leafLabs[j]);
	}
	for(int i = 0; i < nleaf; ++i)
	{
		pairTreeDis[i].resize(nleaf, 0);
		pairLab[i].resize(nleaf);
		const disSearch& curSearch = search.at(leafLabs[i]);
		for(int j = 0 ; j < i; ++j)
		{
			pairTreeDis[i][j] = curSearch.labDis.at(leafLabs[j]);
			pairLab[i][j].resize(nlab, 0);
			for(int k = 0 ; k < nlab ; ++k)
			{
				if(curSearch.nodeCheck.at(leafLabs[j]).count(labels[k]))
					pairLab[i][j][k] = 1;
			}
		}
	}
	rootLab.resize(nleaf);
	for(int i = 0 ; i < nleaf; ++i)
	{
		rootLab[i].resize(nlab, 0);
		rootLab[i][0] = 1;
		for(int j = 1 ; j < nlab ; ++j)
		{
			if(rootSearch.nodeCheck.at(leafLabs[i]).count(labels[j]))
				rootLab[i][j] = 1;
		}
	}
	mu = -1;
}

matchModel::~matchModel()
{
	int nlab = labels.size();
	for(int i = 0 ; i < nlab; ++i)
	{
		delete[] initVal[i];
		delete[] lb[i];
		delete[] ub[i];
	}
	delete[] initVal;
	delete[] lb;
	delete[] ub;
}

int matchModel::segMatchLK(double* ndiff, double length, int recep, double* mat, double* llk) const
{
	int nleaf (leafLabs.size()), nlab(labels.size());
	for(int i = 0 ; i < nleaf; ++i)
		matData.diff[i] = ndiff[i];
	double x;
	double tol(1e-9);
	for(int i = 0 ; i < nlab; ++i)
	{
		nlopt_opt opt;
		opt = nlopt_create(NLOPT_LD_MMA, 1);
		nlopt_set_lower_bounds(opt, lb[i]);
		nlopt_set_upper_bounds(opt, ub[i]);
		curMatFuncData d(matData, length, recep, i);
		nlopt_set_max_objective(opt, matLKFunc, &d);
		x = initVal[i][0];
		nlopt_set_xtol_abs(opt, &tol);
		nlopt_set_maxeval(opt, MAXEVAL);
		double tllk;
		if(nlopt_optimize(opt, &x, &tllk) < 0)
			llk[i] = -HUGE_VAL;
		else
		{
			llk[i] = tllk;
			mat[i] = x;
		}
		nlopt_destroy(opt);
	}
	double maxLK(llk[0]);
	int maxP(0);
	for(int i = 2; i < nlab; ++i)
		if(maxLK < llk[i])
		{
			maxP = i;
			maxLK = llk[i];
		}
	return maxP;
}

void matchModel::modelCalibrate(const std::vector<std::vector<double> >& pairDiff, long len)
{
	int nleaf(leafLabs.size());
	int nlab(labels.size());
	double sum(0), sumLen(0);
	for(int i = 0 ; i < nleaf; ++i)
		for(int j = 0 ; j < i ; ++j)
		{
			sum += pairDiff[i][j] / len;
			sumLen += pairTreeDis[i][j];
		}
	mu = sum / sumLen;
	modelCalibrateData cData(pairDiff, pairLab, rootLab, mu, len);
	double calLen[nlab];
	nlopt_opt opt;
	opt = nlopt_create(NLOPT_LD_MMA, nlab);
	double clb[nlab], cub[nlab], ctol[nlab];
	for(int i = 0 ; i < nlab ; ++i)
	{
		calLen[i] = nodeTestUpbound[i];
		clb[i] = 1e-128;
		cub[i] = HUGE_VAL;
		ctol[i] = 1e-9;
	}
	double llk(0);
	nlopt_set_lower_bounds(opt, clb);
	nlopt_set_upper_bounds(opt, cub);
	nlopt_set_max_objective(opt, modelCalibrateFunc, &cData);
/*	std::vector<modelCalibrateConstraintData> constraintData;
	for(int i = 0 ; i < nleaf ; ++i)
		for(int j = 0 ; j < i ; ++j)
		{
			constraintData.push_back(modelCalibrateConstraintData(rootLab[i], rootLab[j], \
					nodeLen[i] - nodeLen[j]) );
			std::cout << leafLabs[i] << '\t' << leafLabs[j] << '\t' << \
					nodeLen[i] - nodeLen[j] << std::endl;
			nlopt_add_equality_constraint(opt, modelCalibrateConstraintFunc, \
					&constraintData, 1e-9);
		}
		*/
	nlopt_set_xtol_abs(opt, ctol);
	nlopt_set_maxeval(opt, MAXEVAL);
	if(nlopt_optimize(opt, calLen, &llk) < 0)
	{
		std::cout << "Cannot find the opt" << std::endl;
		llk = -HUGE_VAL;
	}

	nlopt_destroy(opt);
	calibratedLen.resize(nlab, 0);
	for(int i = 0 ; i < nlab ; ++i)
		calibratedLen[i] = calLen[i];

	for(int i = 0 ; i < nlab ; ++i)
	{
		access.access.at(labels[i])->dis = calibratedLen[i];
		nodeTestUpbound[i] = calibratedLen[i];
	}
	search.clear();
	for(int i = 0 ; i < nlab ; ++i)
		access.access.at(labels[i])->search(search[labels[i]]);

	const disSearch& rootSearch = search[model.getLab()];
	double maxLen(0);
	for(int i = 0; i < nleaf; ++i)
		nodeLen[i] = rootSearch.labDis.at(leafLabs[i]);
	for(int i = 0; i < nlab ; ++i)
	{
		rootLen[i] = rootSearch.labDis.at(labels[i]);
		if(maxLen < rootLen[i])
			maxLen = rootLen[i];
	}
	nodeTestUpbound[0] = maxLen;
	for(int i = 0; i < nlab; ++i)
	{
		ub[i][0] = nodeTestUpbound[i];
		initVal[i][0] = nodeTestUpbound[i] / 2;
	}
	for(int i = 0; i < nlab; ++i)
	{
		std::vector<double>& curDis = nodeDis[i];
		const disSearch& curSearch = search.at(labels[i]);
		for(int j = 0; j < nleaf; ++j)
			curDis[j] = curSearch.labDis.at(leafLabs[j]);
	}
	for(int i = 0; i < nleaf; ++i)
	{
		const disSearch& curSearch = search.at(leafLabs[i]);
		for(int j = 0 ; j < i; ++j)
		{
			pairTreeDis[i][j] = curSearch.labDis.at(leafLabs[j]);
		}
	}

	for(int i = 0 ; i < nleaf; ++i)
		for(int j = 0 ; j < nlab; ++j)
		{
			double d(nodeLen[i] - rootLen[j]);
			//std::cout << d << std::endl;
			for(int k = 0 ; k < nleaf; ++k)
				matData.dc[i][j][k] = d + nodeDis[j][k];
		}
	matData.mu = mu;
//	std::cout << matData.mu << std::endl;
//	std::cout << mu << std::endl;
}






