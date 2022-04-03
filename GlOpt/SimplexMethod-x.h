#ifndef SIMPLEX_METHODX_H
#define SIMPLEX_METHODX_H

#include <vector>

using Vec    = std::vector<double>;
using Matrix = std::vector<std::vector<double>>;
using Basis  = std::vector<int>;

enum RESSIMP { resOK, resUnlim, systemIsInconsist };

class TSimplexRes
{
public:
	RESSIMP res;
	Vec x;
	Basis xbi;
	double resVal;

	TSimplexRes(const int& sizeX, const int& sizeBasis) {
		xbi.resize(sizeBasis);
		x.resize(sizeX);
		resVal = 0;
	}

	void resize(const int& sizeX, const int& sizeBasis) {
		xbi.resize(sizeBasis);
		x.resize(sizeX);
		resVal = 0;
	}
};

class SimplexMethod {
	int n; // variables
	int m; // constraints

	Vec c;
	Vec c_copy;
	Vec c_copy1;

	Vec b;

	Matrix A;

	Basis basis;
	Basis imit_basis;
	Vec x;

	int count_imit;
	int err_OK;

	double eps;

	TSimplexRes result;
public:
	SimplexMethod(const Vec&, const Vec&, const Matrix&);
	TSimplexRes FindSolution();
	void excludeImitationVariables(const int& s);
	void deleteElemFromImitBasis(const int& s);
	void reBasis(const int& s);
	void DeleteFromB(const int& r);
	void DeleteFromC(const int& s);
	void deleteRow(const int& r);
	void deleteColumn(const int& s);
	void deleteImitColumn();
	double CalcRes(const int& stage);
	int FindBasis();
	int AddImitBasis(const int& cntImit, const Basis& tmp);
	int SM_Algorithm();
	int getLeadingColumn();
	int getLeadingRow(const int& s);
	void gaussOperation(const int& s, const int& r);
	int compareLexVectors(const Vec& v1, const Vec& v2);
	~SimplexMethod() {};
};

#endif