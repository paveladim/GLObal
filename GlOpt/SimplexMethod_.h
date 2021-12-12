#ifndef SIMPLEX_METHOD_H
#define SIMPLEX_METHOD_H

#include <vector>

using Vec = std::vector<double>;
using Basis = std::vector<int>;
using Matrix = std::vector<std::vector<double>>;

enum RESSIMP { resOK, resUnlim, systemIsInconsist };

struct SimplexRes {
    RESSIMP _res;
    Vec _x;
    Basis _xbi;
    double _resVal;

    SimplexRes(int sizeX, int sizeBasis) :
        _xbi(sizeBasis),
        _x(sizeX),
        _resVal(0.0),
        _res(RESSIMP::resUnlim) { }

    void change(int sizeX, int sizeBasis) {
        _xbi.resize(sizeBasis);
        _x.resize(sizeX);
    }
};

class SimplexMethod {
private: // fields
    int _n;
    int _m;
    Vec _c;
    Vec _c_copy;
    Vec _c_copy1;
    Vec _b;
    Matrix _A;
    Basis _basis;
    Vec _x;
    Basis _ImitBasis;
    int _countImit;
    int _err_OK;
    SimplexRes _result;
    double _eps;
private: // methods
    void excludeImitVariables(const int& s);
    void deleteElemFromImitBasis(const int& s);
    void reBasis(const int& s);
    void DeleteFromB(const int& r);
    void DeleteFromC(const int& s);
    void deleteRow(const int& r);
    void deleteColumn(const int& s);
    void deleteImitColumn();
    double CalcRes(const int& stage);
    int FindBasis();
    int AddImitBasis(int cntImit, Basis& tmp);
    int SM_Algoritm();
    int getLeadingColumn();
    int getLeadingRow(int s);
    int compareLexVectors(Vec& v1, Vec& v2);
    void gaussOperation(int s, int r);
public: //methods
    SimplexMethod(const int& nA, const int& mA, const Vec& c, const Vec& b, const Matrix& A);
    SimplexRes solve();
    ~SimplexMethod() {}
};

#endif // SIMPLEX_METHOD_H

