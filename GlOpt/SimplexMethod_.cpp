#include "SimplexMethod_.h"

SimplexMethod::SimplexMethod(const int& nA, const int& mA, const Vec& c, const Vec& b, const Matrix& A) :
	_n(nA),
	_m(mA),
	_c(c),
	_b(b),
	_result(nA, mA),
	_A(mA),
	_countImit(0),
	_eps(1e-6),
    _err_OK(0) {
	for (auto& elem : _A) elem.resize(nA);
}

SimplexRes SimplexMethod::solve() {
    _c_copy.resize(_n);
    for (int i = 0; i < _n; i++)
        _c_copy[i] = _c[i]; // копируем вектор с

    double calc_res;
    int res_OK = -1;

    //Ищем базис
    if (FindBasis() == 0) // если базис найден, то применяем алгоритм 1 (стр.35)
        if (_countImit != 0) // есть искусственные переменные, значит будем решать в 2 этапа
        // ЭТАП 1
        {
            int tmp_1 = 0;
            int tmp_2 = 0;
            for (int k = 0; k < _countImit; k++)
            {
                tmp_1 = _ImitBasis[k];
                for (int i = 0; i < _m; i++)
                    if (_A[i][tmp_1] == 1) tmp_2 = i;
                gaussOperation(tmp_1, tmp_2);
            }

            do
            {
                res_OK = SM_Algoritm(); // применяем алгоритм 1, пока не найдём решение 
                if (res_OK == -1) { 
                    _result._res = RESSIMP::resUnlim; 
                    return _result; 
                }
            } while (res_OK != 0);

            calc_res = CalcRes(1);

            if (calc_res > 0)
            {
                _result._res = RESSIMP::systemIsInconsist;
                return _result;
            }

            if (calc_res < 0) { 
                _result._res = RESSIMP::systemIsInconsist; 
                return _result; 
            }

            if (calc_res == 0)
            {
                bool flag = false;
                for (int ii = 0; ii < _countImit; ii++)
                {
                    if (std::find(_basis.begin(), _basis.end(), _ImitBasis[ii]) != _basis.end())
                    {
                        flag = true;
                        excludeImitVariables(_ImitBasis[ii]);
                        ii--;
                    }
                    flag = false;
                }

                if (!flag)
                    deleteImitColumn();
            }
        }
    // ЭТАП 2
    _c.resize(_n);
    for (int i = 0; i < _n; i++)
        _c[i] = _c_copy[i];

    int tmp1 = 0;
    int tmp2 = 0;

    for (int k = 0; k < _basis.size(); k++)
    {
        tmp1 = _basis[k];
        for (int i = 0; i < _m; i++)
            if (_A[i][tmp1] == 1) tmp2 = i;
        gaussOperation(tmp1, tmp2);
    }
    do
    {
        res_OK = SM_Algoritm(); // применяем алгоритм 1, пока не найдём решение 
        calc_res = CalcRes(2);
    } while (res_OK != 0);

    return _result;
}

void SimplexMethod::excludeImitVariables(const int& s) {
    int r = 0;
    bool isNull = true;
    // находим такое r, для которого q(r,s) = 1
    for (int i = 0; i < _m; i++)
        if (_A[i][s] == 1) r = i;

    for (int j = 0; j < _n - _countImit; j++)
    {
        if (_A[r][j] != 0)
        {
            isNull = false;
            gaussOperation(j, r);//б) выполняем шаг гаусова преобразования
            _basis[r] = j;// устанавливаем новую базисную переменной взамен старой
            deleteColumn(s);//удаляем s-ый столбец
            deleteElemFromImitBasis(s);
            reBasis(s);
            break;
        }

    }
    if (isNull == true) // а) строка r избыточна, нужно её удалить. Тогда столбец s будет полностью нулевым, поэтому его тоже удалим.
    {
        deleteRow(r);
        deleteColumn(s);
        deleteElemFromImitBasis(s);
        reBasis(s);
    }
}

void SimplexMethod::deleteElemFromImitBasis(const int& s) {
    for (int i = 0; i < _countImit; i++)
        if (_ImitBasis[i] == s)
            for (int j = i + 1; j < _countImit; j++)
                _ImitBasis[j - 1] = _ImitBasis[j];
    _countImit--;
    Basis ImitBasis_tmp(_countImit);
    for (int i = 0; i < _countImit; i++)
        ImitBasis_tmp[i] = _ImitBasis[i];

    _ImitBasis = ImitBasis_tmp;
    for (int i = 0; i < _countImit; i++)
    {
        if (_ImitBasis[i] > s) _ImitBasis[i]--;
    }
}

void SimplexMethod::reBasis(const int& s) {
    for (int i = 0; i < _basis.size(); i++)
        if (_basis[i] == s)
            for (int j = i + 1; j < _basis.size(); j++)
                _basis[j - 1] = _basis[j];
    Basis basis_tmp(_m);
    for (int i = 0; i < _m; i++)
        basis_tmp[i] = _basis[i];

    _basis = basis_tmp;
    for (int i = 0; i < _m; i++)
    {
        if (_basis[i] > s) _basis[i]--;
    }
}

void SimplexMethod::DeleteFromB(const int& r) {
    Vec b_tmp(_m - 1);
    for (int j = r + 1; j < _m; j++)
        _b[j - 1] = _b[j];

    for (int i = 0; i < _m - 1; i++)
        b_tmp[i] = _b[i];

    _b = b_tmp;
}

void SimplexMethod::DeleteFromC(const int& s) {
    bool c_copyCont = false;
    bool c_copy1Cont = false;
    if (_c_copy.size() > s) c_copyCont = true;
    if (_c_copy1.size() > s) c_copy1Cont = true;

    Vec C_tmp(_n - 1);
    Vec C_copy_tmp(_n - 1);
    Vec C_copy1_tmp(_n - 1);
    for (int j = s + 1; j < _n; j++)
    {
        _c[j - 1] = _c[j];
        if (c_copyCont) _c_copy[j - 1] = _c_copy[j];
        if (c_copy1Cont) _c_copy1[j - 1] = _c_copy1[j];
    }
    for (int i = 0; i < _n - 1; i++)
    {
        C_tmp[i] = _c[i];
        if (c_copyCont) C_copy_tmp[i] = _c_copy[i];
        if (c_copy1Cont) C_copy1_tmp[i] = _c_copy1[i];
    }
    if (c_copyCont)
    {
        _c_copy = C_copy_tmp;
    }
    if (c_copy1Cont)
    {
        _c_copy1 = C_copy1_tmp;
    }

    _c = C_tmp;
}

void SimplexMethod::deleteRow(const int& r) {
    Matrix A_tmp(_m - 1);
    for (auto& elem : A_tmp)
        elem.resize(_n);

    for (int j = 0; j < _n; j++)
    {
        for (int i = 0; i < r; i++)
            A_tmp[i][j] = _A[i][j];
        for (int i = r; i < _m - 1; i++)
            A_tmp[i][j] = _A[i + 1][j];
    }
    DeleteFromB(r);
    _m--;
    _A = A_tmp;
}

void SimplexMethod::deleteColumn(const int& s) {
    Matrix A_tmp(_m);
    for (auto& elem : A_tmp)
        elem.resize(_n - 1);

    for (int i = 0; i < _m; i++)
    {
        for (int j = 0; j < s; j++)
            A_tmp[i][j] = _A[i][j];
        for (int j = s; j < _n - 1; j++)
            A_tmp[i][j] = _A[i][j + 1];
    }
    DeleteFromC(s);
    _n--;
    _A = A_tmp;
}

void SimplexMethod::deleteImitColumn() {
    Matrix A_tmp(_m);
    for (auto& elem : A_tmp)
        elem.resize(_n - _countImit);

    for (int i = 0; i < _m; i++)
        for (int j = 0; j < _n - _countImit; j++)
            A_tmp[i][j] = _A[i][j];
    for (int im = 0; im < _countImit; im++)
    {
        DeleteFromC(_ImitBasis[im]);
        _n--;
    }
    _countImit = 0;
    _ImitBasis.resize(_countImit);
    _A = A_tmp;
}

double SimplexMethod::CalcRes(const int& stage) {
    double res = 0;
    _x.resize(_n);
    _result.change(_n, _m);

    for (int i = 0; i < _m; i++)
    {
        _x[_basis[i]] = _b[i];
        _result._xbi[i] = _basis[i];
        _result._x[_basis[i]] = _b[i];
    }

    for (int j = 0; j < _n; j++)
    {
        if (stage == 1) { res += _c_copy1[j] * _x[j]; _result._resVal = res; }
        if (stage == 2) { res += _c_copy[j] * _x[j]; _result._resVal = res; }
    }

    return res;
}

int SimplexMethod::FindBasis() {
    int countBasis = 0;
    int k;

    _basis.resize(_m);
    for (int i = 0; i < _m; i++) _basis[i] = -1;

    for (int i = 0; i < _m; i++)
        if (_b[i] < 0)
        {
            _b[i] *= -1;
            for (int j = 0; j < _n; j++)
                _A[i][j] *= -1;
        }

    Basis tmp(_m);
    for (int i = 0; i < _m; i++)
        for (int j = 0; j < _n; j++)
            if (_A[i][j] == 1)
            {
                if (std::find(_basis.begin(), _basis.end(), j) == _basis.end())
                {
                    for (k = 0; k < _m; k++)
                        if (k != i)
                            if (_A[k][j] != 0) break;
                    if (k == _m)
                    {
                        _basis[countBasis] = j; 
                        countBasis++; 
                        tmp[i] = 1;
                    }
                }
            }
    if (countBasis != _m) { _countImit = _m - countBasis; AddImitBasis(_countImit, tmp); }
    return _err_OK;
}

int SimplexMethod::AddImitBasis(int cntImit, Basis& tmp) {
    Matrix copy_A(_m);
    for (auto& elem : copy_A)
        elem.resize(_n);

    copy_A = _A;
    for (auto& elem : _A)
        elem.resize(_n + _countImit);

    Matrix Imit(_m);
    for (auto& elem : Imit)
        elem.resize(_countImit);

    int j = 0;
    while (j != _countImit)
        for (int i = 0; i < _m; i++)
            if (tmp[i] == 0) { Imit[i][j] = 1; j++; }

    for (int i = 0; i < _m; i++)
    {
        for (int k = 0; k < _countImit; k++) _A[i][_n + k] = Imit[i][k];
        for (j = 0; j < _n; j++) _A[i][j] = copy_A[i][j];
    }

    _ImitBasis.resize(_countImit);
    int jj = 0;
    for (int ii = _n; jj < _countImit; ii++)
    {
        _ImitBasis[jj] = ii;
        _basis[_m - _countImit + jj] = ii;
        jj++;
    }

    if (std::find(_basis.begin(), _basis.end(), -1) != _basis.end()) { 
        _err_OK = 1; 
        return _err_OK; 
    }

    _c.resize(_n + _countImit);
    for (int i = 0; i < _countImit; i++)
        _c[_ImitBasis[i]] = 1;

    _n += _countImit;
    _c_copy1.resize(_n);
    for (int i = 0; i < _n; i++)
        _c_copy1[i] = _c[i];

    Basis Basis_copy(_m);
    for (int i = 0; i < _m; i++)
        Basis_copy[i] = _basis[i];

    for (int i = 0; i < _m; i++)
        for (j = 0; j < _m; j++)
        {
            if (_A[i][Basis_copy[j]] == 1) _basis[i] = Basis_copy[j];
        }

    return _err_OK;
}

int SimplexMethod::SM_Algoritm() {
    int s = getLeadingColumn(); // выбираем ведущий столбец s
    if (s == -1) { 
        _result._res = RESSIMP::resOK; 
        return 0; 
    }    

    int r = getLeadingRow(s);
    if (r == -1) { 
        _result._res = RESSIMP::resUnlim; 
        return -1; 
    }

    gaussOperation(s, r);
    _basis[r] = s;

    return 1;
}

int SimplexMethod::getLeadingColumn() {
    int columnNumber = -1;
    int min = 9999999;
    for (int i = 0; i < _n; i++)
        if (_c[i] < 0)
        {
            if (_c[i] < min)
            {
                min = i;
                columnNumber = min;
            }
        }
    return columnNumber;
}

int SimplexMethod::getLeadingRow(int s) {
    int rowNumber = -1;
    int cnt_pos = 0;
    Matrix lexVectors(_m);
    for (auto& elem : lexVectors)
        elem.resize(_n + 1);

    int minlexVector;
    Vec leadRow(_m);
    for (int i = 0; i < _m; i++)
    {
        leadRow[i] = _A[i][s];
        if (leadRow[i] > 0)
        {
            cnt_pos++;
            for (int j = 0; j < _n; j++)
                lexVectors[i][j + 1] = _A[i][j] / leadRow[i];
            lexVectors[i][0] = _b[i] / leadRow[i];
        }
    }

    if (cnt_pos == 0) return rowNumber;

    minlexVector = 0;
    for (int i = 1; i < _m; i++)
        if (compareLexVectors(lexVectors[minlexVector], lexVectors[i]) > 0)
            minlexVector = i;

    rowNumber = minlexVector;
    return rowNumber;
}

int SimplexMethod::compareLexVectors(Vec& v1, Vec& v2) {
    int result = 0;
    int widht = v1.size();
    Vec diff(widht);
    bool flagV1 = true;
    bool flagV2 = true;
    for (int j = 0; j < widht; j++)
    {
        if (v1[j] != 0) flagV1 = false;
        if (v2[j] != 0) flagV2 = false;
        diff[j] = v1[j] - v2[j];
    }
    if (flagV1 && flagV2) return 0;
    if (flagV1) return 1;
    if (flagV2) return -1;


    for (int j = 0; j < widht; j++)
    {
        if (diff[j] == 0) continue;
        else
        {
            if (diff[j] > 0) result = 1;
            if (diff[j] < 0) result = -1;
            break;
        }
    }
    return result;
}

void SimplexMethod::gaussOperation(int s, int r) {
    double delit;
    delit = -_c[s] / _A[r][s];
    for (int j = 0; j < _n; j++)
    {
        double temp1 = _A[r][j] * delit;
        _c[j] += temp1;
    }

    for (int i = 0; i < _m; i++)
    {
        if (i == r) continue;
        double delit2 = -_A[i][s] / _A[r][s];
        for (int j = 0; j < _n; j++)
        {
            double tmp = _A[r][j] * delit2;
            _A[i][j] += tmp;
        }

        _b[i] += _b[r] * delit2;
    }

    double lVal = _A[r][s];
    for (int j = 0; j < _n; j++)
    {
        _A[r][j] = _A[r][j] / lVal;
    }

    _b[r] = _b[r] / lVal;
}
