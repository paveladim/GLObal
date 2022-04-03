#include "SimplexMethod-x.h"

SimplexMethod::SimplexMethod(const Vec& C, const Vec& B, const Matrix& a) :
	c(C), b(B), A(a), n(C.size()), m(B.size()), result(n, m), eps(0.000001), err_OK(0)
{

}

TSimplexRes SimplexMethod::FindSolution() {
    c_copy.resize(n);
    for (int i = 0; i < n; i++)
        c_copy[i] = c[i];

    double calc_res;
    result.resize(n, m);
    int res_OK = -1;

    //»щем базис
    if (FindBasis() == 0) // если базис найден, то примен€ем алгоритм 1 (стр.35)
        if (count_imit != 0) // есть искусственные переменные, значит будем решать в 2 этапа
        // Ё“јѕ 1
        {
            int tmp_1 = 0, tmp_2 = 0;
            for (int k = 0; k < count_imit; k++)
            {
                tmp_1 = imit_basis[k];
                for (int i = 0; i < m; i++)
                    if (A[i][tmp_1] == 1) tmp_2 = i;
                gaussOperation(tmp_1, tmp_2);
            }
            do
            {
                res_OK = SM_Algorithm();
                if (res_OK == -1) { result.res = RESSIMP::resUnlim; return result; }
            } while (res_OK != 0);

            calc_res = CalcRes(1);
            if (calc_res > 0)
            {
                result.res = RESSIMP::systemIsInconsist;
                return result;
            }
            if (calc_res < 0) { result.res = RESSIMP::systemIsInconsist; return result; }
            if (calc_res == 0) 
            {
                bool flag = false; 
                for (int ii = 0; ii < count_imit; ii++)
                {
                    if (std::find(basis.begin(), basis.end(), imit_basis[ii]) != basis.end())
                    {                                  
                        flag = true;
                        excludeImitationVariables(imit_basis[ii]);   
                        ii--;
                    }
                    flag = false;
                }
                if (!flag)
                    deleteImitColumn();
            }
        }
    // Ё“јѕ 2
    c.resize(n);
    for (int i = 0; i < n; i++)
        c[i] = c_copy[i];

    int tmp1 = 0;
    int tmp2 = 0;

    for (int k = 0; k < basis.size(); k++)
    {
        tmp1 = basis[k];
        for (int i = 0; i < m; i++)
            if (A[i][tmp1] == 1) tmp2 = i;
        gaussOperation(tmp1, tmp2);
    }
    do
    {
        res_OK = SM_Algorithm(); // примен€ем алгоритм 1, пока не найдЄм решение 
        calc_res = CalcRes(2);
    } while (res_OK != 0);

    return result;
}

void SimplexMethod::excludeImitationVariables(const int& s) {
    int r = 0;
    bool isNull = true;

    for (int i = 0; i < m; i++)
        if (A[i][s] == 1) r = i;
    for (int j = 0; j < n - count_imit; j++)
    {
        if (A[r][j] != 0)
        {
            isNull = false;
            gaussOperation(j, r);
            basis[r] = j;
            deleteColumn(s);
            deleteElemFromImitBasis(s);
            reBasis(s);
            break;
        }

    }
    if (isNull == true)
    {
        deleteRow(r);
        deleteColumn(s);
        deleteElemFromImitBasis(s);
        reBasis(s);
    }
}

void SimplexMethod::deleteElemFromImitBasis(const int& s) {
    int crit_val = s;
    auto elem = std::find(imit_basis.begin(), imit_basis.end(), s);
    if (elem != imit_basis.end()) {
        imit_basis.erase(elem);
        count_imit--;

        for (int i = 0; i < count_imit; i++)
        {
            if (imit_basis[i] > crit_val) imit_basis[i]--;
        }
    }
}

void SimplexMethod::reBasis(const int& s) {
    int crit_val = s;
    auto elem = std::find(basis.begin(), basis.end(), s);
    if (elem != basis.end()) {
        basis.erase(elem);
        basis.push_back(basis[m - 1]);

        for (int i = 0; i < m; i++)
        {
            if (basis[i] > crit_val) basis[i]--;
        }
    }
}

void SimplexMethod::DeleteFromB(const int& r) {
    Vec b_tmp(m - 1);
    for (int j = r + 1; j < m; j++)
        b[j - 1] = b[j];

    for (int i = 0; i < m - 1; i++)
        b_tmp[i] = b[i];

    b.resize(m - 1);
    b = b_tmp;
}

void SimplexMethod::DeleteFromC(const int& s) {
    bool c_copyCont = false;
    bool c_copy1Cont = false;
    if (c_copy.size() > s) c_copyCont = true;
    if (c_copy1.size() > s) c_copy1Cont = true;

    Vec C_tmp(n - 1);
    Vec C_copy_tmp(n - 1);
    Vec C_copy1_tmp(n - 1);

    for (int j = s + 1; j < n; j++)
    {
        c[j - 1] = c[j];
        if (c_copyCont) c_copy[j - 1] = c_copy[j];
        if (c_copy1Cont) c_copy1[j - 1] = c_copy1[j];
    }

    for (int i = 0; i < n - 1; i++)
    {
        C_tmp[i] = c[i];
        if (c_copyCont) C_copy_tmp[i] = c_copy[i];
        if (c_copy1Cont) C_copy1_tmp[i] = c_copy1[i];
    }

    if (c_copyCont)
    {
        c_copy.resize(n - 1);
        c_copy = C_copy_tmp;
    }

    if (c_copy1Cont)
    {
        c_copy1.resize(n - 1);
        c_copy1 = C_copy1_tmp;
    }

    c.resize(n - 1);
    c = C_tmp;
}

void SimplexMethod::deleteRow(const int& r) {
    Matrix A_tmp(m - 1);
    for (auto& elem : A_tmp)
        elem.resize(n);

    for (int j = 0; j < n; j++)
    {
        for (int i = 0; i < r; i++)
            A_tmp[i][j] = A[i][j];
        for (int i = r; i < m - 1; i++)
            A_tmp[i][j] = A[i + 1][j];
    }

    DeleteFromB(r);
    m--;

    A = A_tmp;
}

void SimplexMethod::deleteColumn(const int& s) {
    Matrix A_tmp(m);
    for (auto& elem : A_tmp)
        elem.resize(n - 1);

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < s; j++)
            A_tmp[i][j] = A[i][j];
        for (int j = s; j < n - 1; j++)
            A_tmp[i][j] = A[i][j + 1];
    }

    DeleteFromC(s);
    n--;
    A = A_tmp;
}

void SimplexMethod::deleteImitColumn() {
    Matrix A_tmp(m);
    for (auto& elem : A_tmp)
        elem.resize(n - count_imit);

    for (int i = 0; i < m; i++)
        for (int j = 0; j < n - count_imit; j++)
            A_tmp[i][j] = A[i][j];

    for (int im = 0; im < count_imit; im++)
    {
        DeleteFromC(imit_basis[im]);
        n--;
    }

    count_imit = 0;

    imit_basis.resize(count_imit);
    A = A_tmp;
}

double SimplexMethod::CalcRes(const int& stage) {
    double res = 0;
    x.resize(n);
    result.resize(n, m);

    for (int i = 0; i < m; i++) {
        x[basis[i]] = b[i];
        result.xbi[i] = basis[i];
        result.x[basis[i]] = b[i];
    }

    for (int j = 0; j < n; j++) {
        if (stage == 1) { res += c_copy1[j] * x[j]; result.resVal = res; }
        if (stage == 2) { res += c_copy[j] * x[j]; result.resVal = res; }
    }

    return res;
}

int SimplexMethod::FindBasis() {
    int countBasis = 0;
    int k;

    basis.resize(m);
    for (int i = 0; i < m; i++) basis[i] = -1;

    for (int i = 0; i < m; i++)
        if (b[i] < 0)
        {
            b[i] *= -1;
            for (int j = 0; j < n; j++)
                A[i][j] *= -1;
        }

    Basis tmp(m);

    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            if (A[i][j] == 1)
            {
                if (!(std::find(basis.begin(), basis.end(), j) != basis.end()))
                {
                    for (k = 0; k < m; k++)
                        if (k != i)
                            if (A[k][j] != 0) break;
                    if (k == m) // если остальные значени€ столбца это нули, то добавл€ем его в набор базисных
                    {
                        basis[countBasis] = j; countBasis++; tmp[i] = 1;
                    }
                }
            }
    if (countBasis != m) { count_imit = m - countBasis; AddImitBasis(count_imit, tmp); }
    return err_OK;
}

int SimplexMethod::AddImitBasis(const int& cntImit, const Basis& tmp) {
    Matrix copy_A { A };
    A.resize(m);
    for (auto& elem : A)
        elem.resize(n + count_imit);

    Matrix Imit(m);
    for (auto& elem : Imit)
        elem.resize(count_imit);

    int j = 0;
    while (j != count_imit)
        for (int i = 0; i < m; i++)
            if (tmp[i] == 0) { Imit[i][j] = 1; j++; }

    for (int i = 0; i < m; i++)
    {
        for (int k = 0; k < count_imit; k++) A[i][n + k] = Imit[i][k];
        for (j = 0; j < n; j++) A[i][j] = copy_A[i][j];
    }

    imit_basis.resize(count_imit); 
    int jj = 0;
    for (int ii = n; jj < count_imit; ii++)
    {
        imit_basis[jj] = ii;
        basis[m - count_imit + jj] = ii;
        jj++;
    }

    if (std::find(basis.begin(), basis.end(), -1) != basis.end()) { 
        err_OK = 1; 
        return err_OK; 
    }

    c.resize(n + count_imit);
    for (int i = 0; i < count_imit; i++)
        c[imit_basis[i]] = 1;

    n += count_imit;
    c_copy1.resize(n);
    for (int i = 0; i < n; i++)
        c_copy1[i] = c[i];

    Basis Basis_copy(m);
    for (int i = 0; i < m; i++)
        Basis_copy[i] = basis[i];

    for (int i = 0; i < m; i++)
        for (j = 0; j < m; j++)
        {
            if (A[i][Basis_copy[j]] == 1) basis[i] = Basis_copy[j];
        }
    return err_OK;
}

int SimplexMethod::SM_Algorithm() {
    int s = getLeadingColumn();
    if (s == -1) { result.res = RESSIMP::resOK; return 0; }
    int r = getLeadingRow(s); 
    if (r == -1) { result.res = RESSIMP::resUnlim; return -1; }
    gaussOperation(s, r);
    basis[r] = s; 

    return 1;
}

int SimplexMethod::getLeadingColumn() {
    int columnNumber = -1;
    int min = 9999999;
    for (int i = 0; i < n; i++)
        if (c[i] < 0) {
            if (c[i] < min) {
                min = i;
                columnNumber = min;
            }
        }
    return columnNumber;
}

int SimplexMethod::getLeadingRow(const int& s) {
    int rowNumber = -1;
    int cnt_pos = 0; 

    Matrix lexVectors(m);
    for (int i = 0; i < m; i++)
        lexVectors[i].resize(n + 1);

    int minlexVector;

    Vec leadRow(m);
    for (int i = 0; i < m; i++)
    {
        leadRow[i] = A[i][s];
        if (leadRow[i] > 0)
        {
            cnt_pos++;
            for (int j = 0; j < n; j++)
                lexVectors[i][j + 1] = A[i][j] / leadRow[i]; 
            lexVectors[i][0] = b[i] / leadRow[i]; 
        }
    }

    if (cnt_pos == 0) return rowNumber;

    minlexVector = 0;
    for (int i = 1; i < m; i++)
        if (compareLexVectors(lexVectors[minlexVector], lexVectors[i]) > 0)
            minlexVector = i;

    rowNumber = minlexVector;

    return rowNumber;
}

void SimplexMethod::gaussOperation(const int& s, const int& r) {
    double delit;
    delit = -c[s] / A[r][s];
    for (int j = 0; j < n; j++)
    {
        double temp1 = A[r][j] * delit;
        c[j] += temp1;
    }

    for (int i = 0; i < m; i++)
    {
        if (i == r) continue;
        double delit2 = -A[i][s] / A[r][s];
        for (int j = 0; j < n; j++)
        {
            double tmp = A[r][j] * delit2;
            A[i][j] += tmp;
        }
        b[i] += b[r] * delit2;
    }

    double lVal = A[r][s];
    for (int j = 0; j < n; j++)
    {
        A[r][j] = A[r][j] / lVal;
    }

    b[r] = b[r] / lVal;
}

int SimplexMethod::compareLexVectors(const Vec& v1, const Vec& v2) {
    int result = 0;
    int width = v1.size();
    Vec razn(width);
    bool flagV1 = true;
    bool flagV2 = true;

    for (int j = 0; j < width; j++)
    {
        if (v1[j] != 0) flagV1 = false; 
        if (v2[j] != 0) flagV2 = false;
        razn[j] = v1[j] - v2[j];
    }
    if (flagV1 && flagV2) return 0;
    if (flagV1) return 1;
    if (flagV2) return -1;


    for (int j = 0; j < width; j++) // проходим по все эллементам и смотрим чтобы первый ненулевой был > 0
    {
        if (razn[j] != 0)
        {
            if (razn[j] > 0) result = 1;
            if (razn[j] < 0) result = -1;
            break;
        }
    }
    return result;
}
