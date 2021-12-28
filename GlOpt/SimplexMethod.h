#ifndef SIMPLEX_METHOD_H
#define SIMPLEX_METHOD_H

#include <vector>

using Vec = std::vector<double>;
using Matrix = std::vector<std::vector<double>>;
using Basis = std::vector<int>;

/*
*		this is linear optimization solver. it uses simplex method.
* 
*		it solves the following task:
*		<c, x> -> min
*		Ax = b
*		x >= 0
* 
*		make sure that you transformed all inequalities to equalities
*/

enum class Solution { optimal, unlimited, inconsistent, notfound };

class SimplexMethod {
private: // Fields
	Matrix _simplex_table;
	// size of objective vector
	int _n;
	// size of constraints
	int _m;
	// numbers of basis columns
	Basis _basis;
	// solution's state
	Solution _state;
	// marker for M-method
	bool _is_m_method;
private: // Methods
	int find_basis();
	int detect_leading_column();
	int detect_leading_row(const int& leading_column);
	void basis_transform(const int& lr, const int& lc);
	void gauss_transform(const int& lr, const int& lc);
	void iterate();
	void solve_imit(const int& quantity_imit);
	void gauss_tranform_imit();
	int check_for_imit(const int& dim);
	void exclude_imit_column(const int& imit_column, const int& dim);
	double get_elem(const int& row, const int& column) const;
	void transform_task(const SimplexMethod& out);
	void transform_for_basis();
	double detect_one(const int& basis);
	Basis get_basis() const;
public:
	SimplexMethod(const Vec& c, const Vec& b, const Matrix& A, const bool& m = false);
	void solve();
	double get_solution() const { return _simplex_table[0][0]; }
	Solution get_state() const { return _state; }
	double get_basis_size() const { return _basis.size(); }
	~SimplexMethod() {}
};

#endif // SIMPLEX_METHOD_H

