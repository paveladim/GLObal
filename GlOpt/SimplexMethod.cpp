#include "SimplexMethod.h"

SimplexMethod::SimplexMethod(const Vec& c, const Vec& b, const Matrix& A, 
	                         const bool& m) :
	_n(c.size()), 
	_m(b.size()),
	_state(Solution::notfound),
	_is_m_method(m),
	x1(c.size() + 1),
	x2(c.size() + 1) {

	// making the simplex table
	_simplex_table.resize(_m + 1);
	for (auto& row : _simplex_table) row.resize(_n + 1);

	_simplex_table[0][0] = 0.0;
	for (int j = 1; j < _n + 1; ++j) _simplex_table[0][j] = -c[j - 1];
	for (int i = 1; i < _m + 1; ++i) _simplex_table[i][0] = b[i - 1];
	for (int i = 1; i < _m + 1; ++i)
		for (int j = 1; j < _n + 1; ++j)
			_simplex_table[i][j] = A[i - 1][j - 1];

	for (int i = 1; i < _m + 1; ++i)
		if (b[i - 1] < 0)
			for (int j = 0; j < _n + 1; ++j)
				_simplex_table[i][j] *= -1;

	// basis is not detected
	_basis.resize(_m);
	for (auto& basis_col_number : _basis) basis_col_number = -1;

	if (m) gauss_tranform_imit();
}

int SimplexMethod::find_basis() {
	int basis_num = 0;
	for (int j = 1; j < _n + 1; ++j) {
		bool is_basis_column = true;
		bool is_one = false;
		for (int i = 0; (i < _m + 1) && (is_basis_column); ++i) {
			if ((_simplex_table[i][j] != 0) && (_simplex_table[i][j] != 1)) 
				is_basis_column = false;

			if (_simplex_table[i][j] == 1) {
				if (!is_one) is_one = true;
				else is_basis_column = false;
			}
		}

		if (is_basis_column) {
			_basis[basis_num] = j;
			++basis_num;
		}
	}

	return _m - basis_num;
}

int SimplexMethod::detect_leading_column() {
	int leading_column{ -1 };
	double max{ 0.0 };

	for (int j = 1; (j < _n + 1) && (leading_column == -1); ++j)
		if (_simplex_table[0][j] > 0) leading_column = j;

	//for (int j = 1; j < _n + 1; ++j)
	//	if (_simplex_table[0][j] > max) {
	//		max = _simplex_table[0][j];
	//		leading_column = j;
	//	} 

	return leading_column;
}

int SimplexMethod::detect_leading_row(const int& leading_column) {
	int leading_row{ -1 };
	int start = -1;
	double temp_ratio = 0.0;

	Vec v1(_n + 1);
	Vec v2(_n + 1);

	for (int i = 1; (i < _m + 1) && (start == -1); ++i)
		if (_simplex_table[i][leading_column] > 0.0) {
			start = i;
			for (int j = 0; j < _n + 1; ++j) 
				v1[j] = _simplex_table[start][j] / _simplex_table[start][leading_column];
		}

	if (start == -1) return leading_row;
	leading_row = start;

	for (int i = start + 1; i < _m + 1; i++) {
		if (_simplex_table[i][leading_column] > 0.0) {
			for (int j = 0; j < _n + 1; ++j)
				v2[j] = _simplex_table[i][j] / _simplex_table[i][leading_column];

			if (lex_compare(v1, v2) == 1) {
				v1 = v2;
				leading_row = i;
			}
		}
	}

	return leading_row;
}

int SimplexMethod::lex_compare(const Vec& a, const Vec& b) {
	for (int i = 0; i < _n + 1; ++i) {
		if (a[i] - b[i] > 0) return 1;
		if (a[i] - b[i] < 0) return -1;
	}

	return 0;
}

// lr means leading row, lc means leading column
void SimplexMethod::basis_transform(const int& lr, const int& lc) {
	bool column_detected = false;
	for (int i = 0; ((i < _basis.size()) && (!column_detected)); ++i)
		if (_simplex_table[lr][_basis[i]] != 0) {
			_basis[i] = lc;
			column_detected = true;
		}
}

// lr means leading row, lc means leading column
void SimplexMethod::gauss_transform(const int& lr, const int& lc) {
	double divider = _simplex_table[lr][lc];
	for (int j = 0; j < _n + 1; ++j) 
		_simplex_table[lr][j] /= divider;

	for (int i = 0; i < lr; ++i) {
		double multiplicator = _simplex_table[i][lc];
		for (int j = 0; j < _n + 1; ++j)
			_simplex_table[i][j] -= _simplex_table[lr][j] * multiplicator;
	}

	for (int i = lr + 1; i < _m + 1; ++i) {
		double multiplicator = _simplex_table[i][lc];
		for (int j = 0; j < _n + 1; ++j)
			_simplex_table[i][j] -= _simplex_table[lr][j] * multiplicator;
	}
}

void SimplexMethod::iterate() {
	int leading_column = detect_leading_column();
	if (leading_column == -1) {
		_state = Solution::optimal; 
		return;
	}

	int leading_row = detect_leading_row(leading_column);
	if (leading_row == -1) {
		_state = Solution::unlimited;
		return;
	}

	basis_transform(leading_row, leading_column);
	gauss_transform(leading_row, leading_column);
}

void SimplexMethod::solve_imit(const int& quantity_imit) {
	Vec c(_n + quantity_imit);
	Vec b(_m);
	Matrix A(_m);
	for (auto& row : A) row.resize(_n + quantity_imit);

	for (int j = 0; j < _n; ++j) c[j] = 0.0;
	for (int j = _n; j < c.size(); ++j) c[j] = 1.0;

	for (int i = 1; i < _m + 1; ++i) b[i - 1] = _simplex_table[i][0];

	for (int j = 1; j < _n + 1; ++j)
		for (int i = 1; i < _m + 1; ++i)
			A[i - 1][j - 1] = _simplex_table[i][j];

	//int k = 0;
	//for (int j = _n; j < _n + quantity_imit; ++j) {
	//	A[k][j] = 1.0;
	//	++k;
	//}

	int for_basis_columns = 0;
	int j = _n;
	while (j < c.size()) {
		bool exist = false;
		for (int i = 0; i < _basis.size(); ++i)
			if (_basis[i] > -1)
				if (_simplex_table[for_basis_columns + 1][_basis[i]] == 1)
					exist = true;

		if (!exist) {
			A[for_basis_columns][j] = 1.0;
			++j;
		}

		++for_basis_columns;
	}

	for (auto& elem : _basis)
		if (elem > -1) c[elem] = 1.0;

	SimplexMethod sm_imit(c, b, A, true);
	sm_imit.solve();
	if (std::abs(sm_imit.get_solution()) > 1e-4) _state = Solution::inconsistent;
	else {
		int imit_row = sm_imit.check_for_imit(_n);
		while (imit_row != -1) {
			sm_imit.exclude_imit_column(imit_row, _n);
			_m = sm_imit.get_basis_size();
			imit_row = sm_imit.check_for_imit(_n);
		}
		transform_task(sm_imit);
	}
}

void SimplexMethod::gauss_tranform_imit() {
	int leading_row{ 0 };
	int leading_column{ _n };

	for (int i = 1; i < _n + 1; ++i)
		if (_simplex_table[0][i] != 0) {
			leading_column = i;
			for (int j = 1; j < _m + 1; ++j)
				if (_simplex_table[j][leading_column] != 0)
					leading_row = j;

			gauss_transform(leading_row, leading_column);
		}

	//int border{ 0 };
	//for (border = 1;
	//	(border < _n + 1) && (_simplex_table[0][border] == 0);
	//	++border);
	//--border;

	//while (leading_column > border) {
	//	for (int i = 1; i < _m + 1; ++i)
	//		if (_simplex_table[i][leading_column] == 1.0)
	//			leading_row = i;

	//	gauss_transform(leading_row, leading_column);
	//	--leading_column;
	//}
}

int SimplexMethod::check_for_imit(const int& dim) {
	for (auto& elem : _basis)
		if (elem > dim) return elem;

	return -1;
}

void SimplexMethod::exclude_imit_column(const int& imit_column,
										const int& dim) {
	int row = -1;
	for (int i = 1; i < _m + 1; ++i)
		if (_simplex_table[i][imit_column] == 1) row = i;

	bool everyone_is_zero = true;
	int leading_column = -1;
	for (int j = 0; j < dim + 1; ++j)
		if ((_simplex_table[row][j] != 0) && (j != imit_column)) {
			everyone_is_zero = false;
			if ((j > 0) && (j < dim + 1) && (leading_column == -1))
				leading_column = j;
		}

	if (everyone_is_zero) {
		_simplex_table.erase(_simplex_table.begin() + row);
		auto for_erase = std::find(_basis.begin(), _basis.end(), imit_column);
		_basis.erase(for_erase);
	}
	else {
		basis_transform(row, leading_column);
		gauss_transform(row, leading_column);
	}
}

double SimplexMethod::get_elem(const int& row, const int& column) const {
	return _simplex_table[row][column];
}

Basis SimplexMethod::get_basis() const {
	return _basis;
}

void SimplexMethod::transform_task(const SimplexMethod& out) {
	Vec copy_c(_n);
	for (int i = 1; i < _n + 1; ++i)
		copy_c[i - 1] = _simplex_table[0][i];

	_simplex_table.resize(_m + 1);
	for (auto& elem : _simplex_table)
		elem.resize(_n + 1);

	for (int i = 1; i < _m + 1; ++i)
		for (int j = 0; j < _n + 1; ++j)
			_simplex_table[i][j] = out.get_elem(i, j);

	for (int i = 1; i < _n + 1; ++i)
		_simplex_table[0][i] = copy_c[i - 1];

	_basis = out.get_basis();
}

void SimplexMethod::transform_for_basis() {
	int leading_row = -1;
	for (int i = 0; i < _basis.size(); ++i) {
		for (int j = 1; j < _m + 1; ++j) {
			if (_simplex_table[j][_basis[i]] == 1.0)
				leading_row = j;
		}

		gauss_transform(leading_row, _basis[i]);
	}
}

double SimplexMethod::detect_one(const int& basis) {
	int row = -1;

	for (int j = 1; j < _m + 1; ++j)
		if (_simplex_table[j][basis] == 1.0) {
			row = j;
			return row;
		}

	return row;
}

void SimplexMethod::solve() {
	int imit_basis_quantity = find_basis();
	if (imit_basis_quantity == 0)
		while (_state == Solution::notfound) iterate();
	else {
		solve_imit(imit_basis_quantity);
		
		if (_state == Solution::inconsistent) return;
		else {
			transform_for_basis();
			while (_state == Solution::notfound) iterate();
		}
	}
}
