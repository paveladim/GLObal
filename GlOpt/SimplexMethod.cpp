#include "SimplexMethod.h"

SimplexMethod::SimplexMethod(const Vec& c, const Vec& b, const Matrix& m) :
	_n(c.size()), 
	_m(b.size()),
	_c(c),
	_c_copy(c), 
	_b(b),
	_matrix(m) {

}

void SimplexMethod::make_canonical() {
	for (uint16_t i = 0; i < _m; ++i)
		if (_b[i] < 0) {
			_b[i] *= (-1);
			for (uint16_t j = 0; j < _n; ++j) _matrix.mul_elem(i, j, -1.0);
		}
}

uint16_t SimplexMethod::find_basis() {
	_basis.resize(_m);
	make_canonical();
	uint16_t pos = 0;
	// запоминаем индексы единичных столбцов
	std::vector<bool> places(_m);
	uint16_t row_num = -1;
	// предполагается, что базис изначально не найден
	for (auto& elem : _basis) elem = -1;
	// ищем базис
	bool basis_column = false;
	for (uint16_t j = 0; j < _n; ++j) {
		// предполагаем, что столбец не единичный
		basis_column = false;
		for (uint16_t i = 0; i < _m; ++i) {
			if (_matrix.get_value(i, j) == 1) {
				basis_column = true;
				for (uint16_t k = 0; (k < _m) && (basis_column); ++k)
					if ((_matrix.get_value(k, j) != 0) && (k != i))
						basis_column = false;

				if (basis_column) row_num = i;
			}
		}

		if (basis_column) {
			_basis[pos] = j;
			places[row_num] = true;
			++pos;
		}
	}

	if (_m - pos != 0) create_imitation_basis(_m - pos, places);
	return (_m - pos);
}

void SimplexMethod::create_imitation_basis(const uint16_t& basis_size,
	const std::vector<bool>& places) {
	Matrix imitation_matrix(_m, basis_size, 'n');
	uint16_t j = 0;
	for (uint16_t i = 0; i < _m; ++i)
		if ((!places[i]) && (j < basis_size)) {
			imitation_matrix.set_value(i, j, 1.0);
			++j;
		}

	_matrix.attach_matrix(imitation_matrix);

	_imit_basis.resize(basis_size);
	uint16_t j = 0;
	for (uint16_t i = _n; j < basis_size; ++i) {
		//добавили в допустимый базис искусственную переменную
		_imit_basis[j] = i;
		_basis[_m - basis_size + j] = i;
		++j;
	}

	// сохраняем исходную задачу и вводим новую
	_c_copy = _c;
	_c.resize(_n + basis_size);
	_n = _n + basis_size;
	for (uint16_t j = 0; j < _n; ++j)
		if (j < basis_size) _c[j] = 1.0;
		else _c[j] = 0.0;

	// упорядочить столбцы единичной матрицы
	std::vector<int16_t> basis_copy(_m);
	for (uint16_t i = 0; i < _m; ++i)
		for (uint16_t j = 0; j < _n; ++j)
			if (std::find(_basis.begin(), _basis.end(), j) != _basis.end()) {
				if (_matrix.get_value(i, j) == 1) basis_copy[i] = j;
			}

	_basis = basis_copy;
}

void SimplexMethod::gauss_transform(const uint16_t& r, const uint16_t& s) {
	double coeff = -_c[s] / _matrix.get_value(r, s);

	for (uint16_t j = 0; j < _n; ++j) {
		double tmp = coeff * _matrix.get_value(r, j);
		_c[j] += tmp;
	}

	// преобразование с матрицей
	for (uint16_t i = 0; i < _m; ++i) {
		if (i != r) {
			coeff = -_matrix.get_value(i, s) / _matrix.get_value(r, s);

			for (int j = 0; j < _n; ++j)
			{
				double tmp = coeff * _matrix.get_value(r, j);
				_matrix.add_elem(i, j, tmp);
			}

			_b[i] += coeff * _b[r];
		}
	}

	// преобразования с ведущей строкой
	coeff = _matrix.get_value(r, s);
	for (uint16_t j = 0; j < _n; ++j)
		_matrix.mul_elem(r, j, 1.0 / coeff);

	// для вектора B
	_b[r] = _b[r] / coeff;
}

void SimplexMethod::find_solution() {
	uint16_t imit_basis_size = find_basis();
	if (imit_basis_size != 0) {
		uint16_t tmp1 = 0;
		uint16_t tmp2 = 0;

		for (uint16_t k = 0; k < imit_basis_size; ++k)
		{
			tmp1 = _imit_basis[k];
			for (uint i = 0; i < _m; ++i)
				if (_matrix.get_value(i, tmp1) == 1) tmp2 = i;

			gauss_transform(tmp1, tmp2);
		}
	}

	_c.resize(_n);
	for (uint16_t i = 0; i < _n; ++i)
		_c[i] = _c_copy[i];

	uint16_t tmp1 = 0;
	uint16_t tmp2 = 0;
	for (uint16_t k = 0; k < _basis.size(); k++)
	{
		tmp1 = _basis[k];
		for (uint16_t i = 0; i < _m; ++i)
			if (_matrix.get_value(i, tmp1) == 1) tmp2 = i;
		gauss_transform(tmp1, tmp2);
	}
}