#include "SimplexMethod.h"

SimplexMethod::SimplexMethod(const Vec& c, const Vec& b, const Matrix& m) :
	_n(c.size()), 
	_m(b.size()),
	_c(c),
	_c_copy(c), 
	_b(b),
	_matrix(m),
	_result(_n, _m) {

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

uint16_t SimplexMethod::get_leading_column() {
	uint16_t column = std::numeric_limits<uint16_t>::max();
	double minimum = std::numeric_limits<double>::max();

	for (uint16_t i = 0; i < _n; ++i)
		if (_c[i] < 0)
		{
			if (_c[i] < minimum)
			{
				minimum = _c[i];
				column = i;
			}
		}
	return column;
}

uint16_t SimplexMethod::get_leading_row(const uint16_t& s) {
	uint16_t rowNumber = std::numeric_limits<uint16_t>::max();
	// кол-во положительных эл-ов в ведущем столбце
	uint16_t cnt_pos = 0;
	Matrix lex_vec(_m, _n + 1, 'o');
	// номер строки лексикографически миним-й вектора
	uint16_t min_lex_vec;
	// вектор ведущего столбца
	Vec lead_row(_m);

	for (uint16_t i = 0; i < _m; ++i)
	{
		lead_row[i] = _matrix.get_value(i, s); // заполняем ведущий столбец
		if (lead_row[i] > 0) //выбираем из ведущего столбца строку с положительным эл-ом
		{
			cnt_pos++;
			// составляем вектор для сравнения
			for (uint16_t j = 0; j < _n; ++j)
				lex_vec.set_value(i, j + 1, _matrix.get_value(i, j) / lead_row[i]);

			lex_vec.set_value(i, 0, _b[i] / lead_row[i]);
		}
	}

	// если нет положительных эл-ов в ведущем столбце, то значение целевой ф-ии не ограничено
	if (cnt_pos == 0) return rowNumber;

	// выбирем лексиграфически минимальный вектор:
	min_lex_vec = 0;
	for (uint16_t i = 1; i < _m; ++i)
		if (lex_vec.compare_lex(min_lex_vec, i) > 0) min_lex_vec = i;

	rowNumber = min_lex_vec; //  номер лексикографически минимального вектора
	return rowNumber;
}

int SimplexMethod::step() {
	//Шаг 1.
	uint16_t s = get_leading_column(); // выбираем ведущий столбец s
	// все эл-ты вектора С неотрицательны, то соответствующее решение оптимально
	if (s == std::numeric_limits<uint16_t>::max()) return 0;
	//Шаг 2.
	uint16_t r = get_leading_row(s); //выбираем ведущую строку
	// значение целевой ф-ии неограничено
	if (r == std::numeric_limits<uint16_t>::max()) return -1;
	//Шаг 3.
	gauss_transform(r, s); //выполняем шаг гауссова преобразования
	//Шаг 4.
	_basis[r] = s; // устанавливаем новую базисную переменной взамен старой

	return 1;
}

int SimplexMethod::calculate_result(const bool& phase) {
	double res = 0;
	// выделяем память под вектор решения
	x.resize(_n);
	_result.change_size(_n, _m);

	for (uint16_t i = 0; i < _m; ++i) {
		x[_basis[i]] = _b[i]; // вектор, соответсвующий текущему базису
		_result._basis[i] = _basis[i];
		_result._x[_basis[i]] = _b[i];
	}

	for (uint16_t j = 0; j < _n; ++j) {
		if (!phase) { res += _c_ccopy[j] * x[j]; _result._value = res; }
		else { res += _c_copy[j] * x[j]; _result._value = res; }
	}

	return res;
}

void SimplexMethod::delete_imit_column(uint16_t& imit_basis) {
	Matrix A_tmp(_m, _n - imit_basis, 'n');

	for (uint16_t i = 0; i < _m; ++i)
		for (uint16_t j = 0; j < _n - imit_basis; ++j)
			A_tmp.set_value(i, j, _matrix.get_value(i, j));

	for (uint16_t i = 0; i < imit_basis; ++i) {
		delete_from_c(_imit_basis[i]);
		--_n;
	}

	imit_basis = 0;
	_imit_basis.resize(imit_basis);
	_matrix = A_tmp;
}

void SimplexMethod::delete_from_c(const uint16_t& s) {
	bool c_copyCont = false;
	bool c_copy1Cont = false;
	if (_c_copy.size() > s) c_copyCont = true;
	if (_c_ccopy.size() > s) c_copy1Cont = true;

	Vec c_tmp(_n - 1);
	Vec c_copy_tmp(_n - 1);
	Vec c_ccopy_tmp(_n - 1);
	
	for (uint16_t j = s + 1; j < _n; ++j) {
		_c[j - 1] = _c[j];
		if (c_copyCont) _c_copy[j - 1] = _c_copy[j];
		if (c_copy1Cont) _c_ccopy[j - 1] = _c_ccopy[j];
	}

	for (uint16_t i = 0; i < _n - 1; ++i) {
		c_tmp[i] = _c[i];
		if (c_copyCont) c_copy_tmp[i] = _c_copy[i];
		if (c_copy1Cont) c_ccopy_tmp[i] = _c_ccopy[i];
	}

	if (c_copyCont) {
		_c_copy.resize(_n - 1);
		_c_copy = c_copy_tmp;
	}

	if (c_copy1Cont) {
		_c_ccopy.resize(_n - 1);
		_c_ccopy = c_ccopy_tmp;
	}

	_c.resize(_n - 1);
	_c = c_tmp;
}

void SimplexMethod::delete_from_b(const uint16_t& r) {
	auto it = _b.begin();
	for (uint16_t i = 0; i < _m; ++i)
		if (i < r) ++it;
	_b.erase(it);
}

void SimplexMethod::exclude_imit_variables(const uint16_t& s,
										   uint16_t& imit_basis) {
	uint16_t r = 0;
	bool isNull = true;
	// находим такое r, для которого q(r,s) = 1
	for (uint16_t i = 0; i < _m; ++i)
		if (_matrix.get_value(i, s) == 1) r = i;

	for (uint16_t j = 0; j < _n - imit_basis; ++j)
	{
		if (_matrix.get_value(r, j) != 0)
		{
			isNull = false;
			gauss_transform(j, r);
			// устанавливаем новую базисную переменной взамен старой
			_basis[r] = j;
			//удаляем s-ый столбец
			_matrix.delete_column(s);
			delete_from_c(s);
			--_n;
			del_from_imit_basis(s, imit_basis);
			rebasis(s);
			break;
		}
	}

	// а) строка r избыточна, нужно её удалить. 
	// Тогда столбец s будет полностью нулевым, поэтому его тоже удалим.
	if (isNull == true) {
		_matrix.delete_row(r);
		delete_from_b(r);
		--_m;
		_matrix.delete_column(s);
		delete_from_c(s);
		--_n;
		del_from_imit_basis(s, imit_basis);
		rebasis(s);
	}
}

void SimplexMethod::del_from_imit_basis(const uint16_t& s, uint16_t& imit_basis) {
	for (uint16_t i = 0; i < imit_basis; ++i)
		if (_imit_basis[i] == s)
			for (uint16_t j = i + 1; j < imit_basis; ++j)
				_imit_basis[j - 1] = _imit_basis[j];

	--imit_basis;
	std::vector<int16_t> imit_basis_tmp(imit_basis);
	for (uint16_t i = 0; i < imit_basis; ++i)
		imit_basis_tmp[i] = _imit_basis[i];

	_imit_basis.resize(imit_basis);
	_imit_basis = imit_basis_tmp;

	for (uint16_t i = 0; i < imit_basis; ++i)
		if (_imit_basis[i] > s) --_imit_basis[i];
}

void SimplexMethod::rebasis(const uint16_t& s) {
	for (uint16_t i = 0; i < _basis.size(); ++i)
		if (_basis[i] == s)
			for (uint16_t j = i + 1; j < _basis.size(); ++j)
				_basis[j - 1] = _basis[j];

	std::vector<int16_t> basis_tmp(_m);
	for (uint16_t i = 0; i < _m; ++i)
		basis_tmp[i] = _basis[i];

	_basis.resize(_m);
	_basis = basis_tmp;

	for (uint16_t i = 0; i < _m; ++i)
		if (_basis[i] > s) --_basis[i];
}

void SimplexMethod::find_solution() {
	uint16_t imit_basis_size = find_basis();
	double calc_res;
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

		int result_state = step();
		while (result_state != 0) {
			result_state = step();
			if (result_state == -1) _state = Solution::unlimited;
		}

		calc_res = calculate_result(false);

		if (calc_res > 0) { _result._state = Solution::inconsistent; return; }
		else if (calc_res < 0) { _result._state = Solution::inconsistent; return; }
		else {
			bool flag = false; // содержит ли базис искусственные переменные
			for (uint16_t i = 0; i < imit_basis_size; ++i) {
				if (std::find(_imit_basis.begin(), _imit_basis.end(), _imit_basis[i]) != 
					_imit_basis.end()) {
					flag = true;
					// исключаем исскуственный столбец i
					exclude_imit_variables(_imit_basis[i], imit_basis_size);
					--i;
				}
				flag = false;
			}

			// случай 2: база не содержит номера искусственных переменных, 
			// значит база допустима для исходной задачи
			// убираем искусственные столбцы
			if (!flag) delete_imit_column(imit_basis_size);
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

	int result_state_ = step();
	while (result_state_ != 0) {
		result_state_ = step();
		if (result_state_ == -1) _state = Solution::unlimited;
		calc_res = calculate_result(true);
	}
}