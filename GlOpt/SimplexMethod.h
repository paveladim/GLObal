#ifndef __SIMPLEX_METHOD__
#define __SIMPLEX_METHOD__

#include "synonymous_types.h"
#include "Matrix.h"

enum class Solution{exists, unlimited, inconsistent, notfound};

struct SimplexResult {
	Solution _state;
	// вектор решения
	Vec _x;
	// номера базисных переменных
	std::vector<uint16_t> _basis;
	// значение целевой функции
	double _value;
	// конструктор
	SimplexResult(const uint16_t& res_size, const uint16_t& basis_size) 
		: _value(std::numeric_limits<double>::max()), _state(Solution::notfound) {
		_x.resize(res_size);
		_basis.resize(basis_size);
	}

	void change_size(const uint16_t& res_size, const uint16_t& basis_size) {
		_x.resize(res_size);
		_basis.resize(basis_size);
	}
};

class SimplexMethod
{
	// состояние решения
	Solution _state;
	// количество переменных
	uint16_t _n;
	// количество ограничений
	uint16_t _m;
	// вектор целевой функции
	Vec _c;
	// копия вектора целевой функции
	Vec _c_copy;
	// копия вектора целевой функции задачи искусственного базиса
	Vec _c_ccopy;
	// правая часть
	Vec _b;
	// матрица ограничений
	Matrix _matrix;
	// номера базисных столбцов
	std::vector<int16_t> _basis;
	// номера столбцов искусственного базиса
	std::vector<int16_t> _imit_basis;
	// значения базисных переменных
	Vec x;
	// результат решения
	SimplexResult _result;
	// привести задачу к канонической
	void make_canonical();
	// отыскать базис
	uint16_t find_basis();
	// создать искусственный базис
	void create_imitation_basis(const uint16_t& basis_size,
							    const std::vector<bool>& places);
	// гауссово преобразование
	void gauss_transform(const uint16_t& s, const uint16_t& r);
	// найти ведущий столбец
	uint16_t get_leading_column();
	// найти ведущую строку
	uint16_t get_leading_row(const uint16_t&);
	// шаг симплекс метода
	int step();
	// посчитать решение
	int calculate_result(const bool& phase);
	void delete_imit_column(uint16_t& imit_basis);
	void delete_from_c(const uint16_t& s);
	void delete_from_b(const uint16_t& r);
	void exclude_imit_variables(const uint16_t& s, uint16_t& imit_basis);
	void del_from_imit_basis(const uint16_t& s, uint16_t& imit_basis);
	void rebasis(const uint16_t& s);
public:
	SimplexMethod() = delete;
	SimplexMethod(const Vec& c, const Vec& b, const Matrix& m);
	void find_solution();
};

#endif // __SIMPLEX_METHOD__