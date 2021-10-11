#ifndef __SIMPLEX_METHOD__
#define __SIMPLEX_METHOD__

#include "synonymous_types.h"
#include "Matrix.h"

enum class Solution{exists, unlimited, inconsistent, notfound};

struct SimplexResult {
	Solution _state;
	// ������ �������
	Vec _x;
	// ������ �������� ����������
	std::vector<uint16_t> _basis;
	// �������� ������� �������
	double _value;
	// �����������
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
	// ��������� �������
	Solution _state;
	// ���������� ����������
	uint16_t _n;
	// ���������� �����������
	uint16_t _m;
	// ������ ������� �������
	Vec _c;
	// ����� ������� ������� �������
	Vec _c_copy;
	// ����� ������� ������� ������� ������ �������������� ������
	Vec _c_ccopy;
	// ������ �����
	Vec _b;
	// ������� �����������
	Matrix _matrix;
	// ������ �������� ��������
	std::vector<int16_t> _basis;
	// ������ �������� �������������� ������
	std::vector<int16_t> _imit_basis;
	// �������� �������� ����������
	Vec x;
	// ��������� �������
	SimplexResult _result;
	// �������� ������ � ������������
	void make_canonical();
	// �������� �����
	uint16_t find_basis();
	// ������� ������������� �����
	void create_imitation_basis(const uint16_t& basis_size,
							    const std::vector<bool>& places);
	// �������� ��������������
	void gauss_transform(const uint16_t& s, const uint16_t& r);
	// ����� ������� �������
	uint16_t get_leading_column();
	// ����� ������� ������
	uint16_t get_leading_row(const uint16_t&);
	// ��� �������� ������
	int step();
	// ��������� �������
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