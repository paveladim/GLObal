#ifndef __SIMPLEX_METHOD__
#define __SIMPLEX_METHOD__

#include "synonymous_types.h"
#include "Matrix.h"

enum class Solution{exists, unlimited, inconsistent};

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
public:
	SimplexMethod() = delete;
	SimplexMethod(const Vec& c, const Vec& b, const Matrix& m);
	void find_solution();
};

#endif // __SIMPLEX_METHOD__