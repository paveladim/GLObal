#ifndef __MATRIX_H__
#define __MATRIX_H__

#include "synonymous_types.h"

class Matrix
{
	uint16_t _rows;
	uint16_t _columns;
	std::vector<Vec> _matrix;
public:
	Matrix() : _rows(0), _columns(0) {}
	Matrix(const uint16_t& rows, const uint16_t& columns, const char& mode);
	Matrix(const Matrix& out);

	Matrix& operator=(const Matrix& out);

	void delete_row(const uint16_t& row);
	double get_value(const uint16_t& i, const uint16_t& j) const { return _matrix[i][j]; }
	void set_value(const uint16_t& i, const uint16_t& j, const double& val) {
		_matrix[i][j] = val;
	}
	void add_elem(const uint16_t& i, const uint16_t& j, const double& elem) {
		_matrix[i][j] += elem;
	}
	void mul_elem(const uint16_t& i, const uint16_t& j, const double& mul) {
		_matrix[i][j] *= mul;
	}

	void delete_column(const uint16_t& col);
	int compare_lex(const uint16_t& i, const uint16_t& j);

	void attach_matrix(const Matrix& out);
};

#endif // __MATRIX_H__
