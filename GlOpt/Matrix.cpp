#include "Matrix.h"

Matrix::Matrix(const uint16_t& rows, const uint16_t& columns, const char& mode) :
	_rows(rows),
	_columns(columns) {
	// �������� ������ ��� �������
	_matrix.resize(_rows);
	for (auto& elem : _matrix) elem.resize(_columns);
	// � ����������� �� mode �������� ����� ������� �������
	switch (mode) {
	case 'n':
		for (uint16_t i = 0; i < _rows; ++i)
			for (uint16_t j = 0; j < _columns; ++j)
				_matrix[i][j] = 0.0;
		break;
	case 'e':
		if (rows == columns) {
			for (uint16_t i = 0; i < _rows; ++i)
				for (uint16_t j = 0; j < _columns; ++j)
					if (i == j) _matrix[i][j] = 1.0;
					else _matrix[i][j] = 0.0;
		}
		break;
	default:
		break;
	}
}

Matrix::Matrix(const Matrix& out) :
	_rows(out._rows),
	_columns(out._columns) {
	// �������� ������ ��� �������
	_matrix.resize(_rows);
	for (auto& elem : _matrix) elem.resize(_columns);

	for (uint16_t i = 0; i < _rows; ++i)
		for (uint16_t j = 0; j < _columns; ++j)
			_matrix[i][j] = out._matrix[i][j];
}

Matrix::Matrix(const std::vector<Vec>& matrix) :
	_rows(matrix.size()),
	_columns(matrix[0].size()) {
	_matrix.resize(_rows);
	for (auto& elem : _matrix) elem.resize(_columns);

	for (uint16_t i = 0; i < _rows; ++i)
		for (uint16_t j = 0; j < _columns; ++j)
			_matrix[i][j] = matrix[i][j];
}

Matrix& Matrix::operator=(const Matrix& out) {
	if (this != &out) {
		_rows = out._rows;
		_columns = out._columns;

		_matrix.clear();
		_matrix.resize(_rows);
		for (auto& elem : _matrix) elem.resize(_columns);

		for (uint16_t i = 0; i < _rows; ++i)
			for (uint16_t j = 0; j < _columns; ++j)
				_matrix[i][j] = out._matrix[i][j];
	}

	return *this;
}

void Matrix::delete_row(const uint16_t& row) {
	auto it = _matrix.begin();
	for (uint16_t i = 0; i < _rows; ++i)
		if (i < row) ++it;

	_matrix.erase(it);
	--_rows;
}

void Matrix::delete_column(const uint16_t& col) {
	for (uint16_t i = 0; i < _rows; ++i) {
		auto it = _matrix[i].begin();
		for (uint16_t j = 0; j < _columns; ++j)
			if (j < col) ++it;
		_matrix[i].erase(it);
	}

	--_columns;
}

int Matrix::compare_lex(const uint16_t& i, const uint16_t& j) {
	int result = 0;
	// �������� ������������ ��������
	Vec diff(_columns);
	// ������� �������� ��� ���� ������ v1 ��� v2 �������, �� ��� ��� ������� ����,
	// ��� ��� ������������
	bool flagV1 = true;
	bool flagV2 = true;
	for (uint16_t k = 0; k < _rows; ++k)
	{
		if (_matrix[k][i] != 0) flagV1 = false; //������ ����� �� �������
		if (_matrix[k][j] != 0) flagV2 = false; //������ ����� �� �������
		diff[k] = _matrix[k][i] - _matrix[k][j];
	}
	if (flagV1 && flagV2) return 0; //��� �������
	if (flagV1) return 1;
	if (flagV2) return -1;

	bool f = true;
	// �������� �� ��� ���������� � ������� ����� ������ ��������� ��� > 0
	for (uint16_t k = 0; ((k < _rows) && (f)); ++k)
		if (diff[k] != 0) {
			if (diff[k] > 0) result = 1;
			if (diff[k] < 0) result = -1;
			f = false;
		}

	return result;
}

void Matrix::attach_matrix(const Matrix& out) {
	if (_rows != out._rows) return;
	std::vector<Vec> tmp(_matrix);
	_matrix.clear();
	_matrix.resize(_rows + out._rows);
	for (auto& elem : _matrix) elem.resize(_columns + out._columns);

	_rows = _rows + out._rows;
	_columns = _columns + out._columns;

	for (uint16_t i = 0; i < out._rows; ++i)
		for (uint16_t j = 0; j < out._columns; ++j)
			_matrix[i][j] = out.get_value(i, j);

	for (uint16_t i = 0; i < _rows; ++i)
		for (uint16_t j = out._columns; j < _columns; ++j)
			_matrix[i][j] = tmp[i][j - out._columns];
}
