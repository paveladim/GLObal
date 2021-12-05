#include "SimplePMwithSM.h"

void SimplePMwithSM::calculate_localLipshConst(const uint& id_hyp) {
	// A, V, U, B
	static EncodedCoordinates points(4 * _dimension);
	static std::vector<double> incs(12);

	Hyperinterval& hyp1 = _intervals[id_hyp];
	Hyperinterval& hyp2 = _intervals[_generated_intervals - 1];
	Hyperinterval& hyp3 = _intervals[_generated_intervals - 2];

	uint pos_a = hyp1.get_coordA();
	uint pos_v = hyp2.get_coordA();
	uint pos_u = hyp2.get_coordB();
	uint pos_b = hyp3.get_coordB();

	// считываем точки
	for (uint i = 0; i < _dimension; ++i) {
		points[i] = _coords[pos_a + i];
		points[i + _dimension] = _coords[pos_v + i];
		points[i + 2 * _dimension] = _coords[pos_u + i];
		points[i + 3 * _dimension] = _coords[pos_b + i];
	}

	calculate_and_project(points, incs, hyp1.get_previous_axis());

	static std::vector<std::vector<double>> st(36);
	for (auto& row : st) row.resize(45, 0);
	static std::vector<double> b(36, 0);
	static std::vector<double> c(45, 0);
	c[8] = 1.0;

	generate_simplex_table(st, incs);

	SimplexMethod sm(c, b, st);
	sm.solve();
}

void SimplePMwithSM::calculate_and_project(const EncodedCoordinates& out,
									       std::vector<double>& incs,
										   const uint& axis) {
	static std::vector<double> increments(6 * _dimension);
	static std::vector<double> transit(_dimension);
	// высчитываем приращения
	for (uint i = 0; i < _dimension; ++i) {
		// 21
		increments[i] = (double)out[i + _dimension] - (double)out[i];
		// 31
		increments[i + _dimension] =
			(double)out[i + 2 * _dimension] - (double)out[i];
		// 41
		increments[i + 2 * _dimension] =
			(double)out[i + 3 * _dimension] - (double)out[i];
		// 32
		increments[i + 3 * _dimension] =
			(double)out[i + 2 * _dimension] - (double)out[i + _dimension];
		// 42
		increments[i + 4 * _dimension] =
			(double)out[i + 3 * _dimension] - (double)out[i + _dimension];
		// 43
		increments[i + 5 * _dimension] =
			(double)out[i + 3 * _dimension] - (double)out[i + 2 * _dimension];
	}

	std::vector<double> e2(_dimension, 0);
	e2[axis] = 1;
	std::vector<double> e1(_dimension, 0);
	double scalar = increments[2 * _dimension + axis];
	for (uint i = 0; i < _dimension; ++i)
		e1[i] = increments[i + 2 * _dimension] - scalar * e1[i];

	/* for (uint j = 0; j < 6; ++j) {
		for (uint i = 0; i < _dimension; ++i)
			transit[i] = increments[i + j * _dimension];

		incs[0 + 2 * j] = scalar_product(transit, e1);
		incs[1 + 2 * j] = scalar_product(transit, e2);
	} */


	for (uint i = 0; i < _dimension; ++i)
		transit[i] = increments[i];

	incs[0] = scalar_product(transit, e1);
	incs[1] = scalar_product(transit, e2);

	for (uint i = 0; i < _dimension; ++i)
		transit[i] = increments[i + _dimension];

	incs[2] = scalar_product(transit, e1);
	incs[3] = scalar_product(transit, e2);

	for (uint i = 0; i < _dimension; ++i)
		transit[i] = increments[i + 2 * _dimension];

	incs[4] = scalar_product(transit, e1);
	incs[5] = scalar_product(transit, e2);

	for (uint i = 0; i < _dimension; ++i)
		transit[i] = increments[i + 3 * _dimension];

	incs[6] = scalar_product(transit, e1);
	incs[7] = scalar_product(transit, e2);

	for (uint i = 0; i < _dimension; ++i)
		transit[i] = increments[i + 4 * _dimension];

	incs[8] = scalar_product(transit, e1);
	incs[9] = scalar_product(transit, e2);

	for (uint i = 0; i < _dimension; ++i)
		transit[i] = increments[i + 5 * _dimension];

	incs[10] = scalar_product(transit, e1);
	incs[11] = scalar_product(transit, e2);
}

double SimplePMwithSM::scalar_product(const std::vector<double>& a,
									  const std::vector<double>& b) {
	double result = 0.0;

	for (uint i = 0; i < _dimension; ++i)
		result += a[i] * b[i];
}

void SimplePMwithSM::generate_simplex_table(std::vector<std::vector<double>>& A,
										    const std::vector<double>& incs) {
	for (uint i = 0; i < 36; ++i) A[i][i + 9] = 1.0;

	static std::vector<double> norms(6);
	for (uint i = 0; i < 6; ++i) {
		norms[i] = incs[2 * i] * incs[2 * i];
		norms[i] += incs[2 * i + 1] * incs[2 * i + 1];
		norms[i] = sqrt(norms[i]);
	}

	for (uint i = 0; i < 6; ++i) {
		A[6 * i][8] = -norms[i];
		A[6 * i + 1][8] = -norms[i];
		A[6 * i + 2][8] = -0.5 * norms[i] * norms[i];
		A[6 * i + 3][8] = -0.5 * norms[i] * norms[i];
		A[6 * i + 4][8] = -0.5 * norms[i] * norms[i];
		A[6 * i + 5][8] = -0.5 * norms[i] * norms[i];
	}

	uint i = 0;
	uint k = 1;
	// check cycle
	while (i < 4)
		for (uint j = 0; j < 6; ++j) {
			A[6 * j][2 * i] = -incs[2 * i];
			A[6 * j][2 * i + 1] = -incs[2 * i + 1];
			A[6 * j][2 * i + 2 * k] = incs[2 * i];
			A[6 * j][2 * i + 2 * k + 1] = incs[2 * i + 1];

			A[6 * j + 1][2 * i] = incs[2 * i];
			A[6 * j + 1][2 * i + 1] = incs[2 * i + 1];
			A[6 * j + 1][2 * i + 2 * k] = -incs[2 * i];
			A[6 * j + 1][2 * i + 2 * k + 1] = -incs[2 * i + 1];

			A[6 * j + 2][2 * i] = -incs[2 * i];
			A[6 * j + 2][2 * i + 1] = -incs[2 * i];

			A[6 * j + 3][2 * i] = incs[2 * i];
			A[6 * j + 3][2 * i + 1] = incs[2 * i];

			A[6 * j + 4][2 * i + 2 * k] = incs[2 * i];
			A[6 * j + 4][2 * i + 2 * k + 1] = incs[2 * i];

			A[6 * j + 5][2 * i + 2 * k] = -incs[2 * i];
			A[6 * j + 5][2 * i + 2 * k + 1] = -incs[2 * i];

			++k;
			if (k + i > 3) {
				++i;
				k = 1;
			}
		}
}

void generate_right_part(std::vector<double>& b,
						 const uint& function,
						 const uint& eval_a,
						 const uint& eval_v,
						 const uint& eval_u,
						 const uint& eval_b) {

}

uint SimplePMwithSM::optimal_to_trisect() {
	return 0;
}

uint SimplePMwithSM::iterate(const uint& id_hyp) {
	trisect_interval(id_hyp);
	calculate_localLipshConst(id_hyp);
	calculate_characteristic(id_hyp);
	calculate_characteristic(_generated_intervals - 1);
	calculate_characteristic(_generated_intervals - 2);
	return optimal_to_trisect();
}