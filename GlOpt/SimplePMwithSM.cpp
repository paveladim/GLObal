#include <iostream>
#include "SimplePMwithSM.h"

SimplePMwithSM::SimplePMwithSM(const uint& dimension,
							   const uint& constraints,
							   Parameters& parameters,
							   Problem& problem) :
	SimplePM(dimension, constraints, parameters, problem) {}

void SimplePMwithSM::calculate_localLipshConst(const uint& id_hyp) {
	// A, V, U, B
	static CoordinatesValues points(4 * _dimension);
	static std::vector<double> incs(12);

	Hyperinterval& hyp1 = _intervals[id_hyp];
	Hyperinterval& hyp2 = _intervals[_generated_intervals - 2];
	Hyperinterval& hyp3 = _intervals[_generated_intervals - 1];

	uint pos_a = hyp1.get_idA();
	uint pos_v = hyp2.get_idA();
	uint pos_u = hyp2.get_idB();
	uint pos_b = hyp3.get_idB();

	uint aa = _coords[pos_a * _dimension + hyp1.get_previous_axis()];
	uint vv = _coords[pos_v * _dimension + hyp2.get_previous_axis()];
	uint uu = _coords[pos_u * _dimension + hyp2.get_previous_axis()];
	uint bb = _coords[pos_b * _dimension + hyp3.get_previous_axis()];

	for (uint i = 0; i < _dimension; ++i)
		_transit1[i] = _coords[pos_a * _dimension + i];
	CoordinatesValues& t = _problem.decode_coordinates(_transit1);
	for (uint i = 0; i < _dimension; ++i)
		points[i] = t[i];

	for (uint i = 0; i < _dimension; ++i)
		_transit1[i] = _coords[pos_v * _dimension + i];
	t = _problem.decode_coordinates(_transit1);
	for (uint i = 0; i < _dimension; ++i)
		points[i + _dimension] = t[i];

	for (uint i = 0; i < _dimension; ++i)
		_transit1[i] = _coords[pos_u * _dimension + i];
	t = _problem.decode_coordinates(_transit1);
	for (uint i = 0; i < _dimension; ++i)
		points[i + 2 * _dimension] = t[i];

	for (uint i = 0; i < _dimension; ++i)
		_transit1[i] = _coords[pos_b * _dimension + i];
	t = _problem.decode_coordinates(_transit1);
	for (uint i = 0; i < _dimension; ++i)
		points[i + 3 * _dimension] = t[i];

	calculate_and_project(points, incs, hyp1.get_previous_axis());

	static std::vector<std::vector<double>> st(36);
	for (auto& row : st) row.resize(45, 0);
	static std::vector<double> b(36, 0);
	static std::vector<double> c(45, 0);
	c[8] = 1.0;

	generate_simplex_table(st, incs);
	
	for (uint i = 0; i < _constraints + 1; ++i) {
		generate_right_part(b, i, pos_a * (_constraints + 1),
								  pos_u * (_constraints + 1),
								  pos_v * (_constraints + 1),
								  pos_b * (_constraints + 1));

		SimplexMethod sm(c, b, st);
		sm.solve();
		_localLipshEval[i] = sm.get_solution() / ((double)MAX_POWER_THREE * (double)MAX_POWER_THREE);
	}

	hyp1.update_localLipQueues(_localLipshEval, 
							   _parameters._delta / 
							   ((double)MAX_POWER_THREE * (double)MAX_POWER_THREE));
	hyp2.update_localLipQueues(_localLipshEval, 
							   _parameters._delta / 
							   ((double)MAX_POWER_THREE * (double)MAX_POWER_THREE));
	hyp3.update_localLipQueues(_localLipshEval, 
							   _parameters._delta / 
							   ((double)MAX_POWER_THREE * (double)MAX_POWER_THREE));
}

void SimplePMwithSM::calculate_and_project(const CoordinatesValues& out,
									       std::vector<double>& incs,
										   const uint& axis) {
	static std::vector<double> increments(6 * _dimension);
	static std::vector<double> transit(_dimension);
	// высчитываем приращения
	for (uint i = 0; i < _dimension; ++i) {
		// 21
		increments[i] = out[i + _dimension] - out[i];
		// 31
		increments[i + _dimension] =
			out[i + 2 * _dimension] - out[i];
		// 41
		increments[i + 2 * _dimension] =
			out[i + 3 * _dimension] - out[i];
		// 32
		increments[i + 3 * _dimension] =
			out[i + 2 * _dimension] - out[i + _dimension];
		// 42
		increments[i + 4 * _dimension] =
			out[i + 3 * _dimension] - out[i + _dimension];
		// 43
		increments[i + 5 * _dimension] =
			out[i + 3 * _dimension] - out[i + 2 * _dimension];
	}

	std::vector<double> e2(_dimension, 0);
	e2[axis] = 1;
	std::vector<double> e1(_dimension, 0);
	double scalar = increments[2 * _dimension + axis];
	for (uint i = 0; i < _dimension; ++i)
		e1[i] = increments[i + 2 * _dimension] - scalar * e2[i];

	scalar = scalar_product(e1, e1);
	for (uint i = 0; i < _dimension; ++i)
		e1[i] /= sqrt(scalar);

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

	return result;
}

void SimplePMwithSM::generate_simplex_table(std::vector<std::vector<double>>& A,
										    const std::vector<double>& incs) {
	A[0][9] = 1.0;
	for (uint i = 1; i < 36; ++i) A[i][i + 9] = -1 * A[i - 1][(i - 1) + 9];

	static std::vector<double> norms(6);
	for (uint i = 0; i < 6; ++i) {
		norms[i] = incs[2 * i] * incs[2 * i];
		norms[i] += incs[2 * i + 1] * incs[2 * i + 1];
	}

	for (uint i = 0; i < 6; ++i) {
		A[6 * i][8] = -norms[i];
		A[6 * i + 1][8] = norms[i];
		A[6 * i + 2][8] = -0.5 * norms[i];
		A[6 * i + 3][8] = 0.5 * norms[i];
		A[6 * i + 4][8] = -0.5 * norms[i];
		A[6 * i + 5][8] = 0.5 * norms[i];
	}

	uint i = 0;
	uint k = 1;
	for (uint j = 0; j < 6; ++j) {
		A[6 * j][2 * i] = -incs[2 * j];
		A[6 * j][2 * i + 1] = -incs[2 * j + 1];
		A[6 * j][2 * k] = incs[2 * j];
		A[6 * j][2 * k + 1] = incs[2 * j + 1];

		A[6 * j + 1][2 * i] = -incs[2 * j];
		A[6 * j + 1][2 * i + 1] = -incs[2 * j + 1];
		A[6 * j + 1][2 * k] = incs[2 * j];
		A[6 * j + 1][2 * k + 1] = incs[2 * j + 1];

		A[6 * j + 2][2 * i] = -incs[2 * j];
		A[6 * j + 2][2 * i + 1] = -incs[2 * j + 1];

		A[6 * j + 3][2 * i] = -incs[2 * j];
		A[6 * j + 3][2 * i + 1] = -incs[2 * j + 1];

		A[6 * j + 4][2 * k] = incs[2 * j];
		A[6 * j + 4][2 * k + 1] = incs[2 * j + 1];

		A[6 * j + 5][2 * k] = incs[2 * j];
		A[6 * j + 5][2 * k + 1] = incs[2 * j + 1];

		++k;
		if (k * 2 == 8) {
			++i;
			k = i + 1;
		}
	}
}

void SimplePMwithSM::generate_right_part(std::vector<double>& b,
						 const uint& function,
						 const uint& eval_a,
						 const uint& eval_v,
						 const uint& eval_u,
						 const uint& eval_b) {
	double fa = _evaluations[eval_a + function];
	double fv = _evaluations[eval_v + function];
	double fu = _evaluations[eval_u + function];
	double fb = _evaluations[eval_b + function];

	b[2] = fa - fv;
	b[3] = fa - fv;
	b[4] = fv - fa;
	b[5] = fv - fa;

	b[8] = fa - fu;
	b[9] = fa - fu;
	b[10] = fu - fa;
	b[11] = fu - fa;

	b[14] = fa - fb;
	b[15] = fa - fb;
	b[16] = fb - fa;
	b[17] = fb - fa;

	b[20] = fv - fu;
	b[21] = fv - fu;
	b[22] = fu - fv;
	b[23] = fu - fv;

	b[26] = fv - fb;
	b[27] = fv - fb;
	b[28] = fb - fv;
	b[29] = fb - fv;

	b[32] = fu - fb;
	b[33] = fu - fb;
	b[34] = fb - fu;
	b[35] = fb - fu;
}

void SimplePMwithSM::calculate_characteristic(const uint& id_hyp) {
	if (_intervals[id_hyp].get_divisions() == 0) return;
	Hyperinterval& hyp = _intervals[id_hyp];

	double mixed_LipshEval = mixedLipEval(hyp, 0);
	double t_min = 0.5 * (_evaluations[hyp.get_evalA()] -
				   _evaluations[hyp.get_evalB()]);
	t_min = t_min / mixed_LipshEval;

	uint index = (hyp.get_divisions() - 1) / _dimension;
	uint j = (hyp.get_divisions() - 1) % _dimension + 1;
	double h1 = (double)HYPER_INTERVAL_SIDE_LENGTHS[20 - index];
	double h2 = (double)HYPER_INTERVAL_SIDE_LENGTHS[20 - index - 1];
	double e = 0.5 * sqrt((_dimension - j) * h1 * h1 + j * h2 * h2);

	t_min = t_min / e;
	double left = -e;
	double right = e;
	give_borders(left, right, hyp);

	if ((left == e) && (right == -e)) {
		hyp.set_charact(std::numeric_limits<double>::max());
	}
	else {
		if ((t_min > left) && (t_min < right)) {
			double charact = -_evaluations[hyp.get_evalA()] * (t_min - e);
			charact = charact + _evaluations[hyp.get_evalB()] * (t_min + e);
			charact = 0.5 * charact / e;
			charact = charact + 0.5 * mixed_LipshEval * (t_min * t_min - e * e);
			hyp.set_charact(charact);
		}
		else {
			double charact1 = -_evaluations[hyp.get_evalA()] * (left - e);
			charact1 = charact1 + _evaluations[hyp.get_evalB()] * (left + e);
			charact1 = 0.5 * charact1 / e;
			charact1 = charact1 + 0.5 * mixed_LipshEval * (left * left - e * e);

			double charact2 = -_evaluations[hyp.get_evalA()] * (right - e);
			charact2 = charact2 + _evaluations[hyp.get_evalB()] * (right + e);
			charact2 = 0.5 * charact2 / e;
			charact2 = charact2 + 0.5 * mixed_LipshEval * (right * right - e * e);

			hyp.set_charact(std::min(charact1, charact2));
		}
	}
}