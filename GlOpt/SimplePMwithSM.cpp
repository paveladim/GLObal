#include <iostream>
#include "SimplePMwithSM.h"

SimplePMwithSM::SimplePMwithSM(const uint& dimension,
							   const uint& constraints,
							   Parameters& parameters,
							   Problem& problem) :
	SimplePM(dimension, constraints, parameters, problem),
	_points(4 * _dimension),
	_incs(6 * 2),
	_non_proj_incs(6 * _dimension) {}

void SimplePMwithSM::decode_and_save(const uint& pos, const uint& order) {
	for (size_t i = 0; i < _dimension; ++i)
		_transit1[i] = _coords[pos * _dimension + i];
	CoordinatesValues& t = _problem.decode_coordinates(_transit1);
	for (size_t i = 0; i < _dimension; ++i)
		_points[i + _dimension * order] = t[i];
}

void SimplePMwithSM::projection(const uint& order, 
								std::vector<double>& e1, 
								std::vector<double>& e2) {
	_incs[2 * order]     = scalar_product(_non_proj_incs.begin() + order * _dimension, 
									      e1.begin());
	_incs[2 * order + 1] = scalar_product(_non_proj_incs.begin() + order * _dimension, 
										  e2.begin());
}

void SimplePMwithSM::calculate_localLipshConst(const uint& id_hyp) {
	// A, V, U, B
	static CoordinatesValues points(4 * _dimension);

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

	decode_and_save(pos_a, 0);
	decode_and_save(pos_v, 1);
	decode_and_save(pos_u, 2);
	decode_and_save(pos_b, 3);

	calculate_and_project(hyp1.get_previous_axis());

	static std::vector<std::vector<double>> st(36);
	for (auto& row : st) row.resize(53, 0);
	static std::vector<double> b(36, 0);
	static std::vector<double> c(53, 0);
	c[0] = 1.0;

	generate_simplex_table(st, _incs);
	
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

void SimplePMwithSM::calculate_and_project(const uint& axis) {
	// высчитываем приращения
	for (uint i = 0; i < _dimension; ++i) {
		// 21
		_non_proj_incs[i] = _points[i + _dimension] - _points[i];
		// 31
		_non_proj_incs[i + _dimension] =
			_points[i + 2 * _dimension] - _points[i];
		// 41
		_non_proj_incs[i + 2 * _dimension] =
			_points[i + 3 * _dimension] - _points[i];
		// 32
		_non_proj_incs[i + 3 * _dimension] =
			_points[i + 2 * _dimension] - _points[i + _dimension];
		// 42
		_non_proj_incs[i + 4 * _dimension] =
			_points[i + 3 * _dimension] - _points[i + _dimension];
		// 43
		_non_proj_incs[i + 5 * _dimension] =
			_points[i + 3 * _dimension] - _points[i + 2 * _dimension];
	}

	std::vector<double> e2(_dimension, 0);
	e2[axis] = 1;
	std::vector<double> e1(_dimension, 0);
	double scalar = _non_proj_incs[2 * _dimension + axis];
	for (uint i = 0; i < _dimension; ++i)
		e1[i] = _non_proj_incs[i + 2 * _dimension] - scalar * e2[i];

	scalar = scalar_product(e1.begin(), e1.begin());
	for (uint i = 0; i < _dimension; ++i)
		e1[i] /= sqrt(scalar);

	projection(0, e1, e2);
	projection(1, e1, e2);
	projection(2, e1, e2);
	projection(3, e1, e2);
	projection(4, e1, e2);
	projection(5, e1, e2);
}

double SimplePMwithSM::scalar_product(const std::vector<double>::iterator& a,
									  const std::vector<double>::iterator& b) {
	double result = 0.0;

	for (size_t i = 0; i < _dimension; ++i)
		result += *(a + i) * *(b + i);

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