#include <iostream>
#include "SimplePMwithSM.h"

SimplePMwithSM::SimplePMwithSM(const uint& dimension,
							   const uint& constraints,
							   Parameters& parameters,
							   Problem& problem) :
	SimplePM(dimension, constraints, parameters, problem),
	_point(4 * _dimension),
	_incs(6 * 2),
	_non_proj_incs(6 * _dimension),
	_st(36),
	_b(36),
	_c(53) {
	for (auto& elem : _st) elem.resize(53);
	_c[0] = 1.0;
}

void SimplePMwithSM::decode_and_save(const uint& pos, const uint& order) {
	for (size_t i = 0; i < _dimension; ++i)
		_transit1[i] = _coords[pos * _dimension + i];
	CoordinatesValues& t = _problem.decode_coordinates(_transit1);
	for (size_t i = 0; i < _dimension; ++i)
		_point[i + _dimension * order] = t[i];
}

void SimplePMwithSM::projection(const uint& order, 
								std::vector<double>& e1, 
								std::vector<double>& e2) {
	_incs[2 * order]     = scalar_product(_non_proj_incs.begin() + order * _dimension, 
									      e2.begin());
	_incs[2 * order + 1] = scalar_product(_non_proj_incs.begin() + order * _dimension, 
										  e1.begin());
}

void SimplePMwithSM::calculate_localLipshConst(const uint& id_hyp) {
	Hyperinterval& hyp1 = _intervals[id_hyp];
	Hyperinterval& hyp2 = _intervals[_generated_intervals - 2];
	Hyperinterval& hyp3 = _intervals[_generated_intervals - 1];

	uint pos_a = hyp1.get_idA();
	uint pos_v = hyp2.get_idA();
	uint pos_u = hyp2.get_idB();
	uint pos_b = hyp3.get_idB();

	decode_and_save(pos_a, 0);
	decode_and_save(pos_v, 1);
	decode_and_save(pos_u, 2);
	decode_and_save(pos_b, 3);

	calculate_and_project(hyp1.get_previous_axis());
	generate_simplex_table();
	
	for (uint i = 0; i < _constraints + 1; ++i) {
		generate_right_part(i, pos_a * (_constraints + 1),
							   pos_v * (_constraints + 1),
							   pos_u * (_constraints + 1),
							   pos_b * (_constraints + 1));

		SimplexMethod sm(_c, _b, _st);
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
		_non_proj_incs[i] = _point[i + _dimension] - _point[i];
		// 31
		_non_proj_incs[i + _dimension] =
			_point[i + 2 * _dimension] - _point[i];
		// 41
		_non_proj_incs[i + 2 * _dimension] =
			_point[i + 3 * _dimension] - _point[i];
		// 32
		_non_proj_incs[i + 3 * _dimension] =
			_point[i + 2 * _dimension] - _point[i + _dimension];
		// 42
		_non_proj_incs[i + 4 * _dimension] =
			_point[i + 3 * _dimension] - _point[i + _dimension];
		// 43
		_non_proj_incs[i + 5 * _dimension] =
			_point[i + 3 * _dimension] - _point[i + 2 * _dimension];
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

size_t SimplePMwithSM::index(const size_t& ind) {
	if (ind == 1) return 1;
	if (ind == 2) return 5;
	if (ind == 3) return 9;
	if (ind == 4) return 13;

	return -1;
}

void SimplePMwithSM::print_table() {
	for (size_t i = 0; i < 17; ++i)
		printf("%10.3lf", _c[i]);

	printf("\n");

	for (size_t j = 0; j < 36; ++j) {
		for (size_t i = 0; i < 17; ++i) {
			printf("%10.3lf ", _st[j][i]);
		}

		printf("%10.3lf ", _b[j]);
		printf("\n");
	}
}

void SimplePMwithSM::generate_simplex_table() {
	size_t k{ 0 };
	for (size_t i = 17; i < 53; ++i) {
		_st[k][i] = 1.0;
		++k;
	}

	for (size_t i = 0; i < 6; ++i) {
		_st[6 * i + 0][0] = -(_incs[2 * i + 0] * _incs[2 * i + 0] + _incs[2 * i + 1] * _incs[2 * i + 1]);
		_st[6 * i + 1][0] = -(_incs[2 * i + 0] * _incs[2 * i + 0] + _incs[2 * i + 1] * _incs[2 * i + 1]);
		_st[6 * i + 2][0] = -0.5 * (_incs[2 * i + 0] * _incs[2 * i + 0] + _incs[2 * i + 1] * _incs[2 * i + 1]);
		_st[6 * i + 3][0] = -0.5 * (_incs[2 * i + 0] * _incs[2 * i + 0] + _incs[2 * i + 1] * _incs[2 * i + 1]);
		_st[6 * i + 4][0] = -0.5 * (_incs[2 * i + 0] * _incs[2 * i + 0] + _incs[2 * i + 1] * _incs[2 * i + 1]);
		_st[6 * i + 5][0] = -0.5 * (_incs[2 * i + 0] * _incs[2 * i + 0] + _incs[2 * i + 1] * _incs[2 * i + 1]);
	}

	size_t i{ 1 };
	size_t j{ 2 };

	for (size_t k = 0; k < 6; ++k) {
		_st[6 * k + 0][index(i) + 0] = -_incs[2 * k + 0];
		_st[6 * k + 0][index(i) + 1] = _incs[2 * k + 0];
		_st[6 * k + 0][index(i) + 2] = -_incs[2 * k + 1];
		_st[6 * k + 0][index(i) + 3] = _incs[2 * k + 1];

		_st[6 * k + 0][index(j) + 0] = _incs[2 * k + 0];
		_st[6 * k + 0][index(j) + 1] = -_incs[2 * k + 0];
		_st[6 * k + 0][index(j) + 2] = _incs[2 * k + 1];
		_st[6 * k + 0][index(j) + 3] = -_incs[2 * k + 1];

		_st[6 * k + 1][index(i) + 0] = _incs[2 * k + 0];
		_st[6 * k + 1][index(i) + 1] = -_incs[2 * k + 0];
		_st[6 * k + 1][index(i) + 2] = _incs[2 * k + 1];
		_st[6 * k + 1][index(i) + 3] = -_incs[2 * k + 1];

		_st[6 * k + 1][index(j) + 0] = -_incs[2 * k + 0];
		_st[6 * k + 1][index(j) + 1] = _incs[2 * k + 0];
		_st[6 * k + 1][index(j) + 2] = -_incs[2 * k + 1];
		_st[6 * k + 1][index(j) + 3] = _incs[2 * k + 1];

		_st[6 * k + 2][index(i) + 0] = -_incs[2 * k + 0];
		_st[6 * k + 2][index(i) + 1] = _incs[2 * k + 0];
		_st[6 * k + 2][index(i) + 2] = -_incs[2 * k + 1];
		_st[6 * k + 2][index(i) + 3] = _incs[2 * k + 1];

		_st[6 * k + 3][index(i) + 0] = _incs[2 * k + 0];
		_st[6 * k + 3][index(i) + 1] = -_incs[2 * k + 0];
		_st[6 * k + 3][index(i) + 2] = _incs[2 * k + 1];
		_st[6 * k + 3][index(i) + 3] = -_incs[2 * k + 1];

		_st[6 * k + 4][index(j) + 0] = -_incs[2 * k + 0];
		_st[6 * k + 4][index(j) + 1] = _incs[2 * k + 0];
		_st[6 * k + 4][index(j) + 2] = -_incs[2 * k + 1];
		_st[6 * k + 4][index(j) + 3] = _incs[2 * k + 1];
		
		_st[6 * k + 5][index(j) + 0] = _incs[2 * k + 0];
		_st[6 * k + 5][index(j) + 1] = -_incs[2 * k + 0];
		_st[6 * k + 5][index(j) + 2] = _incs[2 * k + 1];
		_st[6 * k + 5][index(j) + 3] = -_incs[2 * k + 1];

		++j;
		if (j > 4) {
			++i;
			j = i + 1;
		}
	}
}

void SimplePMwithSM::generate_right_part(const uint& function,
										 const uint& eval_a,
										 const uint& eval_v,
										 const uint& eval_u,
										 const uint& eval_b) {
	double fa = _evaluations[eval_a + function];
	double fv = _evaluations[eval_v + function];
	double fu = _evaluations[eval_u + function];
	double fb = _evaluations[eval_b + function];

	_b[2] = fa - fv;
	_b[3] = fv - fa;
	_b[4] = fa - fv;
	_b[5] = fv - fa;

	_b[8] = fa - fu;
	_b[9] = fu - fa;
	_b[10] = fa - fu;
	_b[11] = fu - fa;

	_b[14] = fa - fb;
	_b[15] = fb - fa;
	_b[16] = fa - fb;
	_b[17] = fb - fa;

	_b[20] = fv - fu;
	_b[21] = fu - fv;
	_b[22] = fv - fu;
	_b[23] = fu - fv;

	_b[26] = fv - fb;
	_b[27] = fb - fv;
	_b[28] = fv - fb;
	_b[29] = fb - fv;

	_b[32] = fu - fb;
	_b[33] = fb - fu;
	_b[34] = fu - fb;
	_b[35] = fb - fu;
}

void SimplePMwithSM::calculate_characteristic(const uint& id_hyp) {
	if (_intervals[id_hyp].get_divisions() == 0) return;
	Hyperinterval& hyp = _intervals[id_hyp];

	double mixed_LipshEval = mixedLipEval(hyp, 0);
	balance(mixed_LipshEval);
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