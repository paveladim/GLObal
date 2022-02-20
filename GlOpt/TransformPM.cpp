#include "TransformPM.h"

TransformPM::TransformPM(const uint& dimension,
						 const uint& constraints,
						 Parameters& parameters,
						 Problem& problem) :
	DivideByThree(dimension, constraints, parameters, problem),
	_areAllCharInfty(false),
	_doesGlobalChange(false),
	_localLipshEval(constraints + 1),
	_globalLipshEval(constraints + 1) {}

void TransformPM::decode_and_save(const uint& pos, const uint& order) {
	for (size_t i = 0; i < _dimension; ++i)
		_transit1[i] = _coords[pos * _dimension + i];
	CoordinatesValues& t = _problem.decode_coordinates(_transit1);
	for (size_t i = 0; i < _dimension; ++i)
		_points[i + _dimension * order] = t[i];
}

void TransformPM::projection(const uint& order,
	std::vector<double>& e1,
	std::vector<double>& e2) {
	_incs[2 * order] = scalar_product(_non_proj_incs.begin() + order * _dimension,
									  e2.begin());
	_incs[2 * order + 1] = scalar_product(_non_proj_incs.begin() + order * _dimension,
										  e1.begin());
}

void TransformPM::calculate_and_project(const uint& axis) {
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

double TransformPM::scalar_product(const std::vector<double>::iterator& a,
	const std::vector<double>::iterator& b) {
	double result = 0.0;

	for (size_t i = 0; i < _dimension; ++i)
		result += *(a + i) * *(b + i);

	return result;
}

void TransformPM::generate_simplex_table() {
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

void TransformPM::generate_right_part(const uint& function,
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

size_t TransformPM::index(const size_t& ind) {
	if (ind == 1) return 1;
	if (ind == 2) return 5;
	if (ind == 3) return 9;
	if (ind == 4) return 13;

	return -1;
}

double TransformPM::calculate_residual(const double& t, const uint& id_hyp) {
	Hyperinterval& hyp = _intervals[id_hyp];

	uint index = (hyp.get_divisions() - 1) / _dimension;
	uint j = (hyp.get_divisions() - 1) % _dimension + 1;
	double h1 = static_cast<double>(HYPER_INTERVAL_SIDE_LENGTHS[20 - index]);
	double h2 = static_cast<double>(HYPER_INTERVAL_SIDE_LENGTHS[20 - index - 1]);
	double e = 0.5 * sqrt((_dimension - j) * h1 * h1 + j * h2 * h2);

	size_t idx_A = hyp.get_evalA();
	size_t idx_B = hyp.get_evalB();

	double result = -0.5 * _evaluations[idx_A] * (t - e) / e;
	result += 0.5 * _evaluations[idx_B] * (t + e) / e;
	result += 0.5 * mixedLipEval(hyp, 0) * (t * t - e * e);

	double candidate = 0.0;
	for (size_t i = 1; i < _constraints; ++i) {
		candidate = -0.5 * _evaluations[idx_A + i] * (t - e) / e;
		candidate += 0.5 * _evaluations[idx_B + i] * (t + e) / e;
		candidate += 0.5 * mixedLipEval(hyp, i) * (t * t - e * e);

		if (result < candidate) result = candidate;
	}

	return candidate;
}

void TransformPM::calculate_localLipshConst(const uint& id_hyp) {
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

void TransformPM::calculate_characteristic(const uint& id_hyp) {
	if (_intervals[id_hyp].get_divisions() == 0) return;
}

FunctionValue TransformPM::give_residual(const FunctionsValues& evals) {
	double max_value = evals[0] - _current_minimum;

	for (size_t i = 1; i < _constraints; ++i)
		if (evals[i] > max_value) max_value = evals[i];

	return max_value;
}

void TransformPM::calculate_globalLipshConst(const uint& id_hyp) {
	Hyperinterval& hyp1 = _intervals[id_hyp];
	Hyperinterval& hyp2 = _intervals[_generated_intervals - 2];
	Hyperinterval& hyp3 = _intervals[_generated_intervals - 1];

	for (uint i = 0; i < _constraints + 1; ++i) {
		if (_globalLipshEval[i] < hyp1.get_maxLipshEvaluations()[i]) {
			_doesGlobalChange = true;
			_globalLipshEval[i] = hyp1.get_maxLipshEvaluations()[i];
		}

		if (_globalLipshEval[i] < hyp2.get_maxLipshEvaluations()[i]) {
			_doesGlobalChange = true;
			_globalLipshEval[i] = hyp2.get_maxLipshEvaluations()[i];
		}

		if (_globalLipshEval[i] < hyp3.get_maxLipshEvaluations()[i]) {
			_doesGlobalChange = true;
			_globalLipshEval[i] = hyp3.get_maxLipshEvaluations()[i];
		}
	}

	if (_doesGlobalChange) {
		update_all_charact();
		_doesGlobalChange = false;
	}
	else {
		calculate_characteristic(id_hyp);
		calculate_characteristic(_generated_intervals - 2);
		calculate_characteristic(_generated_intervals - 1);
	}
}

double TransformPM::mixedLipEval(const Hyperinterval& hyp, const uint& i) {
	double mixed_LipshitzEval = 0.0;
	double ratio = hyp.get_diagonal() / _parameters._criticalSize;

	if (i == 0) {
		if (hyp.get_diagonal() < _parameters._criticalSize) {
			mixed_LipshitzEval = ratio * _parameters._gainGlobalObj *
				_globalLipshEval[i];
			mixed_LipshitzEval += (1 - ratio) * _parameters._gainLocalObj *
				hyp.get_maxLipshEvaluations()[i];
		}
		else mixed_LipshitzEval = _parameters._gainGlobalObj *
			_globalLipshEval[i];
	}
	else {
		if (hyp.get_diagonal() < _parameters._criticalSize) {
			mixed_LipshitzEval = ratio * _parameters._gainGlobalCst *
				_globalLipshEval[i];
			mixed_LipshitzEval += (1 - ratio) * _parameters._gainLocalCst *
				hyp.get_maxLipshEvaluations()[i];
		}
		else mixed_LipshitzEval = _parameters._gainGlobalCst *
			_globalLipshEval[i];
	}

	return mixed_LipshitzEval;
}

void TransformPM::update_all_charact() {
	for (uint i = 0; i < _generated_intervals; ++i)
		calculate_characteristic(i);
}

uint TransformPM::optimal_to_trisect() {
	uint id_optimal_hyp = 0;

	if (_areAllCharInfty) {
		double optimal_charact = _intervals[id_optimal_hyp].get_diagonal();
		double current_charact = 0.0;

		for (uint id_hyp = 1; id_hyp < _generated_intervals; ++id_hyp) {
			current_charact = _intervals[id_hyp].get_diagonal();
			if (std::abs(current_charact - optimal_charact) < std::numeric_limits<double>::epsilon());
			else if (optimal_charact < current_charact) {
				optimal_charact = current_charact;
				id_optimal_hyp = id_hyp;
			}
		}
	}
	else {
		double optimal_charact = _intervals[id_optimal_hyp].get_charact();
		double current_charact = 0.0;

		for (uint id_hyp = 1; id_hyp < _generated_intervals; ++id_hyp) {
			current_charact = _intervals[id_hyp].get_charact();
			if (std::abs(current_charact - optimal_charact) < std::numeric_limits<double>::epsilon());
			else if (optimal_charact > current_charact) {
				optimal_charact = current_charact;
				id_optimal_hyp = id_hyp;
			}
		}
	}

	return id_optimal_hyp;
}

uint TransformPM::iterate(const uint& id_hyp) {
	trisect_interval(id_hyp);
	calculate_localLipshConst(id_hyp);
	calculate_globalLipshConst(id_hyp);
	return optimal_to_trisect();
}