#include "TransformPM.h"

TransformPM::TransformPM(const uint& dimension,
						 const uint& constraints,
						 Parameters& parameters,
						 Problem& problem) :
	DivideByThree(dimension, constraints, parameters, problem),
	_areAllCharInfty(false),
	_doesGlobalChange(false),
	_localLipshEval(constraints + 1),
	_globalLipshEval(constraints + 1),
	_residual_minimum(std::numeric_limits<double>::max()) {}

void TransformPM::update_minimum(const FunctionsValues& evals,
								 const uint& idp) {
	if (evals[0] < _current_minimum) {
		bool flag = true;
		for (uint i = 1; i < _constraints + 1; ++i)
			if (evals[i] > 0) flag = false;

		if (flag) {
			_current_minimum = evals[0];
			_id_minimum = idp;
			update_all_charact();
		}
	}
}

double TransformPM::golden_ratio(double a, 
								 double b, 
								 const double& e, 
								 const uint& id_hyp) {
	if (b < a) std::swap(a, b);
	double x1 = b - phi * (b - a);
	double x2 = a + phi * (b - a);

	while (0.5 * (b - a) > 1e-10)
		if (calculate_residual(x1, e, id_hyp) > calculate_residual(x2, e, id_hyp)) {
			a = x1;
			x1 = b - (x2 - a);
		}
		else {
			b = x2;
			x2 = a + (b - x1);
		}

	return calculate_residual(0.5 * (b - a), e, id_hyp);
}

double TransformPM::calculate_residual(const double& t, const double& e, 
									   const uint& id_hyp) {
	Hyperinterval& hyp = _intervals[id_hyp];

	size_t idx_A = hyp.get_evalA();
	size_t idx_B = hyp.get_evalB();
	double mixed_LipshEval = 0.0;

	double result = -0.5 * (_evaluations[idx_A] - _current_minimum) * (t - e) / e;
	result += 0.5 * (_evaluations[idx_B] - _current_minimum) * (t + e) / e;
	mixed_LipshEval = mixedLipEval(hyp, 0);
	balance(mixed_LipshEval);
	result += mixed_LipshEval * (t * t - e * e);

	double candidate = 0.0;
	for (size_t i = 1; i < _constraints + 1; ++i) {
		candidate = -0.5 * _evaluations[idx_A + i] * (t - e) / e;
		candidate += 0.5 * _evaluations[idx_B + i] * (t + e) / e;
		mixed_LipshEval = mixedLipEval(hyp, i);
		balance(mixed_LipshEval);
		candidate += mixed_LipshEval * (t * t - e * e);

		if (result < candidate) result = candidate;
	}

	return result;
}

void TransformPM::calculate_characteristic(const uint& id_hyp) {
	if (_intervals[id_hyp].get_divisions() == 0) return;

	Hyperinterval& hyp = _intervals[id_hyp];
	uint index = (hyp.get_divisions() - 1) / _dimension;
	uint j = (hyp.get_divisions() - 1) % _dimension + 1;
	double h1 = static_cast<double>(HYPER_INTERVAL_SIDE_LENGTHS[20 - index]);
	double h2 = static_cast<double>(HYPER_INTERVAL_SIDE_LENGTHS[20 - index - 1]);
	double e = 0.5 * sqrt((_dimension - j) * h1 * h1 + j * h2 * h2);

	hyp.set_charact(golden_ratio(-e, e, e, id_hyp));
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

void TransformPM::balance(double& lipshConst) const {
	if (_iteration % 5 == 0) lipshConst *= _parameters._Gain;
	else lipshConst *= _parameters._Reduce;
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

void TransformPM::solve() {
	initialization();
	uint id_current_interval = 0;
	bool flag = true;

	for (uint i = 0; ((i < 500) && (flag)); ++i) {
		id_current_interval = iterate(id_current_interval);
		Hyperinterval& hyp = _intervals[id_current_interval];

		if (_intervals[id_current_interval].get_diagonal() < _parameters._eps)
			flag = false;
	}

	if (!flag) std::cout << "STOPPED BY PRECISION" << std::endl;

	std::cout << "Current minimum: " << _current_minimum << std::endl;
	EncodedCoordinates ec(_dimension);
	Point& point = _points[_id_minimum];
	for (uint i = 0; i < _dimension; ++i)
		ec[i] = _coords[point.get_id_coord() + i];
	CoordinatesValues dc = _problem.decode_coordinates(ec);
	std::cout << "Point: " << std::endl;
	for (uint i = 0; i < _dimension; ++i)
		std::cout << dc[i] << std::endl;
	std::cout << "Num of evaluations: " << _generated_points << std::endl;
}