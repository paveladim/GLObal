#include <iostream>
#include "SimplePMwithLagrange.h"

SimplePMwithLagrange::SimplePMwithLagrange(const uint& dimension,
										   const uint& constraints,
										   Parameters& parameters,
										   Problem& problem) :
	DivideByThree(dimension, constraints, parameters, problem),
	_areAllCharInfty(false),
	_doesGlobalChange(false),
	_localLipshEval(constraints + 1),
	_globalLipshEval(constraints + 1) {}

void SimplePMwithLagrange::balance(double& lipshConst) const {
	if (_iteration % 5 == 0) lipshConst *= _parameters._stochGain;
	else lipshConst *= _parameters._stochReduce;
}

void SimplePMwithLagrange::calculate_localLipshConst(const uint& id_hyp) {
	Hyperinterval& hyp1 = _intervals[id_hyp];
	Hyperinterval& hyp2 = _intervals[_generated_intervals - 2];
	Hyperinterval& hyp3 = _intervals[_generated_intervals - 1];

	double fa, fv, fu, fb;
	double dist = 0.0;
	for (uint i = 0; i < _dimension; ++i)
		dist = dist +
			  ((double)_coords[hyp1.get_coordA() + i] - 
			   (double)_coords[hyp3.get_coordB() + i]) * 
			  ((double)_coords[hyp1.get_coordA() + i] - 
			   (double)_coords[hyp3.get_coordB() + i]);

	for (uint i = 0; i < _constraints + 1; ++i) {
		fa = _evaluations[hyp1.get_evalA() + i];
		fv = _evaluations[hyp2.get_evalA() + i];
		fu = _evaluations[hyp2.get_evalB() + i];
		fb = _evaluations[hyp3.get_evalB() + i];
		_localLipshEval[i] = 4 * (fb - fu - fv + fa) / dist;
		double conv = _localLipshEval[i] * (double)MAX_POWER_THREE * (double)MAX_POWER_THREE;
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

double SimplePMwithLagrange::mixedLipEval(const Hyperinterval& hyp, const uint& i) {
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

void SimplePMwithLagrange::give_borders(double& l,
										double& r,
										Hyperinterval& hyp) {
	double a = 0.0;
	double b = 0.0;
	double c = 0.0;
	double fa = 0.0;
	double fb = 0.0;
	double M = 0.0;
	double e = r;

	for (uint i = 1; i < _constraints + 1; ++i) {
		M = mixedLipEval(hyp, i);
		fa = _evaluations[hyp.get_evalA() + i];
		fb = _evaluations[hyp.get_evalB() + i];

		a = 0.5 * M;
		b = 0.5 * (fb - fa) / e;
		c = 0.5 * (fb + fa - M * e * e);

		if (b * b - 4 * a * c >= 0) {
			fa = (-b - sqrt(b * b - 4 * a * c)) * 0.5 / a;
			fb = (-b + sqrt(b * b - 4 * a * c)) * 0.5 / a;

			if (fa > fb) std::swap(fa, fb);
			if ((fa >= l) && (fa <= r)) l = fa;
			if ((fb <= r) && (fb >= l)) r = fb;
		}
		else {
			l = e;
			r = -e;
			return;
		}
	}
}

void SimplePMwithLagrange::calculate_characteristic(const uint& id_hyp) {
	if (_intervals[id_hyp].get_divisions() == 0) return;
	Hyperinterval& hyp = _intervals[id_hyp];

	double mixed_LipshEval = mixedLipEval(hyp, 0);
	double t_min = 0.5 * (_evaluations[hyp.get_evalA()] -
		_evaluations[hyp.get_evalB()]);
	t_min = t_min / mixed_LipshEval;

	uint index = 20 - (hyp.get_divisions() - 1) / _dimension;
	uint j = hyp.get_previous_axis();
	double h2 = HYPER_INTERVAL_SIDE_LENGTHS[index];
	double h1 = sqrt(j * h2 * h2 / 9.0 + (_dimension - j - 1) * h2 * h2);
	double e = 0.5 * sqrt(h1 * h1 + h2 * h2 / 9.0);

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

void SimplePMwithLagrange::calculate_globalLipshConst(const uint& id_hyp) {
	Hyperinterval& hyp1 = _intervals[id_hyp];
	Hyperinterval& hyp2 = _intervals[_generated_intervals - 1];
	Hyperinterval& hyp3 = _intervals[_generated_intervals - 2];

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

void SimplePMwithLagrange::update_all_charact() {
	for (uint i = 0; i < _generated_intervals; ++i)
		calculate_characteristic(i);
}

uint SimplePMwithLagrange::optimal_to_trisect() {
	uint id_optimal_hyp = 0;

	if (_areAllCharInfty) {
		double optimal_charact = _intervals[id_optimal_hyp].get_charact(_areAllCharInfty);
		double current_charact = 0.0;

		for (uint id_hyp = 1; id_hyp < _generated_intervals; ++id_hyp) {
			current_charact = _intervals[id_hyp].get_charact();
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

uint SimplePMwithLagrange::iterate(const uint& id_hyp) {
	trisect_interval(id_hyp);
	calculate_localLipshConst(id_hyp);
	calculate_globalLipshConst(id_hyp);
	return optimal_to_trisect();
}

void SimplePMwithLagrange::solve() {
	initialization();
	uint id_current_interval = 0;
	bool flag = true;

	_iteration = 0;
	for (; ((_iteration < 500) && (flag)); ++_iteration) {
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