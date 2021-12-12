#include <iostream>
#include <fstream>

#include "SimplePMwithSM.h"

SimplePMwithSM::SimplePMwithSM(const uint& dimension,
							   const uint& constraints,
							   Parameters& parameters,
							   Problem& problem) : 
	DivideByThree(dimension, constraints, parameters, problem),
	_areAllCharInfty(false),
	_doesGlobalChange(false),
	_localLipshEval(constraints + 1),
	_globalLipshEval(constraints + 1) {}

void SimplePMwithSM::calculate_localLipshConst(const uint& id_hyp) {
	// A, V, U, B
	static EncodedCoordinates points(4 * _dimension);
	static std::vector<double> incs(12);

	Hyperinterval& hyp1 = _intervals[id_hyp];
	Hyperinterval& hyp2 = _intervals[_generated_intervals - 2];
	Hyperinterval& hyp3 = _intervals[_generated_intervals - 1];

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
	
	for (uint i = 0; i < _constraints + 1; ++i) {
		generate_right_part(b, i, hyp1.get_idA() * (_constraints + 1),
								  hyp2.get_idA() * (_constraints + 1),
								  hyp2.get_idB() * (_constraints + 1),
								  hyp3.get_idB() * (_constraints + 1));

		SimplexMethod sm(45, 36, c, b, st);
		SimplexRes res = sm.solve();
		_localLipshEval[i] = res._resVal / ((double)MAX_POWER_THREE * (double)MAX_POWER_THREE);
		if (i == 0)
			std::cout << _localLipshEval[i] << std::endl;
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
		norms[i] = sqrt(norms[i]);
	}

	for (uint i = 0; i < 6; ++i) {
		A[6 * i][8] = -norms[i];
		A[6 * i + 1][8] = norms[i];
		A[6 * i + 2][8] = -0.5 * norms[i] * norms[i];
		A[6 * i + 3][8] = 0.5 * norms[i] * norms[i];
		A[6 * i + 4][8] = -0.5 * norms[i] * norms[i];
		A[6 * i + 5][8] = 0.5 * norms[i] * norms[i];
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

		A[6 * j + 2][2 * i] = incs[2 * j];
		A[6 * j + 2][2 * i + 1] = incs[2 * j + 1];

		A[6 * j + 3][2 * i] = incs[2 * j];
		A[6 * j + 3][2 * i + 1] = incs[2 * j + 1];

		A[6 * j + 4][2 * k] = -incs[2 * j];
		A[6 * j + 4][2 * k + 1] = -incs[2 * j + 1];

		A[6 * j + 5][2 * k] = -incs[2 * j];
		A[6 * j + 5][2 * k + 1] = -incs[2 * j + 1];

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

	b[2] = fv - fa;
	b[3] = fv - fa;
	b[4] = fa - fv;
	b[5] = fa - fv;

	b[8] = fu - fa;
	b[9] = fu - fa;
	b[10] = fa - fu;
	b[11] = fa - fu;

	b[14] = fb - fa;
	b[15] = fb - fa;
	b[16] = fa - fb;
	b[17] = fa - fb;

	b[20] = fu - fv;
	b[21] = fu - fv;
	b[22] = fv - fu;
	b[23] = fv - fu;

	b[26] = fb - fv;
	b[27] = fb - fv;
	b[28] = fv - fb;
	b[29] = fv - fb;

	b[32] = fb - fu;
	b[33] = fb - fu;
	b[34] = fu - fb;
	b[35] = fu - fb;
}


double SimplePMwithSM::mixedLipEval(const Hyperinterval& hyp, const uint& i) {
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

void
SimplePMwithSM::give_borders(double& l,
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

void SimplePMwithSM::calculate_globalLipshConst(const uint& id_hyp) {
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

void SimplePMwithSM::update_all_charact() {
	for (uint i = 0; i < _generated_intervals; ++i)
		calculate_characteristic(i);
}

uint SimplePMwithSM::optimal_to_trisect() {
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

uint SimplePMwithSM::iterate(const uint& id_hyp) {
	trisect_interval(id_hyp);
	calculate_localLipshConst(id_hyp);
	calculate_globalLipshConst(id_hyp);
	return optimal_to_trisect();
}

void SimplePMwithSM::solve() {
	initialization();
	uint id_current_interval = 0;
	std::ofstream out;
	out.open("D:\\materials\\projects\\visual_hyperinterval\\minimums.txt");
	if (out.is_open()) {
		bool flag = true;
		for (uint i = 0; ((i < 250) && (flag)); ++i) {
			id_current_interval = iterate(id_current_interval);
			Hyperinterval& hyp = _intervals[id_current_interval];

			/*if ((F_intervals[id_current_interval].get_diagonal() /
				((double)MAX_POWER_THREE * (double)MAX_POWER_THREE)) < eps)
				flag = false; */

			/*if ((_current_minimum < -0.80467) &&
				(std::abs(F_current_minimum + 0.80467) < eps)) flag = false;
			out << F_current_minimum << std::endl; */
		}

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
}