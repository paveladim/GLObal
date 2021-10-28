#include <iostream>
#include <fstream>

#include "TSimplePMwithConstraints.h"

TSimplePMwithConstraints::TSimplePMwithConstraints(const uint& out_dim,
												   const uint& out_constr,
												   const uint& depth,
												   TProblem& out_prob,
												   const GainLipshConstant& out_gainObj,
												   const GainLipshConstant& out_gainCst,
												   const double& beta,
												   const double& _eps
												   ) :
	TMethodDivByThree(out_dim, out_constr, depth, out_prob),
	F_gainObjective(out_gainObj),
	F_gainConstraints(out_gainCst),
	delta(0.0000000001),
	eps(_eps),
	does_LipshConstValue_change(false),
	are_allCharInfty(false) {
	F_criticalSize =
		beta * sqrt(F_dimension * (CoordinateValue)MAX_POWER_THREE * (CoordinateValue)MAX_POWER_THREE);
}

void 
TSimplePMwithConstraints::give_borders(double& l, 
									   double& r, 
									   const THyperinterval& hyp) {
	double a = 0.0;
	double b = 0.0;
	double c = 0.0;
	double fa = 0.0;
	double fb = 0.0;
	double M = 0.0;
	double e = r;

	for (uint i = 1; i < F_constraints + 1; ++i) {
		M = get_mixedLipshitzEval(hyp, i);
		fa = F_evaluations[hyp.get_idEvaluationsA() + i];
		fb = F_evaluations[hyp.get_idEvaluationsB() + i];

		a = 0.5 * M;
		b = 0.5 * (fb - fa) / e;
		c = 0.5 * (fb + fa - M * e * e);

		if (b * b - 4 * a * c >= 0) {
			fa = -0.5 * b / a + sqrt(b * b - 4 * a * c);
			fb = -0.5 * b / a - sqrt(b * b - 4 * a * c);

			if (fa > fb) std::swap(fa, fb);
			if (fa >= l) l = fa;
			if (fb <= r) r = fb;
		}
		else {
			l = e;
			r = -e;
			return;
		}
	}
}

double 
TSimplePMwithConstraints::get_mixedLipshitzEval(const THyperinterval& hyp, 
												const uint& i) {
	double mixed_LipshitzEval = 0.0;

	if (hyp.get_diagonal() < get_critical_size()) {
		double ratio = hyp.get_diagonal() / get_critical_size();
		mixed_LipshitzEval = ratio * F_globalLipshEvaluations[i] +
			(1 - ratio) * hyp.get_maxLipshEvaluations()[i];
	}
	else mixed_LipshitzEval = F_globalLipshEvaluations[i];

	return mixed_LipshitzEval;
}

void 
TSimplePMwithConstraints::compute_characteristic(const uint& id_Hyp) {
	if (F_intervals[id_Hyp].get_divisions() == 0) return;
	THyperinterval& hyp = F_intervals[id_Hyp];
	double mixed_LipshitzEval = get_mixedLipshitzEval(hyp, 0);
	double t_min = 0.5 * (F_evaluations[hyp.get_idEvaluationsA()] -
		F_evaluations[hyp.get_idEvaluationsB()]);
	t_min = t_min / mixed_LipshitzEval;

	uint index = (hyp.get_divisions() - 1) / F_dimension;
	uint j = (hyp.get_divisions() - 1) % F_dimension + 1;
	double h1 = (double)HYPER_INTERVAL_SIDE_LENGTHS[20 - index];
	double h2 = (double)HYPER_INTERVAL_SIDE_LENGTHS[20 - index - 1];
	double e = 0.5 * sqrt((F_dimension - j) * h1 * h1 + j * h2 * h2);

	t_min = t_min / (2 * e);
	double left = -e;
	double right = e;
	give_borders(left, right, hyp);

	if (left == e) {
		hyp.set_characteristic(std::numeric_limits<double>::max());
	}
	else {
		if ((t_min > -e) && (t_min < e)) {
			double charact = -F_evaluations[hyp.get_idEvaluationsA()] * (t_min - e);
			charact = charact + F_evaluations[hyp.get_idEvaluationsB()] * (t_min + e);
			charact = 0.5 * charact / e;
			charact = charact + 0.5 * mixed_LipshitzEval * (t_min * t_min - e * e);
			hyp.set_characteristic(charact);
		}
		else {
			hyp.set_characteristic(
				std::min(F_evaluations[hyp.get_idEvaluationsA()],
					F_evaluations[hyp.get_idEvaluationsB()])
			);
		}
	}
}

void 
TSimplePMwithConstraints::check_ifAllCharAreInfty() {
	are_allCharInfty = true;
	for (uint i = 0; i < F_generated_intervals; ++i)
		if (F_intervals[i].get_characteristic() < std::numeric_limits<double>::max())
			are_allCharInfty = false;

}

void
TSimplePMwithConstraints::compute_localLipshConst(const uint& id_Hyp) {
	static std::vector<LipschitzConstantValue> new_llcv(F_constraints + 1);
	THyperinterval& hyp1 = F_intervals[id_Hyp];
	THyperinterval& hyp2 = F_intervals[F_generated_intervals - 2];
	THyperinterval& hyp3 = F_intervals[F_generated_intervals - 1];

	//std::cout << "HYPERINTERVAL" << std::endl;
	//std::cout << "A=(" << F_coords[hyp1.get_idA()] << ',' << F_coords[hyp1.get_idA() + 1] << ')' << std::endl;
	//std::cout << "B=(" << F_coords[hyp3.get_idB()] << ',' << F_coords[hyp3.get_idB() + 1] << ')' << std::endl;
	//std::cout << "U=(" << F_coords[hyp2.get_idA()] << ',' << F_coords[hyp2.get_idA() + 1] << ')' << std::endl;
	//std::cout << "V=(" << F_coords[hyp2.get_idB()] << ',' << F_coords[hyp2.get_idB() + 1] << ')' << std::endl;

	uint index = (hyp1.get_divisions() - 1) / F_dimension;
	uint j = (hyp1.get_divisions() - 1) % F_dimension + 1;

	//std::cout << "index=" << 20 - index << " j = " << j << std::endl << std::endl;

	double h1 = (double)HYPER_INTERVAL_SIDE_LENGTHS[20 - index];
	double h2 = (double)HYPER_INTERVAL_SIDE_LENGTHS[20 - index - 1];

	double e1 = 0.5 * sqrt(((double)(F_dimension - j) + 1) * h1 * h1 +
		(j - 1) * h2 * h2);
	double e2 = 0.5 * sqrt((F_dimension - j) * h1 * h1 + j * h2 * h2);

	for (uint i = 0; i < F_constraints + 1; ++i) {
		new_llcv[i] = F_evaluations[hyp1.get_idEvaluationsA() + i]
			+ F_evaluations[hyp3.get_idEvaluationsB() + i];
		new_llcv[i] -= (F_evaluations[hyp2.get_idEvaluationsA() + i]
			+ F_evaluations[hyp2.get_idEvaluationsB() + i]);

		if (i == 0)
			new_llcv[i] = F_gainObjective * std::abs(new_llcv[i] / (e2 * e2 - e1 * e1));
		else new_llcv[i] = F_gainConstraints * std::abs(new_llcv[i] / (e2 * e2 - e1 * e1));
	}

	hyp1.update_queuesLipshEvaluations(new_llcv, delta / ((double)MAX_POWER_THREE * (double)MAX_POWER_THREE));
	hyp2.update_queuesLipshEvaluations(new_llcv, delta / ((double)MAX_POWER_THREE * (double)MAX_POWER_THREE));
	hyp3.update_queuesLipshEvaluations(new_llcv, delta / ((double)MAX_POWER_THREE * (double)MAX_POWER_THREE));
}

void
TSimplePMwithConstraints::update_globalLipshEval(const uint& id_Hyp) {
	THyperinterval& hyp1 = F_intervals[id_Hyp];
	THyperinterval& hyp2 = F_intervals[F_generated_intervals - 1];
	THyperinterval& hyp3 = F_intervals[F_generated_intervals - 2];

	for (uint i = 0; i < F_constraints + 1; ++i) {
		if (F_globalLipshEvaluations[i] < hyp1.get_maxLipshEvaluations()[i]) {
			does_LipshConstValue_change = true;
			F_globalLipshEvaluations[i] = hyp1.get_maxLipshEvaluations()[i];
		}

		if (F_globalLipshEvaluations[i] < hyp2.get_maxLipshEvaluations()[i]) {
			does_LipshConstValue_change = true;
			F_globalLipshEvaluations[i] = hyp2.get_maxLipshEvaluations()[i];
		}

		if (F_globalLipshEvaluations[i] < hyp3.get_maxLipshEvaluations()[i]) {
			does_LipshConstValue_change = true;
			F_globalLipshEvaluations[i] = hyp3.get_maxLipshEvaluations()[i];
		}
	}

	if (does_LipshConstValue_change) {
		update_all_characteristics();
		does_LipshConstValue_change = false;
	}
	else {
		compute_characteristic(id_Hyp);
		compute_characteristic(F_generated_intervals - 2);
		compute_characteristic(F_generated_intervals - 1);
	}
}

void
TSimplePMwithConstraints::update_all_characteristics() {
	for (uint i = 0; i < F_generated_intervals; ++i)
		compute_characteristic(i);
}

uint
TSimplePMwithConstraints::choose_optimal_to_trisect() {
	uint id_optimal_hyp = 0;

	if (are_allCharInfty) {
		double optimal_charact = F_intervals[id_optimal_hyp].get_diagonal();
		double current_charact = 0.0;

		for (uint id_hyp = 1; id_hyp < F_generated_intervals; ++id_hyp) {
			current_charact = F_intervals[id_hyp].get_diagonal();
			if (std::abs(current_charact - optimal_charact) < std::numeric_limits<double>::epsilon());
			else if (optimal_charact < current_charact) {
				optimal_charact = current_charact;
				id_optimal_hyp = id_hyp;
			}
		}
	}
	else {
		double optimal_charact = F_intervals[id_optimal_hyp].get_characteristic();
		double current_charact = 0.0;

		for (uint id_hyp = 1; id_hyp < F_generated_intervals; ++id_hyp) {
			current_charact = F_intervals[id_hyp].get_characteristic();
			if (std::abs(current_charact - optimal_charact) < std::numeric_limits<double>::epsilon());
			else if (optimal_charact > current_charact) {
				optimal_charact = current_charact;
				id_optimal_hyp = id_hyp;
			}
		}
	}

	THyperinterval& hyp = F_intervals[id_optimal_hyp];

	EncodedCoordinates a(F_dimension);
	EncodedCoordinates b(F_dimension);

	for (uint i = 0; i < F_dimension; ++i) {
		a[i] = F_coords[hyp.get_idA() + i];
		b[i] = F_coords[hyp.get_idB() + i];
	}

	CoordinatesValues dec_a = Fp.decode_coordinates(a);
	CoordinatesValues dec_b = Fp.decode_coordinates(b);

	std::cout << "(" << std::setprecision(4)
		<< 0.5 * (dec_a[0] + dec_b[0]) << ';'
		<< 0.5 * (dec_a[1] + dec_b[1]) << ')';

	std::cout << " D = " << hyp.get_diagonal() / MAX_POWER_THREE;
	std::cout << " Char: " << hyp.get_characteristic();

	std::vector<double> evals = hyp.get_maxLipshEvaluations();
	std::cout << " O: " << evals[0] * (double)MAX_POWER_THREE * (double)MAX_POWER_THREE
		<< " C: " << evals[1] * (double)MAX_POWER_THREE * (double)MAX_POWER_THREE;

	std::cout << " (" << F_evaluations[hyp.get_idEvaluationsA()] << ';'
		<< F_evaluations[hyp.get_idEvaluationsB()] << ')' << std::endl;

	return id_optimal_hyp;
}

uint
TSimplePMwithConstraints::do_step(const uint& id_divHyp) {
	trisect_interval(id_divHyp);
	compute_localLipshConst(id_divHyp);
	update_globalLipshEval(id_divHyp);
	check_ifAllCharAreInfty();
	return choose_optimal_to_trisect();
}

void
TSimplePMwithConstraints::launch_method() {
	initialization();
	uint id_current_interval = 0;
	std::ofstream out;
	out.open("D:\\materials\\projects\\visual_hyperinterval\\minimums.txt");
	if (out.is_open()) {
		bool flag = true;
		for (uint i = 0; (i < 100) && (flag); ++i) {
			id_current_interval = do_step(id_current_interval);
			THyperinterval& hyp = F_intervals[id_current_interval];

			if (F_intervals[id_current_interval].get_diagonal() < eps * F_criticalSize) flag = false;
			if (std::min(F_evaluations[hyp.get_idEvaluationsA()], F_evaluations[hyp.get_idEvaluationsB()]) < F_current_minimum)
				F_current_minimum = std::min(F_evaluations[hyp.get_idEvaluationsA()], F_evaluations[hyp.get_idEvaluationsB()]);
			out << F_current_minimum << std::endl;
		}

		std::cout << "Current minimum: " << F_current_minimum << std::endl;
	}
}