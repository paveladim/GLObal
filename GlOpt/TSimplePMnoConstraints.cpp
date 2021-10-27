#include <iostream>
#include <fstream>

#include "TSimplePMnoConstraints.h"

TSimplePMnoConstraints::TSimplePMnoConstraints(const uint& out_dim,
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
	does_LipshConstValue_change(false) {
	F_criticalSize =
		beta * sqrt(F_dimension * (CoordinateValue)MAX_POWER_THREE * (CoordinateValue)MAX_POWER_THREE);
}

void 
TSimplePMnoConstraints::compute_characteristic(const uint& id_Hyp) {
	if (F_intervals[id_Hyp].get_divisions() == 0) return;
	THyperinterval& hyp = F_intervals[id_Hyp];
	double mixed_LipshitzEval = get_mixedLipshitzEval(hyp);
	double t_min = 0.5 * (F_evaluations[hyp.get_idEvaluationsA()] - 
				   F_evaluations[hyp.get_idEvaluationsB()]);
	t_min = t_min / mixed_LipshitzEval;

	uint index = (hyp.get_divisions() - 1) / F_dimension;
	uint j = F_dimension - hyp.get_div_axis();
	double h1 = (double)HYPER_INTERVAL_SIDE_LENGTHS[20 - index];
	double h2 = (double)HYPER_INTERVAL_SIDE_LENGTHS[20 - index - 1];
	double e = 0.5 * sqrt((F_dimension - j) * h1 * h1 + j * h2 * h2);

	t_min = t_min / (2 * e);

	if ((t_min >= -e) && (t_min <= e)) {
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

double 
TSimplePMnoConstraints::get_mixedLipshitzEval(const THyperinterval& hyp) {
	double mixed_LipshitzEval = 0.0;

	if (hyp.get_diagonal() < F_criticalSize) {
		double ratio = hyp.get_diagonal() / F_criticalSize;
		mixed_LipshitzEval = ratio * F_globalLipshEvaluations[0] +
							 (1 - ratio) * hyp.get_maxLipshEvaluations()[0];
	}
	else mixed_LipshitzEval = F_globalLipshEvaluations[0];

	return mixed_LipshitzEval;
}

void 
TSimplePMnoConstraints::compute_localLipshConst(const uint& id_Hyp) {
	static std::vector<LipschitzConstantValue> new_llcv(F_constraints + 1);
	THyperinterval& hyp1 = F_intervals[id_Hyp];
	THyperinterval& hyp2 = F_intervals[F_generated_intervals - 2];
	THyperinterval& hyp3 = F_intervals[F_generated_intervals - 1];

	uint index = (hyp1.get_divisions() - 1) / F_dimension;
	uint j = F_dimension - hyp1.get_div_axis();

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
		new_llcv[i] = std::abs(new_llcv[i] / (e2 * e2 - e1 * e1));
	}

	hyp1.update_queuesLipshEvaluations(new_llcv, delta / ((double)MAX_POWER_THREE * (double)MAX_POWER_THREE));
	hyp2.update_queuesLipshEvaluations(new_llcv, delta / ((double)MAX_POWER_THREE * (double)MAX_POWER_THREE));
	hyp3.update_queuesLipshEvaluations(new_llcv, delta / ((double)MAX_POWER_THREE * (double)MAX_POWER_THREE));
}

void 
TSimplePMnoConstraints::update_globalLipshEval(const uint& id_Hyp) {
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
TSimplePMnoConstraints::update_all_characteristics() {
	for (uint i = 0; i < F_generated_intervals; ++i)
		compute_characteristic(i);
}

uint 
TSimplePMnoConstraints::choose_optimal_to_trisect() {
	uint id_optimal_hyp = 0;
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
TSimplePMnoConstraints::do_step(const uint& id_divHyp) {
	trisect_interval(id_divHyp);
	compute_localLipshConst(id_divHyp);
	update_globalLipshEval(id_divHyp);
	return choose_optimal_to_trisect();
}

void 
TSimplePMnoConstraints::launch_method() {
	initialization();
	uint id_current_interval = 0;
	std::ofstream out;
	out.open("D:\\materials\\projects\\visual_hyperinterval\\minimums.txt");
	if (out.is_open()) {
		bool flag = true;
		for (uint i = 0; (i < 500) && (flag); ++i) {
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