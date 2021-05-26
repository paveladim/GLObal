#include "TPiyavskiiMethod.h"

TPiyavskiiMethod::TPiyavskiiMethod(const uint& out_dim, const uint& out_constr, const uint& depth,
	TProblem& out_prob, const GainLipshConstant& out_gainObj, const GainLipshConstant& out_gainCst, const double& beta) : 
	TMethodDivByThree(out_dim, out_constr, depth, out_prob), F_gainObjective(out_gainObj),
	F_gainConstraints(out_gainCst), delta(0.0000000001) {
	F_criticalSize = beta * sqrt(F_dimension * (CoordinateValue)MAX_POWER_THREE * (CoordinateValue)MAX_POWER_THREE);
	does_LipshConstValue_change = false;
}

void TPiyavskiiMethod::compute_characteristic(const uint& id_Hyp) {
	THyperinterval& hyp = F_intervals[id_Hyp];
	double charact = 0.5 * (F_evaluations[hyp.get_idEvaluationsA()] + F_evaluations[hyp.get_idEvaluationsB()]);

	if (hyp.get_diagonal() < F_criticalSize) {
		double ratio = hyp.get_diagonal() / F_criticalSize;
		charact = charact - 0.5 * hyp.get_diagonal() * (ratio * F_gainObjective *
			F_globalLipshEvaluations[0] + (1 - ratio) * hyp.get_maxLipshEvaluations()[0]);
	}
	else {
		charact = charact - 0.5 * hyp.get_diagonal() * F_globalLipshEvaluations[0];
	}

	hyp.set_characteristic(charact);
}

uint TPiyavskiiMethod::choose_optimal_to_trisect() {
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

	return id_optimal_hyp;
}

void TPiyavskiiMethod::compute_localLipshConst(const uint& id_Hyp) {
	static std::vector<LipschitzConstantValue> new_llcv(F_constraints + 1);
	THyperinterval& hyp = F_intervals[id_Hyp];

	for (uint i = 0; i < F_constraints + 1; ++i) {
		new_llcv[i] = F_evaluations[hyp.get_idEvaluationsA() + i] - F_evaluations[hyp.get_idEvaluationsB() + i];
		new_llcv[i] = std::abs(new_llcv[i]) / hyp.get_diagonal();
	}

	hyp.update_queuesLipshEvaluations(new_llcv, delta);
}

void TPiyavskiiMethod::update_globalLipshEval(const uint& id_Hyp) {
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
}

void TPiyavskiiMethod::update_all_characteristics() {
	for (uint i = 0; i < F_generated_intervals; ++i)
		compute_characteristic(i);
}

uint TPiyavskiiMethod::do_step(const uint& id_divHyp) {
	trisect_interval(id_divHyp);
	compute_localLipshConst(id_divHyp);
	compute_localLipshConst(F_generated_intervals - 1);
	compute_localLipshConst(F_generated_intervals - 2);
	update_globalLipshEval(id_divHyp);
	return choose_optimal_to_trisect();
}

void TPiyavskiiMethod::launch_method() {
	initialization();
	compute_localLipshConst(0);
	uint id_current_interval = 0;
	for (uint i = 0; i < 100; ++i) {
		id_current_interval = do_step(id_current_interval);
		THyperinterval& hyp = F_intervals[id_current_interval];
		F_current_minimum = std::min(F_evaluations[hyp.get_idEvaluationsA()], F_evaluations[hyp.get_idEvaluationsB()]);
	}
}