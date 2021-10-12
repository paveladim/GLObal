#include <iostream>
#include <fstream>
#include "TSimplePMwithoutSM.h"

TSimplePMwithoutSM::TSimplePMwithoutSM(const uint& out_dim,
									   const uint& out_constr,
									   const uint& depth,
									   TProblem& out_prob,
									   const GainLipshConstant& out_gainObj,
									   const GainLipshConstant& out_gainCst,
									   const double& beta,
									   const double& _eps) :
	TMethodDivByThree(out_dim, out_constr, depth, out_prob),
	F_gainObjective(out_gainObj),
	F_gainConstraints(out_gainCst),
	delta(0.0000000001),
	eps(_eps) {
	F_criticalSize =
		beta * sqrt(F_dimension * (CoordinateValue)MAX_POWER_THREE * (CoordinateValue)MAX_POWER_THREE);
	does_LipshConstValue_change = false;
}

void TSimplePMwithoutSM::initialization() {
	THyperinterval::init_static(F_dimension, F_constraints, F_queueDepth);
	resize_points_deque();
	resize_coords_deque();

	F_points[0].F_idThis = get_new_id();
	F_points[1].F_idThis = get_new_id();

	TPoint& a = F_points[0];
	TPoint& b = F_points[1];

	for (uint i = 0; i < F_dimension; ++i)
		F_coords[a.F_idThis * F_dimension + i] = 0;

	for (uint i = 0; i < F_dimension; ++i)
		F_coords[b.F_idThis * F_dimension + i] = MAX_POWER_THREE;

	compute_evaluations(a.F_idThis);
	compute_evaluations(b.F_idThis);

	resize_intervals_deque();
	F_intervals[0].set_idPointA(a.F_idThis);
	F_intervals[0].set_idPointB(b.F_idThis);
	F_intervals[0].set_idThis(get_new_interval());
	F_intervals[0].init_queues();
	compute_diagonal(F_intervals[0].get_idThis());
	//compute_localLipshConst(F_intervals[0].get_idThis());

	//for (uint i = 0; i < F_constraints + 1; ++i)
		//if (F_globalLipshEvaluations[i] < F_intervals[0].get_maxLipshEvaluations()[i]) {
			//F_globalLipshEvaluations[i] = F_intervals[0].get_maxLipshEvaluations()[i];
		//}

	compute_characteristic(F_intervals[0].get_idThis());
}

void TSimplePMwithoutSM::compute_characteristic(const uint& id_Hyp) {
	THyperinterval& hyp = F_intervals[id_Hyp];
	double charact = 0;

	if (hyp.get_diagonal() < F_criticalSize) {
		double ratio = hyp.get_diagonal() / F_criticalSize;
		charact = F_gainObjective * F_globalLipshEvaluations[0];
		charact += (1 - ratio) * hyp.get_maxLipshEvaluations()[0];
	}
	else {
		charact = F_globalLipshEvaluations[0];
	}

	hyp.set_characteristic(charact);
}

uint TSimplePMwithoutSM::choose_optimal_to_trisect() {
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

void TSimplePMwithoutSM::compute_localLipshConst(const uint& id_Hyp1, 
												 const uint& id_Hyp2, 
												 const uint& id_Hyp3) {
	static std::vector<LipschitzConstantValue> new_llcv(F_constraints + 1);
	THyperinterval& hyp1 = F_intervals[id_Hyp1];
	THyperinterval& hyp2 = F_intervals[id_Hyp2];
	THyperinterval& hyp3 = F_intervals[id_Hyp3];

	uint n = F_iter / F_dimension;
	uint j = F_iter % F_dimension;

	double h2 = double(MAX_POWER_THREE) / double(HYPER_INTERVAL_SIDE_LENGTHS[n]);
	double h1 = sqrt(j * h2 * h2 / 9 + (F_dimension - j - 1) * h2 * h2);

	double e1 = sqrt(0.25 * h1 * h1 + h2 * h2 / 36.0);
	double e2 = 0.5 * sqrt(h1 * h1 + h2 * h2);

	for (uint i = 0; i < F_constraints + 1; ++i) {
		new_llcv[i] = F_evaluations[hyp1.get_idEvaluationsA() + i] + F_evaluations[hyp3.get_idEvaluationsB() + i];
		new_llcv[i] -= F_evaluations[hyp2.get_idEvaluationsA() + i] + F_evaluations[hyp2.get_idEvaluationsB() + i];
		new_llcv[i] = new_llcv[i] / (e2 * e2 - e1 * e1);
	}

	hyp1.update_queuesLipshEvaluations(new_llcv, delta);
	hyp2.update_queuesLipshEvaluations(new_llcv, delta);
	hyp3.update_queuesLipshEvaluations(new_llcv, delta);
}

void TSimplePMwithoutSM::update_globalLipshEval(const uint& id_Hyp) {
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

void TSimplePMwithoutSM::update_all_characteristics() {
	for (uint i = 0; i < F_generated_intervals; ++i)
		compute_characteristic(i);
}

uint TSimplePMwithoutSM::do_step(const uint& id_divHyp) {
	trisect_interval(id_divHyp);
	compute_localLipshConst(id_divHyp,
							F_generated_intervals - 2,
							F_generated_intervals - 1);
	update_globalLipshEval(id_divHyp);
	return choose_optimal_to_trisect();
}

void TSimplePMwithoutSM::launch_method() {
	initialization();
	uint id_current_interval = 0;
	std::ofstream out;
	out.open("D:\\materials\\projects\\visual_hyperinterval\\minimums.txt");
	if (out.is_open()) {
		for (F_iter = 0;
			(F_iter < 200) && (std::abs(F_current_minimum - F_intervals[id_current_interval].get_characteristic()) >= eps);
			++F_iter) {
			id_current_interval = do_step(id_current_interval);
			THyperinterval& hyp = F_intervals[id_current_interval];

			if (std::min(F_evaluations[hyp.get_idEvaluationsA()], F_evaluations[hyp.get_idEvaluationsB()]) < F_current_minimum)
				F_current_minimum = std::min(F_evaluations[hyp.get_idEvaluationsA()], F_evaluations[hyp.get_idEvaluationsB()]);
			out << F_current_minimum << std::endl;
		}

		std::cout << "Current minimum: " << F_current_minimum << std::endl;
	}
}

void TSimplePMwithoutSM::write_generated_points_to_file() {
	std::ofstream out;
	out.open("D:\\materials\\projects\\visual_hyperinterval\\points.txt");
	if (out.is_open())
	{
		for (uint i = 0; i < F_generated_points * F_dimension; ++i)
			out << F_coords[i] << std::endl;
	}
}

void TSimplePMwithoutSM::write_intervals_to_file() {
	std::ofstream out;
	out.open("D:\\materials\\projects\\visual_hyperinterval\\hyp.txt");
	if (out.is_open())
	{
		for (uint id_hyp = 0; id_hyp < F_generated_intervals; ++id_hyp) {
			for (uint i = 0; i < F_dimension; ++i)
				out << F_coords[F_intervals[id_hyp].get_idA() + i] << std::endl;
			for (uint i = 0; i < F_dimension; ++i)
				out << F_coords[F_intervals[id_hyp].get_idB() + i] << std::endl;
		}
	}
}
