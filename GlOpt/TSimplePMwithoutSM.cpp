#include <iostream>
#include <fstream>
#include "TSimplePMwithoutSM.h"

EncodedCoordinates TSimplePMwithoutSM::F_decodeA;
EncodedCoordinates TSimplePMwithoutSM::F_decodeB;
CoordinatesValues TSimplePMwithoutSM::F_decodedA;
CoordinatesValues TSimplePMwithoutSM::F_decodedB;

TSimplePMwithoutSM::TSimplePMwithoutSM(const uint& out_dim,
									   const uint& out_constr,
									   const uint& depth,
									   TProblem& out_prob,
									   const GainLipshConstant& out_gainObj,
									   const GainLipshConstant& out_gainCst,
									   const double& beta,
									   const double& _eps,
									   const double& _nu) :
	TMethodDivByThree(out_dim, out_constr, depth, out_prob),
	F_gainObjective(out_gainObj),
	F_gainConstraints(out_gainCst),
	delta(0.0000000001),
	eps(_eps),
	nu(_nu),
	F_iter(0),
	do_AllCharAreInfty(false) {
	F_criticalSize =
		beta * sqrt(F_dimension * (CoordinateValue)MAX_POWER_THREE * (CoordinateValue)MAX_POWER_THREE);
	does_LipshConstValue_change = false;

	F_decodeA.resize(F_dimension);
	F_decodeB.resize(F_dimension);
	F_decodedA.resize(F_dimension);
	F_decodedB.resize(F_dimension);
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
}

void TSimplePMwithoutSM::compute_characteristic(const uint& id_Hyp) {
	THyperinterval& hyp = F_intervals[id_Hyp];
	static std::vector<LipschitzConstantValue> mixed_lipshEvaluations(F_constraints + 1);
	double charact = 0.0;

	if (hyp.get_diagonal() < F_criticalSize) {
		double ratio = hyp.get_diagonal() / F_criticalSize;
		mixed_lipshEvaluations[0] = F_gainObjective * F_globalLipshEvaluations[0];
		mixed_lipshEvaluations[0] += (1 - ratio) * hyp.get_maxLipshEvaluations()[0];
		for (uint i = 1; i < F_constraints + 1; ++i) {
			mixed_lipshEvaluations[i] = F_gainConstraints * F_globalLipshEvaluations[i];
			mixed_lipshEvaluations[i] += (1 - ratio) * hyp.get_maxLipshEvaluations()[i];
		}
	}
	else {
		for (uint i = 0; i < F_constraints + 1; ++i)
			mixed_lipshEvaluations[i] = F_globalLipshEvaluations[i];
	}

	double e = 0.0;
	uint id_A = F_intervals[id_Hyp].get_idA();
	uint id_B = F_intervals[id_Hyp].get_idB();

	for (uint i = 0; i < F_dimension; ++i) {
		F_decodeA[i] = F_coords[id_A + i];
		F_decodeB[i] = F_coords[id_B + i];
	}

	F_decodedA = Fp.decode_coordinates(F_decodeA);
	F_decodedB = Fp.decode_coordinates(F_decodeB);

	for (uint i = 0; i < F_dimension; ++i)
		e += (((double)F_decodeB[i] - (double)F_decodeA[i]) * ((double)F_decodeB[i] - (double)F_decodeA[i]));

	/*for (uint i = 0; i < F_dimension; ++i)
		e += ((F_decodedB[i] - F_decodedA[i]) * (F_decodedB[i] - F_decodedA[i]));*/

	e = sqrt(e) * 0.5;

	double left_border_MPOS = -e;
	double right_border_MPOS = e;
	double left_border_MNEG = e;
	double right_border_MNEG = -e;
	bool does_empty = false;

	id_A = F_intervals[id_Hyp].get_idEvaluationsA();
	id_B = F_intervals[id_Hyp].get_idEvaluationsB();

	static double a = 0.0;
	static double b = 0.0;
	static double c = 0.0;
	for (uint i = 1; (i < F_constraints + 1) && (!does_empty); ++i) {
		a = 0.5 * mixed_lipshEvaluations[i];
		b = 0.5 * (F_evaluations[id_B + i] - F_evaluations[id_A + i]) / e;
		c = 0.5 * (F_evaluations[id_A + i] + F_evaluations[id_B + i]);
		c = c - mixed_lipshEvaluations[i] * e * e * 0.5;

		if ((b * b - 4 * a * c) < 0) {
			if (mixed_lipshEvaluations[i] > 0) does_empty = true;
		}
		else {
			double left_border_tmp = -0.5 * b / a -
				sqrt(0.25 * b * b / a / a - a * c);
			double right_border_tmp = -0.5 * b / a +
				sqrt(0.25 * b * b / a / a - a * c);
			if (right_border_tmp < left_border_tmp)
				std::swap(right_border_tmp, left_border_tmp);

			if (mixed_lipshEvaluations[i] > 0) {
				if ((left_border_tmp >= left_border_MPOS) &&
					(left_border_tmp <= right_border_MPOS))
					left_border_MPOS = left_border_tmp;
				else if (left_border_tmp > right_border_MPOS)
					does_empty = true;

				if ((right_border_tmp >= left_border_MPOS) &&
					(right_border_tmp <= right_border_MPOS))
					right_border_MPOS = right_border_tmp;
				else if (right_border_tmp < left_border_MPOS)
					does_empty = true;
			}
			else {
				if ((left_border_tmp >= -e) &&
					left_border_tmp <= left_border_MNEG)
					left_border_MNEG = left_border_tmp;
				else if (left_border_tmp < -e) left_border_MNEG = -e;

				if ((right_border_tmp <= e) &&
					(right_border_tmp >= right_border_MNEG))
					right_border_MNEG = right_border_tmp;
				else if (right_border_tmp > e) right_border_MNEG = e;

				if ((left_border_MNEG == 0) && (right_border_MPOS))
					does_empty = true;
			}
		}
	}

	if ((left_border_MNEG < left_border_MPOS) && (right_border_MPOS < right_border_MNEG))
		does_empty = true;

	if (does_empty)
		hyp.set_characteristic(std::numeric_limits<double>::max());
	else {
		if (mixed_lipshEvaluations[0] > 0) {
			double apex = F_evaluations[id_A] - F_evaluations[id_B];
			apex = 0.5 * apex / (mixed_lipshEvaluations[0] * e);

			if ((left_border_MNEG == e) &&
				(left_border_MPOS == -e) &&
				(right_border_MNEG == -e) &&
				(right_border_MPOS == e)) {

				charact = -F_evaluations[id_A] * (apex - e) * 0.5 / e;
				charact = charact + F_evaluations[id_B] * (apex + e) * 0.5 / e;
				charact = charact + 0.5 * mixed_lipshEvaluations[0] * (apex * apex - e * e);
			}
			else if ((apex >= left_border_MPOS) && (apex <= left_border_MNEG) ||
					 (apex >= right_border_MNEG) && (apex <= right_border_MPOS))
				apex = apex;
			else
				apex = std::min({ 
					left_border_MPOS + nu * (left_border_MNEG - left_border_MPOS),
					left_border_MNEG - nu * (left_border_MNEG - left_border_MPOS),
					right_border_MNEG + nu * (right_border_MPOS - right_border_MNEG),
					right_border_MPOS - nu * (right_border_MPOS - right_border_MNEG) 
					});

			charact = -F_evaluations[id_A] * (apex - e) * 0.5 / e;
			charact = charact + F_evaluations[id_B] * (apex + e) * 0.5 / e;
			charact = charact + 0.5 * mixed_lipshEvaluations[0] * (apex * apex - e * e);
		}
		else {
			double apex = 0.0;
			if ((left_border_MNEG == e) &&
				(left_border_MPOS == -e) &&
				(right_border_MNEG == -e) &&
				(right_border_MPOS == e)) {
				apex = 
					std::min({
					left_border_MPOS + nu * (left_border_MNEG - left_border_MPOS),
					right_border_MPOS - nu * (right_border_MPOS - right_border_MNEG)
					});
			} 
			else apex = 
					std::min({
					left_border_MPOS + nu * (left_border_MNEG - left_border_MPOS),
					left_border_MNEG - nu * (left_border_MNEG - left_border_MPOS),
					right_border_MNEG + nu * (right_border_MPOS - right_border_MNEG),
					right_border_MPOS - nu * (right_border_MPOS - right_border_MNEG)
					});

			charact = -F_evaluations[id_A] * (apex - e) * 0.5 / e;
			charact = charact + F_evaluations[id_B] * (apex + e) * 0.5 / e;
			charact = charact + 0.5 * mixed_lipshEvaluations[0] * (apex * apex - e * e);
		}
		hyp.set_characteristic(charact);
	}
}

uint TSimplePMwithoutSM::choose_optimal_to_trisect() {
	uint id_optimal_hyp = 0;

	check_ifAllCharAreInfty();
	if (!do_AllCharAreInfty) {
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
	else {
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
	std::cout << " O: " << evals[0]
			  << " C: " << evals[1];

	std::cout << " (" << F_evaluations[hyp.get_idEvaluationsA()] << ';'
					  << F_evaluations[hyp.get_idEvaluationsB()] << ')' << std::endl;

	return id_optimal_hyp;
}

void TSimplePMwithoutSM::check_ifAllCharAreInfty() {
	do_AllCharAreInfty = true;
	for (uint i = 0; i < F_generated_intervals; ++i)
		if (F_intervals[i].get_characteristic() < std::numeric_limits<double>::max())
			do_AllCharAreInfty = false;
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

	/*double e1 = 0.0;
	double e2 = 0.0;

	for (uint i = 0; i < F_dimension; ++i) {
		F_decodeA[i] = F_coords[hyp2.get_idA() + i];
		F_decodeB[i] = F_coords[hyp2.get_idB() + i];
	}

	F_decodedA = Fp.decode_coordinates(F_decodeA);
	F_decodedB = Fp.decode_coordinates(F_decodeB);

	for (uint i = 0; i < F_dimension; ++i)
		e1 += ((F_decodedB[i] - F_decodedA[i]) * (F_decodedB[i] - F_decodedA[i]));

	e1 = sqrt(e1);

	for (uint i = 0; i < F_dimension; ++i) {
		F_decodeA[i] = F_coords[hyp1.get_idA() + i];
		F_decodeB[i] = F_coords[hyp3.get_idB() + i];
	}

	F_decodedA = Fp.decode_coordinates(F_decodeA);
	F_decodedB = Fp.decode_coordinates(F_decodeB);

	for (uint i = 0; i < F_dimension; ++i)
		e2 += ((F_decodedB[i] - F_decodedA[i]) * (F_decodedB[i] - F_decodedA[i]));

	e2 = sqrt(e2); */

	for (uint i = 0; i < F_constraints + 1; ++i) {
		new_llcv[i] = F_evaluations[hyp1.get_idEvaluationsA() + i] + F_evaluations[hyp3.get_idEvaluationsB() + i];
		new_llcv[i] -= F_evaluations[hyp2.get_idEvaluationsA() + i] + F_evaluations[hyp2.get_idEvaluationsB() + i];
		new_llcv[i] = 0.5 * std::abs(new_llcv[i] / (e2 * e2 - e1 * e1));
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
	else {
		compute_characteristic(id_Hyp);
		compute_characteristic(F_generated_intervals - 2);
		compute_characteristic(F_generated_intervals - 1);
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
			(F_iter < 300) && (F_intervals[id_current_interval].get_diagonal() > eps * F_intervals[0].get_diagonal());
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
