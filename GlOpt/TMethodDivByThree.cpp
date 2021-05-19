#include <iostream>
#include "TMethodDivByThree.h"

uint TPoint::F_dimension;
uint TPoint::F_constraints;

TMethodDivByThree::TMethodDivByThree(const uint& out_dim, const uint& out_constr, const uint& depth,
	TProblem& out_prob, const GainLipshConstant& out_gainObj, const GainLipshConstant& out_gainCst, const double& beta) :
	F_dimension(out_dim), F_queueDepth(depth), F_constraints(out_constr), Fp(out_prob), F_gainObjective(out_gainObj),
	F_gainConstraints(out_gainCst), delta(0.0000000001) {
	F_generated_points = 0;
	F_generated_intervals = 0;
	F_divide_axis = 0;
	F_current_minimum = 0.0;
	F_id_current_minimum = 0;
	transit_coord_1.resize(F_dimension);
	transit_coord_2.resize(F_dimension);
	F_globalLipshEvaluations.resize(F_constraints + 1);
	for (auto& elem : F_globalLipshEvaluations) elem = 0.0;

	F_criticalSize = beta * sqrt(F_dimension * (CoordinateValue)MAX_POWER_THREE * (CoordinateValue)MAX_POWER_THREE);

	TPoint::F_dimension = out_dim;
	TPoint::F_constraints = out_constr;
}

void TMethodDivByThree::initialization() {
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
	compute_diagonal(F_intervals[0].get_idThis());
	compute_characteristic(F_intervals[0].get_idThis());
	F_intervals[0].init_queues();
	compute_localLipshConst(F_intervals[0].get_idThis());
}

bool TMethodDivByThree::trisect_interval(const uint& id_divHyp) {
	resize_points_deque();
	resize_coords_deque();
	THyperinterval& div = F_intervals[id_divHyp];
	TPoint& point_a = F_points[div.get_idPointA()];
	TPoint& point_b = F_points[div.get_idPointB()];
	uint pos_a = div.get_idA();
	uint pos_b = div.get_idB();

	F_divide_axis = div.get_div_axis(); // по какой оси будем разбивать гиперинтервал?

	// считываем координаты из базы
	for (uint i = 0; i < F_dimension; ++i) {
		transit_coord_1[i] = F_coords[pos_a + i];
		transit_coord_2[i] = F_coords[pos_b + i];
	}

	uint pos = div.get_div_tag();
	EncodedCoordinate new_coord_u = transit_coord_2[F_divide_axis] - HYPER_INTERVAL_DOUBLE_SIDE_LENGTHS[pos]; // новая координата порождённой точки u по делимой оси
	TPoint::direction direction_u = TPoint::direction::DECREASE;
	EncodedCoordinate new_coord_v = transit_coord_1[F_divide_axis] + HYPER_INTERVAL_DOUBLE_SIDE_LENGTHS[pos]; // новая координата порождённой точки v по делимой оси
	TPoint::direction direction_v = TPoint::direction::INCREASE;

	if (transit_coord_1[F_divide_axis] > transit_coord_2[F_divide_axis]) {
		new_coord_u = transit_coord_2[F_divide_axis] + HYPER_INTERVAL_DOUBLE_SIDE_LENGTHS[pos];
		direction_u = TPoint::direction::INCREASE;
		new_coord_v = transit_coord_1[F_divide_axis] - HYPER_INTERVAL_DOUBLE_SIDE_LENGTHS[pos];
		direction_v = TPoint::direction::DECREASE;
	}

	// порождаем точку u от точки b
	uint new_id_u = point_b.does_point_exist(new_coord_u, direction_u, F_divide_axis, F_coords);
	// если точка не нашлась, то порождаем новую
	if (new_id_u == 0) {
		new_id_u = get_new_id();
		F_points[new_id_u].F_idThis = new_id_u;
		for (uint i = 0; i < F_dimension; ++i)
			F_coords[F_points[new_id_u].get_id_coord() + i] = transit_coord_2[i];

		F_coords[F_points[new_id_u].get_id_coord() + F_divide_axis] = new_coord_u;
		compute_evaluations(new_id_u);
		point_b.connect_points(new_id_u, F_divide_axis, direction_u);
	}

	// порождаем точку v от точки a
	uint new_id_v = point_a.does_point_exist(new_coord_v, direction_v, F_divide_axis, F_coords);
	// если точка не нашлась, то порождаем новую
	if (new_id_v == 0) {
		new_id_v = get_new_id();
		F_points[new_id_v].F_idThis = new_id_v;
		for (uint i = 0; i < F_dimension; ++i)
			F_coords[F_points[new_id_v].get_id_coord() + i] = transit_coord_1[i];

		F_coords[F_points[new_id_v].get_id_coord() + F_divide_axis] = new_coord_v;
		compute_evaluations(new_id_v);
		point_a.connect_points(new_id_v, F_divide_axis, direction_v);
	}

	div.increase_division();
	fill_intervals(div, new_id_u, new_id_v);
	return true;
}

void TMethodDivByThree::fill_intervals(THyperinterval& parent, const uint& id_u, const uint& id_v) {
	TPoint& point_u = F_points[id_u];
	TPoint& point_v = F_points[id_v];
	resize_intervals_deque();

	uint pos_hyp_2 = get_new_interval();
	uint pos_hyp_3 = get_new_interval();

	F_intervals[pos_hyp_2] = parent;
	F_intervals[pos_hyp_3] = parent;
	THyperinterval& new_hyp_2 = F_intervals[pos_hyp_2];
	THyperinterval& new_hyp_3 = F_intervals[pos_hyp_3];

	parent.set_idPointB(id_u);
	compute_diagonal(parent.get_idThis());
	compute_localLipshConst(parent.get_idThis());
	compute_characteristic(parent.get_idThis());

	new_hyp_2.set_idThis(pos_hyp_2);
	new_hyp_2.set_idPointA(id_u);
	new_hyp_2.set_idPointB(id_v);
	compute_diagonal(new_hyp_2.get_idThis());
	compute_localLipshConst(new_hyp_2.get_idThis());
	compute_characteristic(new_hyp_2.get_idThis());

	new_hyp_3.set_idThis(pos_hyp_3);
	new_hyp_3.set_idPointA(id_v);
	compute_diagonal(new_hyp_3.get_idThis());
	compute_localLipshConst(new_hyp_3.get_idThis());
	compute_characteristic(new_hyp_3.get_idThis());
}

void TMethodDivByThree::compute_diagonal(const uint& id_Hyp) {
	THyperinterval& hyp = F_intervals[id_Hyp];
	uint pos_coord_a = hyp.get_idA();
	uint pos_coord_b = hyp.get_idB();

	double diagonal = 0.0;
	double temp = 0.0;

	for (uint i = 0; i < F_dimension; ++i) {
		transit_coord_1[i] = F_coords[pos_coord_a + i];
		transit_coord_2[i] = F_coords[pos_coord_b + i];
	}

	CoordinatesValues decoded_coord_a = Fp.decode_coordinates(transit_coord_1);
	CoordinatesValues decoded_coord_b = Fp.decode_coordinates(transit_coord_2);

	for (uint i = 0; i < F_dimension; ++i) {
		temp = decoded_coord_b[i] - decoded_coord_a[i];
		diagonal = diagonal + temp * temp;
	}

	diagonal = sqrt(diagonal);
	hyp.set_diagonal(diagonal);
}

void TMethodDivByThree::compute_characteristic(const uint& id_Hyp) {
	THyperinterval& hyp = F_intervals[id_Hyp];
	hyp.set_characteristic(hyp.get_diagonal());
}

void TMethodDivByThree::compute_evaluations(const uint& out_idPoint) {
	resize_evaluations_deque();
	TPoint& point = F_points[out_idPoint];
	uint pos = point.get_id_coord();
	for (uint i = 0; i < F_dimension; ++i)
		transit_coord_2[i] = F_coords[pos + i];

	FunctionsValues& evals = Fp(transit_coord_2);

	for (uint i = 0; i < F_constraints + 1; ++i)
		F_evaluations[point.F_idThis * (F_constraints + 1) + i] = evals[i];
} 

void TMethodDivByThree::compute_localLipshConst(const uint& id_Hyp) {
	std::vector<LipschitzConstantValue> new_llcv;
	new_llcv.resize(F_constraints + 1);
	THyperinterval& hyp = F_intervals[id_Hyp];

	for (uint i = 0; i < F_constraints + 1; ++i) {
		new_llcv[i] = F_evaluations[hyp.get_idEvaluationsA() + i] - F_evaluations[hyp.get_idEvaluationsB() + i];
		new_llcv[i] = std::abs(new_llcv[i]) / hyp.get_diagonal();
	}

	hyp.update_queuesLipshEvaluations(new_llcv, delta);
}

void TMethodDivByThree::compute_globalLipshConst() {
	for (uint i = 0; i < F_generated_intervals; ++i) {
		const std::vector<double>& max_i = F_intervals[i].get_maxLipshEvaluations();
		for (uint j = 0; j < F_constraints + 1; ++j)
			if (max_i[j] > F_globalLipshEvaluations[j])
				F_globalLipshEvaluations[j] = max_i[j];
	}
}

uint TMethodDivByThree::choose_optimal_to_trisect() {
	uint id_optimal_hyp = 0;
	double optimal_charact = F_intervals[id_optimal_hyp].get_characteristic();
	double current_charact = 0.0;

	for (uint id_hyp = 1; id_hyp < F_generated_intervals; ++id_hyp) {
		current_charact = F_intervals[id_hyp].get_characteristic();
		if (std::abs(current_charact - optimal_charact) < std::numeric_limits<double>::epsilon());
		else if (optimal_charact < current_charact) {
			optimal_charact = current_charact;
			id_optimal_hyp = id_hyp;
		}
	}

	return id_optimal_hyp;
}

uint TMethodDivByThree::do_step(const uint& id_divHyp) {
	trisect_interval(id_divHyp);
	return choose_optimal_to_trisect();
}

void TMethodDivByThree::launch_method() {
	initialization();
	uint id_current_interval = 0;
	for (uint i = 0; i < 500; ++i) {
		id_current_interval = do_step(id_current_interval);
		std::cout << id_current_interval << std::endl;
	}
}