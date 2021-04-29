#include "TMethodDivByThree.h"

uint TPoint::F_dimension;
uint TPoint::F_constraints;

TMethodDivByThree::TMethodDivByThree(const uint& out_dim, const uint& out_constr, const uint& depth, const double& out_eps, TProblem& out_prob) : 
	F_dimension(out_dim), 
	F_queueDepth(depth), F_constraints(out_constr), F_eps(out_eps), Fp(out_prob) {
	F_generated_points = 0;
	F_generated_intervals = 0;
	F_divide_axis = 0;
	F_current_minimum = 0.0;
	new_coord_u.resize(F_dimension);
	new_coord_v.resize(F_dimension);

	TPoint::F_dimension = out_dim;
	TPoint::F_constraints = out_constr;
}

void TMethodDivByThree::initialization() {
	THyperinterval::init_static(F_dimension, F_constraints, F_queueDepth);
	TPoint point_a;
	TPoint point_b;
	point_a.F_idThis = get_new_id();
	point_b.F_idThis = get_new_id();

	for (uint i = 0; i < F_dimension; ++i)
		F_coords.push_back(0);

	for (uint i = 0; i < F_dimension; ++i)
		F_coords.push_back(MAX_POWER_THREE);

	F_points.push_back(point_a);
	F_points.push_back(point_b);

	THyperinterval first;
	first.set_idPointA(point_a.F_idThis);
	first.set_idPointB(point_b.F_idThis);
	first.set_idA(point_a.get_id_coord());
	first.set_idB(point_b.get_id_coord());
	first.set_idThis(get_new_interval());
	F_intervals.push_back(first);
}

bool TMethodDivByThree::divideInterval(const uint& id_divHyp) {
	THyperinterval& div = F_intervals[id_divHyp];
	TPoint& point_a = F_points[div.get_idPointA()];
	TPoint& point_b = F_points[div.get_idPointB()];
	uint pos_a = div.get_idA();
	uint pos_b = div.get_idB();

	// считываем координаты из базы
	for (uint i = 0; i < F_dimension; ++i) {
		new_coord_u[i] = F_coords[pos_a + i];
		new_coord_v[i] = F_coords[pos_b + i];
	}

	uint pos = div.get_div_tag();

	// порождаем точку u от точки b
	uint new_id_u = point_b.does_point_exist(HYPER_INTERVAL_SIDE_LENGTHS[pos], TPoint::direction::BACKWARD, F_divide_axis, F_coords);
	// если точка не нашлась, то порождаем новую
	if (new_id_u == 0) new_id_u = get_new_id();
	if (new_id_u == F_points.size()) {
		TPoint point_u;
		point_u.F_idThis = new_id_u;
		for (uint i = 0; i < F_dimension; ++i)
			F_coords.push_back(new_coord_u[i]);
		F_coords[point_u.get_id_coord() + F_divide_axis] = HYPER_INTERVAL_SIDE_LENGTHS[pos];
		F_points.push_back(point_u);
	}

	// порождаем точку v от точки a
	uint new_id_v = point_a.does_point_exist(HYPER_INTERVAL_SIDE_LENGTHS[pos], TPoint::direction::FORWARD, F_divide_axis, F_coords);
	if (new_id_v == 0) new_id_v = get_new_id();
	if (new_id_v == F_points.size()) {
		TPoint point_v;
		point_v.F_idThis = new_id_v;
		for (uint i = 0; i < F_dimension; ++i)
			F_coords.push_back(new_coord_v[i]);
		F_coords[point_v.get_id_coord() + F_divide_axis] = HYPER_INTERVAL_DOUBLE_SIDE_LENGTHS[pos];
		F_points.push_back(point_v);
	}

	div.increase_division();
	fillIntervals(div, new_id_u, new_id_v);
	F_divide_axis = (++F_divide_axis) % F_dimension;
	return true;
}

void TMethodDivByThree::fillIntervals(THyperinterval& parent, const uint& id_u, const uint& id_v) {
	TPoint& point_u = F_points[id_u];
	TPoint& point_v = F_points[id_v];

	THyperinterval new_hyp_2 = parent;
	THyperinterval new_hyp_3 = parent;

	parent.set_idPointB(id_u);
	parent.set_idB(point_u.get_id_coord());
	parent.set_idEvaluationsB(point_u.get_id_evaluations());
	compute_diagonal(parent.get_idThis());

	new_hyp_2.set_idThis(get_new_interval());
	new_hyp_2.set_idPointA(id_u);
	new_hyp_2.set_idPointB(id_v);
	new_hyp_2.set_idA(point_u.get_id_coord());
	new_hyp_2.set_idB(point_v.get_id_coord());
	new_hyp_2.set_idEvaluationsA(point_u.get_id_evaluations());
	new_hyp_2.set_idEvaluationsB(point_v.get_id_evaluations());
	F_intervals.push_back(new_hyp_2);
	compute_diagonal(new_hyp_2.get_idThis());

	new_hyp_3.set_idThis(get_new_interval());
	new_hyp_3.set_idPointA(id_v);
	new_hyp_3.set_idA(point_v.get_id_coord());
	new_hyp_3.set_idEvaluationsA(point_v.get_id_evaluations());
	F_intervals.push_back(new_hyp_3);
	compute_diagonal(new_hyp_3.get_idThis());
}

void TMethodDivByThree::compute_diagonal(const uint& id_Hyp) {
	THyperinterval& hyp = F_intervals[id_Hyp];
	uint pos_coord_a = hyp.get_idA();
	uint pos_coord_b = hyp.get_idB();

	double diagonal = 0.0;
	double temp;
	for (uint i = 0; i < F_dimension; ++i) {
		temp = (double)F_coords[pos_coord_b + i] - (double)F_coords[pos_coord_a + i];
		temp = fabs(temp);
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
	TPoint point = F_points[out_idPoint];
	uint pos = point.get_id_coord();
	for (uint i = 0; i < F_dimension; ++i)
		new_coord_u[i] = F_coords[pos + i];

	FunctionsValues& evals = Fp(new_coord_u);

	for (uint i = 0; i < F_constraints + 1; ++i)
		F_evaluations.push_back(evals[i]);
} 

/* uint TMethodDivByThree::does_point_exist(TPoint& parent, const uint& des_val, const direction& dir) {
	if (dir == direction::FORWARD) {
		auto beg = parent.inc_coords[F_divide_axis].begin();
		auto end = parent.inc_coords[F_divide_axis].end();
		for (auto it = beg; it != end; ++it)
		{
			TPoint& point = F_points[*it];
			uint pos = point.get_id_coord() + F_divide_axis;
			if (F_coords[pos] == des_val) return *it;
			if (F_coords[pos] < des_val) {
				uint new_id = get_new_id();
				parent.inc_coords[F_divide_axis].insert(it, new_id);
				return new_id;
			}
		}

		if (beg == end) {
			uint new_id = get_new_id();
			parent.inc_coords[F_divide_axis].push_back(new_id);
			return new_id;
		}
	}
	else {
		auto beg = parent.dec_coords[F_divide_axis].begin();
		auto end = parent.dec_coords[F_divide_axis].end();
		for (auto it = beg; it != end; ++it)
		{
			TPoint& point = F_points[*it];
			uint pos = point.get_id_coord() + F_divide_axis;
			if (F_coords[pos] == des_val) return *it;
			if (F_coords[pos] < des_val) {
				uint new_id = get_new_id();
				parent.inc_coords[F_divide_axis].insert(it, new_id);
				return new_id;
			}
		}

		if (beg == end) {
			uint new_id = get_new_id();
			parent.inc_coords[F_divide_axis].push_back(new_id);
			return new_id;
		}
	}
} */