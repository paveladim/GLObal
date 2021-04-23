#include "TMethodDivByThree.h"

TMethodDivByThree::TMethodDivByThree(const uint& out_dim, const uint& out_constr, const uint& depth, const double& out_eps,
	const uint& max_gen_points, const uint& max_gen_interv, TProblem& out_prob) : F_dimension(out_dim), 
	F_queueDepth(depth), F_constraints(out_constr), F_eps(out_eps), max_generated_points(max_gen_points), 
	max_generated_intervals(max_gen_interv), Fp(out_prob) {
	F_generated_points = 0;
	F_generated_intervals = 0;
	F_divide_axis = 0;
	F_current_minimum = 0.0;
	coord_a.resize(F_dimension);
	coord_b.resize(F_dimension);
}

void TMethodDivByThree::createFirstInterval() {
	THyperinterval::init_static(F_dimension, F_constraints, F_queueDepth);
	TPoint point_a;
	TPoint point_b;
	point_a.dec_coords.resize(F_dimension);
	point_a.inc_coords.resize(F_dimension);
	point_b.dec_coords.resize(F_dimension);
	point_b.inc_coords.resize(F_dimension);
	point_a.F_idThis = get_new_id();
	point_b.F_idThis = get_new_id();
	point_a.F_idCoords = point_a.F_idThis * F_dimension;
	point_b.F_idCoords = point_b.F_idThis * F_dimension;


	for (uint i = 0; i < F_dimension; ++i)
		F_coords.push_back(0);

	for (uint i = 0; i < F_dimension; ++i)
		F_coords.push_back(MAX_POWER_THREE);

	F_points.push_back(point_a);
	F_points.push_back(point_b);

	THyperinterval first;
	first.set_idPointA(point_a.F_idThis);
	first.set_idPointB(point_b.F_idThis);
	first.set_idA(point_a.F_idCoords);
	first.set_idB(point_b.F_idCoords);
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
		coord_a[i] = F_coords[pos_a + i];
		coord_b[i] = F_coords[pos_b + i];
	}

	if (div.increase_div_tag(F_divide_axis) == false) return false;
	uint pos = div.get_div_tag(F_divide_axis);
	TPoint point_u;
	point_u.F_idThis = get_new_id();
	point_u.F_idCoords = get_id_coord();
	for (uint i = 0; i < F_dimension; ++i)
		F_coords.push_back(coord_b[i]);
	F_coords[point_u.F_idCoords + F_divide_axis] = HYPER_INTERVAL_SIDE_LENGTHS[MAX_EXPONENT_THREE - pos];

	TPoint point_v;
	point_v.F_idThis = get_new_id();
	point_v.F_idCoords = get_id_coord();
	for (uint i = 0; i < F_dimension; ++i)
		F_coords.push_back(coord_a[i]);
	F_coords[point_v.F_idCoords + F_divide_axis] = HYPER_INTERVAL_DOUBLE_SIDE_LENGTHS[MAX_EXPONENT_THREE - pos];

	F_divide_axis = (++F_divide_axis) % F_dimension;
	return true;
}

uint TMethodDivByThree::get_new_id() {
	if (F_generated_points < max_generated_points)
		return F_generated_points++;
	else
		return -1;
}

uint TMethodDivByThree::get_id_coord() {
	return F_coords.size();
}

uint TMethodDivByThree::get_new_interval() {
	if (F_generated_intervals < max_generated_intervals)
		return F_generated_intervals++;
	else
		return -1;
}

void TMethodDivByThree::delete_point() {
	if (F_generated_points > 0) --F_generated_points;
}

void TMethodDivByThree::delete_hyperinterval() {
	if (F_generated_intervals > 0) --F_generated_intervals;
}

uint TMethodDivByThree::does_point_exist(TPoint& parent, const uint& des_val, bool f) {
	if (f) {
		auto beg = parent.inc_coords[F_divide_axis].begin();
		auto end = parent.inc_coords[F_divide_axis].end();
		for (auto it = beg; it != end; ++it)
		{
			TPoint& point = F_points[*it];
			uint pos = point.F_idCoords + F_divide_axis;
			if (F_coords[pos] == des_val) return *it;
			if (F_coords[pos] < des_val) {
				parent.inc_coords[F_divide_axis].insert(it, get_new_id());
				return -1;
			}
		}
	}
	else {
		auto beg = parent.dec_coords[F_divide_axis].begin();
		auto end = parent.dec_coords[F_divide_axis].end();
		for (auto it = beg; it != end; ++it)
		{
			TPoint& point = F_points[*it];
			uint pos = point.F_idCoords + F_divide_axis;
			if (F_coords[pos] == des_val) return *it;
		}
	}
	return get_new_id();
}