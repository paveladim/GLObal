#include <iostream>
#include <fstream>

#include "DivideByThree.h"

DivideByThree::DivideByThree(const uint& dimension,
			  const uint& constraints,
			  Parameters& parameters,
			  Problem& problem) :
	_dimension(dimension),
	_constraints(constraints),
	_parameters(parameters),
	_problem(problem),
	_generated_points(0),
	_generated_intervals(0),
	_transit1(dimension),
	_transit2(dimension),
	_current_minimum(std::numeric_limits<double>::max()),
	_id_minimum(0),
	_divided_axis(0),
	_iteration(0) {
	Point::_dimension = _dimension;
	Point::_constraints = _constraints;
}

void DivideByThree::initialization() {
	Hyperinterval::init_static(_dimension, _constraints, _parameters._queueDepth);
	resize_points_deque();
	resize_coords_deque();

	_points[0].set_id(generate_id());
	_points[1].set_id(generate_id());

	Point& a = _points[0];
	Point& b = _points[1];

	for (uint i = 0; i < _dimension; ++i)
		_coords[a.get_id_coord() + i] = 0;

	for (uint i = 0; i < _dimension; ++i)
		_coords[b.get_id_coord() + i] = MAX_POWER_THREE;

	compute_evaluations(a.get_id());
	compute_evaluations(b.get_id());

	resize_intervals_deque();
	_intervals[0].set_idA(a.get_id());
	_intervals[0].set_idB(b.get_id());
	_intervals[0].set_id(generate_hyp());
	compute_diagonal(_intervals[0].get_id());
	calculate_characteristic(0);
}

void DivideByThree::trisect_interval(const uint& id_hyp) {
	Hyperinterval& div = _intervals[id_hyp];
	uint pos_a = div.get_coordA();
	uint pos_b = div.get_coordB();
	Point& point_a = _points[div.get_idA()];
	Point& point_b = _points[div.get_idB()];

	// по какой оси будем разбивать гиперинтервал?
	_divided_axis = div.get_axis();

	// считываем координаты из базы
	for (uint i = 0; i < _dimension; ++i) {
		_transit1[i] = _coords[pos_a + i];
		_transit2[i] = _coords[pos_b + i];
	}

	uint pos = div.get_shift();
	// новая координата порождённой точки u по делимой оси
	EncodedCoordinate new_coord_u = 
		_transit1[_divided_axis] + HYPER_INTERVAL_DOUBLE_SIDE_LENGTHS[pos];
	Point::Direction direction_u = Point::Direction::INCREASE;
	// новая координата порождённой точки v по делимой оси
	EncodedCoordinate new_coord_v = 
		_transit2[_divided_axis] - HYPER_INTERVAL_DOUBLE_SIDE_LENGTHS[pos];
	Point::Direction direction_v = Point::Direction::DECREASE;

	if (_transit1[_divided_axis] > _transit2[_divided_axis]) {
		new_coord_u = _transit1[_divided_axis] - 
					  HYPER_INTERVAL_DOUBLE_SIDE_LENGTHS[pos];
		direction_u = Point::Direction::DECREASE;
		new_coord_v = _transit2[_divided_axis] + 
					  HYPER_INTERVAL_DOUBLE_SIDE_LENGTHS[pos];
		direction_v = Point::Direction::INCREASE;
	}

	// порождаем точку v от точки b
	uint new_id_v = point_b.does_point_exist(new_coord_v, 
											 direction_v, 
											 _divided_axis, 
											 _coords);
	// если точка не нашлась, то порождаем новую
	if (new_id_v == 0) {
		new_id_v = generate_id();
		resize_points_deque();
		resize_coords_deque();
		_points[new_id_v].set_id(new_id_v);
		for (uint i = 0; i < _dimension; ++i)
			_coords[_points[new_id_v].get_id_coord() + i] = _transit2[i];

		_coords[_points[new_id_v].get_id_coord() + _divided_axis] = new_coord_v;
		compute_evaluations(new_id_v);
		point_b.connect_points(new_id_v, _divided_axis, direction_v);
	}

	// порождаем точку u от точки a
	uint new_id_u = point_a.does_point_exist(new_coord_u, 
											 direction_u, 
											 _divided_axis, 
											 _coords);
	// если точка не нашлась, то порождаем новую
	if (new_id_u == 0) {
		new_id_u = generate_id();
		resize_points_deque();
		resize_coords_deque();
		_points[new_id_u].set_id(new_id_u);
		for (uint i = 0; i < _dimension; ++i)
			_coords[_points[new_id_u].get_id_coord() + i] = _transit1[i];

		_coords[_points[new_id_u].get_id_coord() + _divided_axis] = new_coord_u;
		compute_evaluations(new_id_u);
		point_a.connect_points(new_id_u, _divided_axis, direction_u);
	}

	div.increase_divisions();
	fill_intervals(div, new_id_u, new_id_v);
}

void DivideByThree::fill_intervals(Hyperinterval& parent, 
								   const uint& id_u, 
								   const uint& id_v) {
	Point& point_u = _points[id_u];
	Point& point_v = _points[id_v];
	resize_intervals_deque();

	uint pos_hyp_2 = generate_hyp();
	uint pos_hyp_3 = generate_hyp();

	_intervals[pos_hyp_2] = parent;
	_intervals[pos_hyp_3] = parent;
	Hyperinterval& new_hyp_2 = _intervals[pos_hyp_2];
	Hyperinterval& new_hyp_3 = _intervals[pos_hyp_3];

	parent.set_idB(id_v);
	compute_diagonal(parent.get_id());

	new_hyp_2.set_id(pos_hyp_2);
	new_hyp_2.set_idA(id_v);
	new_hyp_2.set_idB(id_u);
	compute_diagonal(new_hyp_2.get_id());

	new_hyp_3.set_id(pos_hyp_3);
	new_hyp_3.set_idA(id_u);
	compute_diagonal(new_hyp_3.get_id());
}

uint DivideByThree::generate_id() {
	return _generated_points++;
}

uint DivideByThree::generate_hyp() {
	return _generated_intervals++;
}

void DivideByThree::compute_evaluations(const uint& id_point) {
	resize_evaluations_deque();
	Point& point = _points[id_point];
	uint pos = point.get_id_coord();
	for (uint i = 0; i < _dimension; ++i)
		_transit2[i] = _coords[pos + i];

	FunctionsValues& evals = _problem(_transit2);

	// обновление оценки минимума
	if (evals[0] < _current_minimum) {
		bool flag = true;
		for (uint i = 1; i < _constraints + 1; ++i)
			if (evals[i] > 0) flag = false;

		if (flag) {
			_current_minimum = evals[0];
			_id_minimum = id_point;
		}
	}

	for (uint i = 0; i < _constraints + 1; ++i)
		_evaluations[point.get_id_evaluations() + i] = evals[i];
}

void DivideByThree::compute_diagonal(const uint& id_hyp) {
	Hyperinterval& hyp = _intervals[id_hyp];
	uint pos_a = hyp.get_idA() * _dimension;
	uint pos_b = hyp.get_idB() * _dimension;

	double diagonal = 0.0;
	double temp = 0.0;

	for (uint i = 0; i < _dimension; ++i) {
		_transit1[i] = _coords[pos_a + i];
		_transit2[i] = _coords[pos_b + i];
	}

	for (uint i = 0; i < _dimension; ++i) {
		temp = (CoordinateValue)_transit1[i] - (CoordinateValue)_transit2[i];
		diagonal = diagonal + temp * temp;
	}

	diagonal = sqrt(diagonal);
	hyp.set_diagonal(diagonal);
}

void DivideByThree::resize_points_deque() {
    if (_points.size() - _generated_points < 1)
        _points.resize(_points.size() + 100);
}

void DivideByThree::resize_coords_deque() {
    if (_coords.size() - _generated_points * _dimension < _dimension)
        _coords.resize(_coords.size() + 100 * _dimension);
}

void DivideByThree::resize_evaluations_deque() {
    if ((_evaluations.size() - _generated_points * (_constraints + 1) < 
        (_constraints + 1)) || _evaluations.size() == 0)
        _evaluations.resize(_evaluations.size() + 100 * (_constraints + 1));
}

void DivideByThree::resize_intervals_deque() {
    if (_intervals.size() < _generated_intervals + 2)
        _intervals.resize(_intervals.size() + 100);
}

void DivideByThree::write_generated_points() {
	std::ofstream out;
	out.open("D:\\materials\\projects\\visual_hyperinterval\\points.txt");
	if (out.is_open())
	{
		for (uint i = 0; i < _generated_points * _dimension; ++i)
			out << _coords[i] << std::endl;
	}
}

void DivideByThree::write_generated_intervals() {
	std::ofstream out;
	out.open("D:\\materials\\projects\\visual_hyperinterval\\hyp.txt");
	if (out.is_open())
	{
		for (uint id_hyp = 0; id_hyp < _generated_intervals; ++id_hyp) {
			for (uint i = 0; i < _dimension; ++i)
				out << _coords[_intervals[id_hyp].get_coordA() + i] << std::endl;
			for (uint i = 0; i < _dimension; ++i)
				out << _coords[_intervals[id_hyp].get_coordB() + i] << std::endl;
		}
	}
}