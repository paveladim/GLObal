#include "Point.h"

Point::Point() : _id(0),
				 _inc_coord(_dimension),
				 _dec_coord(_dimension) {}

uint Point::get_id() const {
	return _id;
}

uint Point::get_id_coord() const {
	return _id * _dimension;
}

uint Point::get_id_evaluations() const {
	return _id * (_constraints + 1);
}

void Point::set_id(const uint& id) {
	_id = id;
}

uint Point::does_point_exist(const uint& value,
							 const Direction& dir,
							 const uint& div_axis,
							 const std::deque<uint>& coords) {

}

void Point::connect_points(const uint& child_id,
						   const uint& div_axis,
						   const Direction& dir) {
	if (dir == Direction::DECREASE)
		_dec_coord[div_axis].push_back(child_id);
	else
		_inc_coord[div_axis].push_back(child_id);
}