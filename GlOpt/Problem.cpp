#include "Problem.h"
#include "dicretization.h"

Problem::Problem(const uint& dimension,
				 const uint& constraints,
				 const CoordinatesValues& left_brd,
				 const CoordinatesValues& right_brd,
				 const FunctionsCalculator& problem) :
	_dimension(dimension),
	_constraints(constraints),
	_left_borders(left_brd),
	_right_borders(right_brd),
	_problem(problem),
	_decoded_point(dimension),
	_function_values(constraints + 1) {}

FunctionsValues& Problem::f(const CoordinatesValues& out) {
	return _problem(_function_values, out);
}

CoordinatesValues& Problem::decode_coordinates(const EncodedCoordinates& out) {
	for (uint i = 0; i < _dimension; ++i) {
		_decoded_point[i] = (double)out[i] / (double)MAX_POWER_THREE;
		_decoded_point[i] = _decoded_point[i] * 
							(_right_borders[i] - _left_borders[i]) + 
							_left_borders[i];
	}

	return _decoded_point;
}

FunctionsValues& Problem::operator()(const EncodedCoordinates& out) {
	return f(decode_coordinates(out));
}