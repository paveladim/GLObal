#ifndef PROBLEM_H
#define PROBLEM_H

#include "synonymous_types.h"

class Problem {
	// размерность задачи и число ограничений
	uint _dimension;
	uint _constraints;
	// область поиска
	CoordinatesValues _left_borders;
	CoordinatesValues _right_borders;
	// для вычисления значений целевой функции и ограничений
	FunctionsCalculator _problem;
	// вектор для хранения расшифрованных координат
	CoordinatesValues _decoded_point;
	// 
	FunctionsValues _function_values;
private:
	// вычислить функцию в точке для вещественных координат
	FunctionsValues& f(const CoordinatesValues& out);
public:
	Problem() {}
	Problem(const uint& dimension, 
			const uint& constraints,
			const CoordinatesValues& left_brd,
			const CoordinatesValues& right_brd,
			const FunctionsCalculator& problem);
	// расшифровать координаты
	CoordinatesValues& decode_coordinates(const EncodedCoordinates& out);
	// вычислить значений целевой функции и ограничений
	FunctionsValues& operator()(const EncodedCoordinates& out);
	~Problem() {}
};

#endif // PROBLEM_H
