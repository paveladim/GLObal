#include "Uniform.h"

Uniform::Uniform(const uint& dimension,
				 const uint& constraints,
				 Parameters& parameters,
				 Problem& problem) : 
	DivideByThree(dimension, constraints, parameters, problem) {}

void Uniform::calculate_characteristic(const uint& id_hyp){
	_intervals[id_hyp].set_charact(_intervals[id_hyp].get_diagonal());
}

uint Uniform::optimal_to_trisect() {
	uint optimal = 0;
	for (uint i = 1; i < _generated_intervals; ++i)
		if (_intervals[i].get_diagonal() > _intervals[optimal].get_diagonal())
			optimal = i;

	return optimal;
}
uint Uniform::iterate(const uint& id_hyp) {
	trisect_interval(id_hyp);
	calculate_characteristic(id_hyp);
	return optimal_to_trisect();
}

void Uniform::solve() {
	initialization();
	uint id_current_interval = 0;
	for (uint i = 0; i < 5; ++i)
		id_current_interval = iterate(id_current_interval);
}