#include "SimplePMwithSM.h"

void SimplePMwithSM::calculate_localLipshConst(const uint& id_hyp) {
	Hyperinterval& hyp1 = _intervals[id_hyp];
	Hyperinterval& hyp2 = _intervals[_generated_intervals - 1];
	Hyperinterval& hyp3 = _intervals[_generated_intervals - 2];
}

uint SimplePMwithSM::optimal_to_trisect() {
	return 0;
}

uint SimplePMwithSM::iterate(const uint& id_hyp) {
	trisect_interval(id_hyp);
	calculate_localLipshConst(id_hyp);
	calculate_characteristic(id_hyp);
	calculate_characteristic(_generated_intervals - 1);
	calculate_characteristic(_generated_intervals - 2);
	return optimal_to_trisect();
}