#include "TransformPMwithConj.h"

TransformPMwithConj::TransformPMwithConj(const uint& dimension,
										 const uint& constraints,
										 Parameters& parameters,
										 Problem& problem) :
	TransformPM(dimension, constraints, parameters, problem) {}

void TransformPMwithConj::calculate_localLipshConst(const uint& id_hyp) {
	Hyperinterval& hyp1 = _intervals[id_hyp];
	Hyperinterval& hyp2 = _intervals[_generated_intervals - 2];
	Hyperinterval& hyp3 = _intervals[_generated_intervals - 1];

	uint index = 20 - (hyp1.get_divisions() - 1) / _dimension;
	uint j = hyp1.get_previous_axis();
	double h2 = HYPER_INTERVAL_SIDE_LENGTHS[index];
	double h1 = sqrt(j * h2 * h2 / 9.0 + (_dimension - j - 1) * h2 * h2);

	double e1 = 0.5 * sqrt(h1 * h1 + h2 * h2 / 9.0);
	double e2 = 0.5 * sqrt(h1 * h1 + h2 * h2);

	for (uint i = 0; i < _constraints + 1; ++i) {
		_localLipshEval[i] = _evaluations[hyp1.get_evalA() + i]
			+ _evaluations[hyp3.get_evalB() + i];
		_localLipshEval[i] -= (_evaluations[hyp2.get_evalA() + i]
			+ _evaluations[hyp2.get_evalB() + i]);

		if (i == 0)
			_localLipshEval[i] = std::abs(_localLipshEval[i] / (e2 * e2 - e1 * e1));
		else
			_localLipshEval[i] = std::abs(_localLipshEval[i] / (e2 * e2 - e1 * e1));
	}

	hyp1.update_localLipQueues(_localLipshEval,
		_parameters._delta /
		((double)MAX_POWER_THREE * (double)MAX_POWER_THREE));
	hyp2.update_localLipQueues(_localLipshEval,
		_parameters._delta /
		((double)MAX_POWER_THREE * (double)MAX_POWER_THREE));
	hyp3.update_localLipQueues(_localLipshEval,
		_parameters._delta /
		((double)MAX_POWER_THREE * (double)MAX_POWER_THREE));
}

void TransformPMwithConj::calculate_characteristic(const uint& id_hyp) {
	if (_intervals[id_hyp].get_divisions() == 0) return;
	Hyperinterval& hyp = _intervals[id_hyp];

	uint index = 20 - (hyp.get_divisions() - 1) / _dimension;
	uint j = hyp.get_previous_axis();
	double h2 = HYPER_INTERVAL_SIDE_LENGTHS[index];
	double h1 = sqrt(j * h2 * h2 / 9.0 + (_dimension - j - 1) * h2 * h2);
	double e = 0.5 * sqrt(h1 * h1 + h2 * h2 / 9.0);

	hyp.set_charact(golden_ratio(-e, e, e, id_hyp));
}