#include <iostream>
#include "SimplePMwithConj.h"

SimplePMwithConj::SimplePMwithConj(const uint& dimension,
								   const uint& constraints,
								   Parameters& parameters,
							 	   Problem& problem) :
	SimplePM(dimension, constraints, parameters, problem) {}

void SimplePMwithConj::calculate_localLipshConst(const uint& id_hyp) {
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

void SimplePMwithConj::calculate_characteristic(const uint& id_hyp) {
	if (_intervals[id_hyp].get_divisions() == 0) return;
	Hyperinterval& hyp = _intervals[id_hyp];

	double mixed_LipshEval = mixedLipEval(hyp, 0);
	balance(mixed_LipshEval);
	double t_min = 0.5 * (_evaluations[hyp.get_evalA()] -
						  _evaluations[hyp.get_evalB()]);
	t_min = t_min / mixed_LipshEval;

	uint index = 20 - (hyp.get_divisions() - 1) / _dimension;
	uint j = hyp.get_previous_axis();
	double h2 = HYPER_INTERVAL_SIDE_LENGTHS[index];
	double h1 = sqrt(j * h2 * h2 / 9.0 + (_dimension - j - 1) * h2 * h2);
	double e = 0.5 * sqrt(h1 * h1 + h2 * h2 / 9.0);

	t_min = t_min / e;
	double left = -e;
	double right = e;
	give_borders(left, right, hyp);

	if ((left == e) && (right == -e)) {
		hyp.set_charact(std::numeric_limits<double>::max());
	}
	else {
		if ((t_min > left) && (t_min < right)) {
			double charact = -_evaluations[hyp.get_evalA()] * (t_min - e);
			charact = charact + _evaluations[hyp.get_evalB()] * (t_min + e);
			charact = 0.5 * charact / e;
			charact = charact + 0.5 * mixed_LipshEval * (t_min * t_min - e * e);
			hyp.set_charact(charact);
		}
		else {
			double charact1 = -_evaluations[hyp.get_evalA()] * (left - e);
			charact1 = charact1 + _evaluations[hyp.get_evalB()] * (left + e);
			charact1 = 0.5 * charact1 / e;
			charact1 = charact1 + 0.5 * mixed_LipshEval * (left * left - e * e);

			double charact2 = -_evaluations[hyp.get_evalA()] * (right - e);
			charact2 = charact2 + _evaluations[hyp.get_evalB()] * (right + e);
			charact2 = 0.5 * charact2 / e;
			charact2 = charact2 + 0.5 * mixed_LipshEval * (right * right - e * e);

			hyp.set_charact(std::min(charact1, charact2));
		}
	}
}