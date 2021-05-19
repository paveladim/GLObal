#include "THyperinterval.h"

uint THyperinterval::F_dimension;
uint THyperinterval::F_constraints;
uint THyperinterval::F_queue_depth;

THyperinterval::THyperinterval() : F_idThis(0), F_idPointA(0), F_idPointB(0), F_divisions(0),
									F_characteristic(0), F_diagonal(0) {
	F_localsLipshEvaluations.resize(F_constraints + 1);
	F_maxLocalLipshEvaluations.resize(F_constraints + 1);
	for (auto& elem : F_maxLocalLipshEvaluations) elem = 0.0;
}

std::vector<LipschitzConstantValue>& THyperinterval::getLocalLipschitzConstantValues(std::vector<LipschitzConstantValue>& llcv) {
	uint i = 0;
	for (auto& elem : F_localsLipshEvaluations) {
		llcv[i] = elem.back();
		++i;
	}

	return llcv;
}

void THyperinterval::init_queues() {
	for (auto& elem : F_localsLipshEvaluations)
		for (uint i = 0; i < F_queue_depth; ++i) elem.push(0.0);
}

void THyperinterval::update_queuesLipshEvaluations(std::vector<LipschitzConstantValue>& new_llcv, const double& _delta) {
	for (uint i = 0; i < F_constraints + 1; ++i) {
		F_localsLipshEvaluations[i].pop();
		if (new_llcv[i] > _delta) {
			if (new_llcv[i] > F_maxLocalLipshEvaluations[i])
				F_maxLocalLipshEvaluations[i] = new_llcv[i];
			F_localsLipshEvaluations[i].push(new_llcv[i]);
		}
		else {
			if (_delta > F_maxLocalLipshEvaluations[i])
				F_maxLocalLipshEvaluations[i] = _delta;
			F_localsLipshEvaluations[i].push(_delta);
		}
	}
}
