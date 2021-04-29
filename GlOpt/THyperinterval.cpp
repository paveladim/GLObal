#include "THyperinterval.h"

uint THyperinterval::F_dimension;
uint THyperinterval::F_constraints;
uint THyperinterval::F_queue_depth;

THyperinterval::THyperinterval() : F_idThis(0), F_idPointA(0), F_idPointB(0), F_idA(0), F_idB(0),
									F_idEvaluationsA(0), F_idEvaluationsB(0), F_divisions(0),
									F_characteristic(0), F_diagonal(0) {
	F_localsLipshEvaluations.resize(F_constraints + 1);
}

std::vector<LipschitzConstantValue>& THyperinterval::getLocalLipschitzConstantValues(std::vector<LipschitzConstantValue>& llcv) {
	uint i = 0;
	for (auto it = F_localsLipshEvaluations.begin(); it != F_localsLipshEvaluations.end(); ++it) {
		llcv[i] = (*it).back();
		++i;
	}

	return llcv;
}
