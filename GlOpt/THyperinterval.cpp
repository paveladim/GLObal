#include "THyperinterval.h"

uint THyperinterval::F_dimension;
uint THyperinterval::F_constraints;
uint THyperinterval::F_queue_depth;

THyperinterval::THyperinterval() : 
				F_idThis(0), 
				F_idPointA(0), 
				F_idPointB(0), 
				F_divisions(0),
				F_characteristic(0.0), 
				F_diagonal(0.0) {
	F_localsLipshEvaluations.resize(F_constraints + 1);
	F_maxLocalLipshEvaluations.resize(F_constraints + 1);
	for (auto& elem : F_maxLocalLipshEvaluations) elem = 0.0;
}

void THyperinterval::init_queues() {
	for (auto& elem : F_localsLipshEvaluations)
		for (uint i = 0; i < F_queue_depth; ++i) elem.push(0.0);
}

void THyperinterval::find_max_localLipshEvaluations() {
	LipschitzConstantValue potential_max = 0;
	for (uint i = 0; i < F_constraints + 1; ++i) {
		for (uint j = 0; j < F_queue_depth; ++j) {
			potential_max = F_localsLipshEvaluations[i].front();
			F_localsLipshEvaluations[i].pop();
			if (F_maxLocalLipshEvaluations[i] < potential_max)
				F_maxLocalLipshEvaluations[i] = potential_max;
			F_localsLipshEvaluations[i].push(potential_max);
		}
	}
}

void THyperinterval::update_queuesLipshEvaluations(std::vector<LipschitzConstantValue>& new_llcv, 
												   const double& _delta) {
	for (uint i = 0; i < F_constraints + 1; ++i) {
		F_localsLipshEvaluations[i].pop();
		if (new_llcv[i] > _delta)
			F_localsLipshEvaluations[i].push(new_llcv[i]);
		else
			F_localsLipshEvaluations[i].push(_delta);
	}
	
	find_max_localLipshEvaluations();
}
