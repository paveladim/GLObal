#include "THyperinterval.h"

uint THyperinterval::F_dimension;
uint THyperinterval::F_constraints;
uint THyperinterval::F_queue_depth;

THyperinterval::THyperinterval() : F_idThis(0), F_idPointA(0), F_idPointB(0), F_idA(0), F_idB(0),
									F_idEvaluationsA(0), F_idEvaluationsB(0), F_divisions(0),
									F_characteristic(0), F_diagonal(0) {
	F_division_tags.resize(F_dimension);
	for (auto it = F_division_tags.begin(); it != F_division_tags.end(); ++it)
		*it = 0;
	F_evaluations.resize(F_constraints + 1);
}

THyperinterval& THyperinterval::operator=(const THyperinterval& out_interval) {
	if (this == &out_interval) return *this;
	uint F_idThis; // итендификатор гиперинтервала

	F_idPointA = out_interval.F_idPointA; 
	F_idPointB = out_interval.F_idPointB;

	F_idA = out_interval.F_idA; 
	F_idB = out_interval.F_idB; 

	F_idEvaluationsA = out_interval.F_idEvaluationsA; 
	F_idEvaluationsB = out_interval.F_idEvaluationsB; 

	F_divisions = out_interval.F_divisions; 
	F_characteristic = out_interval.F_characteristic; 
	F_diagonal = out_interval.F_divisions; 
	F_evaluations = out_interval.F_evaluations;
	F_division_tags = out_interval.F_division_tags;

	return *this;
}

void THyperinterval::increase_division() {
	++F_divisions;
}

uint THyperinterval::get_idThis() {
	return F_idThis;
}

uint THyperinterval::get_idPointA() {
	return F_idPointA;
}

uint THyperinterval::get_idPointB() {
	return F_idPointB;
}

uint THyperinterval::get_idA() {
	return F_idA;
}

uint THyperinterval::get_idB() {
	return F_idB;
}

uint THyperinterval::get_idEvaluationsA() {
	return F_idEvaluationsA;
}

uint THyperinterval::get_idEvaluationsB() {
	return F_idEvaluationsB;
}

double THyperinterval::get_diagonal() {
	return F_diagonal;
}

void THyperinterval::set_idThis(const uint& out_idThis) {
	F_idThis = out_idThis;
}

void THyperinterval::set_idPointA(const uint& out_idPA) {
	F_idPointA = out_idPA;
}

void THyperinterval::set_idPointB(const uint& out_idPB) {
	F_idPointB = out_idPB;
}

void THyperinterval::set_idA(const uint& out_idA) {
	F_idA = out_idA;
}

void THyperinterval::set_idB(const uint& out_idB) {
	F_idB = out_idB;
}

void THyperinterval::set_idEvaluationsA(const uint& out_idEvalA) {
	F_idEvaluationsA = out_idEvalA;
}

void THyperinterval::set_idEvaluationsB(const uint& out_idEvalB) {
	F_idEvaluationsB = out_idEvalB;
}

void THyperinterval::set_characteristic(const double& out_charact) {
	F_characteristic = out_charact;
}

void THyperinterval::set_diagonal(const double& value) {
	F_diagonal = value;
}

double THyperinterval::get_characteristic() {
	return F_characteristic;
}

uint THyperinterval::get_div_tag(const uint& out_index) {
	return F_division_tags[out_index];
}

bool THyperinterval::increase_div_tag(const uint& out_index) {
	if (F_division_tags[out_index] == 20) return false;
	++F_division_tags[out_index];
	return true;
}

std::vector<LipschitzConstantValue>& THyperinterval::getLocalLipschitzConstantValues(std::vector<LipschitzConstantValue>& llcv) {
	uint i = 0;
	for (auto it = F_evaluations.begin(); it != F_evaluations.end(); ++it) {
		llcv[i] = (*it).back();
		++i;
	}

	return llcv;
}
