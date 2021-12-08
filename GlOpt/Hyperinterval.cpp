#include "Hyperinterval.h"

uint Hyperinterval::_dimension;
uint Hyperinterval::_constraints;
uint Hyperinterval::_queueDepth;

Hyperinterval::Hyperinterval() : _id(0),
								 _idA(0),
								 _idB(0),
								 _divisions(0),
								 _diagonal(0.0),
								 _charact(0.0),
								 _localLipEvaluations(_constraints + 1),
								 _maxLipEvaluations(_constraints + 1, 0) {
	for (auto& elem : _localLipEvaluations)
		for (uint i = 0; i < _queueDepth; ++i) elem.push(0.0);
}

Hyperinterval::Hyperinterval(const Hyperinterval& hyperinterval) : 
	_id(hyperinterval._id),
	_idA(hyperinterval._idA),
	_idB(hyperinterval._idB),
	_divisions(hyperinterval._divisions),
	_diagonal(hyperinterval._diagonal),
	_charact(hyperinterval._charact),
	_localLipEvaluations(hyperinterval._localLipEvaluations),
	_maxLipEvaluations(hyperinterval._maxLipEvaluations) {}

Hyperinterval& Hyperinterval::operator=(const Hyperinterval& hyperinterval) {
	if (this == &hyperinterval) return *this;
	_id = hyperinterval._id;
	_idA = hyperinterval._idA;
	_idB = hyperinterval._idB;
	_divisions = hyperinterval._divisions;
	_diagonal = hyperinterval._diagonal;
	_charact = hyperinterval._charact;
	_localLipEvaluations = hyperinterval._localLipEvaluations;
	_maxLipEvaluations = hyperinterval._maxLipEvaluations;
	return *this;
}

void Hyperinterval::init_static(const uint& dimension, 
								const uint& constraints, 
								const uint& queueDepth) {
	_dimension = dimension;
	_constraints = constraints;
	_queueDepth = queueDepth;
}

uint Hyperinterval::get_id() const {
	return _id;
}

uint Hyperinterval::get_idA() const {
	return _idA;
}

uint Hyperinterval::get_idB() const {
	return _idB;
}

uint Hyperinterval::get_coordA() const {
	return _idA * _dimension;
}

uint Hyperinterval::get_coordB() const {
	return _idB * _dimension;
}

uint Hyperinterval::get_evalA() const {
	return _idA * (_constraints + 1);
}

uint Hyperinterval::get_evalB() const {
	return _idB * (_constraints + 1);
}

double Hyperinterval::get_charact(const bool& mode) const {
	if (mode) return _diagonal;
	return _charact;
}

void Hyperinterval::set_id(const uint& id) {
	_id = id;
}

void Hyperinterval::set_idA(const uint& idA) {
	_idA = idA;
}

void Hyperinterval::set_idB(const uint& idB) {
	_idB = idB;
}

void Hyperinterval::set_diagonal(const double& diagonal) {
	_diagonal = diagonal;
}

void Hyperinterval::set_charact(const double& charact) {
	_charact = charact;
}

void Hyperinterval::increase_divisions() {
	++_divisions;
}

uint Hyperinterval::get_divisions() const {
	return _divisions;
}

uint Hyperinterval::get_axis() const {
	return _divisions % _dimension;
}

uint Hyperinterval::get_previous_axis() const {
	return (_divisions - 1) % _dimension;
}

uint Hyperinterval::get_shift() const {
	return MAX_EXPONENT_THREE - 1 - _divisions / _dimension;
}

double Hyperinterval::get_diagonal() const {
	return _diagonal;
}

void Hyperinterval::update_localLipQueues(std::vector<LipschitzConstantValue>& new_llcv,
										  const double& delta) {
	for (uint i = 0; i < _constraints + 1; ++i) {
		_localLipEvaluations[i].pop();
		if (new_llcv[i] > delta)
			_localLipEvaluations[i].push(new_llcv[i]);
		else
			_localLipEvaluations[i].push(delta);
	}

	find_maxLipEval();
}

void Hyperinterval::find_maxLipEval() {
	LipschitzConstantValue potential_max = 0;
	for (uint i = 0; i < _constraints + 1; ++i) {
		for (uint j = 0; j < _queueDepth; ++j) {
			potential_max = _localLipEvaluations[i].front();
			_localLipEvaluations[i].pop();
			if (_maxLipEvaluations[i] < potential_max)
				_maxLipEvaluations[i] = potential_max;
			_localLipEvaluations[i].push(potential_max);
		}
	}
}

const std::vector<LipschitzConstantValue>& 
Hyperinterval::get_maxLipshEvaluations() const {
	return _maxLipEvaluations;
}