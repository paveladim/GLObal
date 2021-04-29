#ifndef THYPERINTERVAL_H
#define THYPERINTERVAL_H

#include "dicretization.h"
#include "synonymous_types.h"

class THyperinterval
{
	// для связи с хранилищем
	uint F_idThis; // итендификатор гиперинтервала

	uint F_idPointA; // итендификатор точки активной диагонали
	uint F_idPointB; // итендификатор точки активной диагонали

	uint F_idA; // позиция начала размещения блока координат точки А в хранилище
	uint F_idB; // позиция начала размещения блока координат точки B в хранилище

	uint F_idEvaluationsA; // позиция начала размещения блока значений функций точки А в хранилище
	uint F_idEvaluationsB; // позиция начала размещения блока значений функций точки B в хранилище

	// собственные поля
	uint F_divisions; // число уже выполненных делений гиперинтервала
	double F_characteristic; // характеристика гиперинтервала
	double F_diagonal; // длина главной диагонали гиперинтервала
	std::vector<std::queue<LipschitzConstantValue>> F_localsLipshEvaluations; // очередь локальных оценок констант 

	// статические поля
	static uint F_dimension; // размерность задачи
	static uint F_constraints; // количество ограничений задачи
	static uint F_queue_depth; // глубина очереди
public:
	THyperinterval();
	THyperinterval(const THyperinterval& out_interval) = default;
	THyperinterval(THyperinterval&& out_interval) = default;
	static void init_static(const uint& out_dim, const uint& out_constr, const uint& out_depth) {
		F_dimension = out_dim;
		F_constraints = out_constr;
		F_queue_depth = out_depth;
	} 

	THyperinterval& operator=(const THyperinterval& out_interval) = default;

	void increase_division() { ++F_divisions; } // увеличить число делений
	uint get_idThis() const { return F_idThis; }
	uint get_idPointA() const { return F_idPointA; }
	uint get_idPointB() const { return F_idPointB; }
	uint get_idA() const { return F_idA; }
	uint get_idB() const { return F_idB; }
	uint get_idEvaluationsA() const { return F_idEvaluationsA; }
	uint get_idEvaluationsB() const { return F_idEvaluationsB; }
	uint get_div_tag() const { return MAX_EXPONENT_THREE - 1 - F_divisions / F_dimension; }
	std::vector<LipschitzConstantValue>& const getLocalLipschitzConstantValues(std::vector<LipschitzConstantValue>& llcv);
	double get_diagonal() const { return F_diagonal; }
	double get_characteristic() const { return F_characteristic; }
	void set_idThis(const uint& out_idThis) { F_idThis = out_idThis; }
	void set_idPointA(const uint& out_idPA) { F_idPointA = out_idPA; }
	void set_idPointB(const uint& out_idPB) { F_idPointB = out_idPB; }
	void set_idA(const uint& out_idA) { F_idA = out_idA; }
	void set_idB(const uint& out_idB) { F_idB = out_idB; }
	void set_idEvaluationsA(const uint& out_idEvalA) { F_idEvaluationsA = out_idEvalA; }
	void set_idEvaluationsB(const uint& out_idEvalB) { F_idEvaluationsB = out_idEvalB; }
	void set_diagonal(const double& diagonal) { F_diagonal = diagonal; }
	void set_characteristic(const double& out_charact) { F_characteristic = out_charact; }
};

#endif
