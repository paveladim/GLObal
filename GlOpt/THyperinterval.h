#ifndef THYPERINTERVAL_H
#define THYPERINTERVAL_H

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
	std::vector<std::queue<LipschitzConstantValue>> F_evaluations; // очередь локальных оценок констант 
	std::vector<uint> F_division_tags; // ??

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

	void increase_division(); // увеличить число делений
	uint get_idThis();
	uint get_idPointA();
	uint get_idPointB();
	uint get_idA();
	uint get_idB();
	uint get_idEvaluationsA();
	uint get_idEvaluationsB();
	void set_idThis(const uint& out_idThis);
	void set_idPointA(const uint& out_idPA);
	void set_idPointB(const uint& out_idPB);
	void set_idA(const uint& out_idA);
	void set_idB(const uint& out_idB);
	void set_idEvaluationsA(const uint& out_idEvalA);
	void set_idEvaluationsB(const uint& out_idEvalB);
	double get_characteristic();
	void set_characteristic(const double& out_charact);
	uint get_div_tag(const uint& out_index);
	bool increase_div_tag(const uint& out_index);
	std::vector<LipschitzConstantValue>& getLocalLipschitzConstantValues(std::vector<LipschitzConstantValue>& llcv);
};

#endif
