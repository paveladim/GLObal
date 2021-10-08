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

	// собственные поля
	uint F_divisions; // число уже выполненных делений гиперинтервала
	double F_characteristic; // характеристика гиперинтервала
	double F_diagonal; // длина главной диагонали гиперинтервала
	std::vector<std::queue<LipschitzConstantValue>> F_localsLipshEvaluations; // очередь локальных оценок констант 
	std::vector<LipschitzConstantValue> F_maxLocalLipshEvaluations; // максимальные значения оценок констант

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
	} // инициализация статических полей

	THyperinterval& operator=(const THyperinterval& out_interval) = default;

	void increase_division() { ++F_divisions; } // увеличить число делений
	uint get_idThis() const { return F_idThis; } // получить итендификатор гиперинтервала
	uint get_idPointA() const { return F_idPointA; } // получить итендификатор точки активной диагонали
	uint get_idPointB() const { return F_idPointB; } // получить итендификатор точки активной диагонали
	// получить итендификатор места начала хранения координат точки активной диагонали
	uint get_idA() const { return F_idPointA * F_dimension; }
	// получить итендификатор места начала хранения координат точки активной диагонали
	uint get_idB() const { return F_idPointB * F_dimension; }
	// получить итендификатор значения целевой функции и ограничений в точке активной диагонали
	uint get_idEvaluationsA() const { return F_idPointA * (F_constraints + 1); }
	// получить итендификатор значения целевой функции и ограничений в точке активной диагонали
	uint get_idEvaluationsB() const { return F_idPointB * (F_constraints + 1); }
	uint get_div_tag() const { return MAX_EXPONENT_THREE - 1 - F_divisions / F_dimension; }
	uint get_div_axis() const { return F_divisions % F_dimension; }
	double get_diagonal() const { return F_diagonal; } // получить размер диагонали гиперинтервала
	double get_characteristic() const { return F_characteristic; } // получить значение характеристики гиперинтервала
	void set_idThis(const uint& out_idThis) { F_idThis = out_idThis; } // установить гиперинтервалу его итендификатор
	void set_idPointA(const uint& out_idPA) { F_idPointA = out_idPA; } // установить точке главной диагонали её итендификатор
	void set_idPointB(const uint& out_idPB) { F_idPointB = out_idPB; } // установить точке главной диагонали её итендификатор
	void set_diagonal(const double& diagonal) { F_diagonal = diagonal; } // установить размер диагонали гиперинтервала
	void set_characteristic(const double& out_charact) { F_characteristic = out_charact; }
	void init_queues(); // инициализировать очереди оценок констант Липшица нулями
	void find_max_localLipshEvaluations();
	// обновление оценок констант Липшица, поддержка длины очереди
	void update_queuesLipshEvaluations(std::vector<LipschitzConstantValue>& new_llcv, const double& _delta);
	// получить вектор максимальных локальных оценок констант
	const std::vector<LipschitzConstantValue>& 
	get_maxLipshEvaluations() const { return F_maxLocalLipshEvaluations; }
};

#endif
