#ifndef HYPERINTERVAL_H
#define HYPERINTERVAL_H

#include "dicretization.h"
#include "synonymous_types.h"

class Hyperinterval {
	// идентификатор гиперинтервала
	uint _id;

	// идентификатор точек главной диагонали
	uint _idA;
	uint _idB;

	// количество выполненных делений гиперинтервала
	uint _divisions;

	// характеристика гиперинтервала
	double _charact;
	// длина активной диагонали гиперинтервала
	double _diagonal;

	// очередь локальных оценок констант Липшица
	std::vector<std::queue<LipschitzConstantValue>> _localLipEvaluations;
	// максимальное значение локальных оценок констант Липшица
	std::vector<LipschitzConstantValue> _maxLipEvaluations;

	// размерность задачи
	static uint _dimension;
	// количество ограничений
	static uint _constraints;
	// глубина очереди
	static uint _queueDepth;

public:
	Hyperinterval();
	Hyperinterval(const Hyperinterval& hyperinterval);
	Hyperinterval& operator=(const Hyperinterval& hyperinterval);

	// инициализировать статические поля
	static void init_static(const uint&, const uint&, const uint&);

	// получить идентификатор гиперинтервала
	uint get_id() const;
	// получить идентификатор точки главной диагонали
	uint get_idA() const;
	uint get_idB() const;
	// получить начало хранения координат точки главной диагонали
	uint get_coordA() const;
	uint get_coordB() const;
	// получить значение характеристики гиперинтервала
	double get_charact(const bool& = false) const;

	// установить идентификатор гиперинтервалу
	void set_id(const uint&);
	// установить идентификатор точки главной диагонали
	void set_idA(const uint&);
	void set_idB(const uint&);
	// установить значение главной диагонали гиперинтервала
	void set_diagonal(const double&);
	// установить значение характеристики гиперинтервала
	void set_charact(const double&);

	// увеличить число делений гиперинтервала
	void increase_divisions();
	// получить число делений
	uint get_divisions() const;
	// получить номер координатной оси, по которой будет делить гиперинтервал
	uint get_axis() const;
	// получить сдвиг вдоль координатной оси для деления гиперинтервала
	uint get_shift() const;

	// обновление оценок констант Липшица, поддержка длины очереди, 
	// вычисление максимальных оценок
	void update_localLipQueues(std::vector<LipschitzConstantValue>&, 
							   const double&);

	// получить значения максимальных оценок констант Липшица
	const std::vector<LipschitzConstantValue>& get_maxLipshEvaluations() const;
};

#endif HYPERINTERVAL_H
