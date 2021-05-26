#ifndef TPIYAVSKIIMETHOD_H
#define TPIYAVSKIIMETHOD_H

#include "TMethodDivByThree.h"

class TPiyavskiiMethod : protected TMethodDivByThree
{
	double F_criticalSize; // пороговый размер гиперинтервала
	GainLipshConstant F_gainObjective; // константа завышения константы Липшица для целевой функции
	GainLipshConstant F_gainConstraints; // константа завышения констант Липшица для ограничений
	double delta; // параметр для осторожных локальных оценок
	bool does_LipshConstValue_change; // изменилась ли какая-то из глобальных оценок констант Липшица
public:
	TPiyavskiiMethod() = delete;
	TPiyavskiiMethod(const uint& out_dim, const uint& out_constr, const uint& depth,
		TProblem& out_prob, const GainLipshConstant& out_gainObj, const GainLipshConstant& out_gainCst, const double& beta);

	virtual void compute_characteristic(const uint& id_Hyp); // вычислить характеристику
	virtual uint choose_optimal_to_trisect(); // найти оптимальный для деления на три
	void compute_localLipshConst(const uint& id_Hyp); // вычислить локальную оценку констант в гиперинтервале
	void update_globalLipshEval(const uint& id_Hyp); // обновить глобальные оценки констант
	void update_all_characteristics(); // обновить характеристики всех гиперинтервалов
	virtual uint do_step(const uint& if_divHyp); // совершить шаг метода
	virtual void launch_method(); // запуск работы метода
};

#endif