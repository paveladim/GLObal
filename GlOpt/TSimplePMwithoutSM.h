#ifndef __TSIMPLEPM_WITHOUT_SM__
#define __TSIMPLEPM_WITHOUT_SM__

#include "TMethodDivByThree.h"

class TSimplePMwithoutSM : protected TMethodDivByThree {
	double F_criticalSize; // пороговый размер гиперинтервала
	GainLipshConstant F_gainObjective; // константа завышения константы Липшица для целевой функции
	GainLipshConstant F_gainConstraints; // константа завышения констант Липшица для ограничений
	double delta; // параметр для осторожных локальных оценок
	bool does_LipshConstValue_change; // изменилась ли какая-то из глобальных оценок констант Липшица
	double eps; // для остановки по точности
	uint F_iter;
public:
	TSimplePMwithoutSM() = delete;
	TSimplePMwithoutSM(const uint& out_dim,
					   const uint& out_constr,
					   const uint& depth,
					   TProblem& out_prob,
					   const GainLipshConstant& out_gainObj,
					   const GainLipshConstant& out_gainCst,
					   const double& beta,
					   const double& _eps);

	virtual void initialization(); // инициализация гиперинтервала
	virtual void compute_characteristic(const uint& id_Hyp); // вычислить характеристику
	virtual uint choose_optimal_to_trisect(); // найти оптимальный для деления на три
	// вычислить локальную оценку констант в гиперинтервале
	void compute_localLipshConst(const uint& id_Hyp1, const uint& id_Hyp2, const uint& id_Hyp3);
	void update_globalLipshEval(const uint& id_Hyp); // обновить глобальные оценки констант
	void update_all_characteristics(); // обновить характеристики всех гиперинтервалов
	virtual uint do_step(const uint& id_divHyp); // совершить шаг метода
	virtual void launch_method(); // запуск работы метода
	void write_generated_points_to_file(); // записать в файл все порождённые методом точки испытания
	void write_intervals_to_file(); // записать в файл точки активной диагонали всех гиперинтервалов
};

#endif // __TSIMPLEPM_WITHOUT_SM__

