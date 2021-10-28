#ifndef TSIMPLEPM_WITH_CONSTRAINTS_H
#define TSIMPLEPM_WITH_CONSTRAINTS_H

#include "TMethodDivByThree.h"

class TSimplePMwithConstraints : public TMethodDivByThree {
	// пороговый размер гиперинтервала
	double F_criticalSize;
	// константа завышения константы Липшица для целевой функции
	GainLipshConstant F_gainObjective;
	// константа завышения констант Липшица для ограничений
	GainLipshConstant F_gainConstraints;
	// параметр для осторожных локальных оценок
	double delta;
	// изменилась ли какая-то из глобальных оценок констант Липшица?
	bool does_LipshConstValue_change;
	// для остановки по точности
	double eps;
	// флаг, показывающий все ли характеристики бесконечно большие
	bool are_allCharInfty;
private:
	virtual void compute_characteristic(const uint& id_Hyp) override;
	void compute_localLipshConst(const uint& id_Hyp);
	void update_globalLipshEval(const uint& id_Hyp);
	double get_mixedLipshitzEval(const THyperinterval& hyp, const uint& i);
	void give_borders(double& l, double& r, const THyperinterval& hyp);
	void check_ifAllCharAreInfty();
	// обновить характеристики всех гиперинтервалов
	void update_all_characteristics();
	// найти оптимальный для деления на три
	uint choose_optimal_to_trisect() override;
	uint do_step(const uint& id_divHyp) override;
public:
	TSimplePMwithConstraints() = delete;
	TSimplePMwithConstraints(const uint& out_dim,
							 const uint& out_constr,
							 const uint& depth,
							 TProblem& out_prob,
							 const GainLipshConstant& out_gainObj,
							 const GainLipshConstant& out_gainCst,
							 const double& beta,
							 const double& _eps
	                         );
	virtual double get_critical_size() { return F_criticalSize; }
	void launch_method() override;
};

#endif // TSIMPLEPM_WITH_CONSTRAINTS_H
