#ifndef TSIMPLEPM_WITH_CONSTRAINTS_H
#define TSIMPLEPM_WITH_CONSTRAINTS_H

#include "TMethodDivByThree.h"

class TSimplePMwithConstraints : public TMethodDivByThree {
	// ��������� ������ ��������������
	double F_criticalSize;
	// ��������� ��������� ��������� ������� ��� ������� �������
	GainLipshConstant F_gainObjective;
	// ��������� ��������� �������� ������� ��� �����������
	GainLipshConstant F_gainConstraints;
	// �������� ��� ���������� ��������� ������
	double delta;
	// ���������� �� �����-�� �� ���������� ������ �������� �������?
	bool does_LipshConstValue_change;
	// ��� ��������� �� ��������
	double eps;
	// ����, ������������ ��� �� �������������� ���������� �������
	bool are_allCharInfty;
private:
	virtual void compute_characteristic(const uint& id_Hyp) override;
	void compute_localLipshConst(const uint& id_Hyp);
	void update_globalLipshEval(const uint& id_Hyp);
	double get_mixedLipshitzEval(const THyperinterval& hyp, const uint& i);
	void give_borders(double& l, double& r, const THyperinterval& hyp);
	void check_ifAllCharAreInfty();
	// �������� �������������� ���� ���������������
	void update_all_characteristics();
	// ����� ����������� ��� ������� �� ���
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
