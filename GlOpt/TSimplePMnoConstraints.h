#ifndef TSIMPLEPM_NOCONSTRAINTS_H
#define TSIMPLEPM_NOCONSTRAINTS_H

#include "TMethodDivByThree.h"

class TSimplePMnoConstraints : public TMethodDivByThree {
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
private:
	void compute_characteristic(const uint& id_Hyp) override;
	double get_mixedLipshitzEval(const THyperinterval& hyp);
	void compute_localLipshConst(const uint& id_Hyp);
	void update_globalLipshEval(const uint& id_Hyp);
	// �������� �������������� ���� ���������������
	void update_all_characteristics();
	// ����� ����������� ��� ������� �� ���
	uint choose_optimal_to_trisect() override;
	uint do_step(const uint& id_divHyp) override;
public:
	TSimplePMnoConstraints() = delete;
	TSimplePMnoConstraints(const uint& out_dim,
						   const uint& out_constr,
						   const uint& depth,
						   TProblem& out_prob,
						   const GainLipshConstant& out_gainObj,
					       const GainLipshConstant& out_gainCst,
						   const double& beta,
						   const double& _eps
						   );
	void launch_method() override;
};

#endif // TSIMPLEPM_NOCONSTRAINTS_H