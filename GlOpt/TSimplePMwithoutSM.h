#ifndef __TSIMPLEPM_WITHOUT_SM__
#define __TSIMPLEPM_WITHOUT_SM__

#include "TMethodDivByThree.h"

class TSimplePMwithoutSM : protected TMethodDivByThree {
	double F_criticalSize; // ��������� ������ ��������������
	GainLipshConstant F_gainObjective; // ��������� ��������� ��������� ������� ��� ������� �������
	GainLipshConstant F_gainConstraints; // ��������� ��������� �������� ������� ��� �����������
	double delta; // �������� ��� ���������� ��������� ������
	bool does_LipshConstValue_change; // ���������� �� �����-�� �� ���������� ������ �������� �������
	double eps; // ��� ��������� �� ��������
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

	virtual void initialization(); // ������������� ��������������
	virtual void compute_characteristic(const uint& id_Hyp); // ��������� ��������������
	virtual uint choose_optimal_to_trisect(); // ����� ����������� ��� ������� �� ���
	// ��������� ��������� ������ �������� � ��������������
	void compute_localLipshConst(const uint& id_Hyp1, const uint& id_Hyp2, const uint& id_Hyp3);
	void update_globalLipshEval(const uint& id_Hyp); // �������� ���������� ������ ��������
	void update_all_characteristics(); // �������� �������������� ���� ���������������
	virtual uint do_step(const uint& id_divHyp); // ��������� ��� ������
	virtual void launch_method(); // ������ ������ ������
	void write_generated_points_to_file(); // �������� � ���� ��� ���������� ������� ����� ���������
	void write_intervals_to_file(); // �������� � ���� ����� �������� ��������� ���� ���������������
};

#endif // __TSIMPLEPM_WITHOUT_SM__

