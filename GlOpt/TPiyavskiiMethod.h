#ifndef TPIYAVSKIIMETHOD_H
#define TPIYAVSKIIMETHOD_H

#include "TMethodDivByThree.h"

class TPiyavskiiMethod : protected TMethodDivByThree
{
	double F_criticalSize; // ��������� ������ ��������������
	GainLipshConstant F_gainObjective; // ��������� ��������� ��������� ������� ��� ������� �������
	GainLipshConstant F_gainConstraints; // ��������� ��������� �������� ������� ��� �����������
	double delta; // �������� ��� ���������� ��������� ������
	bool does_LipshConstValue_change; // ���������� �� �����-�� �� ���������� ������ �������� �������
public:
	TPiyavskiiMethod() = delete;
	TPiyavskiiMethod(const uint& out_dim, const uint& out_constr, const uint& depth,
		TProblem& out_prob, const GainLipshConstant& out_gainObj, const GainLipshConstant& out_gainCst, const double& beta);

	virtual void compute_characteristic(const uint& id_Hyp); // ��������� ��������������
	virtual uint choose_optimal_to_trisect(); // ����� ����������� ��� ������� �� ���
	void compute_localLipshConst(const uint& id_Hyp); // ��������� ��������� ������ �������� � ��������������
	void update_globalLipshEval(const uint& id_Hyp); // �������� ���������� ������ ��������
	void update_all_characteristics(); // �������� �������������� ���� ���������������
	virtual uint do_step(const uint& if_divHyp); // ��������� ��� ������
	virtual void launch_method(); // ������ ������ ������
};

#endif