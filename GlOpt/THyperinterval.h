#ifndef THYPERINTERVAL_H
#define THYPERINTERVAL_H

#include "synonymous_types.h"

class THyperinterval
{
	// ��� ����� � ����������
	uint F_idThis; // ������������� ��������������

	uint F_idPointA; // ������������� ����� �������� ���������
	uint F_idPointB; // ������������� ����� �������� ���������

	uint F_idA; // ������� ������ ���������� ����� ��������� ����� � � ���������
	uint F_idB; // ������� ������ ���������� ����� ��������� ����� B � ���������

	uint F_idEvaluationsA; // ������� ������ ���������� ����� �������� ������� ����� � � ���������
	uint F_idEvaluationsB; // ������� ������ ���������� ����� �������� ������� ����� B � ���������

	// ����������� ����
	uint F_divisions; // ����� ��� ����������� ������� ��������������
	double F_characteristic; // �������������� ��������������
	double F_diagonal; // ����� ������� ��������� ��������������
	std::vector<std::queue<LipschitzConstantValue>> F_evaluations; // ������� ��������� ������ �������� 
	std::vector<uint> F_division_tags; // ??

	// ����������� ����
	static uint F_dimension; // ����������� ������
	static uint F_constraints; // ���������� ����������� ������
	static uint F_queue_depth; // ������� �������

public:
	THyperinterval();
	THyperinterval(const THyperinterval& out_interval) = default;
	THyperinterval(THyperinterval&& out_interval) = default;
	static void init_static(const uint& out_dim, const uint& out_constr, const uint& out_depth) {
		F_dimension = out_dim;
		F_constraints = out_constr;
		F_queue_depth = out_depth;
	} 

	void increase_division(); // ��������� ����� �������
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
