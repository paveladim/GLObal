#ifndef THYPERINTERVAL_H
#define THYPERINTERVAL_H

#include "dicretization.h"
#include "synonymous_types.h"

class THyperinterval
{
	// ��� ����� � ����������
	uint F_idThis; // ������������� ��������������

	uint F_idPointA; // ������������� ����� �������� ���������
	uint F_idPointB; // ������������� ����� �������� ���������

	// ����������� ����
	uint F_divisions; // ����� ��� ����������� ������� ��������������
	double F_characteristic; // �������������� ��������������
	double F_diagonal; // ����� ������� ��������� ��������������
	std::vector<std::queue<LipschitzConstantValue>> F_localsLipshEvaluations; // ������� ��������� ������ �������� 
	std::vector<LipschitzConstantValue> F_maxLocalLipshEvaluations; // ������������ �������� ������ ��������

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
	} // ������������� ����������� �����

	THyperinterval& operator=(const THyperinterval& out_interval) = default;

	void increase_division() { ++F_divisions; } // ��������� ����� �������
	uint get_idThis() const { return F_idThis; } // �������� ������������� ��������������
	uint get_idPointA() const { return F_idPointA; } // �������� ������������� ����� �������� ���������
	uint get_idPointB() const { return F_idPointB; } // �������� ������������� ����� �������� ���������
	// �������� ������������� ����� ������ �������� ��������� ����� �������� ���������
	uint get_idA() const { return F_idPointA * F_dimension; }
	// �������� ������������� ����� ������ �������� ��������� ����� �������� ���������
	uint get_idB() const { return F_idPointB * F_dimension; }
	// �������� ������������� �������� ������� ������� � ����������� � ����� �������� ���������
	uint get_idEvaluationsA() const { return F_idPointA * (F_constraints + 1); }
	// �������� ������������� �������� ������� ������� � ����������� � ����� �������� ���������
	uint get_idEvaluationsB() const { return F_idPointB * (F_constraints + 1); }
	uint get_div_tag() const { return MAX_EXPONENT_THREE - 1 - F_divisions / F_dimension; }
	uint get_div_axis() const { return F_divisions % F_dimension; }
	double get_diagonal() const { return F_diagonal; } // �������� ������ ��������� ��������������
	double get_characteristic() const { return F_characteristic; } // �������� �������� �������������� ��������������
	void set_idThis(const uint& out_idThis) { F_idThis = out_idThis; } // ���������� �������������� ��� �������������
	void set_idPointA(const uint& out_idPA) { F_idPointA = out_idPA; } // ���������� ����� ������� ��������� � �������������
	void set_idPointB(const uint& out_idPB) { F_idPointB = out_idPB; } // ���������� ����� ������� ��������� � �������������
	void set_diagonal(const double& diagonal) { F_diagonal = diagonal; } // ���������� ������ ��������� ��������������
	void set_characteristic(const double& out_charact) { F_characteristic = out_charact; }
	void init_queues(); // ���������������� ������� ������ �������� ������� ������
	void find_max_localLipshEvaluations();
	// ���������� ������ �������� �������, ��������� ����� �������
	void update_queuesLipshEvaluations(std::vector<LipschitzConstantValue>& new_llcv, const double& _delta);
	// �������� ������ ������������ ��������� ������ ��������
	const std::vector<LipschitzConstantValue>& 
	get_maxLipshEvaluations() const { return F_maxLocalLipshEvaluations; }
};

#endif
