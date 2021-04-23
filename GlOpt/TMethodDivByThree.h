#ifndef TMETHODDIVBYTHREE_H
#define TMETHODDIVBYTHREE_H

#include <deque>
#include "synonymous_types.h"
#include "TProblem.h"
#include "THyperinterval.h"
#include "TPoint.h"

class TMethodDivByThree
{
    std::deque<uint> F_coords; // ��� ��������� ���������� �����
    std::deque<FunctionValue> F_evaluations; // ��� ��������� � ���������� ������
    std::deque<THyperinterval> F_intervals; // ��� ���������� ���������������
    std::deque<TPoint> F_points; // ��� ���������� �����

    CoordinatesValues coord_a;
    CoordinatesValues coord_b;

    uint F_dimension; // ����������� ������
    uint F_constraints; // ����� ����������� � ������
    uint F_queueDepth; // ������� ������� ������ �������� �������

    uint F_divide_axis; // ����� ��� �����

    uint F_generated_points; // ����� ��� ��������������� ������� �����
    uint max_generated_points; // ������������ ���������� ����� �����, ����� �������������

    double F_eps; // �������� ���������
    double F_current_minimum; // ������� ����������� �������� ������� �������
    uint F_generated_intervals; // ����� ��� ��������������� ������� ���������������
    uint max_generated_intervals; // ������������ ���������� ����� ���������������, ������� ����� �������������

    TProblem& Fp; // ������ ������� ������� � �������-�����������
    std::vector<std::queue<double>> global_evaluations; // ��� �������� ��� ���������� ������ �������� �������
public:
    TMethodDivByThree() = delete;
    TMethodDivByThree(const uint& out_dim, const uint& out_constr, const uint& depth, const double& out_eps,
        const uint& max_gen_points, const uint& max_gen_interv, TProblem& out_prob);
    ~TMethodDivByThree() {};

    virtual void createFirstInterval(); // ������� ����� ������ �������������
    virtual bool divideInterval(const uint& id_divHyp); // �������� �������������
    // virtual uint compute_new_coord(const uint& out_id, bool f); // ������������� ����� ����������
    uint get_new_id(); // ������ �������������� �����
    uint get_id_coord(); // �������� ������������� ��� ���������
    uint get_new_interval(); // ������ �������������� ��������������
    void delete_point(); // ������� ������������� �����
    void delete_hyperinterval(); // ������� ������������� ��������������
    // bool get_coord(const uint& out_id_a); // ���������� ����� �� ���� ������
    // bool get_coords(const uint& out_id_a, const uint& out_id_b); // ���������� ����� �� ���� ������

    uint does_point_exist(TPoint& parent, const uint& des_val, bool f); // ���������� �� ����� � ���� ������ ���
};

#endif