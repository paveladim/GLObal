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

    EncodedCoordinates new_coord_u;
    EncodedCoordinates new_coord_v;

    uint F_dimension; // ����������� ������
    uint F_constraints; // ����� ����������� � ������
    uint F_queueDepth; // ������� ������� ������ �������� �������

    uint F_divide_axis; // ����� ��� ������ ������?

    uint F_generated_points; // ����� ��� ��������������� ������� �����
    uint F_generated_intervals; // ����� ��� ��������������� ������� ���������������

    double F_eps; // �������� ���������
    double F_current_minimum; // ������� ����������� �������� ������� �������

    TProblem& Fp; // ����������� ������� ������� � �������-�����������
    std::vector<double> F_globalLipshEvaluations; // ��� �������� ��� ���������� ������ �������� �������
public:
    TMethodDivByThree() = delete;
    TMethodDivByThree(const uint& out_dim, const uint& out_constr, const uint& depth, const double& out_eps, TProblem& out_prob);
    ~TMethodDivByThree() {};

    virtual void initialization(); // ������� ����� ������ �������������
    virtual bool divideInterval(const uint& id_divHyp); // �������� �������������
    virtual void fillIntervals(THyperinterval& parent, const uint& id_u, const uint& id_v);
    void compute_diagonal(const uint& id_Hyp);
    virtual void compute_characteristic(const uint& id_Hyp);
    void compute_evaluations(const uint& out_idPoint);
    uint get_new_id() { return F_generated_points++; } // ������ �������������� �����
    uint get_id_coord() const { return F_generated_points * F_dimension; } // �������� ������������� ��� ���������
    uint get_id_evaluations() const { return F_generated_points * (F_constraints + 1); }
    uint get_new_interval() { return F_generated_intervals++; } // ������ �������������� ��������������

    // uint does_point_exist(TPoint& parent, const uint& des_val, const TPoint::direction& dir); // ���������� �� ����� � ���� ������ ���
};

#endif