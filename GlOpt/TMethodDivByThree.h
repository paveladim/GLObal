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

    // ���� ��� ������ ������������
    EncodedCoordinates transit_coord_1;
    EncodedCoordinates transit_coord_2;

    uint F_dimension; // ����������� ������
    uint F_constraints; // ����� ����������� � ������
    uint F_queueDepth; // ������� ������� ������ �������� �������

    uint F_divide_axis; // ����� ��� ������ ������?

    uint F_generated_points; // ����� ��� ��������������� ������� �����
    uint F_generated_intervals; // ����� ��� ��������������� ������� ���������������

    double F_criticalSize; // ��������� ������ ��������������
    double F_current_minimum; // ������� ����������� �������� ������� �������
    uint F_id_current_minimum; // ������������� ����������� ��������

    TProblem& Fp; // ����������� ������� ������� � �������-�����������
    GainLipshConstant F_gainObjective; // ��������� ��������� ��������� ������� ��� ������� �������
    GainLipshConstant F_gainConstraints; // ��������� ��������� �������� ������� ��� �����������
    double delta; // �������� ��� ���������� ��������� ������
    std::vector<double> F_globalLipshEvaluations; // ������ ���������� ������ �������� ������� ������� ������� � �����������
public:
    TMethodDivByThree() = delete;
    TMethodDivByThree(const uint& out_dim, const uint& out_constr, const uint& depth, 
        TProblem& out_prob, const GainLipshConstant& out_gainObj, const GainLipshConstant& out_gainCst, const double& beta);
    ~TMethodDivByThree() {};

    virtual void initialization(); // ������� ����� ������ �������������
    virtual bool trisect_interval(const uint& id_divHyp); // �������� �������������
    virtual void fill_intervals(THyperinterval& parent, const uint& id_u, const uint& id_v); // ��������� ������������� �������
    void compute_diagonal(const uint& id_Hyp); // ��������� ������ �������� ��������� ��������������
    virtual void compute_characteristic(const uint& id_Hyp); // ��������� �������������� ��������������
    void compute_evaluations(const uint& out_idPoint); // ��������� �������� ������� � �����
    uint get_new_id() { return F_generated_points++; } // ������ �������������� �����
    uint get_id_coord() const { return F_generated_points * F_dimension; } // �������� ������������� ��� ���������
    // �������� ������������� �� ������ ���������� �������� ������� ������� � ����������� � �����
    uint get_id_evaluations() const { return F_generated_points * (F_constraints + 1); }
    uint get_new_interval() { return F_generated_intervals++; } // ������ �������������� ��������������
    void compute_localLipshConst(const uint& id_Hyp); // ��������� ��������� ������ �������� �������
    void compute_globalLipshConst(); // ��������� ���������� ������ �������� �������
    // ������� "������" ������������� ��� �������
    uint choose_optimal_to_trisect(); // ����� ����������� ��� ������� �� ���
    uint do_step(const uint& if_divHyp); // ������� ��� ������
    void launch_method(); // ��� ������������
    // ������ ���������� �����
    void resize_points_deque() {
        if (F_points.size() - F_generated_points < 1)
            F_points.resize(F_points.size() + 100);
    }

    void resize_coords_deque() { 
        if (F_coords.size() - F_generated_points * F_dimension < F_dimension) 
            F_coords.resize(F_coords.size() + 100); 
    }

    void resize_evaluations_deque() {
        if ((F_evaluations.size() - F_generated_points * (F_constraints + 1) < (F_constraints + 1)) || F_evaluations.size() == 0)
            F_evaluations.resize(F_evaluations.size() + 100);
    }

    void resize_intervals_deque() {
        if (F_intervals.size() < F_generated_intervals + 2)
            F_intervals.resize(F_intervals.size() + 100);
    }
};

#endif