#ifndef TMETHODDIVBYTHREE_H
#define TMETHODDIVBYTHREE_H

#include <deque>
#include "synonymous_types.h"
#include "TProblem.h"
#include "THyperinterval.h"
#include "TPoint.h"

class TMethodDivByThree
{
    std::deque<uint> F_coords; // дек координат порождённых точек
    std::deque<FunctionValue> F_evaluations; // дек измерений в порождённых точках
    std::deque<THyperinterval> F_intervals; // дек порождённых гиперинтервалов
    std::deque<TPoint> F_points; // дек порождённых точек

    EncodedCoordinates new_coord_u;
    EncodedCoordinates new_coord_v;

    uint F_dimension; // размерность задачи
    uint F_constraints; // число ограничений в задаче
    uint F_queueDepth; // глубина очереди оценок констант Липшица

    uint F_divide_axis; // какую ось должны делить?

    uint F_generated_points; // число уже сгенерированных методом точек
    uint F_generated_intervals; // число уже сгенерированных методом гиперинтервалов

    double F_eps; // параметр остановки
    double F_current_minimum; // текущее минимальное значение целевой функции

    TProblem& Fp; // вычислитель целевой функции и функций-ограничений
    std::vector<double> F_globalLipshEvaluations; // дек очередей для глобальных оценок констант Липшица
public:
    TMethodDivByThree() = delete;
    TMethodDivByThree(const uint& out_dim, const uint& out_constr, const uint& depth, const double& out_eps, TProblem& out_prob);
    ~TMethodDivByThree() {};

    virtual void initialization(); // создать самый первый гиперинтервал
    virtual bool divideInterval(const uint& id_divHyp); // поделить гиперинтервал
    virtual void fillIntervals(THyperinterval& parent, const uint& id_u, const uint& id_v);
    void compute_diagonal(const uint& id_Hyp);
    virtual void compute_characteristic(const uint& id_Hyp);
    void compute_evaluations(const uint& out_idPoint);
    uint get_new_id() { return F_generated_points++; } // выдача идентификатора точки
    uint get_id_coord() const { return F_generated_points * F_dimension; } // получить итендификатор для координат
    uint get_id_evaluations() const { return F_generated_points * (F_constraints + 1); }
    uint get_new_interval() { return F_generated_intervals++; } // выдача идентификатора гиперинтервала

    // uint does_point_exist(TPoint& parent, const uint& des_val, const TPoint::direction& dir); // существует ли точка в базе данных уже
};

#endif