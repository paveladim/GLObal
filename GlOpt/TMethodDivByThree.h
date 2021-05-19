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

    // поля для обмена координатами
    EncodedCoordinates transit_coord_1;
    EncodedCoordinates transit_coord_2;

    uint F_dimension; // размерность задачи
    uint F_constraints; // число ограничений в задаче
    uint F_queueDepth; // глубина очереди оценок констант Липшица

    uint F_divide_axis; // какую ось должны делить?

    uint F_generated_points; // число уже сгенерированных методом точек
    uint F_generated_intervals; // число уже сгенерированных методом гиперинтервалов

    double F_criticalSize; // пороговый размер гиперинтервала
    double F_current_minimum; // текущее минимальное значение целевой функции
    uint F_id_current_minimum; // итендификатор глобального минимума

    TProblem& Fp; // вычислитель целевой функции и функций-ограничений
    GainLipshConstant F_gainObjective; // константа завышения константы Липшица для целевой функции
    GainLipshConstant F_gainConstraints; // константа завышения констант Липшица для ограничений
    double delta; // параметр для осторожных локальных оценок
    std::vector<double> F_globalLipshEvaluations; // вектор глобальных оценок констант Липшица целевой функции и ограничений
public:
    TMethodDivByThree() = delete;
    TMethodDivByThree(const uint& out_dim, const uint& out_constr, const uint& depth, 
        TProblem& out_prob, const GainLipshConstant& out_gainObj, const GainLipshConstant& out_gainCst, const double& beta);
    ~TMethodDivByThree() {};

    virtual void initialization(); // создать самый первый гиперинтервал
    virtual bool trisect_interval(const uint& id_divHyp); // поделить гиперинтервал
    virtual void fill_intervals(THyperinterval& parent, const uint& id_u, const uint& id_v); // заполнить гиперинтервал данными
    void compute_diagonal(const uint& id_Hyp); // вычислить размер активной диагонали гиперинтервала
    virtual void compute_characteristic(const uint& id_Hyp); // вычислить характеристику гиперинтервала
    void compute_evaluations(const uint& out_idPoint); // вычислить значение функции в точке
    uint get_new_id() { return F_generated_points++; } // выдача идентификатора точки
    uint get_id_coord() const { return F_generated_points * F_dimension; } // получить итендификатор для координат
    // получить итендификатор на начало размещения значений целевой функции и ограничений в точке
    uint get_id_evaluations() const { return F_generated_points * (F_constraints + 1); }
    uint get_new_interval() { return F_generated_intervals++; } // выдача идентификатора гиперинтервала
    void compute_localLipshConst(const uint& id_Hyp); // вычислить локальные оценки констант Липшица
    void compute_globalLipshConst(); // вычислить глобальные оценки констант Липшица
    // выбрать "лучший" гиперинтервал для деления
    uint choose_optimal_to_trisect(); // найти оптимальный для деления на три
    uint do_step(const uint& if_divHyp); // сделать шаг метода
    void launch_method(); // для тестирования
    // методы расширения деков
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