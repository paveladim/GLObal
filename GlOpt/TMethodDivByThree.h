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

    CoordinatesValues coord_a;
    CoordinatesValues coord_b;

    uint F_dimension; // размерность задачи
    uint F_constraints; // число ограничений в задаче
    uint F_queueDepth; // глубина очереди оценок констант Липшица

    uint F_divide_axis; // какую ось делим

    uint F_generated_points; // число уже сгенерированных методом точек
    uint max_generated_points; // максимальное допустимое число точек, можно сгенерировать

    double F_eps; // параметр остановки
    double F_current_minimum; // текущее минимальное значение целевой функции
    uint F_generated_intervals; // число уже сгенерированных методом гиперинтервалов
    uint max_generated_intervals; // максимальное допустимое число гиперинтервалов, которое можно сгенерировать

    TProblem& Fp; // вектор целевой функции и функций-ограничений
    std::vector<std::queue<double>> global_evaluations; // дек очередей для глобальных оценок констант Липшица
public:
    TMethodDivByThree() = delete;
    TMethodDivByThree(const uint& out_dim, const uint& out_constr, const uint& depth, const double& out_eps,
        const uint& max_gen_points, const uint& max_gen_interv, TProblem& out_prob);
    ~TMethodDivByThree() {};

    virtual void createFirstInterval(); // создать самый первый гиперинтервал
    virtual bool divideInterval(const uint& id_divHyp); // поделить гиперинтервал
    virtual void fillIntervals(THyperinterval& parent, const uint& id_u, const uint& id_v);
    void compute_diagonal(const uint& id_Hyp);
    virtual void compute_characteristic(const uint& id_Hyp);
    void compute_evaluations(const uint& out_idPoint);
    // virtual uint compute_new_coord(const uint& out_id, bool f); // сгенерировать новые координаты
    uint get_new_id(); // выдача идентификатора точки
    uint get_id_coord(); // получить итендификатор для координат
    uint get_new_interval(); // выдача идентификатора гиперинтервала
    void delete_point(); // удалить идентификатор точки
    void delete_hyperinterval(); // удалить итендификатор гиперинтервала

    uint does_point_exist(TPoint& parent, const uint& des_val, bool f); // существует ли точка в базе данных уже
};

#endif