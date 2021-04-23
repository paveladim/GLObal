#ifndef TPROBLEM_H
#define TPROBLEM_H

#include "synonymous_types.h"
#include "dicretization.h"

class TProblem
{
public:
    enum class TProblemTypeError { OK };
private:
    uint F_dimension; // размерность задачи
    uint F_constraints_count; // число ограничений задачи
    CoordinatesValues F_left_borders; // вектор левых границ функции
    CoordinatesValues F_right_borders; // вектор правых границ функции
    FunctionsCalculator Fp; // для вычислений функции
    CoordinatesValues Fx_tmp;
    FunctionsValues Ffuncs_tmp;
    TProblemTypeError Ferror{ TProblemTypeError::OK };
public:
    TProblem() = delete;
    TProblem(uint dim, uint constrCount, const CoordinatesValues& left_brd,
        const CoordinatesValues& right_brd, const FunctionsCalculator& prob);
    ~TProblem() {}

    CoordinatesValues& decode_coordinates(const EncodedCoordinates& out); // расшифровать координаты
    FunctionsValues& F(const CoordinatesValues& out); // вычислить функцию в точке
    FunctionsValues& operator()(const EncodedCoordinates& out);
};

#endif