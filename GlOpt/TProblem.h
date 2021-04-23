#ifndef TPROBLEM_H
#define TPROBLEM_H

#include "synonymous_types.h"
#include "dicretization.h"

class TProblem
{
public:
    enum class TProblemTypeError { OK };
private:
    uint F_dimension; // ����������� ������
    uint F_constraints_count; // ����� ����������� ������
    CoordinatesValues F_left_borders; // ������ ����� ������ �������
    CoordinatesValues F_right_borders; // ������ ������ ������ �������
    FunctionsCalculator Fp; // ��� ���������� �������
    CoordinatesValues Fx_tmp;
    FunctionsValues Ffuncs_tmp;
    TProblemTypeError Ferror{ TProblemTypeError::OK };
public:
    TProblem() = delete;
    TProblem(uint dim, uint constrCount, const CoordinatesValues& left_brd,
        const CoordinatesValues& right_brd, const FunctionsCalculator& prob);
    ~TProblem() {}

    CoordinatesValues& decode_coordinates(const EncodedCoordinates& out); // ������������ ����������
    FunctionsValues& F(const CoordinatesValues& out); // ��������� ������� � �����
    FunctionsValues& operator()(const EncodedCoordinates& out);
};

#endif