#include "TProblem.h"

TProblem::TProblem(uint dim, uint constrCount, const CoordinatesValues& left_brd,
    const CoordinatesValues& right_brd, const FunctionsCalculator& prob)
    :F_dimension(dim), F_constraints_count(constrCount),
    F_left_borders(left_brd), F_right_borders(right_brd), Fp(prob)
{
    Fx_tmp.resize(F_dimension);
    Ffuncs_tmp.resize(F_constraints_count + 1);
}

CoordinatesValues& TProblem::decode_coordinates(const EncodedCoordinates& out)
{
    CoordinatesValues& left = F_left_borders;
    CoordinatesValues& right = F_right_borders;

    for (uint i = 0; i < F_dimension; ++i) {
        Fx_tmp[i] = (double)out[i] / (double)MAX_POWER_THREE;
        Fx_tmp[i] = Fx_tmp[i] * (right[i] - left[i]) + left[i];
    }

    return Fx_tmp;
}

FunctionsValues& TProblem::F(const CoordinatesValues& out)
{
    return Fp(Ffuncs_tmp, out);
}

FunctionsValues& TProblem::operator()(const EncodedCoordinates& out)
{
    Fx_tmp = decode_coordinates(out);
    return F(Fx_tmp);
}