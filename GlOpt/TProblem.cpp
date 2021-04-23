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
    size_t i{ 0 };
    CoordinatesValues& left = F_left_borders;
    CoordinatesValues& right = F_right_borders;

    auto doOperation = [&i, left, right](uint ix)
    {
        return left[i] + (right[i] - left[i]) *
            (((CoordinateValue)ix) / MAX_POWER_THREE);
    };

    std::transform(out.begin(), out.end(), Ffuncs_tmp.begin(), doOperation);
    return Ffuncs_tmp;
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