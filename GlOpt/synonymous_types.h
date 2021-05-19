#ifndef SYNONYMOUS_TYPES_H
#define SYNONYMOUS_TYPES_H

#include <vector>
#include <deque>
#include <queue>
#include <list>
#include <functional>
#include <cmath>

using uint = unsigned int;
using EncodedCoordinate = uint;
using CoordinateValue = double;
using FunctionValue = double;
using CoordinatesValues = std::vector<CoordinateValue>;
using EncodedCoordinates = std::vector<EncodedCoordinate>;
using FunctionsValues = std::vector<FunctionValue>;
using FunctionsCalculator =
std::function< FunctionsValues& (FunctionsValues&, const CoordinatesValues&)>;
using GainLipshConstant = double;
using FunctionsCalculators = std::vector<FunctionsCalculator>;
using LipschitzConstantValue = double;
using FeatureValue = double;

#endif