#ifndef SYNONYMOUS_TYPES_H
#define SYNONYMOUS_TYPES_H

#include <vector>
#include <deque>
#include <queue>
#include <list>
#include <functional>
#include <cmath>
#include <iomanip>

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
using Vec = std::vector<double>;

const double phi = 0.6180339887498948;

#endif