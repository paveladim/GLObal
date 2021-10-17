#define _USE_MATH_DEFINES
#include <cmath>

#include <iostream>
#include <vector>
#include "TMethodDivByThree.h"
#include "TPiyavskiiMethod.h"
#include "TSimplePMwithoutSM.h"
#include "Matrix.h"
#include "SimplexMethod.h"

FunctionsValues& f(FunctionsValues& res, const CoordinatesValues& x)
{
	res[0] = (x[0] - 0.7) * (x[0] - 0.7) + (x[1] - 0.2) * (x[1] - 0.2);
	res[1] = (x[0] - 0.1) * (x[0] - 0.1) + (x[1] - 0.8) * (x[1] - 0.8);
	return res;
}

FunctionsValues& g(FunctionsValues& res, const CoordinatesValues& x)
{
	res[0] = 0.01 * (x[0] * x[1] + (x[0] - M_PI) * (x[0] - M_PI) + 3 * (x[1] - M_PI) * (x[1] - M_PI));
	res[0] = res[0] - sin(x[0]) * sin(x[0]) * sin(2 * x[1]) * sin(2 * x[1]);
	res[1] = (x[0] - 0.1) * (x[0] - 0.1) + (x[1] - 0.8) * (x[1] - 0.8);
	return res;
}

int main() {
	THyperinterval test;
	uint dim{ 2 };
	uint cst{ 1 };
	uint dep{ 3 };
	TProblem testProblem1(dim, cst, { 0, 0 }, { 1, 1 }, &g);
	GainLipshConstant obj_const{ 1 };
	GainLipshConstant cst_const{ 1 };
	double beta = 1e-8;
	double eps = 1e-3;
	double nu = 1e-8;

	//TMethodDivByThree testMethod(dim, cst, dep, testProblem1);
	//TPiyavskiiMethod testMethod(dim, cst, dep, testProblem1, obj_const, cst_const, beta, eps);
	TSimplePMwithoutSM testMethod(dim, cst, dep, testProblem1, obj_const, cst_const, beta, eps, nu);
	testMethod.launch_method();
	testMethod.write_generated_points_to_file();
	testMethod.write_intervals_to_file();

	return 0;
}