#define _USE_MATH_DEFINES
#include <cmath>

#include <iostream>
#include <vector>
#include "TMethodDivByThree.h"
#include "TPiyavskiiMethod.h"
#include "TSimplePMnoConstraints.h"
#include "TSimplePMwithConstraints.h"
#include "Matrix.h"
#include "SimplexMethod.h"

FunctionsValues& f(FunctionsValues& res, const CoordinatesValues& x)
{
	res[0] = (x[0] - 0.7) * (x[0] - 0.7) + (x[1] - 0.2) * (x[1] - 0.2);
	// res[1] = -(x[0] - 0.1) * (x[0] - 0.1) - (x[1] - 0.8) * (x[1] - 0.8);
	// res[1] = 1;
	res[1] = 5 * (x[0] - 0.3) * (x[0] - 0.3) * (x[0] - 0.3);
	return res;
}

FunctionsValues& g(FunctionsValues& res, const CoordinatesValues& x) {
	res[0] = 0.01 * (x[0] * x[1] + (x[0] - M_PI) * (x[0] - M_PI) + 3 * (x[1] - M_PI) * (x[1] - M_PI));
	res[0] = res[0] - sin(x[0]) * sin(x[0]) * sin(2 * x[1]) * sin(2 * x[1]);
	res[1] = (x[0] - 0.1) * (x[0] - 0.1) + (x[1] - 0.8) * (x[1] - 0.8);
	return res;
}

FunctionsValues& task1(FunctionsValues& res, const CoordinatesValues& x) {
	res[0] = -1.5 * x[0] * x[0] * exp(1 - x[0] * x[0] - 20.25 * (x[0] - x[1]) * (x[0] - x[1]));
	res[0] = res[0] - 0.5 * (x[0] - 1) * (x[1] - 1) * 0.5 * (x[0] - 1) * (x[1] - 1) *
			 0.5 * (x[0] - 1) * (x[1] - 1) * 0.5 * (x[0] - 1) * (x[1] - 1) *
			 exp(2 - 0.5 * (x[0] - 1) * 0.5 * (x[0] - 1) * 0.5 * (x[0] - 1) * 0.5 * (x[0] - 1) -
			 (x[1] - 1) * (x[1] - 1) * (x[1] - 1) * (x[1] - 1));
	res[1] = 0.001 * ((x[0] - 2.2) * (x[0] - 2.2) + 
					  (x[1] - 1.2) * (x[1] - 1.2) - 2.25);
	res[2] = 100 * (1 - ((x[0] - 2.0) / 1.2) * ((x[0] - 2.0) / 1.2) - 0.25 * x[1] * x[1]);
	res[3] = 10 * (x[1] - 1.5 - 1.5 * sin(2 * M_PI * (x[0] - 1.75)));
	return res;
}

FunctionsValues& task2(FunctionsValues& res, const CoordinatesValues& x) {
	res[0] = (4 - 2.1 * x[0] * x[0] + x[0] * x[0] * x[0] * x[0] / 3) * x[0] * x[0];
	res[0] = res[0] + x[0] * x[1] + (4 * x[1] * x[1] - 4) * x[1] * x[1];
	res[1] = -(1.5 * x[0] - x[1] - 0.2) * (1.5 * x[0] - x[1] - 0.2);
	res[1] = res[1] - (2 * sin(2 * x[1]) + 0.2) * (2 * sin(2 * x[1]) + 0.2) + 7;
	res[2] = -14 + abs(x[0] + 0.1) * abs(x[0] + 0.1) * abs(x[0] + 0.1);
	res[2] = res[2] + 2 * abs(x[1] - 0.2) * abs(x[1] - 0.2) * abs(x[1] - 0.2);

	return res;
}

int main() {
	uint dim_t1{ 2 };
	uint cst_t1{ 3 };
	uint dim_t2{ 2 };
	uint cst_t2{ 2 };

	uint dep{ 3 };
	TProblem testProblem1(dim_t1, cst_t1, { 0.0, -1.0 }, { 4.0, 3.0 }, &task1);
	TProblem testProblem2(dim_t2, cst_t2, { -2.5, -1.5 }, { 2.5, 1.5 }, &task2);

	GainLipshConstant obj_const{ 3.0 };
	GainLipshConstant cst_const{ 2.5 };
	double beta = 0.25;
	double eps = 1e-5;

	TSimplePMwithConstraints solvetask1(dim_t1, cst_t1, dep, testProblem1, obj_const, cst_const, beta, eps);
	TSimplePMwithConstraints solvetask2(dim_t2, cst_t2, dep, testProblem2, obj_const, cst_const, beta, eps);

	solvetask1.launch_method();
	solvetask1.write_generated_points_to_file();
	solvetask1.write_intervals_to_file();

	return 0;
}