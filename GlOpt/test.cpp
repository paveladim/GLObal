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
	res[1] = 5 * x[0] + x[1] * x[1] - 5;
	return res;
}

FunctionsValues& g(FunctionsValues& res, const CoordinatesValues& x)
{
	res[0] = 0.01 * (x[0] * x[1] + (x[0] - M_PI) * (x[0] - M_PI) + 3 * (x[1] - M_PI) * (x[1] - M_PI));
	res[0] = res[0] - sin(x[0]) * sin(x[0]) * sin(2 * x[1]) * sin(2 * x[1]);
	res[1] = 5 * x[0] + x[1] * x[1] - 5;
	return res;
}

int main() {
	/*std::vector<Vec> vec{{3.0, 4.0, 1.0, 0.0, 0.0}, {3.0, 1.0, 0.0, 1.0, 0.0}, {1.0, 0.0, 0.0, 0.0, 1.0}};
	Vec c{ -2.0, -1.0, 0.0, 0.0, 0.0 };
	Vec b{ 32.0, 17.0, 5.0 };
	Matrix m(vec);
	SimplexMethod sm(c, b, m);
	sm.find_solution(); */

	THyperinterval test;
	uint dim{ 2 };
	uint cst{ 1 };
	uint dep{ 3 };
	TProblem testProblem1(dim, cst, { 0, 0 }, { 1, 1 }, &g);
	GainLipshConstant obj_const{ 1 };
	GainLipshConstant cst_const{ 1 };
	double beta = 1e-8;
	double eps = 1e-8;
	double nu = 1e-4;

	//TMethodDivByThree testMethod(dim, cst, dep, testProblem1);
	//TPiyavskiiMethod testMethod(dim, cst, dep, testProblem1, obj_const, cst_const, beta, eps);
	TSimplePMwithoutSM testMethod(dim, cst, dep, testProblem1, obj_const, cst_const, beta, eps, nu);
	testMethod.launch_method();
	testMethod.write_generated_points_to_file();
	testMethod.write_intervals_to_file();

	return 0;
}