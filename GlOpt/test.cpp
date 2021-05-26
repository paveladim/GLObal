#include <iostream>
#include <vector>
#include "TMethodDivByThree.h"
#include "TPiyavskiiMethod.h"

FunctionsValues& f(FunctionsValues& res, const CoordinatesValues& x)
{
	res[0] = x[0] * x[0] + x[1] * x[1];
	res[1] = 5 * x[0] + x[1] * x[1] - 5;
	return res;
}

int main() {
	THyperinterval test;
	uint dim{ 2 };
	uint cst{ 1 };
	uint dep{ 3 };
	TProblem testProblem1(dim, cst, { 0, 0 }, { 1, 1 }, &f);
	GainLipshConstant obj_const{ 1 };
	GainLipshConstant cst_const{ 1 };
	double beta = 0.5;

	TPiyavskiiMethod testMethod(dim, cst, dep, testProblem1, obj_const, cst_const, beta);
	testMethod.launch_method();

	return 0;
}