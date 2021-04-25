#include <iostream>
#include <vector>
#include "TMethodDivByThree.h"

FunctionsValues& f(FunctionsValues& res, const CoordinatesValues& x)
{
	res[0] = x[0] * x[0] + 10 * x[1] * x[1];
	res[1] = 5 * x[0] + x[1] * x[1] - 5;
	return res;
}

FunctionsValues& f1(FunctionsValues& res, const CoordinatesValues& x)
{
	res[0] = x[0] * x[0] + 10 * x[1] * x[1];
	res[1] = 5 * x[0] + x[1] * x[1] - 5;
	return res;
}

int main() {
	THyperinterval test;
	uint dim{ 2 };
	uint cst{ 1 };
	uint dep{ 3 };
	TProblem testProblem1(dim, cst, { 0, 0 }, { 1, 1 }, &f);
	//TProblem testProblem2(dim, cst, { 0, 0 }, { 1, 1 }, &f1);

	TMethodDivByThree testMethod(dim, cst, dep, 0.001, 100, 100, testProblem1);
	testMethod.createFirstInterval();
	testMethod.divideInterval(0);
	testMethod.compute_evaluations(1);

	return 0;
}