#include <iostream>
#include <vector>
#include "TMethodDivByThree.h"

FunctionsValues& f(FunctionsValues& res, const CoordinatesValues& x)
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
	TProblem testProblem(dim, cst, { 0, 0 }, { 1, 1 }, &f);
	TMethodDivByThree testMethod(dim, cst, dep, 0.001, 100, 100, testProblem);
	testMethod.createFirstInterval();
	testMethod.divideInterval(0);

	return 0;
}