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
	TProblem testProblem1(dim, cst, { 0, 0 }, { 1, 1 }, &f);
	GainConstants gc{ 1, 1 };

	TMethodDivByThree testMethod(dim, cst, dep, 0.001, testProblem1, gc);
	testMethod.initialization();
	std::cout << testMethod.do_step(0) << std::endl;

	return 0;
}