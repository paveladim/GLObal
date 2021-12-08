#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#define _USE_MATH_DEFINES
#include <cmath>

#include "synonymous_types.h"

FunctionsValues& f(FunctionsValues& res, const CoordinatesValues& x)
{
	res[0] = (x[0] - 0.7) * (x[0] - 0.7) + (x[1] - 0.2) * (x[1] - 0.2);
	res[1] = -(x[0] - 0.1) * (x[0] - 0.1) - (x[1] - 0.8) * (x[1] - 0.8);
	// res[1] = 1;
	// res[1] = 5 * (x[0] - 0.3) * (x[0] - 0.3) * (x[0] - 0.3);
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

FunctionsValues& task3(FunctionsValues& res, const CoordinatesValues& x) {
	res[0] = 0.01 * (x[0] * x[1] + (x[0] - M_PI) * (x[0] - M_PI) + 3 * (x[1] - M_PI) * (x[1] - M_PI));
	res[0] = res[0] - sin(x[0]) * sin(x[0]) * sin(2 * x[1]) * sin(2 * x[1]);
	res[1] = 6 - (x[0] - M_PI + 0.1) * (x[0] - M_PI + 0.1);
	res[1] = res[1] - (2 * sin(x[1]) + 0.2) * (2 * sin(x[1]) + 0.2);

	return res;
}

FunctionsValues& task4(FunctionsValues& res, const CoordinatesValues& x) {
	res[0] = 0;
	for (uint i = 0; i < x.size(); ++i)
		res[0] = res[0] + (x[i] * x[i] - cos(18 * x[i] * x[i]));

	res[1] = 0.1 - x[0];
	for (uint i = 1; i < x.size(); ++i)
		res[1] = res[1] + std::abs(x[i]);

	return res;
}

FunctionsValues& task5(FunctionsValues& res, const CoordinatesValues& x) {
	double prev = 0;
	double next = 0;

	next = 1 + 0.25 * (x[0] - 1);
	res[0] = 10 * sin(M_PI * next) * sin(M_PI * next);
	for (uint i = 0; i < x.size() - 2; ++i) {
		prev = next;
		next = 1 + 0.25 * (x[i + 1] - 1);
		res[0] = res[0] + (prev - 1) * (prev - 1) * (1 + 10 * sin(M_PI * next) * sin(M_PI * next));
	}
	res[0] = res[0] + (next - 1) * (next - 1);
	res[0] = res[0] * M_PI / x.size();

	res[1] = 0.4 - 10 * (x[0] - 1.1) * (x[0] - 1.1);
	for (uint i = 1; i < x.size(); ++i)
		res[1] = res[1] + (i + 1.0) * (x[i] - 1) * (x[i] - 1);

	return res;
}

FunctionsValues& task6(FunctionsValues& res, const CoordinatesValues& x) {
	res[0] = 0.01 * (x[0] * x[1] + (x[0] - M_PI) * (x[0] - M_PI) + 3 * (x[1] - M_PI) * (x[1] - M_PI));
	res[0] = res[0] - sin(x[0]) * sin(x[0]) * sin(2 * x[1]) * sin(2 * x[1]);
	res[1] = 6 - (x[0] - M_PI + 0.1) * (x[0] - M_PI + 0.1);
	res[1] = res[1] - (2 * sin(x[1]) + 0.2) * (2 * sin(x[1]) + 0.2);

	if (res[1] > 0) res[0] = res[0] - 1;

	return res;
}

#endif // FUNCTIONS_H