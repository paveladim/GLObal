#define _USE_MATH_DEFINES
#include <cmath>

#include <iostream>
#include <vector>
#include <memory>

#include "synonymous_types.h"
#include "Functions.h"
#include "Problem.h"
#include "Parameters.h"
#include "DivideByThree.h"
#include "SolverFactory.h"

int main() {
	uint dim_t1{ 2 };
	uint cst_t1{ 3 };
	uint dim_t2{ 2 };
	uint cst_t2{ 2 };
	uint dim_t3{ 2 };
	uint cst_t3{ 1 };
	uint dim_t4{ 2 };
	uint cst_t4{ 1 };
	uint dim_t5{ 2 };
	uint cst_t5{ 1 };

	uint dep{ 5 };
	Problem testProblem(2, 1, { 0.0, 0.0 }, { 1.0, 1.0 }, &f);
	Problem testProblem1(dim_t1, cst_t1, { 0.0, -1.0 }, { 4.0, 3.0 }, &task1);
	Problem testProblem2(dim_t2, cst_t2, { -2.5, -1.5 }, { 2.5, 1.5 }, &task2);
	Problem testProblem3(dim_t3, cst_t3, { 0, 0 }, { 2 * M_PI, 2 * M_PI }, &task3);
	Problem testProblem4(dim_t4, cst_t4, { -1.0, -1.0 }, { 1.0, 1.0 }, &task4);
	Problem testProblem5(dim_t5, cst_t5, { -10, -10 }, { 10, 10 }, &task5);
	Problem testProblem6(dim_t3, cst_t3, { 0, 0 }, { 2 * M_PI, 2 * M_PI }, &task6);

	uint dim{ dim_t5 };
	uint cst{ cst_t5 };
	double globalObj{ 4.0 };
	double globalCst{ 4.0 };
	double localObj{ 2.5 };
	double localCst{ 2.5 };
	double delta{ 1e-10 };
	double beta{ 0.4 };
	double eps{ 1e-7 };
	double diag{ (double)MAX_POWER_THREE };
	eps = eps * sqrt((double)dim) * diag;

	/*while (globalObj < 6.0) {
		while (globalCst < 6.0) {
			Parameters param{ dim, cst, dep, localObj, localCst, globalObj, globalCst, delta, beta * diag, eps, 1, 1 };
			std::shared_ptr<DivideByThree> solver(create_solver("SimplexMethod", dim, cst, param, testProblem6));
			solver->solve();
			globalCst += 0.4;
		}
		globalObj += 0.4;
		globalCst = 3.0;
	} */

	Parameters param{ dim, cst, dep, localObj, localCst, globalObj, globalCst, delta, beta * diag, eps, 1.5, 0.5 };
	std::shared_ptr<DivideByThree> solver(create_solver("SimplexMethod", 2, 1, param, testProblem));
	solver->solve();

	solver->write_generated_points();
	solver->write_generated_intervals();

	return 0;
}