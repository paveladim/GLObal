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

	uint dep{ 3 };
	Problem testProblem(2, 1, { 0.0, 0.0 }, { 1.0, 1.0 }, &f);
	Problem testProblem1(dim_t1, cst_t1, { 0.0, -1.0 }, { 4.0, 3.0 }, &task1);
	Problem testProblem2(dim_t2, cst_t2, { -2.5, -1.5 }, { 2.5, 1.5 }, &task2);
	Problem testProblem3(dim_t3, cst_t3, { 0, 0 }, { 2 * M_PI, 2 * M_PI }, &task3);
	Problem testProblem4(dim_t4, cst_t4, { -1.0, -1.0 }, { 1.0, 1.0 }, &task4);
	Problem testProblem5(dim_t5, cst_t5, { -10, -10 }, { 10, 10 }, &task5);
	Problem testProblem6(dim_t3, cst_t3, { 0, 0 }, { 2 * M_PI, 2 * M_PI }, &task6);

	uint dim{ dim_t1 };
	uint cst{ cst_t1 };
	double localObj{ 2.5 };
	double globalObj{ 2.5 };
	double localCst{ 2.5 };
	double globalCst{ 2.5 };
	double delta{ 1e-10 };
	double beta{ 0.3 };
	double eps{ 1e-6 };
	double diag{ (double)MAX_POWER_THREE };
	eps = eps * sqrt((double)dim) * diag;

	double best_min{ 10.0 };
	uint best_chal{ 500 };
	double bestlocalObj{ 1.5 };
	double bestglobalObj{ 1.5 };
	double bestlocalCst{ 1.5 };
	double bestglobalCst{ 1.5 };
	
	/*while (globalObj < 4.0) {
		while (globalCst < 4.0) {
			while (localObj < 4.0) {
				while (localCst < 4.0) {
					std::cout << localObj << " " << localCst << " " << globalObj << " " << globalCst << std::endl;
					Parameters param{ dim, cst, dep, localObj, localCst, globalObj, globalCst, delta, beta * diag, eps };
					std::shared_ptr<DivideByThree> solver(create_solver("Lagrange", dim, cst, param, testProblem1));
					solver->solve();

					if (solver->get_min() < best_min) {
						best_min = solver->get_min();
						best_chal = solver->get_gen();
						bestlocalObj = localObj;
						bestglobalObj = globalObj;
						bestlocalCst = localCst;
						bestglobalCst = globalCst;
					}

					localCst += 0.2;
				}
				localCst = 1.5;
				localObj += 0.2;
			}
			localObj = 1.5;
			localCst = 1.5;
			globalCst += 0.2;
		}
		localObj = 1.5;
		localCst = 1.5;
		globalCst = 1.5;
		globalObj += 0.2;
	} 

	std::cout << best_chal << std::endl;
	std::cout << best_min << std::endl;
	std::cout << bestlocalObj << std::endl;
	std::cout << bestlocalCst << std::endl;
	std::cout << bestglobalObj << std::endl;
	std::cout << bestglobalCst << std::endl; */

	Parameters param{ dim, cst, dep, localObj, localCst, globalObj, globalCst, delta, beta * diag, eps };
	std::shared_ptr<DivideByThree> solver(create_solver("SimplexMethod", dim, cst, param, testProblem1));
	solver->solve();

	solver->write_generated_points();
	solver->write_generated_intervals();

	return 0;
}