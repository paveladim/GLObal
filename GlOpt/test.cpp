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

	uint dim{ dim_t3 };
	uint cst{ cst_t3 };
	double globalObj{ 2.5 };
	double globalCst{ 2.5 };
	double localObj{ 1.5 };
	double localCst{ 1.5 };
	double delta{ 1e-10 };
	double beta{ 0.4 };
	double eps{ 1e-8 };
	double diag{ static_cast<double>(MAX_POWER_THREE) };
	eps = eps * sqrt((double)dim) * diag;

	//while (globalObj < 5.0) {
	//	while (globalCst < 5.0) {
	//		while (localObj < 3.0) {
	//			while (localCst < 3.0) {
	//				std::cout << "TASK\t" << globalObj << " " << globalCst << " " << localObj << " " << localCst << endl;
	//				Parameters param{ dim, cst, dep, localObj, localCst, globalObj, globalCst, delta, beta * diag, eps, 1.0, 1.0 };
	//				std::shared_ptr<DivideByThree> solver(create_solver("SimplexMethodT", 2, 1, param, testProblem3));
	//				solver->solve();
	//				cout << "\n\n";
	//				
	//				localCst += 0.5;
	//			}

	//			localCst = 1.5;
	//			localObj += 0.5;
	//		}

	//		localCst = 1.5;
	//		localObj = 1.5;
	//		globalCst += 0.5;
	//	}

	//	localCst = 1.5;
	//	localObj = 1.5;
	//	globalCst = 2.5;
	//	globalObj += 0.5;
	//}

	Parameters param{ dim, cst, dep, localObj, localCst, globalObj, globalCst, delta, beta * diag, eps, 1.0, 1.0 };
	std::shared_ptr<DivideByThree> solver(create_solver("SimplexMethodT", 2, 1, param, testProblem3));

	solver->solve();
	solver->write_generated_points();
	solver->write_generated_intervals();

	return 0;
}