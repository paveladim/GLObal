#define _USE_MATH_DEFINES
#include <cmath>

#include <iostream>
#include <vector>
#include "synonymous_types.h"
#include "Functions.h"
#include "Problem.h"

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
	Problem testProblem1(dim_t1, cst_t1, { 0.0, -1.0 }, { 4.0, 3.0 }, &task1);
	Problem testProblem2(dim_t2, cst_t2, { -2.5, -1.5 }, { 2.5, 1.5 }, &task2);
	Problem testProblem3(dim_t3, cst_t3, { 0, 0 }, { 2 * M_PI, 2 * M_PI }, &task3);
	Problem testProblem4(dim_t4, cst_t4, { -1.0, -1.0 }, { 1.0, 1.0 }, &task4);
	Problem testProblem5(dim_t5, cst_t5, { -10, -10 }, { 10, 10 }, &task5);
	Problem testProblem6(dim_t3, cst_t3, { 0, 0 }, { 2 * M_PI, 2 * M_PI }, &task6);

	GainLipshConstant obj_const{ 1.5 };
	GainLipshConstant cst_const{ 1.5 };
	double beta = 0.4;
	double eps = 1e-5;

	//TSimplePMwithConstraints solvetask1(dim_t1, cst_t1, dep, testProblem1, obj_const, cst_const, beta, eps);
	//TSimplePMwithConstraints solvetask2(dim_t2, cst_t2, dep, testProblem2, obj_const, cst_const, beta, eps);
	TSimplePMwithConstraints solvetask3(dim_t3, cst_t3, dep, testProblem3, obj_const, cst_const, beta, eps);
	TSimplePMwithConstraints solvetask4(dim_t4, cst_t4, dep, testProblem4, obj_const, cst_const, beta, eps);
	TSimplePMwithConstraints solvetask5(dim_t5, cst_t5, dep, testProblem5, obj_const, cst_const, beta, eps);
	TSimplePMwithConstraints solvetask6(dim_t3, cst_t3, dep, testProblem6, obj_const, cst_const, beta, eps);

	while (obj_const < 4.0) {
		while (cst_const < 4.0) {
			TSimplePMwithConstraints solvetask2(dim_t2, cst_t2, dep, testProblem2, obj_const, cst_const, beta, eps);
			std::cout << obj_const << std::endl;
			std::cout << cst_const << std::endl;
			solvetask2.launch_method();
			//solvetask5.write_generated_points_to_file();
			//solvetask5.write_intervals_to_file();
			cst_const += 0.1;
		}
		cst_const = 1.5;
		obj_const += 0.1;
	}

	return 0;
}