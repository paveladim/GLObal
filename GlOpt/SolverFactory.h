#ifndef SOLVER_FACTORY_H
#define SOLVER_FACTORY_H

#include <string>

#include "DivideByThree.h"
#include "Uniform.h"
#include "SimplePMwithConj.h"
#include "SimplePMwithSM.h"
#include "SimplePMwithLagrange.h"

static DivideByThree* create_solver(const std::string& name,
									const uint& dim,
									const uint& cst,
									Parameters& param,
									Problem& problem) {
	if (name == "Conjugate") return new SimplePMwithConj(dim, cst, param, problem);
	if (name == "SimplexMethod") return new SimplePMwithSM(dim, cst, param, problem);
	if (name == "Lagrange") return new SimplePMwithLagrange(dim, cst, param, problem);
	if (name == "Uniform") return new Uniform(dim, cst, param, problem);

	return nullptr;
}

#endif // SOLVER_FACTORY_H
