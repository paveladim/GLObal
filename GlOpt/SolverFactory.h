#ifndef SOLVER_FACTORY_H
#define SOLVER_FACTORY_H

#include <string>
#include <memory>

#include "DivideByThree.h"
#include "Uniform.h"
#include "SimplePM.h"
#include "SimplePMwithConj.h"
#include "SimplePMwithSM.h"
#include "SimplePMwithLagrange.h"

using std::shared_ptr;

static shared_ptr<DivideByThree> create_solver(const std::string& name,
											   const uint& dim,
											   const uint& cst,
											   Parameters& param,
											   Problem& problem) {
	if (name == "Conjugate") 
		return std::make_shared<DivideByThree>(SimplePMwithConj(dim, cst, param, problem));
	if (name == "SimplexMethod") 
		return std::make_shared<DivideByThree>(SimplePMwithSM(dim, cst, param, problem));
	if (name == "Lagrange") 
		return std::make_shared<DivideByThree>(SimplePMwithLagrange(dim, cst, param, problem));
	if (name == "Uniform") 
		return std::make_shared<DivideByThree>(Uniform(dim, cst, param, problem));

	return nullptr;
}

#endif // SOLVER_FACTORY_H
