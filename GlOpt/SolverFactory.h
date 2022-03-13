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
#include "TransformPM.h"
#include "TransformPMwithConj.h"
#include "TransformPMwithSM.h"

using std::shared_ptr;

static shared_ptr<DivideByThree> create_solver(const std::string& name,
											   const uint& dim,
											   const uint& cst,
											   Parameters& param,
											   Problem& problem) {
	if (name == "ConjugateS")
		return std::make_shared<SimplePMwithConj>(dim, cst, param, problem);
	if (name == "SimplexMethodS")
		return std::make_shared<SimplePMwithSM>(dim, cst, param, problem);
	if (name == "LagrangeS")
		return std::make_shared<SimplePMwithLagrange>(dim, cst, param, problem);
	if (name == "Uniform")
		return std::make_shared<Uniform>(dim, cst, param, problem);
	if (name == "ConjugateT")
		return std::make_shared<TransformPMwithConj>(dim, cst, param, problem);
	if (name == "SimplexMethodT")
		return std::make_shared<TransformPMwithSM>(dim, cst, param, problem);

	return nullptr;
}

#endif // SOLVER_FACTORY_H
