#ifndef SIMPLEPMWITHSM_H
#define SIMPLEPMWITHSM_H

#include "DivideByThree.h"
#include "SimplexMethod.h"

class SimplePMwithSM : protected DivideByThree {
	bool _areAllCharInfty;
	bool _doesGlobalChange;
	FunctionsValues _localLipshEval;
	FunctionsValues _globalLipshEval;
private:
	void calculate_localLipshConst(const uint& id_hyp);
	void calculate_globalLipshConst(const uint& id_hyp);
	void calculate_characteristic(const uint& id_hyp) override;
	double mixedLipEval(const Hyperinterval& hyp, const uint& i);
	void update_all_charact();
	uint optimal_to_trisect() override;
	uint iterate(const uint& id_hyp) override;
	void solve();
private:
	void calculate_and_project(const EncodedCoordinates& out,
							   std::vector<double>& incs,
							   const uint& axis);
	double scalar_product(const std::vector<double>& a,
						  const std::vector<double>& b);
	void generate_simplex_table(std::vector<std::vector<double>>& A,
								const std::vector<double>& incs);
	void generate_right_part(std::vector<double>& b,
							 const uint& function,
							 const uint& eval_a,
							 const uint& eval_v,
							 const uint& eval_u,
							 const uint& eval_b);
	void give_borders(double& l, double& r, Hyperinterval& hyp);
public:
};

#endif // SIMPLEPMWITHSM_H
