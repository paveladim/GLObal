#ifndef SIMPLEPMWITHCONJ_H
#define SIMPLEPMWITHCONJ_H

#include "DivideByThree.h"

class SimplePMwithConj : public DivideByThree {
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
private:
	void give_borders(double& l, double& r, Hyperinterval& hyp);
	void balance(double& _lipshConst) const;
public:
	SimplePMwithConj(const uint& dimension,
					 const uint& constraints,
					 Parameters& parameters,
					 Problem& problem);
	void solve() override;
	~SimplePMwithConj() {}
};

#endif // SIMPLEPMWITHCONJ_H
