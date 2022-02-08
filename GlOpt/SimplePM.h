#ifndef SIMPLEPM_H
#define SIMPLEPM_H

#include "DivideByThree.h"

class SimplePM : public DivideByThree {
protected:
	bool _areAllCharInfty;
	bool _doesGlobalChange;
	FunctionsValues _localLipshEval;
	FunctionsValues _globalLipshEval;
protected:
	virtual void calculate_localLipshConst(const uint& id_hyp) = 0;
	void calculate_globalLipshConst(const uint& id_hyp);
	double mixedLipEval(const Hyperinterval& hyp, const uint& i);
	void update_all_charact();
	uint optimal_to_trisect() override;
	uint iterate(const uint& id_hyp) override;
protected:
	void give_borders(double& l, double& r, Hyperinterval& hyp);
	void balance(double& _lipshConst) const;
public:
	SimplePM(const uint& dimension,
			 const uint& constraints,
			 Parameters& parameters,
			 Problem& problem);
	void solve() override;
	~SimplePM() {}
};

#endif SIMPLEPM_H
