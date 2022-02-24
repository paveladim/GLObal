#ifndef TRANSFORMPM_H
#define TRANSFORMPM_H

#include <iostream>

#include "DivideByThree.h"

using std::cout;
using std::endl;

class TransformPM : public DivideByThree {
protected:
    bool _areAllCharInfty;
    bool _doesGlobalChange;
    FunctionsValues _localLipshEval;
    FunctionsValues _globalLipshEval;
protected:
	double golden_ratio(double a, double b, const double& e, const uint& id_hyp);
	double calculate_residual(const double& t, const double& e, const uint& id_hyp);
    virtual void calculate_localLipshConst(const uint& id_hyp) = 0;
	void calculate_characteristic(const uint& id_hyp) override;
    void calculate_globalLipshConst(const uint& id_hyp);
    double mixedLipEval(const Hyperinterval& hyp, const uint& i);
    void update_all_charact();
    uint optimal_to_trisect() override;
    uint iterate(const uint& id_hyp) override;
public:
	TransformPM(const uint& dimension,
                const uint& constraints,
                Parameters& parameters,
                Problem& problem);
	void solve();
    ~TransformPM() {}
};

#endif // TRANSFORMPM_H

