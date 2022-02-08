#ifndef SIMPLEPMWITHLAGRANGE_H
#define SIMPLEPMWITHLAGRANGE_H

#include "SimplePM.h"

class SimplePMwithLagrange : public SimplePM {
private:
	void calculate_localLipshConst(const uint& id_hyp) override;
	void calculate_characteristic(const uint& id_hyp) override;
public:
	SimplePMwithLagrange(const uint& dimension,
						 const uint& constraints,
						 Parameters& parameters,
						 Problem& problem);
	~SimplePMwithLagrange() {}
};

#endif // SIMPLEPMWITHLAGRANGE_H
