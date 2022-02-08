#ifndef SIMPLEPMWITHCONJ_H
#define SIMPLEPMWITHCONJ_H

#include "SimplePM.h"

class SimplePMwithConj : public SimplePM {
private:
	void calculate_localLipshConst(const uint& id_hyp) override;
	void calculate_characteristic(const uint& id_hyp) override;
public:
	SimplePMwithConj(const uint& dimension,
					 const uint& constraints,
					 Parameters& parameters,
					 Problem& problem);
	~SimplePMwithConj() {}
};

#endif // SIMPLEPMWITHCONJ_H
