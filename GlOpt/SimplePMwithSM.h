#ifndef SIMPLEPMWITHSM_H
#define SIMPLEPMWITHSM_H

#include "DivideByThree.h"
#include "SimplexMethod.h"

class SimplePMwithSM : protected DivideByThree {
	bool _areAllCharInfty;
private:
	void calculate_localLipshConst(const uint& id_hyp);
	void calculate_characteristic(const uint& id_hyp) override;
	uint optimal_to_trisect() override;
	uint iterate(const uint& id_hyp) override;
public:
};

#endif // SIMPLEPMWITHSM_H
