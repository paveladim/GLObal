#ifndef TRANSFORMPMWITHCONJ_H
#define TRANSFORMPMWITHCONJ_H

#include "TransformPM.h"

class TransformPMwithConj : public TransformPM {
private:
	void calculate_localLipshConst(const uint& id_hyp) override;
	void calculate_characteristic(const uint& id_hyp) override;
public:
	TransformPMwithConj(const uint& dimension,
						const uint& constraints,
						Parameters& parameters,
						Problem& problem);
	~TransformPMwithConj() {}
};

#endif // TRANSFORMPMWITHCONJ_H