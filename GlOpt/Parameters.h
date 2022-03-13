#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "synonymous_types.h"

struct Parameters {
	uint _dimension;
	uint _constraints;
	uint _queueDepth;
	GainLipshConstant _gainLocalObj;
	GainLipshConstant _gainLocalCst;
	GainLipshConstant _gainGlobalObj;
	GainLipshConstant _gainGlobalCst;
	double _delta;
	double _criticalSize;
	double _eps;

	GainLipshConstant _Gain;
	GainLipshConstant _Reduce;
};

#endif // PARAMETERS_H