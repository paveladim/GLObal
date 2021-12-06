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
	double delta;
	double _criticalSize;
};

#endif // PARAMETERS_H