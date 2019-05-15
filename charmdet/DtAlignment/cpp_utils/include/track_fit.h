#ifndef DTALIGNMENT_TRACKFIT_H
#define DTALIGNMENT_TRACKFIT_H

#include "TObject.h"
#include <iostream>

namespace Dt_alignment
{
class Linear_Fitter : public TObject
{
public:
	Linear_Fitter();
	virtual ~Linear_Fitter();

	ClassDef(Linear_Fitter,1);
};
}

#endif //DTALIGNMENT_TRACKFIT_H
