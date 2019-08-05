#ifndef CHARMDET_MILLEPEDECALLER_H_
#define CHARMDET_MILLEPEDECALLER_H_

#include "TObject.h"
#include "Mille.h"
//includes for GBL fitter from genfit
#include <vector>
#include "GblPoint.h"
#include "GblTrajectory.h"
#include "MilleBinary.h"
#include "Track.h"
#include "TMatrixD.h"
#include <map>
#include "TVector3.h"

/**
 * A class for wrapping the millepede function call such that it can be called from
 * within a python script
 *
 * @author Stefan Bieschke
 * @date Apr. 9, 2019
 */
class MillepedeCaller//: public TObject
{
public:
	MillepedeCaller(const char *outFileName, bool asBinary = true, bool writeZero = false);
	~MillepedeCaller();

	void call_mille(int n_local_derivatives,
					const float *local_derivatives,
					int n_global_derivatives,
					const float *global_derivatives,
					const int *label,
					float measured_residual,
					float sigma);

	std::vector<gbl::GblPoint> list_hits(genfit::Track* track) const;
	const int* labels() const;

	ClassDef(MillepedeCaller,2);

private:
	Mille mille;

	//helper methods
	TMatrixD* calc_jacobian(genfit::Track* track, const unsigned int hit_id_1, const unsigned int hit_id_2) const;
	std::multimap<double,TMatrixD*> jacobians_with_arclength(genfit::Track* track) const;
};

#endif /* CHARMDET_MILLEPEDECALLER_H_ */
