/*
 * GblSeedTrack.h
 *
 *  Created on: Feb 11, 2020
 *      Author: bieschke
 */
#include <vector>
#include "TVector3.h"
#include "Track.h"

#ifndef CHARMDET_GBLSEEDTRACK_H_
#define CHARMDET_GBLSEEDTRACK_H_

class GblSeedTrack
{
public:
	GblSeedTrack(const genfit::Track& genfit_seed);
	GblSeedTrack(const std::vector<TVector3>& MC_seed);
	~GblSeedTrack();

private:
	std::vector<std::pair<TVector3,TVector3>> hit_list;
};

#endif /* CHARMDET_GBLSEEDTRACK_H_ */
