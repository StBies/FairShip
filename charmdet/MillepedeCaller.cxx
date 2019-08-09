#include "MillepedeCaller.h"
#include <iostream>

using namespace std;


/**
 * Constructor. Initializes the wrapper to the Mille class. All paramaters of this constructor
 * are passed to the constructor of the Mille class.
 *
 * @brief Constructor
 *
 * @author Stefan Bieschke
 * @date Apr 9, 2019
 * @version 1.0
 *
 * @param outFileName string containing the filename written by mille function.
 * @param asBinary Flag that states if file should be written as binary
 * @param writeZero Flag that states if zero values should be kept or not
 */
MillepedeCaller::MillepedeCaller(const char *outFileName, bool asBinary, bool writeZero)
: mille(outFileName, asBinary, writeZero)
{

}

/**
 * Default destructor. Tears down the object
 *
 * @brief Destructor
 */
MillepedeCaller::~MillepedeCaller()
{

}

/**
 * Call the mille function of the Mille object contained in this class. For a documentation of the parameters, see the Millepede documentation at:
 * http://www.desy.de/~kleinwrt/MP2/doc/html/draftman_page.html
 *
 * @brief Mille caller
 *
 * @author Stefan Bieschke
 * @date May 9, 2019
 * @version 1.0
 *
 * @param n_local_derivatives Number of local derivatives.
 * @param local_derivatives Pointer to values of the local derivatives. Array of floats as long as n_local_derivatives states.
 * @param n_global_derivatives Number of global derivatives.
 * @param global_derivatives Pointer to values of the global derivatives. Array of floats as long as n_global_derivatives states.
 * @param label Array of unique labels to identify global parameters. This must be as long as n_global_derivatives states.
 * @param measured_residual Residual of the hit for that mille is called.
 * @param sigma Sigma of this hit.
 */
void MillepedeCaller::call_mille(int n_local_derivatives,
					const float *local_derivatives,
					int n_global_derivatives,
					const float *global_derivatives,
					const int *label,
					float measured_residual,
					float sigma)
{
	mille.mille(n_local_derivatives,local_derivatives,n_global_derivatives,global_derivatives,label,measured_residual,sigma);
}

/**
 * List the hits used to fit the seed track from a fit predecessing the GBL fit as a @c std::vector<gbl::GblPoint>. The GblPoint objects will
 * also contain measurements after the call of this method. The GBL points are ordered by arclength, which is the distance on the track between
 * two consecutive hits.
 *
 * @brief List hits from seed track as vector<GblPoint>
 *
 * @author Stefan Bieschke
 * @date Aug. 06, 2019
 * @version 1.0
 *
 * @param track Seed track, in this case genfit::Track from Kalman fitter
 *
 * @return std::vector<gbl::GblPoint> containing the hits ordered by arclen with measurement added
 */
vector<gbl::GblPoint> MillepedeCaller::list_hits(const genfit::Track* track) const
{
	std::vector<gbl::GblPoint> result = {};

	vector<genfit::TrackPoint* > points = track->getPointsWithMeasurement();
	size_t n_points = points.size();

	//define a struct to handle track parameters at a certain point as well as measurement and residual
	struct hit_info
	{
		TMatrixD* jacobian;
		double rt_measurement;
		TVector3 closest_approach;
	};

	multimap<double,struct hit_info,less<double>> jacobians_with_arclen;

	//#pragma omp parallel for
	for(size_t i = 1; i < n_points; i++)
	{
		struct hit_info hit;
		genfit::TrackPoint* point = points[i];
		//get coords of upper and lower end of hit tube
		genfit::AbsMeasurement* raw_measurement = point->getRawMeasurement();
		TVectorD raw = raw_measurement->getRawHitCoords();
		TVector3 vbot(raw[0],raw[1],raw[2]);
		TVector3 vtop(raw[3],raw[4],raw[5]);
		double measurement = raw[6]; //rt distance [cm]

		TVector3 fit_pos = track->getFittedState(i).getPos();
		TVector3 fit_mom = track->getFittedState(i).getMom();
		TVector3 closest_approach = calc_shortest_distance(vtop,vbot,fit_pos,fit_mom);
		pair<double,TMatrixD*> jacobian_with_arclen = single_jacobian_with_arclength(*track,i);

		//hit struct seems to be copied correctly into multimap
		hit.jacobian = jacobian_with_arclen.second;
		hit.closest_approach = closest_approach;
		hit.rt_measurement = measurement;
		jacobians_with_arclen.insert(make_pair(jacobian_with_arclen.first,hit));
	}
	//TODO calculate GBLpoint for hit no. 0

	for(auto it = jacobians_with_arclen.begin(); it != jacobians_with_arclen.end(); it++)
	{
		TMatrixD* jacobian = it->second.jacobian;
		//TODO test if the GblPoint constructor stores a copy so that original jacobian can be deleted
		result.push_back(gbl::GblPoint(*jacobian));
		TRotation rot = calc_rotation_of_vector(it->second.closest_approach);
		TMatrixD rot_mat = rot_to_matrix(rot);
		//TODO fix, fit not working
		TVectorD rotated_residual(3);
		rotated_residual[0] = it->second.closest_approach.Mag() - it->second.rt_measurement;
		rotated_residual[1] = 0;
		rotated_residual[2] = 0;
		TVectorD precision(rotated_residual);
		precision[0] = 250 * 1e-4; //250 um in cm
		result.back().addMeasurement(rot_mat,rotated_residual,precision);
	}

	return result;
}


//TODO implement
const int* MillepedeCaller::labels() const
{
	return new int[100];
}


//TODO add mathematical description of matrix parameter determination to doc comment
/**
 * Calculate the jacobian between two consecutive hits on a predetermined track. This track (in this case) is the result of a predecessing
 * fit with a KalmanFitter from genfit. It can be rewritten to be any track, even a "guessed" one from any kind of pattern recognition.
 *
 * @brief Calculate jacobian between two consecutive hits
 *
 * @author Stefan Bieschke
 * @version 1.0
 * @date July 30, 2019
 *
 * @param track genfit::Track object, of whom the jacobian for two consecutive hits is calculated
 * @param hit_id_1 ID of the first hit on the track for that the jacobian should be calculated
 * @param hit_id_2 ID of the second hit on the track for that the jacobian should be calculated
 *
 * @return Pointer to a heap object of TMatrixD type with dimensions 5x5
 *
 * @warning Requires hit_id_1 = hit_id_2 - 1
 * @warning Requires hit_id_1 < track.getNumPointsWithMeasurement() && hit_id_2 < track.getNumPointsWithMeasurement()
 * @warning Heap object without auto deletion
 */
TMatrixD* MillepedeCaller::calc_jacobian(const genfit::Track* track, const unsigned int hit_id_1, const unsigned int hit_id_2) const
{
	TMatrixD* jacobian = new TMatrixD(5,5);

	// 1.) init unity matrix
	for(uint8_t row = 0; row < jacobian->GetNrows(); row++)
	{
		for(uint8_t col = 0; col < jacobian->GetNcols(); col++)
		{
			if(row == col)
			{
				(*jacobian)[row][col] = 1;
			}
			else
			{
				(*jacobian)[row][col] = 0;
			}
		}

	}

	//2.) enter non-zero partial differentials
	//2.1) get the two points on track where reconstruction happened
	genfit::MeasuredStateOnPlane state_at_id_1 = track->getFittedState(hit_id_1);
	genfit::MeasuredStateOnPlane state_at_id_2 = track->getFittedState(hit_id_2);

	TVector3 pos1 = state_at_id_1.getPos();
	TVector3 pos2 = state_at_id_2.getPos();

	double dx = pos2.X() - pos1.X();
	double dy = pos2.Y() - pos1.Y();

	//2.2) enter dx and dy to jacobian
	(*jacobian)[3][1] = dx;
	(*jacobian)[4][2] = dy;


	return jacobian;
}


//TODO add doc comment
/**
 *
 */
multimap<double,TMatrixD*> MillepedeCaller::jacobians_with_arclength(const genfit::Track* track) const
{
	multimap<double,TMatrixD*,less<double>> result;

	unsigned int n_hits = track->getNumPointsWithMeasurement();


	//TODO add unity matrix as first entry
	for (unsigned int hit_id = 1; hit_id < n_hits; hit_id++)
	{
		result.insert(single_jacobian_with_arclength(*track,hit_id));
	}

	//return a copy of the map. Since it's small it is probably better than handling memory and matrices are stored as pointers to heap
	return result;
}

/**
 * Calculate the jacobian matrix for the transport of track parameters from the previous hit to a hit labelled by @c hit_id. Alongside
 * the jacobian matrix, an arclength is calculated, which is the distance on the track between the previous hit and the hit specified by
 * @c hit_id. Note that as the transport of parameters from the previous hit is calculated, the parameter @c hit_id is required to be
 * greater than zero. See the requirements below for more details.
 * The result is returned as std::pair<double,TMatrixD*>. The TMatrixD* pointer points to heap memory holding the 5x5 jacobian matrix.
 * For more info on the jacobian see the documentation of the method
 * @c calc_jacobian(const genfit::Track* track, const unsigned int hit_id_1, const unsigned int hit_id_2) const.
 *
 * @brief Calculate a single pair of jacobian with arclength for a given @c hit_id
 *
 * @author Stefan Bieschke
 * @date Aug. 06, 2019
 * @version 1.0
 *
 * @param track Seed track from a fit predecessing the GBL refit, e.g the genfit Kalman fitter
 * @param hit_id ID (number) of the hit for that the pair shall be computed. See requirements below
 *
 * @require 0 < hit_id < track.getNumPointsWithMeasurement()
 *
 * @return std::pair<double,TMatrixD*> containing the the arclength and the jacobian matrix
 *
 * @warning TMatrixD* is heap memory and undeleted. Must be deleted by caller when no longer needed.
 */
pair<double,TMatrixD*> MillepedeCaller::single_jacobian_with_arclength(const genfit::Track& track, const unsigned int hit_id) const
{
	TVector3 fitted_pos_1 = track.getFittedState(hit_id - 1).getPos();
	TVector3 fitted_pos_2 = track.getFittedState(hit_id).getPos();
	TVector3 between_hits = fitted_pos_2 - fitted_pos_1;

	double distance = between_hits.Mag();

	TMatrixD* jacobian = calc_jacobian(&track, hit_id - 1, hit_id);

	return make_pair(distance, jacobian);
}

//TODO document
/**
 *
 */
double MillepedeCaller::perform_GBL_refit(const genfit::Track& track) const
{
	vector<gbl::GblPoint> points = list_hits(&track);
	gbl::GblTrajectory traj(points);

	int rc, ndf;
	double chi2, lostWeight;

	cout << "------------performing refit--------------" << endl;
	cout << "Seed track chi2: " << track.getFitStatus()->getChi2() << " Ndf: " << track.getFitStatus()->getNdf() << endl;

	rc = traj.fit(chi2,ndf,lostWeight);
	cout << "Refit chi2: " << chi2 << " Ndf: " << ndf << endl;

	return chi2;
}

//TODO document
//Reimplementation of python function
/**
 *
 */
TVector3 MillepedeCaller::calc_shortest_distance(const TVector3& wire_top, const TVector3& wire_bot, const TVector3& track_pos, const TVector3& track_mom) const
{
	TVector3 wire_dir = wire_top - wire_bot;

	TVector3 plane_pos = track_pos - wire_bot;
	TVector3 plane_dir_1(track_mom);
	TVector3 plane_dir_2(-1 * wire_dir);

	TVectorD const_vector(2);
	TMatrixD coeff_matrix(2,2);

	const_vector[0] = -(plane_pos.Dot(track_mom));
	const_vector[1] = -(plane_pos.Dot(wire_dir));

	coeff_matrix[0][0] = plane_dir_1.Dot(track_mom);
	coeff_matrix[0][1] = plane_dir_2.Dot(track_mom);
	coeff_matrix[1][0] = plane_dir_1.Dot(wire_dir);
	coeff_matrix[1][1] = plane_dir_2.Dot(wire_dir);

	TDecompLU solvable_matrix(coeff_matrix);
	TVectorD result(const_vector);
	int rc = solvable_matrix.Solve(result);

	TVector3 PCA_on_track(track_pos + result[0] * track_mom);
	TVector3 PCA_on_wire(wire_bot + result[1] * wire_dir);

	return TVector3(PCA_on_track - PCA_on_wire);
}

/**
 * Calculate the rotation matrix in a way such that a given vector v is the x-axis of the rotated coordinate frame.
 *
 * @brief Rotation of a given vector in lab frame so that it is new x axis
 *
 * @author Stefan Bieschke
 * @date Aug. 09, 2019
 * @version 1.0
 *
 * @param v TVector3 that is meant to be the new x axis of the rotated coordinate frame
 *
 * @return TRotation containing the rotation matrix
 */
TRotation MillepedeCaller::calc_rotation_of_vector(const TVector3& v) const
{
	TRotation rot;
	rot.SetXAxis(v);

	return rot;
}

TMatrixD MillepedeCaller::rot_to_matrix(const TRotation& rot) const
{
	TMatrixD result(3,3);
	for(uint8_t i = 0; i < 3; i++)
	{
		for(uint8_t j = 0; j < 3; j++)
		{
			result[i][j] = rot[i][j];
		}
	}
	return result;
}
