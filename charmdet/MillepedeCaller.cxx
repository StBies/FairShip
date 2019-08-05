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

//TODO document
vector<gbl::GblPoint> MillepedeCaller::list_hits(const genfit::Track* track) const
{
	std::vector<gbl::GblPoint> result = {};

	size_t n_points = track->getNumPointsWithMeasurement();
	vector<genfit::TrackPoint* > points = track->getPointsWithMeasurement();

	//TODO calculate measurements and add before sorting
	multimap<double,TMatrixD*> jacobians_with_arclen = jacobians_with_arclength(track);

	for(auto it = jacobians_with_arclen.begin(); it != jacobians_with_arclen.end(); it++)
	{
		TMatrixD* jacobian = it->second;
		result.push_back(gbl::GblPoint(*jacobian));
	}

	return result;
}



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
	for(unsigned int row = 0; row < jacobian->GetNrows(); row++)
	{
		for(unsigned int col = 0; col < jacobian->GetNcols(); col++)
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

//TODO document - especially 0 < hit_id < max trackpoints
pair<double,TMatrixD*> MillepedeCaller::single_jacobian_with_arclength(const genfit::Track& track, const unsigned int hit_id) const
{
	TVector3 fitted_pos_1 = track.getFittedState(hit_id - 1).getPos();
	TVector3 fitted_pos_2 = track.getFittedState(hit_id).getPos();
	TVector3 between_hits = fitted_pos_2 - fitted_pos_1;

	double distance = between_hits.Mag();

	TMatrixD* jacobian = calc_jacobian(track, hit_id - 1, hit_id);

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
	TVector3 wire_dir(wire_top - wire_bot);

	TVector3 plane_pos(track_pos - wire_bot);
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
