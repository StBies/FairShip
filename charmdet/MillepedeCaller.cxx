#include "MillepedeCaller.h"
#include <iostream>


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
: m_mille(outFileName, asBinary, writeZero)
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
	m_mille.mille(n_local_derivatives,local_derivatives,n_global_derivatives,global_derivatives,label,measured_residual,sigma);
}


/**
 * Call the end function of Mille. This stores a set of local parameters to disk.
 * Local parameters are specific to one (straight) track for our context. So after
 * the processing of one track this method sould be called before beginnung to process
 * the next track.
 *
 * @brief Wrapper for Mille's method end() call.
 *
 * @author Stefan Bieschke
 * @date May 21,2019
 * @version 1.0
 */
void MillepedeCaller::call_end()
{
	m_mille.end();
}

/**
 * Getter for the Mille object stored in this call wrapper. This way, you can manilpulate
 * the original Mille object.
 *
 * @brief Raw Mille obj getter
 *
 * @author Stefan Bieschke
 * @date May 21, 2019
 * @version 1.0
 *
 * @returns const reference to stored Mille object
 */
const Mille& MillepedeCaller::get_mille() const
{
	return m_mille;
}
