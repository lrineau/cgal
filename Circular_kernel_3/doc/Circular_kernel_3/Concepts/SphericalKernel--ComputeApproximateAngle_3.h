
/*!
\ingroup PkgCircularKernel3GeometricConcepts
\cgalConcept
*/

class SphericalKernel::ComputeApproximateAngle_3 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
Computes an approximation of the angle of the arc in radian 
 and between \f$\left[0,2\pi\right] \f$. 
*/ 
double 
operator()(const SphericalKernel::Circular_arc_3 & a); 

///@}

}; /* end SphericalKernel::ComputeApproximateAngle_3 */

