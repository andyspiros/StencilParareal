#ifndef CONVECTION_H_
#define CONVECTION_H_

#include "SharedInfrastructure.h"

#ifndef NDEBUG
#  include "MatFile.h"
#endif

class Stencil;

// Define data field
static const int convectionBoundaryLines = 2;
typedef DataFieldIJBoundary<-convectionBoundaryLines, convectionBoundaryLines, -convectionBoundaryLines, convectionBoundaryLines> ConvectionIJBoundary;

#ifdef __CUDA_BACKEND__
typedef DataFieldCUDA<Real, DataFieldStorageFormat<ConvectionIJBoundary, StorageOrder::KJI, CUDARealAlignment> > ConvectionField;
#else
typedef DataFieldOpenMP<Real, DataFieldStorageFormat<ConvectionIJBoundary, StorageOrder::JIK, OpenMPAlignment> > ConvectionField;
#endif

class Convection
{
public:
    Convection(int isize, int jsize, int ksize, Real dx, Real nu0, Real nufreq, Real cx, Real cy, Real cz);
    ~Convection();

    /**
     * Performs RK4 timesteps
     *
     * Performs the integration of the equation from tstart to tend with
     * the given timestep size dt.  The input field is preserved and the result
     * is stored into the given output field.
     *
     * \param inputField The field whence the condition at tstart is read
     * \param outputField The field where the solution at ttend is written
     * \param dt The timestep size
     * \param tstart The initial time
     * \param tend The end time
     *
     * \return The actual end time ttend is returned.  This is the first number
     *         ttend = min(tstart + k * dt),  for integer k
     */
    void DoRK4(ConvectionField& inputField, ConvectionField& outputField,
                 double tStart, double dt, int timesteps);

    /**
     * Performs Euler timesteps
     *
     * Performs the integration of the equation from tstart to tend with
     * the given timestep size dt.  The input field is preserved and the result
     * is stored into the given output field.
     *
     * \param inputField The field whence the condition at tstart is read
     * \param outputField The field where the solution at ttend is written
     * \param dt The timestep size
     * \param tstart The initial time
     * \param tend The end time
     *
     * \return The actual end time ttend is returned.  This is the first number
     *         ttend = min(tstart + k * dt),  for integer k
     */
    void DoEuler(ConvectionField& inputField, ConvectionField& outputField,
                   double tStart, double dt, int timesteps);

    /**
     * Performs a single RK4 timestep
     */
    void DoRK4Timestep(ConvectionField& inputField, ConvectionField& outputField,
                       double t, double dt);
    /**
     * Performs a single Euler timestep
     */
    void DoEulerTimestep(ConvectionField& inputField, ConvectionField& outputField,
                         double t, double dt);


private:
    typedef JokerDataField<ConvectionField> JokerField; 

    inline double computenu(double t) const
    {
        return nu0_ * (1. + .5 * std::sin(nufreq_ * t));
    }

    void InitializeStencils(const IJKSize& domain);

    // Scalars
    Real nu0_, nufreq_, nu_, cx_, cy_, cz_, dx_, dx2_;
    Real nuparam_;
    Real dtparam_;

    // Internal data fields
    ConvectionField qInternal_;
    ConvectionField k1_, k2_, k3_, k4_;

    // Joker for fields
    JokerField qMainIn_, qMainOut_;
    JokerField qRHSIn_, qEulerOut_;
    JokerField k_;

    // Stencils
    Stencil *rhs2Stencil_, *rhs4Stencil_, *eulerStencil_, *rkStencil_;

#ifndef NDEBUG
    MatFile matfile;
#endif
};

#endif // CONVECTION_H

