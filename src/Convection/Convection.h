#ifndef CONVECTION_H_
#define CONVECTION_H_

#include "SharedInfrastructure.h"

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
    Convection(int isize, int jsize, int ksize, Real dx, Real nu, Real cx, Real cy, Real cz);
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
    double DoRK4(ConvectionField& inputField, ConvectionField& outputField,
                 double dt, unsigned timesteps);

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
    double DoEuler(ConvectionField& inputField, ConvectionField& outputField,
                   double dt, unsigned timesteps);

    /**
     * Performs a single RK4 timestep
     */
    void DoRK4Timestep(ConvectionField& inputField, ConvectionField& outputField,
                       double dt);
    /**
     * Performs a single Euler timestep
     */
    void DoEulerTimestep(ConvectionField& inputField, ConvectionField& outputField,
                         double dt);


private:
    typedef JokerDataField<ConvectionField> JokerField; 

    void InitializeStencils(const IJKSize& domain);

    // Scalars
    Real nu_, cx_, cy_, cz_, dx_, dx2_;
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
};

#endif // CONVECTION_H

