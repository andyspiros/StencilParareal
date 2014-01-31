#ifndef RUNTIMECONFIGURATION_H_
#define RUNTIMECONFIGURATION_H_

#include <iostream>
#include <cstdlib>
#include <cmath>

enum RunMode
{
    ModeCompare, ModeSerial, ModeParallel
};

class RuntimeConfiguration
{
public:
    RuntimeConfiguration(int argc, char **argv);

    // Getters
    double nu() const { return nu_; }
    double cx() const { return cx_; }
    double cy() const { return cy_; }
    double cz() const { return cz_; }
    int timeSlices() const { return timeSlices_; }
    int timeStepsFine() const { return timeStepsFine_; }
    int gridSize() const { return gridSize_; }
    double endTime() const { return endTime_; }
    double kmax() const { return kmax_; }
    bool mat() const { return mat_; }
    bool async() const {  return async_; }
    RunMode mode() const { return mode_; }

    // Setters
    void set_nu(double x) { nu_ = x; }
    void set_cx(double x) { cx_ = x; }
    void set_cy(double x) { cy_ = x; }
    void set_cz(double x) { cz_ = x; }
    void set_gridSize(int x) { gridSize_ = x; }
    void set_endTime(double x) { endTime_ = x; }
    void set_timeStepsFine(int x) { timeStepsFine_ = x; }
    void set_timeStepsCoarse(int x) { timeStepsCoarse_ = x; }
    void set_kmax(int x) { kmax_ = x; }
    void set_mat(bool x) { mat_ = x; }
    void set_async(bool x) { async_ = x; }
    void set_mode(RunMode x) { mode_ = x; }


    // Non-stored info

    /**
     * \return The spatial step size is returned
     */
    double dx() const
    {
        return 1. / gridSize_;
    }

    /**
     * \return The timestep used by the fine propagator is returned
     */
    double dtFine() const
    {
        return endTime_ / timeStepsFine_;
    }

    /**
     * \return The timestep used by the coarse propagator is returned
     */
    double dtCoarse() const
    {
        return endTime_ / timeStepsCoarse_;
    }

    int timeStepsFinePerTimeSlice() const {
        return timeStepsFine_ / timeSlices_;
    }

    int timeStepsCoarsePerTimeSlice() const
    {
        return timeStepsCoarse_ / timeSlices_;
    }


    /**
     * \return The size of a time slice is returned
     */
    double timeSliceSize() const
    {
        return endTime_ / timeSlices_;
    }

    /**
     * \return The CFL number used by the fine propagator is returned
     */
    double cflFine() const { return cfl(dtFine()); }

    /**
     * \return The CFL number used by the coarse propagator is returned
     */
    double cflCoarse() const { return cfl(dtCoarse()); }

private:

    inline double cfl(double dt) const
    {
        const double c = std::max(std::max(cx_, cy_), cz_);
        const double dt_dx = dt / dx();
        const double dt_dx2 = dt_dx / dx();
        return std::max(nu_*dt_dx2, c*dt_dx);
    }

    double nu_;
    double cx_, cy_, cz_;
    int gridSize_;
    double endTime_;
    int timeSlices_;
    int timeStepsFine_, timeStepsCoarse_;
    int kmax_;
    bool mat_;
    bool async_;
    RunMode mode_;
};

#endif // RUNTIMECONFIGURATION_H_

