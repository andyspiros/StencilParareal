#ifndef RUNTIMECONFIGURATION_H_
#define RUNTIMECONFIGURATION_H_

#include <iostream>
#include <cstdlib>
#include <cmath>

class RuntimeConfiguration
{
public:
    RuntimeConfiguration(int argc, char **argv);

    // Getters
    double nu() const { return nu_; }
    double cx() const { return cx_; }
    double cy() const { return cy_; }
    double cz() const { return cz_; }
    int gridSize() const { return gridSize_; }
    double endTime() const { return endTime_; }
    double cflFine() const { return cflFine_; }
    double cflCoarse() const { return cflCoarse_; }
    double kmax() const { return kmax_; }

    // Setters
    void set_nu(double x) { nu_ = x; }
    void set_cx(double x) { cx_ = x; }
    void set_cy(double x) { cy_ = x; }
    void set_cz(double x) { cz_ = x; }
    void set_gridSize(int x) { gridSize_ = x; }
    void set_endTime(double x) { endTime_ = x; }
    void set_cflFine(double x) { cflFine_ = x; }
    void set_cflCoarse(double x) { cflCoarse_ = x; }
    void set_kmax(int x) { kmax_ = x; }


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
        const double dx = this->dx();
        return dx*dx * cflFine(); 
    }

    /**
     * \return The timestep used by the coarse propagator is returned
     */
    double dtCoarse() const
    {
        const double dx = this->dx();
        return dx*dx * cflCoarse(); 
    }

private:
    double nu_;
    double cx_, cy_, cz_;
    int gridSize_;
    double endTime_;
    double cflFine_, cflCoarse_;
    int kmax_;
};

#endif // RUNTIMECONFIGURATION_H_

