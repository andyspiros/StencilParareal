#ifndef PARAREAL_H_
#define PARAREAL_H_

#include "SharedInfrastructure.h"
#include "StencilFramework.h"
#include "mpi.h"
#include "RuntimeConfiguration.h"

/*
 * Update stage and parameters
 */


enum
{
    _q, _qfine, _qcoarsenew, _qcoarseold
};

template<typename TEnv>
struct UpdateStage
{
    STENCIL_STAGE(TEnv)

    STAGE_PARAMETER(FullDomain, _q)
    STAGE_PARAMETER(FullDomain, _qfine)
    STAGE_PARAMETER(FullDomain, _qcoarseold)
    STAGE_PARAMETER(FullDomain, _qcoarsenew)

    __ACC__
    static void Do(Context ctx, FullDomain)
    {
        ctx[_q::Center()] =
            ctx[_qcoarsenew::Center()] + ctx[_qfine::Center()] - ctx[_qcoarseold::Center()];
    }
};


template<typename PropagatorT, typename FieldType>
class Parareal
{
#ifdef __CUDA_BACKEND__
    typedef BlockSize<32, 4> BSize;
#else
    typedef BlockSize<8, 8> BSize;
#endif

public:
    Parareal(PropagatorT& propagator, FieldType& initial, FieldType& solution,
            double tStart,
             const RuntimeConfiguration& conf, MPI_Comm comm
            )
        : propagator_(propagator), kmax_(conf.kmax())
        , tsFine_(conf.timeStepsFinePerTimeSlice()), tsCoarse_(conf.timeStepsCoarsePerTimeSlice())
        , tStart_(tStart)
        , dtFine_(conf.dtFine()), dtCoarse_(conf.dtCoarse())
        , comm_(comm)
        , async_(conf.async())
        , qinitial_(initial), q_(solution)
    {
        MPI_Comm_size(comm, &commsize_);
        MPI_Comm_rank(comm, &commrank_);
        isFirst_ = commrank_ == 0;
        isLast_ = commrank_ == commsize_-1;

        // Get size and boundary
        const IJKSize& domain = initial.calculationDomain();
        const IJKBoundary& boundary = initial.boundary();
        KBoundary kboundary;
        kboundary.Init(boundary.kMinusOffset(), boundary.kPlusOffset());

        // Initialize internal fields
        qcoarseold_.Init("qcoarseold", domain, kboundary);
        qcoarsenew_.Init("qcoarsenew", domain, kboundary);
        qfine_.Init("qfine", domain, kboundary);

        // Initialize pointers
        pRecv_ = qinitial_.storage().pStorageBase();
        pSend_ = q_.storage().pStorageBase();

        // Data size
        const IJKSize& psize = q_.storage().paddedSize();
        dataSize_ = psize.iSize() * psize.jSize() * psize.kSize();

        // Initialize the stencil
        StencilCompiler::Build(
            updateStencil_,
            "UpdateStencil",
            domain,
            StencilConfiguration<Real, BlockSize<32, 4> >(),
            pack_parameters(
                // Input fields
                Param<_qfine, cIn>(qfine_),
                Param<_qcoarseold, cIn>(qcoarseold_),
                Param<_qcoarsenew, cIn>(qcoarsenew_),
                // Output fields
                Param<_q, cInOut>(q_)
            ),
            define_loops(
                define_sweep<cKIncrement>(
                    define_stages(
                        StencilStage<UpdateStage, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                    )
                )
            )
        );

        tStartFine_ = tStart_ + dtFine_;
    }

    void DoParallel()
    {
        // Initialize
        propagator_.DoEuler(qinitial_, qinitial_, 0., dtCoarse_, tsCoarse_*commrank_);
        propagator_.DoEuler(qinitial_, qcoarseold_, tStart_, dtCoarse_, tsCoarse_);

        for (int k = 0; k < kmax_; ++k)
        {
            // Fine propagation
            if (async_ && !isFirst_)
            {
                propagator_.DoRK4(qinitial_, qfine_, tStart_, dtFine_, 1);
                MPI_Irecv(pRecv_, dataSize_, MPI_DOUBLE, commrank_-1, k, comm_, &reqRecv_);
                propagator_.DoRK4(qfine_, qfine_, tStartFine_, dtFine_, tsFine_-1);
            }
            else
            {
                propagator_.DoRK4(qinitial_, qfine_, tStart_, dtFine_, tsFine_);
            }


            // Receive data
            if (async_ && !isFirst_)
                MPI_Wait(&reqRecv_, &status_);
            else if (!isFirst_)
            {
                SynchronizeCUDA();
                MPI_Recv(pRecv_, dataSize_, MPI_DOUBLE, commrank_-1, k, comm_, &status_);
                SynchronizeCUDA();
            }


            // Coarse propagation
            propagator_.DoEuler(qinitial_, qcoarsenew_, tStart_, dtCoarse_, tsCoarse_);


            // Update solution
            if (async_ && !isLast_ && k > 0)
                MPI_Wait(&reqSend_, &status_);
            updateStencil_.Apply();


            // Send solution
            if (async_ && !isLast_)
            {
                MPI_Isend(pSend_, dataSize_, MPI_DOUBLE, commrank_+1, k, comm_, &reqSend_);
            }
            else if(!isLast_)
            {
                SynchronizeCUDA();
                MPI_Send(pSend_, dataSize_, MPI_DOUBLE, commrank_+1, k, comm_);
                SynchronizeCUDA();
            }

            qcoarseold_.SwapWith(qcoarsenew_);
        }


        // Wait for last send
        if (async_ && !isLast_)
            MPI_Wait(&reqSend_, &status_);
    }

    void DoTimedParallel(std::vector<double>& times)
    {
        times.resize(5*kmax_ + 2);

        double e = MPI_Wtime();

        // Initialize
        propagator_.DoEuler(qinitial_, qinitial_, 0., dtCoarse_, tsCoarse_*commrank_);
        propagator_.DoEuler(qinitial_, qcoarseold_, tStart_, dtCoarse_, tsCoarse_);

        times[0] = MPI_Wtime() - e;

        for (int k = 0; k < kmax_; ++k)
        {
            // Fine propagation
            if (async_ && !isFirst_)
            {
                propagator_.DoRK4(qinitial_, qfine_, tStart_, dtFine_, 1);
                MPI_Irecv(pRecv_, dataSize_, MPI_DOUBLE, commrank_-1, k, comm_, &reqRecv_);
                propagator_.DoRK4(qfine_, qfine_, tStartFine_, dtFine_, tsFine_-1);
            }
            else
            {
                propagator_.DoRK4(qinitial_, qfine_, tStart_, dtFine_, tsFine_);
            }

            times[5*k + 1] = MPI_Wtime() - e;

            // Receive data
            if (async_ && !isFirst_)
                MPI_Wait(&reqRecv_, &status_);
            else if (!isFirst_)
            {
                SynchronizeCUDA();
                MPI_Recv(pRecv_, dataSize_, MPI_DOUBLE, commrank_-1, k, comm_, &status_);
                SynchronizeCUDA();
            }

            times[5*k + 2] = MPI_Wtime() - e;


            // Coarse propagation
            propagator_.DoEuler(qinitial_, qcoarsenew_, tStart_, dtCoarse_, tsCoarse_);
            times[5*k + 3] = MPI_Wtime() - e;


            // Update solution
            if (async_ && !isLast_ && k > 0)
                MPI_Wait(&reqSend_, &status_);
            updateStencil_.Apply();

            times[5*k + 4] = MPI_Wtime() - e;

            // Send solution
            if (async_ && !isLast_)
            {
                MPI_Isend(pSend_, dataSize_, MPI_DOUBLE, commrank_+1, k, comm_, &reqSend_);
            }
            else if(!isLast_)
            {
                SynchronizeCUDA();
                MPI_Send(pSend_, dataSize_, MPI_DOUBLE, commrank_+1, k, comm_);
                SynchronizeCUDA();
            }
            times[5*k + 5] = MPI_Wtime() - e;

            qcoarseold_.SwapWith(qcoarsenew_);
        }


        // Wait for last send
        if (async_ && !isLast_)
            MPI_Wait(&reqSend_, &status_);

        times[5*kmax_ + 1] = MPI_Wtime() - e;
    }

    void DoSerial()
    {
        propagator_.DoRK4(qinitial_, q_, tStart_, dtFine_, tsFine_*commsize_);
    }

private:
    PropagatorT& propagator_;
    int kmax_;
    int tsFine_, tsCoarse_;
    double tStart_, tStartFine_;
    double dtFine_, dtCoarse_;
    MPI_Comm comm_;
    int commsize_, commrank_;
    bool isFirst_, isLast_;
    bool async_;

    double *pRecv_, *pSend_;
    int dataSize_;
    MPI_Request reqRecv_, reqSend_;
    MPI_Status status_;

    FieldType& qinitial_;
    FieldType& q_;
    FieldType qcoarseold_, qcoarsenew_, qfine_;
    Stencil updateStencil_;



    // Private methods

    /**
     * Synchronizes the CUDA device in the case of a GPU build
     */
    inline void SynchronizeCUDA()
    {
#ifdef __CUDA_BACKEND__
        cudaDeviceSynchronize();
#endif
    }
};

#endif // PARAREAL_H_

