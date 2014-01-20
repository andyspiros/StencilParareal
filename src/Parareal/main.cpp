#include <iostream>
#include "RuntimeConfiguration.h"
#include "ConvectionSolution.h"
#include "mpi.h"
#include "Convection.h"
#include "SharedInfrastructure.h"
#include "StencilFramework.h"

enum
{
    q_, qfine_, qcoarsenew_, qcoarseold_
};

template<typename TEnv>
struct UpdateStage
{
    STENCIL_STAGE(TEnv)

    STAGE_PARAMETER(FullDomain, q_)
    STAGE_PARAMETER(FullDomain, qfine_)
    STAGE_PARAMETER(FullDomain, qcoarseold_)
    STAGE_PARAMETER(FullDomain, qcoarsenew_)

    __ACC__
    static void Do(Context ctx, FullDomain)
    {
        ctx[q_::Center()] = 
            ctx[qcoarsenew_::Center()] + ctx[qfine_::Center()] - ctx[qcoarseold_::Center()];
    }
};


int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int commsize, commrank;
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &commrank);
    const bool isRoot = commrank == 0;
    const bool isLast = commrank == commsize - 1;

    RuntimeConfiguration conf(argc, argv);

    // Compute teimsteps
    const double dt = conf.endTime() / commsize;
    double dtFine = conf.dtFine();
    double dtCoarse = conf.dtCoarse();
    int timestepsFine = static_cast<int>(dt / dtFine + .5);
    int timestepsCoarse = static_cast<int>(dt / dtCoarse + .5);
    dtFine = dt / timestepsFine;
    dtCoarse = dt / timestepsCoarse;
    double cflFine = conf.dx()*conf.dx() / dtFine;
    double cflCoarse = conf.dx()*conf.dx() / dtCoarse;

    if (isRoot);
    std::cout << "Running with:\n"
        << " - spatial discretization step: " << conf.dx() << "\n"
        << " - endtime: " << conf.endTime() << "\n"
        << " - heat coefficient: " << conf.nu() << "\n"
        << " - advection velocity in x: " << conf.cx() << "\n"
        << " - advection velocity in y: " << conf.cy() << "\n"
        << " - advection velocity in z: " << conf.cz() << "\n"
        << " - CFL fine: " << cflFine << "\n"
        << " - CFL coarse: " << cflCoarse << "\n"
        << " - timestep size fine: " << dtFine << "\n"
        << " - timestep size coarse: " << dtCoarse << "\n"
        << "\n";

    // Reduce timestepsFine by one
    --timestepsFine;

    // Calculation domain and boundaries
    IJKSize domain; domain.Init(conf.gridSize(), conf.gridSize(), conf.gridSize()); 
    KBoundary kboundary; kboundary.Init(-convectionBoundaryLines, convectionBoundaryLines);

    // Instantiate fields
    ConvectionField q, qfine, qcoarseold, qcoarsenew;
    q.Init("q", domain, kboundary);
    qfine.Init("qfine", domain, kboundary);
    qcoarseold.Init("qcoarseold", domain, kboundary);
    qcoarsenew.Init("qcoarsenew", domain, kboundary);

    // Field with initial conditions for root only,
    // otherwise, compute initial conditions in q
    ConvectionField* qinitial;
    if (isRoot)
    {
        qinitial = new ConvectionField;
        qinitial->Init("qinitial", domain, kboundary);
        fillQ(*qinitial, conf.nu(), conf.cx(), conf.cy(), conf.cz(), 0, 0., 1., 0., 1., 0., 1.);
        qinitial->SynchronizeDeviceStorage();
        q = *qinitial;
    }
    else
        fillQ(q, conf.nu(), conf.cx(), conf.cy(), conf.cz(), 0, 0., 1., 0., 1., 0., 1.);

    q.SynchronizeDeviceStorage();

    // Convectionj object
    Convection convection(conf.gridSize(), conf.gridSize(), conf.gridSize(), conf.dx(), conf.nu(), conf.cx(), conf.cy(), conf.cz());

    // MPI requests for asynchronous communication
    MPI_Request req;
    MPI_Status status;

    // Base storage for q and size of data
    const IJKSize& psize = q.storage().paddedSize();
    const int datasize = psize.iSize()*psize.jSize()*psize.kSize();
    double *qBuffer = q.storage().pStorageBase();

    // Initialize update stencil
    Stencil updateStencil;
    StencilCompiler::Build(
        updateStencil,
        "UpdateStencil",
        domain,
        StencilConfiguration<Real, BlockSize<32, 4> >(),
        pack_parameters(
            // Input fields
            Param<qfine_, cIn>(qfine),
            Param<qcoarseold_, cIn>(qcoarseold),
            Param<qcoarsenew_, cIn>(qcoarsenew),
            // Output fields
            Param<q_, cInOut>(q)
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<UpdateStage, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                )
            )
        )
    );





    /**********************************
     * PARAREAL ALGORITHM STARTS HERE *
     **********************************/


    /*
     * 1. Initialization
     */
    convection.DoEuler(q, q, dtCoarse, timestepsCoarse);
    qcoarseold = q;


    /*
     * Begin of iteration
     */
    for (int k = 0; k < conf.kmax(); ++k)
    {
        /*
         * 2. Fine propagation
         */
        convection.DoRK4Timestep(q, qfine, dtFine);

        if (!isRoot)
            MPI_Irecv(qBuffer, datasize, MPI_DOUBLE, commrank-1, k, MPI_COMM_WORLD, &req);

        convection.DoRK4(qfine, qfine, dtFine, timestepsFine);


        /**
         * 3. Coarse propagation
         */

        // Wait for receive to complete
        if (!isRoot)
        {
            MPI_Wait(&req, &status);
            if (status.MPI_ERROR)
            {
                std::cerr << "MPI error on proc " << commrank << ": " << status.MPI_ERROR << std::endl;
            }
        }

        // Integrate with coarse propagator
        // TODO useless computation for root
        if (!isRoot)
            convection.DoEuler(q, qcoarsenew, dtCoarse, timestepsCoarse);
        else
            convection.DoEuler(*qinitial, qcoarsenew, dtCoarse, timestepsCoarse);

        // Update solution
        updateStencil.Apply();

        // Send to next process
        if (!isLast)
            MPI_Send(qBuffer, datasize, MPI_DOUBLE, commrank+1, k+1, MPI_COMM_WORLD);
    }

    // Delete root field
    if (isRoot)
        delete qinitial;
}

