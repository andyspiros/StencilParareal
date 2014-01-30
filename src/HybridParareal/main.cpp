int main()
{
    ConvectionCPU convectionCPU;
    ConvectionGPU convectionGPU;

    double *pInit = convectionCPU.pInit();
    double *pFine = convectionCPU.pFine();
    double *pQ    = convectionCPU.pQ();
    int dataSize = convectionCPU.dataSize();

    MPI_Request req;
    MPI_Status status;

    /*
     * Part 1: initialization
     */

    convectionCPU.DoInit();
    convectionGPU.UploadAndStartFine();
    convectionCPU.DoCoarse();

    for (int k = 0; k < kmax; ++k)
    {
        convectionGPU.StartDownload();
        convectionCPU.SwapCoarse();
        convectionGPU.EnsureUpload();
        convectionCPU.Receive();

        if (k < kmax-1)
            convectionGPU.UploadAndStartFine();
            
        convectionCPU.DoCoarse();
        convectionGPU.EnsureDownload();

        if (k > 0)
            convectionCPU.EnsureSend();

        convectionCPU.DoUpdate();
        convectionCPU.StartSend();
    }

    convectionCPU.EnsureSend();
}
