function [ output_args ] = timingParareal(nprocs, kmax)

fnameFine = sprintf('result_%d_%d_fine.dat', nprocs, kmax);
fnameCoarse = sprintf('result_%d_%d_coarse.dat', nprocs, kmax);
fnameSend = sprintf('result_%d_%d_send.dat', nprocs, kmax);
fnameRecv = sprintf('result_%d_%d_recv.dat', nprocs, kmax);

Fine = load(fnameFine);
Coarse = load(fnameCoarse);
Send = load(fnameSend);
Recv = load(fnameRecv);

%% First figure: mean and std of tasks
figure

subplot(1,4,1)
hist(Fine);
title 'Fine propagation'
xlabel 'msec'

subplot(1,4,2)
hist(Coarse);
title 'Coarse propagation'
xlabel 'msec'

subplot(1,4,3)
hist(Send);
title 'MPI\_Wait for Send'
xlabel 'msec'

subplot(1,4,4)
hist(Recv);
title 'MPI\_Wait for Recv'
xlabel 'msec'

%% Reshape to matrices
Fine = reshape(Fine, kmax, nprocs-2);
Coarse = reshape(Coarse, kmax, nprocs-2);
Send = reshape(Send, kmax-1, nprocs-2);
Recv = reshape(Recv, kmax, nprocs-2);


%% Second figure: relevance of task for each process
figure

FineP = sum(Fine, 1);
CoarseP = sum(Coarse, 1);
SendP = sum(Send, 1);
RecvP = sum(Recv, 1);

PP = 1:(nprocs-2);

%subplot(1,4,1)
bar(PP, [FineP' CoarseP' SendP' RecvP'], 'stacked');
legend({
    'Fine propagation',
    'Coarse propagation',
    'Wait for Send',
    'Wait for Recv'
    });
title 'Relevance of task for rank'
xlabel 'Rank'
ylabel 'msec'


%% Third figure: relevance of task per iteration
figure

FineP = sum(Fine, 2);
CoarseP = sum(Coarse, 2);
SendP = [0;sum(Send, 2)];
RecvP = sum(Recv, 2);

PP = 1:kmax;

%subplot(1,4,1)
bar(PP', [FineP CoarseP SendP RecvP], 'stacked');
legend({
    'Fine propagation',
    'Coarse propagation',
    'Wait for Send',
    'Wait for Recv'
    });
title 'Relevance of task for iteration'
xlabel 'Iteration'
ylabel 'msec'
end