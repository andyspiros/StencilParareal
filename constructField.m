function qq = constructField(directory)

gettuple = @(i, Px)[idivide(i,int32(Px^2)), idivide(mod(i,int32(Px^2)),int32(Px)), mod(i,int32(Px))];

listdir = what(directory);
listdir = listdir(end);
results = 0;

for i = 1:length(listdir.mat)
    fname = listdir.mat{i};
    if (length(fname) < 10)
        continue;
    end
    
    if (fname(1:6) == 'result')
        results = results +1;
    end
end

if (results < 1)
    fprintf('No results found');
    return
end

Px = int32(results ^ (1/3));

% Load first result
qq = load(sprintf('%s/result_0_0_0.mat', directory), 'q');
qq = qq.q;
localgrid = size(qq, 1) - 6;
gridsize = localgrid * Px;
qq = zeros(gridsize, gridsize, gridsize);

r = 4:(localgrid+3);

for p = 0:(results-1)
    I = gettuple(p, Px);
    i = I(1);
    j = I(2);
    k = I(3);
    fname = sprintf('%s/result_%d_%d_%d.mat', directory, i, j, k);
    q = load(fname, 'q');
    q = q.q;
    iStart = i*localgrid + 1;
    jStart = j*localgrid + 1;
    kStart = k*localgrid + 1;
    iEnd = (i+1)*localgrid;
    jEnd = (j+1)*localgrid;
    kEnd = (k+1)*localgrid;
    qq(iStart:iEnd, jStart:jEnd, kStart:kEnd) = q(r, r, r);
end
