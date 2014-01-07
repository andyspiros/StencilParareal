function qq = constructField(directory, varargin)

gettuple = @(i, Px)[idivide(i,int32(Px^2)), idivide(mod(i,int32(Px^2)),int32(Px)), mod(i,int32(Px))];
getfname = @(i,j,k) sprintf('%s/result_%d_%d_%d.mat', directory, i-1, j-1, k-1);

listdir = what(directory);
listdir = listdir(end);
results = 0;

expr = 'result_(?<PI>\d+)_(?<PJ>\d+)_(?<PK>\d+)\.mat';
matches = regexp(listdir.mat, expr, 'names');

matfiles = cell(0,0,0);

toti = 0;
totj = 0;
totk = 0;

for i = 1:length(matches)
    m = matches{i};
    if (~isempty(m))
        myPI = str2double(m.PI)+1;
        myPJ = str2double(m.PJ)+1;
        myPK = str2double(m.PK)+1;
        
        m = matfile(getfname(myPI, myPJ, myPK));
        matfiles{myPI, myPJ, myPK} = m;
        s = size(m.q);
        
        if (myPJ == 1 && myPK == 1)
            toti = toti + s(1) - 6;
        end
        if (myPI == 1 && myPK == 1)
            totj = totj + s(2) - 6;
        end
        if (myPJ == 1 && myPI == 1)
            totk = totk + s(3) - 6;
        end
    end
end

qq = zeros(toti, totj, totk);

istart = 1;
jstart = 1;
kstart = 1;
for pi = 1:size(matfiles,1)
    jstart = 1;
    for pj = 1:size(matfiles,2)
        kstart = 1;
        for pk = 1:size(matfiles,3)
            m = matfiles{pi, pj, pk};
            s = size(m.q);
            
            isize = s(1)-6;
            jsize = s(2)-6;
            ksize = s(3)-6;
            
            iend = istart + isize - 1;
            jend = jstart + jsize - 1;
            kend = kstart + ksize - 1;
            
            % Load field into resulting field
            is = 4:(isize+3);
            js = 4:(jsize+3);
            ks = 4:(ksize+3);
            q = m.q;
            qq(istart:iend, jstart:jend, kstart:kend) = q(is, js, ks);
            
            kstart = kend + 1;
        end
        jstart = jend + 1;
    end
    istart = iend + 1;
end

