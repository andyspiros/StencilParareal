function q = convectionQ(q, istart, iend, t, nu, cx, cy, cz)
size = iend-istart + 1;
dx = 1 / size;
for i = istart:iend
    for j = istart:iend
        for k = istart:iend
            x = (i-istart+1/2)*dx;
            y = (j-istart+1/2)*dx;
            z = (k-istart+1/2)*dx;
            
            q(i, j, k) = ...
                sin(2*pi*(x-cx*t)) * ...
                sin(2*pi*(y-cy*t)) * ...
                sin(2*pi*(z-cz*t)) * ...
                exp(-12 * pi*pi * nu * t);
        end
    end
end
end
