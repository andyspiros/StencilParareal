function [q, exact, error] = heatEquation(size, nu, cx, cy, cz)
dx = 1 / (size+1);
endtime = 0.05;
cfl = 0.1;
timesteps = round(endtime * (size+1)*(size+1) / cfl);
dt = endtime / timesteps;

cfl = dt / dx / dx;
%fprintf('Actual CFL: %f\n', cfl);

q = zeros(size+4, size+4, size+4);
istart = 3;
iend = size+2;
q = exactQ(q, istart, iend, 0);

for t = 1:timesteps
    %fprintf('Timestep %d/%d\n', t, timesteps);
    q = rk(q, istart, iend, dt, dx);
end

exact = exactQ(zeros(size+4, size+4, size+4), istart, iend, endtime);
errorfield = q-exact;
error = norm(errorfield(:), 'inf') / norm(exact(:), 'inf');

end

function q = exactQ(q, istart, iend, t, nu, cx, cy, cz)
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

function tens = laplace(q, istart, iend, dt, dx)
tens = zeros(size(q));
dx2 = dx*dx;
for i = istart:iend
    for j = istart:iend
        for k = istart:iend
            tens(i, j, k) = -24 * q(i, j, k) + ...
              2 * ( ...
                q(i-1, j  , k  ) + ...
                q(i+1, j  , k  ) + ...
                q(i  , j-1, k  ) + ...
                q(i  , j+1, k  ) + ...
                q(i  , j  , k-1) + ...
                q(i  , j  , k+1)   ...
              ) + ...
                q(i+0, j-1, k-1) + ...
                q(i+0, j-1, k+1) + ...
                q(i+0, j+1, k-1) + ...
                q(i+0, j+1, k+1) + ...
                q(i-1, j+0, k-1) + ...
                q(i-1, j+0, k+1) + ...
                q(i+1, j+0, k-1) + ...
                q(i+1, j+0, k+1) + ...
                q(i-1, j-1, k+0) + ...
                q(i-1, j+1, k+0) + ...
                q(i+1, j-1, k+0) + ...
                q(i+1, j+1, k+0)   ...
            ;
            tens(i, j, k) = tens(i,j,k) / 6 / dx2;
        end
    end
end
end

function q = rk(q, istart, iend, dt, dx)
k1 = laplace(q          , istart, iend, dt, dx);
k2 = laplace(q + dt/2*k1, istart, iend, dt, dx);
k3 = laplace(q + dt/2*k2, istart, iend, dt, dx);
k4 = laplace(q + dt  *k3, istart, iend, dt, dx);

q = q + dt/6 * (k1 + 2*k2 + 2*k3 + k4);
end

