function qout = integrateConvection(qin, istart, iend, nu, cx, cy, cz, dt, timesteps)

    size = iend-istart + 1;
    dx = 1 / size;

    if (timesteps > 0)
        qout = integrate4(qin, istart, iend, nu, cx, cy, cz, dt, dx);
    end
    timesteps = timesteps - 1;
    while timesteps > 0
        %fprintf('%d timesteps to go\n', timesteps)
        qout = integrate4(qout, istart, iend, nu, cx, cy, cz, dt, dx);
        timesteps = timesteps - 1;
    end
    
end

function qnew = integrate4(q, istart, iend, nu, cx, cy, cz, dt, dx)
    k1 = rhs4(q          , istart, iend, nu, cx, cy, cz, dx);
    k2 = rhs4(q + dt/2*k1, istart, iend, nu, cx, cy, cz, dx);
    k3 = rhs4(q + dt/2*k2, istart, iend, nu, cx, cy, cz, dx);
    k4 = rhs4(q + dt  *k3, istart, iend, nu, cx, cy, cz, dx);
    
    qnew = q + dt/6 * (k1 + 2*k2 + 2*k3 + k4);
end

function qnew = integrate2(q, istart, iend, nu, cx, cy, cz, dt, dx)
    k = rhs4(q, istart, iend, nu, cx, cy, cz, dx);
    qnew = q + dt*k;
end

function v = rhs4(q, istart, iend, nu, cx, cy, cz, dx)
    v = zeros(size(q));
    dx2 = dx*dx;
    r = istart:iend;
    for i = r
        for j = r
            for k = r
                v(i, j, k) = ...
                    advection4(q, i, j, k, cx, cy, cz, dx) + ...
                    nu * laplace4(q, i, j, k, dx2);
            end
        end
    end
end

function v = rhs2(q, istart, iend, nu, cx, cy, cz, dx)
    v = zeros(size(q));
    dx2 = dx*dx;
    r = istart:iend;
    for i = r
        for j = r
            for k = r
                v(i, j, k) = ...
                    advection2(q, i, j, k, cx, cy, cz, dx) + ...
                    nu * laplace2(q, i, j, k, dx2);
            end
        end
    end
end

function v = advection4(q, i, j, k, cx, cy, cz, dx)
    dq_dx = - 1 * q(i-2, j, k) + ...
            + 8 * q(i-1, j, k) + ...
            - 8 * q(i+1, j, k) + ...
            + 1 * q(i+2, j, k);
    dq_dy = - 1 * q(i, j-2, k) + ...
            + 8 * q(i, j-1, k) + ...
            - 8 * q(i, j+1, k) + ...
            + 1 * q(i, j+2, k);
    dq_dz = - 1 * q(i, j, k-2) + ...
            + 8 * q(i, j, k-1) + ...
            - 8 * q(i, j, k+1) + ...
            + 1 * q(i, j, k+2);
        
    v = (cx*dq_dx + cy*dq_dy + cz*dq_dz) / (12 * dx);
end

function v = advection2(q, i, j, k, cx, cy, cz, dx)
    dq_dx = + q(i-1, j, k) + ...
            - q(i+1, j, k); 
    dq_dy = + q(i, j-1, k) + ...
            - q(i, j+1, k); 
    dq_dz = + q(i, j, k-1) + ...
            - q(i, j, k+1); 
        
    v = (cx*dq_dx + cy*dq_dy + cz*dq_dz) / (2 * dx);
end

function v = laplace4(q, i, j, k, dx2)
    d2q_dx2 = -  1 * q(i-2, j, k) ...
              + 16 * q(i-1, j, k) ...
              - 30 * q(i  , j, k) ...
              + 16 * q(i+1, j, k) ...
              -  1 * q(i+2, j, k);
    d2q_dy2 = -  1 * q(i, j-2, k) ...
              + 16 * q(i, j-1, k) ...
              - 30 * q(i, j  , k) ...
              + 16 * q(i, j+1, k) ...
              -  1 * q(i, j+2, k);
    d2q_dz2 = -  1 * q(i, j, k-2) ...
              + 16 * q(i, j, k-1) ...
              - 30 * q(i, j, k  ) ...
              + 16 * q(i, j, k+1) ...
              -  1 * q(i, j, k+2);
    v = (d2q_dx2 + d2q_dy2 + d2q_dz2) / (12 * dx2);
end

function v = laplace2(q, i, j, k, dx2)
    d2q_dx2 = + 1 * q(i-1, j, k) ...
              - 2 * q(i  , j, k) ...
              + 1 * q(i+1, j, k);
    d2q_dy2 = + 1 * q(i, j-1, k) ...
              - 2 * q(i, j  , k) ...
              + 1 * q(i, j+1, k);
    d2q_dz2 = + 1 * q(i, j, k-1) ...
              - 2 * q(i, j, k  ) ...
              + 1 * q(i, j, k+1);
    v = (d2q_dx2 + d2q_dy2 + d2q_dz2) / dx2;
end

