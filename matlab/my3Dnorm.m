function n = my3Dnorm(field)
    l = size(field, 3);
    n = zeros(1, l);
    
    n = norm(field(:), 'inf');
end

