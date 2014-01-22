function n = my3Dnorm(field)
    l = size(field, 3);
    n = zeros(1, l);
    
    for i = 1:l
        n(i) = norm(field(:,:,i));
    end
    
    n = norm(n, 'inf');

end

