function daintEnergy(data)

sep = find(diff(data(:,1)) < 0);

for i = [0 1]
    clear data3d
    
    d = data((sep*i+1):(sep*(i+1)), :);
    data3d(1, :, :) = repmat(mean(d(:,6:8)), 3, 1);
    
    groupLabels = {'Serial'};
    
    g = 2;
    for n = unique(d(:,1))'
        data3d(g, :, :) = d(d(:,1)==n, 3:5);
        groupLabels{g} = sprintf('%d nodes', n);
        g = g + 1;
    end
    
    plotBarStackGroups(data3d, groupLabels);
    if (i == 0)
        t = 'CPU';
    else
        t = 'GPU';
    end
    title(sprintf('%s power consumption', t));
    ylabel 'Power consumption (J)'
    legend({'Node', 'Communication', 'Blower'})
    xlabel 'kmax'
    grid on
    
end

end

