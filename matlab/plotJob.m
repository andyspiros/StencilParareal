function plotJob(fname)

% Get data
fid = fopen(fname);
data = fscanf(fid, ' parareal%cPU %d %d %e %e %e %e %e %e %e');
fclose(fid);
rows = length(data) / 10;
data = reshape(data, 10, rows)';

% Split CPU / GPU
targets = unique(data(:, 1));

for target = targets(:)'
    d = data(data(:, 1) == target, 2:end);
    
    % Extract columns
    nodes       = unique(d(:, 1));
    kmax        = unique(d(:, 2));
    
    % Choose colors based on kmax
    kcolors = hsv(length(kmax));
    
    % First plot: runtime
    figure
    hold on
    for i = 1:length(kmax)
        k = kmax(i);
        plot(nodes, d(d(:,2) == k, 7), ...
            'Color', kcolors(i, :), ...
            'Marker', '^', ...
            'DisplayName', sprintf('k = %d', k) ...
          );
    end
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
    xlabel 'Nodes'
    ylabel 'Runtime (sec)'
    grid on
    set(gca, 'XMinorGrid', 'off');
    set(gca, 'YMinorGrid', 'off');
    set(gca, 'XTick', nodes);
    set(gca, 'XLim', [3/4*nodes(1) nodes(end)*4/3]);
    title(sprintf('Parareal runtime on %cPU', target));
    legend show
    set(gcf, 'Units', 'centimeters');
    pos = get(gcf, 'Position');
    pos(3:4) = [12 10];
    set(gcf, 'Position', pos);
    set(gcf,'PaperPositionMode','auto');
    print(gcf, '-depsc2', sprintf('1%c.eps', target));
    
    % Second plot: error
    figure
    hold on
    for i = 1:length(kmax)
        k = kmax(i);
        plot(nodes, d(d(:,2) == k, 5), ...
            'Color', kcolors(i, :), ...
            'Marker', '^', ...
            'DisplayName', sprintf('k = %d', k) ...
          );
    end
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
    xlabel 'Nodes'
    ylabel 'Relative error'
    grid on
    set(gca, 'XMinorGrid', 'off');
    set(gca, 'YMinorGrid', 'off');
    set(gca, 'XTick', nodes);
    set(gca, 'XLim', [3/4*nodes(1) nodes(end)*4/3]);
    title(sprintf('Parareal error on %cPU', target));
    legend show
    set(gcf, 'Units', 'centimeters');
    pos = get(gcf, 'Position');
    pos(3:4) = [12 10];
    set(gcf, 'Position', pos);
    set(gcf,'PaperPositionMode','auto');
    print(gcf, '-depsc2', sprintf('2%c.eps', target));
    
    % Third plot: energy
    figure
    hold on
    for i = 1:length(kmax)
        k = kmax(i);
        plot(nodes, d(d(:,2) == k, 9) ./ d(d(:,2) == k, 8), ...
            'Color', kcolors(i, :), ...
            'Marker', '^', ...
            'DisplayName', sprintf('k = %d', k) ...
          );
    end
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
    xlabel 'Nodes'
    ylabel 'E(P)/E(S)'
    grid on
    set(gca, 'XMinorGrid', 'off');
    set(gca, 'YMinorGrid', 'off');
    set(gca, 'XTick', nodes);
    set(gca, 'XLim', [3/4*nodes(1) nodes(end)*4/3]);
    title(sprintf('Energy consumption ratio on %cPU', target));
    legend show
    set(gcf, 'Units', 'centimeters');
    pos = get(gcf, 'Position');
    pos(3:4) = [12 10];
    set(gcf, 'Position', pos);
    set(gcf,'PaperPositionMode','auto');
    print(gcf, '-depsc2', sprintf('3%c.eps', target));
end

end

