function PlotQ(varargin)
    
    % Load file
    if (nargin >= 1)
        fname = varargin{1};
    else
        fname = 'heat.mat';
    end
    load(fname, 'q');
    
    % Choose level
    if (nargin >= 2)
        level = varargin{2};
    else
        level = size(q{1}, 3) / 2;
    end
    
    % Plotting
    fprintf('Plotting level %d\n', level);
    fprintf('%d frames found\n', length(q));
    imagesc(q{1}(:,:,level));
    set(gca, 'NextPlot', 'replaceChildren');
    caxis([0 1]);
    colorbar
    for t = 1:length(q)
        qq = q{t};
        imagesc(qq(:,:,level));
        title(sprintf('Frame %d', t));
        caxis([0 1]);
        pause(0.1);
    end

end

