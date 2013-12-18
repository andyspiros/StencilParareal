function plotFile(fname, varargin)

if (nargin > 1)
    removeZeros = varargin{1};
else
    removeZeros = 1;
end

% Get domain sizes
cmd = ['grep domain ' fname ' | cut -d "x" -f 3'];
[~, out] = system(cmd);
domain = str2num(out);

% Get number of timesteps
cmd = ['grep timesteps ' fname ' | cut -d " " -f 4'];
[~, out] = system(cmd);
timesteps = str2num(out);

% Get CFL
cmd = ['grep CFL ' fname ' | cut -d " " -f 4'];
[~, out] = system(cmd);
cfl = str2num(out);

% Get error
cmd = ['grep "Error at end" ' fname ' | cut -d " " -f 4'];
[~, out] = system(cmd);
error = str2num(out);

% Remove data with too big error or zero error
data = [domain timesteps cfl error];
idx = find(data(:,4) < 1);
if (removeZeros > 0)
    idx = intersect(idx, find(data(:,4) > 0));
end
data = data(idx, :);

% Plot CFL -> error (space / time)
figure
loglog(data(:,3), max(eps, data(:,4)), '*')
xlabel 'CFL value'
ylabel 'Infinity norm of relative error'

% Plot resolution -> error
figure
hold on
sizes = unique(data(:,1));
h = hsv(length(sizes));
for i = 1:length(sizes)
    s = sizes(i);
    idx = find(data(:,1) == s);
    plot(data(idx, 2), max(eps, data(idx, 4)), ...
         'Marker', '+', ...
         'DisplayName', sprintf('%d points', s), ...
         'Color', h(i, :) ...
    );
end
grid on
xlabel 'Number of timesteps'
ylabel 'Infinity norm of relative error'
legend('Location', 'Best')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

% Plot resolution -> error (time / space)
figure
hold on
steps = unique(data(:,2));
h = hsv(length(steps));
for i = 1:length(steps)
    s = steps(i);
    idx = find(data(:,2) == s);
    plot(data(idx, 1), max(eps, data(idx, 4)), ...
         'Marker', '+', ...
         'DisplayName', sprintf('%d timesteps', s), ...
         'Color', h(i, :) ...
    );
end
grid on
xlabel 'Grid size'
ylabel 'Infinity norm of relative error'
legend('Location', 'Best')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

