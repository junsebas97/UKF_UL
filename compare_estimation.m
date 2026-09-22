function compare_estimation(x, y, t, mx, Smat, names)
%{
this function compares the given estimations. For each variable at each
DOF, it plots all the given estimation with the the truth and the measurements.
It creates an interactive plot showing the time series.

INPUTS:
x:     true state/loads                         [m], [m/s], [m/s2], [m] or [kN]
y:     measurements                             [m], [m/s], [m/s2]
t:     analysed times                           [s]
mx:    state/load mean of different estimations [m], [m/s], [m/s2], [m] or [kN]
Smat:  selection matrixes of the measurements   [-]
names: names or labels of the estimations       [-]

MADE BY: junsebas97
%}

Nx     = size(x, 1);          % number of state components
N_DOFs = size(Smat{1}, 2);    % number of DOFs
Nx_DOF = Nx/N_DOFs;           % state components per DOF

% create UI figure and set its layout
fig            = uifigure();
gl             = uigridlayout(fig, [2, 1]);
gl.RowHeight   = {'5x', '1x'};
gl.ColumnWidth = {'1x', '1x', '4x'};

% allocate the subplots
plotPanel  = uipanel(gl);
plotLayout = tiledlayout(plotPanel, Nx_DOF, 1);
ax         = cell(Nx_DOF, 1);
for i = 1:Nx_DOF
    ax{i} = nexttile(plotLayout);
end

% create the slider
sld        = uislider(gl, 'Limits', [1, max(N_DOFs, 1.1)], ...
                          'Value',                      1, ...
                          'MajorTicks', round(1:5:N_DOFs));

btnForward = uibutton(gl, 'Text', 'Forward', ...
                          'ButtonPushedFcn', @(btn, event) stepSlider(sld,  1, ax, x, y, t, mx, Smat, names));

btnBack    = uibutton(gl, 'Text','Back', ...
                          'ButtonPushedFcn', @(btn, event) stepSlider(sld, -1, ax, x, y, t, mx, Smat, names));

% organize the components
plotPanel.Layout.Column  = [1, 3];
sld.Layout.Row           = 2;
sld.Layout.Column        = 3;
btnBack.Layout.Row       = 2;
btnBack.Layout.Column    = 1;
btnForward.Layout.Row    = 2;
btnForward.Layout.Column = 2;

% initialize the plot
updatePlots(sld.Value, ax, x, y, t, mx, Smat, names);

% update the plot to the slider value
sld.ValueChangedFcn = @(src, event) updatePlots(src.Value, ax, x, y, t, mx, ...
                                                Smat, names);
end

function updatePlots(i, ax, x, y, t, mx, Smat, names)
i = round(i);

colors = ["blue", "red", "green", "magenta", "cyan"];
N_est  = size(names, 2);      % number of estimations
Nx     = size(x, 1);          % number of state components
N_DOFs = size(Smat{1}, 2);    % number of DOFs
Nx_DOF = Nx/N_DOFs;           % state components per DOF

for k = 1:Nx_DOF
    x_idx   = i + N_DOFs*(k - 1);                  % index current state
    Ny_prev = size(cell2mat(Smat(1:k - 1)), 1);    % number of previous
                                                   % measurements
    cla( ax{k});
    hold(ax{k}, 'on');

    if size(Smat{k}, 1) ~= 0
        if any(Smat{k}(:, i))
            y_idx = Ny_prev + find(Smat{k}(:, i));
            plot(ax{k}, t, y(y_idx, :), 'k.', 'DisplayName', 'Measurements')
        end
    end
    for j = 1:N_est
        plot(ax{k}, t, mx{j}(x_idx, :), '--', 'Color', colors(j), ...
             'DisplayName', names(j))
    end

    plot(  ax{k}, t, x(x_idx, :), 'k-', 'DisplayName', 'True')
    hold(  ax{k},   'off');
    xlabel(ax{k}, 't [s]');
    axis(  ax{k}, 'tight');
    grid(  ax{k},    'on');
    legend(ax{k});
end
end

function stepSlider(sld, step, ax, x, y, t, mx, Smat, names)
%{
this function to moves the slider by step
%}

% compute the new slider value withint the limits
newValue = round(sld.Value + step);
newValue = max(min(newValue, sld.Limits(2)), sld.Limits(1));

% update the slider and the plot
sld.Value = newValue;
updatePlots(newValue, ax, x, y, t, mx, Smat, names);
end