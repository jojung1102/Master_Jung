% Set scriptPath current directory
directory = 'C:\Users\lenovo\Documents\Masterarbeit\Programming';
cd(directory);
% Add MatCont files
matContPath = fullfile('C:\Users\lenovo\Documents\MATLAB\MatCont7p4');
addpath(genpath(matContPath));

% Create folder to save outputs
outputFolder = [directory '\output'];
savePath = fullfile(outputFolder);

% System of ODEs
MySystem = Dynamic_system;

%% Initial conditions
% Spherical condition and Arc length derivative weighed by max(A_c)/max(A_I)
h = @(x) x - x.^2;
dh = @(x) sqrt(1 + (2 - 4*x).^2);

% Define the x range
x_values = linspace(0, 1, 10000);  % Fine grid for integration

% Compute arc length at each x' int(dh, 0, x')
arc_length = zeros(size(x_values));
for i = 2:length(x_values)
    arc_length(i) = integral(dh, 0, x_values(i));
end

% Get equally spaced arc lengths
spacings = 21;
s_values = linspace(0, arc_length(end), spacings);  % 21 equal spacings (reduced later for unintersting conditions)

% Interpolate to find x values at these arc lengths
A_I_vals = interp1(arc_length, x_values, s_values, 'linear');
A_cont_vals = h(A_I_vals);

% Parameter values t_n, t_a, zeta, kappa are rescaled by a/A
n = 0.1;
C0 = 0.4;
c = 0;
epsilon = -1;
t_a = 2;
kappa = 2;

% Event function for stopping the ODE
OPTIONS = odeset('Events', @eventFunction, 'RelTol', 1e-09, 'AbsTol', 1e-09);

%% t_n and zeta
combis = [1, 1; 1, 2.5; 1, 5.2; 2, 5.2; 1, 7];

for idx = 1:size(combis, 1)
    t_n = combis(idx, 1);
    zeta = combis(idx, 2);
    j = 1;
    
    % Loop through the grid of initial points
    for i = 1:length(A_I_vals)
        A_I_Init = A_I_vals(i);
        A_cont_Init = A_cont_vals(i);
    
        % Check for allowed geometry
        if min(A_I_Init, 1 - A_I_Init) > 0.2
            % Integrate in time
            [t, X] = ode45(@(t, y) MySystem{2}(t, y, C0, c, zeta, kappa, n/2, n/2, epsilon, t_n, t_a), ...
                [0 2.5], [A_I_Init, A_cont_Init, (n/2)*A_I_Init, (n/2)*A_I_Init], OPTIONS);
    
            % Filter out rows with imaginary values
            valid_rows = all(imag(X) == 0, 2);
            t_real = t(valid_rows);
            X_real = X(valid_rows, :);
            new_column = X_real(:,1) .* X_real(:,3) .* (X_real(:,4) - n/2) / (n/2)^2;
            X_real = [X_real, X_real(:,1) .* X_real(:,3) .* (n/2-X_real(:,4)) / (n/2)^2]; % calculcate the score
    
            % Save trajectory for this initial condition
            filename = fullfile(savePath, 'trajectories_t_zeta', sprintf('trajectory_%.f_%.f_%d.csv', t_n, zeta, j));
            csvwrite(filename, [X_real, t_real]);
            j = j+1;
        end
    end
end

%% t_a and zeta
combis = [0.7, 5; 4, 5];

t_n = 0.5;

for idx = 1:size(combis, 1)
    t_a = combis(idx, 1);
    zeta = combis(idx, 2);
    j = 1;
    
    % Loop through the grid of initial points
    for i = 1:length(A_I_vals)
        A_I_Init = A_I_vals(i);
        A_cont_Init = A_cont_vals(i);
    
        % Check for allowed geometry
        if min(A_I_Init, 1 - A_I_Init) > 0.2
            % Integrate in time
            [t, X] = ode45(@(t, y) MySystem{2}(t, y, C0, c, zeta, kappa, n/2, n/2, epsilon, t_n, t_a), ...
                [0 2.5], [A_I_Init, A_cont_Init, (n/2)*A_I_Init, (n/2)*A_I_Init], OPTIONS);
    
            % Filter out rows with imaginary values
            valid_rows = all(imag(X) == 0, 2);
            t_real = t(valid_rows);
            X_real = X(valid_rows, :);
            new_column = X_real(:,1) .* X_real(:,3) .* (X_real(:,4) - n/2) / (n/2)^2;
            X_real = [X_real, X_real(:,1) .* X_real(:,3) .* (n/2-X_real(:,4)) / (n/2)^2]; % calculcate the score
    
            % Save trajectory for this initial condition
            filename = fullfile(savePath, 'trajectories_t_zeta', sprintf('trajectory_ta_%.1f_%.f_%d.csv', t_a, zeta, j));
            csvwrite(filename, [X_real, t_real]);
            j = j+1;
        end
    end
end

%% Different curvatures
combis = [0.5, 5; 1.5, 6; 4.5, 2];

t_n = 0.5;
t_a = 2;

for idx = 1:size(combis, 1)
    c = combis(idx, 1);
    zeta = combis(idx, 2);
    j = 1;
    
    % Loop through the grid of initial points
    for i = 1:length(A_I_vals)
        A_I_Init = A_I_vals(i);
        A_cont_Init = A_cont_vals(i);
    
        % Check for allowed geometry
        if min(A_I_Init, 1 - A_I_Init) > 0.2
            % Integrate in time
            [t, X] = ode45(@(t, y) MySystem{2}(t, y, C0, c, zeta, kappa, n/2, n/2, epsilon, t_n, t_a), ...
                [0 2.5], [A_I_Init, A_cont_Init, (n/2)*A_I_Init, (n/2)*A_I_Init], OPTIONS);
    
            % Filter out rows with imaginary values
            valid_rows = all(imag(X) == 0, 2);
            t_real = t(valid_rows);
            X_real = X(valid_rows, :);
            new_column = X_real(:,1) .* X_real(:,3) .* (X_real(:,4) - n/2) / (n/2)^2;
            X_real = [X_real, X_real(:,1) .* X_real(:,3) .* (n/2-X_real(:,4)) / (n/2)^2]; % calculcate the score
    
            % Save trajectory for this initial condition
            filename = fullfile(savePath, 'trajectories_c_zeta', sprintf('trajectory_c_%.1f_%.f_%d.csv', c, zeta, j));
            csvwrite(filename, [X_real, t_real]);
            j = j+1;
        end
    end
end

%% Event function for stopping ODE
function [value, isterminal, direction] = eventFunction(t, y)
    A_I = y(1);
    A_cont = y(2);
    n_fI = y(3);
    n_dI = y(4);

    % Check for topological boundary reached
    [vesicleFormation, invalid] = final_state(A_I, A_cont, n_fI, n_dI);

    % Stop integration if any condition met (we stop when vesicle formation or invalid)
    value = ~(vesicleFormation || invalid);  % 0 if any condition met
    isterminal = 1;
    direction = 0;
end

function [vesicleFormation, invalid] = final_state(A_I, A_cont, n_fI, n_dI)

    % Initialize output
    vesicleFormation = 0;
    invalid = 0;

    % Check for vesicle formation condition
    if A_cont < 1e-04
        vesicleFormation = 1;

    % Check invalid conditions
    elseif any(imag([A_I, A_cont, n_fI, n_dI]) ~= 0)
        invalid = 1;
    end
end
