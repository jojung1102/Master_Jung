
% set scriptPath current directory
directory = 'C:\Users\lenovo\Documents\Masterarbeit\Programming';
cd(directory);
% add MatCont files
matContPath = fullfile('C:\Users\lenovo\Documents\MATLAB\MatCont7p4');
addpath(genpath(matContPath));

% Create folder to save outputs
outputFolder = [directory '\output'];
savePath = fullfile(outputFolder);

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
spacings = 31;
s_values = linspace(0, arc_length(end), spacings);  % 30 equal spacings

% Interpolate to find x values at these arc lengths
A_I_vals = interp1(arc_length, x_values, s_values, 'linear');
A_cont_vals = h(A_I_vals);

%% Calculate number of trajectories and score of initial states
% Time scales, rescaled by t_cont, since topological transition is in 
% direction of t_cont and changing t_cont would alter running time till
% transition significantly)
running_time = 3;

% Calculate number of trajectories considered and score expected from (relatively large zeta)
traj_number = 0;
abs_score = 0;
% Loop over the grid of initial points
for i = 1:length(A_I_vals)
    A_I_Init = A_I_vals(i);
    A_cont_Init = A_cont_vals(i);

    % Exclude small vesicles
    if min(A_I_Init, 1-A_I_Init) > 0.2
        traj_number = traj_number+1;
        abs_score = abs_score + A_I_Init^2*(1-A_I_Init); % A_I * n_fI/n_f * n_dII/n_d at uniform distribution
    end
end
disp(abs_score);
disp(traj_number/spacings);
disp(traj_number);

%% Calculate score array

% Parameter ranges, zeta, kappa are rescaled by a/A
zeta_list = linspace(0, 8, 21); % zeta' = [0,80], a/A = 1/10

t_n_list = linspace(1, 10, 21);
t_a_list = linspace(0.5, 5, 21);
c_list = linspace(0, 5, 21);
n_list = linspace(0.001, 0.2, 21);
C0_list = linspace(0.4, 1.2, 21);
epsilon_list = linspace(-1, -3, 21);
kappa_list = linspace(1, 4, 21); %kappa' = [10-40], a/A =1/10

param_lists = {n_list, t_a_list, t_n_list, c_list, C0_list, epsilon_list, kappa_list};
param_names = {'n', 't_a', 't_n', 'c', 'C0', 'epsilon', 'kappa'};


for p = 1:length(param_lists)

    % Default parameters
    kappa = 2;
    n = 0.1;
    epsilon = -1;
    t_n = 5;
    t_a = 2;
    C0 = 0.4;
    c = 0;

    param_values = param_lists{p};
    param_name = param_names{p};

    result = [];
    
    for i_param = 1:length(param_values)
        param_value = param_values(i_param);
        eval([param_name ' = param_value;']);

        for zeta = zeta_list
            [score_array, error, slow_mode_count] = calculate_score(A_I_vals, A_cont_vals, n/2, n/2, epsilon, kappa, zeta, C0, c, t_n, t_a, running_time);
            
            if size(score_array, 1) > 0
                score_sum = sum(score_array(:,1));
            else
                score_sum = 0;
            end
            
            % Store results including param_value
            result = [result; [param_value, zeta, score_sum/abs_score, error/length(A_I_vals), slow_mode_count/length(A_I_vals)]];
        end
    end
    % Define filename with param_name and param_value
    filename = sprintf('Score_%s.csv', param_name);
    csvwrite(fullfile(savePath, 'Score_main', filename), result);
end

%% Interplay of nu and c
% Parameter ranges
c_list = linspace(0, 8, 31);
n_list = linspace(0.05, 0.2, 31);

kappa = 2; %a/A = 1/10
epsilon = -1;
t_a = 2;
t_n = 0.5;
C0 = 0.4;
zeta = 4; %a/A = 1/10
result = [];

for c = c_list
    for n = n_list
        [score_array, error, slow_mode_count] = calculate_score(A_I_vals, A_cont_vals, n/2, n/2, epsilon, kappa, zeta, C0, c, t_n, t_a, running_time);

        if size(score_array, 1) > 0
            score_sum = sum(score_array(:,1));
        else
            score_sum = 0;
        end
        
        % Store results including param_value
        result = [result; [n, c, score_sum/abs_score, error/length(A_I_vals), slow_mode_count/length(A_I_vals)]];
    end
end
% Define filename with param_name and param_value
filename = sprintf('Score_c_n.csv');
csvwrite(fullfile(savePath, 'Score_main', filename), result);

%% Interplay of nu and c for fast equilibration?
% Parameter ranges
c_list = linspace(0, 8, 31);
n_list = linspace(0.05, 0.2, 31);

kappa = 2; %a/A = 1/10
epsilon = -1;
t_a = 1;
t_n = 1;
C0 = 0.4;
zeta = 6; %a/A = 1/10
result = [];

for c = c_list
    for n = n_list
        [score_array, error, slow_mode_count] = calculate_score(A_I_vals, A_cont_vals, n/2, n/2, epsilon, kappa, zeta, C0, c, t_n, t_a, running_time);

        if size(score_array, 1) > 0
            score_sum = sum(score_array(:,1));
        else
            score_sum = 0;
        end
        
        % Store results including param_value
        result = [result; [n, c, score_sum/abs_score, error/length(A_I_vals), slow_mode_count/length(A_I_vals)]];
    end
end
% Define filename with param_name and param_value
filename = sprintf('Score_c_n_large_zeta.csv');
csvwrite(fullfile(savePath, 'Score_main', filename), result);



%% Calculation of the score
function [score_array, error, slow_mode_count] = calculate_score(A_I_vals, A_cont_vals, n_f, n_d, epsilon, kappa, zeta, C0, c, t_n, t_A, running_time)
    
    % System for ODEs
    MySystem = Dynamic_system;

    % Initalize score  and slow mode and error count
    score_array = [];
    error = 0;
    slow_mode_count = 0;

    % Loop over the grid of initial points
    for i = 1:length(A_I_vals)
        A_I_Init = A_I_vals(i);
        A_cont_Init = A_cont_vals(i);

        % Check for allowed geometry (not too far in the edges)
        if min(A_I_Init, 1-A_I_Init) > 0.2

            % Set resolution of  ODE and initialize timer for slow modes
            t_start = tic;
            resolution = 1e-09;

            OPTIONS = odeset('Events', @(t, y) eventFunction(t, y, t_start, n_f, n_d), 'RelTol', resolution, 'AbsTol', resolution);

            % Integrate in time, stop if integration takes too long (slow
            % mode)
            [t, X] = ode45(@(t, y) MySystem{2}(t, y, C0, c, zeta, kappa, n_f, n_d, epsilon, t_n, t_A), ...
                [0 running_time], [A_I_Init, A_cont_Init, n_f*A_I_Init, n_d*A_I_Init], OPTIONS);

            final_point = X(end, :);

            [vesicleFormation, invalid, slow_mode] = final_state(t_start, final_point(1), final_point(2), final_point(3), final_point(4), n_f, n_d);

            if vesicleFormation
                score_val = final_point(1) * final_point(3)/n_f * (n_d-final_point(4))/n_d;
                score_array = [score_array; score_val];
            elseif invalid
                if real(min(final_point(1), 1-final_point(1))-final_point(2)) < 0.0001
                    slow_mode_count = slow_mode_count+1;
                else
                    disp(X([end-5:end], :));
                    disp([n_f, c, zeta])
                    error=error+1;
                end
            elseif slow_mode
                slow_mode_count = slow_mode_count+1;
            end
        end
    end
end

function [value, isterminal, direction] = eventFunction(t, y, t_start, n_f, n_d)
    A_I = y(1);
    A_cont = y(2);
    n_fI = y(3);
    n_dI = y(4);

    % Check for topological boundary reached
    [vesicleFormation, invalid, slow_mode] = final_state(t_start, A_I, A_cont, n_fI, n_dI, n_f, n_d);

    % Stop integration if any condition met
    value = ~(vesicleFormation || invalid || slow_mode);  % 0 if any condition met
    isterminal = 1;
    direction = 0;
end

function [vesicleFormation, invalid, slow_mode] = final_state(t_start, A_I, A_cont, n_fI, n_dI, n_f, n_d)

    % Initialize output
    vesicleFormation = 0;
    invalid = 0;
    slow_mode = 0;

    % Check for vesicle formation condition
    if A_cont < 0.0001
        vesicleFormation = 1;

    % Check invalid conditions (imaginary values, invalid geometries,
    % maximum densities)
    elseif any(imag([A_I, A_cont, n_fI, n_dI]) ~= 0) ...
        || min(A_I, 1-A_I)-A_cont < 0  ...
        || min([n_fI, n_dI, 1-n_fI/A_I, 1-n_dI/A_I, ...
        1-(n_f - n_fI)/(1-A_I), 1-(n_d - n_dI)/(1-A_I)]) < 0

        invalid = 1;

    % Integration time takes too long
    elseif toc(t_start) > 1 %in seconds
        slow_mode = 1;
    end
end