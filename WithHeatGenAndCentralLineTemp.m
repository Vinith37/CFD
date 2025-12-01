clc; clear; close all;

%Defining the geometry of slab
L = 1;
W = 1;

%Defining boundary conditions
T_left = 50;
T_right = 50;
T_top = 100;
T_bottom = 100;

%Tolerance for grid independent test
tol = 1e-6;

%Maximum iteration allowed
Max_Iteration = 100000;

%Defining grid size
grid_size = [81];

% Storage for center temperature
center_T = 75*ones(size(grid_size));

% --- New storage for central line profiles ---
centerline_X_profiles = {}; % temperature along mid-y (vary x) for each grid
centerline_Y_profiles = {}; % temperature along mid-x (vary y) for each grid
grid_x_coords = {};         % x coordinates for each grid
grid_y_coords = {};         % y coordinates for each grid

for g = 1:length(grid_size)
    
    % Grid definition
    n_x = grid_size(g);
    n_y = grid_size(g);
    
    dx = L/(n_x-1);
    dy = W/(n_y-1);
    
    % Create coordinate vectors (useful for plotting)
    x = linspace(0, L, n_x);
    y = linspace(0, W, n_y);
    
    % Initialize temperature matrix
    T = 75*ones(n_y, n_x);
    S = 0;
    T(:,1) = T_left; %assigning temperature to left wall 
    T(:,end) = T_right; %assigning temperature to right wall 
    T(end,:) = T_top; %assigning temperature to top wall 
    T(1,:) = T_bottom; %assigning temperature to bottom wall
    
    %This loop repeat until convergence.
    for iteration = 1:Max_Iteration 
        T_old = T; %storing old value
        
        %sweeping temperature values in interior nodes
        for j = 2:n_y-1 %j=1&ny is already assigned with boundary values
            for i = 2:n_x-1 %i=1&nx is already assigned with boundary
                T(j,i) = 0.25 * ( T(j+1,i) + T(j-1,i) + T(j,i+1) + T(j,i-1) + (S*(dx^2)) );
            end
        end
        
        % Error check (solution convergence)
        error = max(max(abs(T - T_old)));
        if error < tol
            fprintf('Grid %dx%d converged in %d iterations (error=%.5f)\n', n_x, n_y, iteration, error);
            break;
        end
    end
    
        % Assume T is your temperature grid
    % and x, y are your coordinates

    % Find middle row and column
    mid_row = round(size(T,1)/2);   % middle in y
    mid_col = round(size(T,2)/2);   % middle in x

    % Get centerlines
    centerline_X = T(mid_row, :);   % across x
    centerline_Y = T(:, mid_col);   % across y

    % Plot centerline along x (middle of y)
    figure;
    plot(x, centerline_X);
    xlabel('x');
    ylabel('Temperature');
    title('Centerline along x');

    % Plot centerline along y (middle of x)
    figure;
    plot(y, centerline_Y);
    xlabel('y');
    ylabel('Temperature');
    title('Centerline along y');

% Make contour plot of temperature
figure;

[Xg, Yg] = meshgrid(x, y);   % create grid
contourf(Xg, Yg, T);         % filled contour plot

colorbar;                    % show color scale
xlabel('x');
ylabel('y');
title('Temperature Distribution');

end



