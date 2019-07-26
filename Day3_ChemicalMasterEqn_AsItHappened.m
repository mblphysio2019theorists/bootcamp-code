% Day 3 - Chemical Master Equation chassis!  (Spread the Butter = Melt the Butter??)

% Goal: Develop a code "chassis" for stochastic simulations of processes
% with on/off rates, left/right hops, macromolecule production/decay,
% cytoskeletal growth/depolymerization, etc.

%% Part I. Diffusion in 1D.
% We're going to simulate the diffusion of molecules, which start in the
% middle of an E. coli cell.

clear variables
close all

% Define the model parameters:
D =   10;  % diffusion coefficient, in um^2/s
L = 2; % size of our system (i.e. a E coli cell), in um
nTotal = 100; % total number of molecules

% Simulation parameters:
dx = 0.01; % box size (position step) in um
numBoxes = L/dx; % number of boxes in the cell
k = D/dx^2; % jump rate of a single molecule, in 1/s
% The rate of jumping comes from D = L^2/t

dt = (1/k)/10; % time step in sec
numTimeSteps = 20000; % number of time steps in our simulation

% Initialize the simulation. We create a matrix N where each column is a
% box of size dx, and each row corresponds to a time point.
N = zeros(numTimeSteps, numBoxes);

% Define our initial conditions to "seed" the simulation;
N(1, numBoxes/2) = nTotal; % we put ALL our molecules at the cell's midpoint

% Do the simulation!  Spread the butter!  Melt the butter?!
for i = 2:numTimeSteps
    
    % First we deal with the non-boundary boxes
    for j = 2:numBoxes-1
        % number now = number before + in from left + in from right
        % - out to left - out to right
        N(i,j) = N(i-1,j) + N(i-1, j-1)*k*dt + N(i-1, j+1)*k*dt...
            - N(i-1,j)*k*dt - N(i-1,j)*k*dt;
    end
    
    % Specially account for our cell boundaries
    % First box
    N(i,1) = N(i-1,1) + N(i-1,2)*k*dt - N(i-1,1)*k*dt;
    
    % Last box
    N(i,numBoxes) = N(i-1,numBoxes) + N(i-1,numBoxes-1)*k*dt - N(i-1,numBoxes)*k*dt;
    
end

% We'll use a 3D bar plot to visualize our results. position, time, # of
% molecules

bar3(N(1:50:end, : ) )
xlabel('position')
ylabel('times')
zlabel('number of molecules')
xlim([numBoxes/2 - 10, numBoxes/2 + 10  ]) % zoom in to center of cell
ylim([0,10 ])

% Let's plot in 2D the distribution at a couple of time points
positionVector = (0:dx:L-dx); % vector of positions in um to plot against

figure
plot(positionVector, N(100,:), '-k') % distribution at time 100, which is 1e-4 sec
hold on
plot(positionVector, N(1000,:), '-r') % distribution at time 1000, which is 1e-3 sec
plot(positionVector, N(10000,:), '-b') % distribution at time 100, which is 0.01 sec
hold off
xlabel('position (um)')
ylabel('number of molecules')
legend('0.0001 sec', '0.001 sec', '0.01 sec')

totalCheck = sum(N(400,:))

%% Diffusion in 1D: FRAP!

% Let's copy our code from above, and just modulate the initial condition.

clear variables
close all

% Define the model parameters:
D =   10;  % diffusion coefficient, in um^2/s
L = 2; % size of our system (i.e. a E coli cell), in um
nTotal = 100; % total number of molecules

% Simulation parameters:
dx = 0.01; % box size (position step) in um
numBoxes = L/dx; % number of boxes in the cell
k = D/dx^2; % jump rate of a single molecule, in 1/s
% The rate of jumping comes from D = L^2/t

dt = (1/k)/10; % time step in sec
numTimeSteps = 20000; % number of time steps in our simulation

% Initialize the simulation. We create a matrix N where each column is a
% box of size dx, and each row corresponds to a time point.
N = zeros(numTimeSteps, numBoxes);


%%%%% Define our initial conditions to "seed" the simulation; %%%%
N(1,:) = 1; % some molecules evenly throughout the cell (AU fluorescence of 1)
N(1, numBoxes/2 - 10 : numBoxes/2 + 10) = 0; % we FRAP the center of the cell


% Do the simulation!  Spread the butter!  Melt the butter?!
for i = 2:numTimeSteps
    
    % First we deal with the non-boundary boxes
    for j = 2:numBoxes-1
        % number now = number before + in from left + in from right
        % - out to left - out to right
        N(i,j) = N(i-1,j) + N(i-1, j-1)*k*dt + N(i-1, j+1)*k*dt...
            - N(i-1,j)*k*dt - N(i-1,j)*k*dt;
    end
    
    % Specially account for our cell boundaries
    % First box
    N(i,1) = N(i-1,1) + N(i-1,2)*k*dt - N(i-1,1)*k*dt;
    
    % Last box
    N(i,numBoxes) = N(i-1,numBoxes) + N(i-1,numBoxes-1)*k*dt - N(i-1,numBoxes)*k*dt;
    
end

% We'll use a 3D bar plot to visualize our results. position, time, # of
% molecules

bar3(N(1:50:end, : ) )
xlabel('position')
ylabel('times')
zlabel('number of molecules')


% Let's plot in 2D the distribution at a couple of time points
positionVector = (0:dx:L-dx); % vector of positions in um to plot against

figure
plot(positionVector, N(100,:), '-k') % distribution at time 100, which is 1e-4 sec
hold on
plot(positionVector, N(1000,:), '-r') % distribution at time 1000, which is 1e-3 sec
plot(positionVector, N(10000,:), '-b') % distribution at time 100, which is 0.01 sec
hold off
xlabel('position (um)')
ylabel('number of molecules')
legend('0.0001 sec', '0.001 sec', '0.01 sec')

totalCheck1 = sum(N(1,:))
totalCheck400 = sum(N(400,:))

























