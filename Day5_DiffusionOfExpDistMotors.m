
%% Diffusion of exponentially-distributed motors by spread the butter

clear variables
close all

% Define the model parameters:
D =   10;  % diffusion coefficient, in um^2/s
L = 200; % size of our system - we'll make it 20um (small camera chip) % CHANGE
%nTotal = 100; % total number of molecules

% Simulation parameters:
dx = 1; % box size (position step) in um %%% CHANGE FROM 0.01
numBoxes = L/dx; % number of boxes in the cell
k = D/dx^2; % jump rate of a single molecule, in 1/s
% The rate of jumping comes from D = L^2/t

dt = (1/k)/10; % time step in sec
numTimeSteps = 20000; % number of time steps in our simulation

% Initialize the simulation. We create a matrix N where each column is a
% box of size dx, and each row corresponds to a time point.
N = zeros(numTimeSteps, numBoxes);



%%%%%%%%% Define the initial conditions to "seed" the simulation %%%%%%%%%%
positionVector = (0:dx:L-dx); % first, we need a vector of positions to plot against

positionVector = positionVector - round(L/2); % center at zero

A = 10; % scale height of exponential peak
lambda = 10; % scale steepness of exponential. bigger lambda = less steep
B = 1; % scale height of reservoir

% Set up motors to have exponential decay
N(1,:) = A * exp(-1 * abs(positionVector)/lambda) + B;


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
% 
% % We'll use a 3D bar plot to visualize our results. position, time, # of
% % molecules
% 
% bar3(N(1:50:end, : ) )
% xlabel('position')
% ylabel('times')
% zlabel('number of molecules')
% xlim([numBoxes/2 - 10, numBoxes/2 + 10  ]) % zoom in to center of cell
% ylim([0,10 ])

% REMEMBER TO DELETE OLD POSITION VECTOR

figure(1)
plot(positionVector, N(1,:), '-k') % distribution at time 1, which is 0 sec
hold on
plot(positionVector, N(1000,:), '-r') % distribution at time 1000, which is 1e-3 sec
plot(positionVector, N(20000,:), '-b') % distribution at time 100, which is 0.02 sec
hold off
xlabel('position (um)')
ylabel('number of molecules')
legend('0 sec', '0.001 sec', '0.02 sec')
title('Diffusion of exponentially-distributed motors - continuum equation')

totalCheck = sum(N(400,:));

%% Diffusion of exponentially-distributed motors by continuum diffusion equation

% Define the model parameters:
D = 10;  % diffusion coefficient, in um^2/s
L = 200; % size of our system - we'll make it 20um (small camera chip) % CHANGE
nTotal = 100; % total number of molecules

% Simulation parameters:
dx = 1; % box size (position step) in um %%% CHANGE FROM 0.01
numBoxes = L/dx; % number of boxes in the cell
k = D/dx^2; % jump rate of a single molecule, in 1/s
% The rate of jumping comes from D = L^2/t

dt = (1/k)/10; % time step in sec
numTimeSteps = 20000; % number of time steps in our simulation

% Initialize the simulation. We create a matrix N where each column is a
% box of size dx, and each row corresponds to a time point.
N = zeros(numTimeSteps, numBoxes);



%%%%%%%%% Define the initial conditions to "seed" the simulation %%%%%%%%%%
positionVector = (0:dx:L-dx); % first, we need a vector of positions to plot against

positionVector = positionVector - round(L/2); % center at zero

A = 10; % scale height of exponential peak
lambda = 10; % scale steepness of exponential. bigger lambda = less steep
B = 1; % scale height of reservoir

% Set up motors to have exponential decay
N(1,:) = A * exp(-1 * abs(positionVector)/lambda) + B;


% Do the simulation!  Spread the butter!  Melt the butter?!
for i = 2:numTimeSteps
    
    % First we deal with the non-boundary boxes
    for j = 2:numBoxes-1
  
        N(i,j) = N(i-1,j) + D*( (N(i-1,j+1)-N(i-1,j)) - (N(i-1,j)-N(i-1,j-1)) / dx^2) * dt;
        
    end
    
    % Specially account for our cell boundaries
%     % First box
%       N(i,1) = N(i-1,1) + D*( (N(i-1,3)-N(i-1,2)) - (N(i-1,2)-N(i-1,1)) / dx^2) * dt;
    N(i,1) = N(i-1,1) + N(i-1,2)*k*dt - N(i-1,1)*k*dt;
%     
%     % Last box
%      N(i,1) = N(i-1,numBoxes) + D*( (N(i-1,numBoxes)-N(i-1,numBoxes-1)) - (N(i-1,numBoxes-1)-N(i-1,numBoxes-2)) / dx^2) * dt;
    N(i,numBoxes) = N(i-1,numBoxes) + N(i-1,numBoxes-1)*k*dt - N(i-1,numBoxes)*k*dt;
    
end
% 
% % We'll use a 3D bar plot to visualize our results. position, time, # of
% % molecules
% 
% bar3(N(1:50:end, : ) )
% xlabel('position')
% ylabel('times')
% zlabel('number of molecules')
% xlim([numBoxes/2 - 10, numBoxes/2 + 10  ]) % zoom in to center of cell
% ylim([0,10 ])

% REMEMBER TO DELETE OLD POSITION VECTOR

figure(1)
pause
hold on
plot(positionVector, N(1,:), '.k') % distribution at time 1, which is 0 sec
plot(positionVector, N(1000,:), '.r') % distribution at time 1000, which is 1e-3 sec
plot(positionVector, N(20000,:), '.b') % distribution at time 100, which is 0.02 sec
hold off
xlabel('position (um)')
ylabel('number of molecules')
legend('0 sec', '0.001 sec', '0.02 sec')
title('Diffusion of exponentially-distributed motors - continuum equation vs spread the butter')

totalCheck1 = sum(N(1,:));
totalCheck4000 = sum(N(4000,:));