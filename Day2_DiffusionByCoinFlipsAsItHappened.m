% Day 2 - Diffusion by coin flips

% Goals:
% 1) Simulate a random walker/diffuser by coin flips. For example:
% diffusion of a particle over time.
% 2) By looking at many trajectories and doing statistical averages over
% all the outcomes... we can explore implications of our model and also
% learn things about diffusion.

clear variables
close all

% Define simulation parameters
nSteps = 100; %number of coin flips / steps

% Create an array to store all our positions over time
pos = zeros(nSteps,1);

% Define the initial condition
pos(1) = 0;

% Define a time vector
time = 1:nSteps;

% Loop over all time steps

for i = 2:nSteps

    coinFlip = rand(); % generate a number between 0 and 1
    
    if coinFlip > 0.5
        pos(i) = pos(i-1)+1; % move to the right
    else
        pos(i) = pos(i-1)-1; % move to the left
    end

end

% Plot position over our particle over time
figure
plot(time,pos)
xlabel('time step')
ylabel('position')
ylim([-nSteps/2, nSteps/2]); % force matlab to plot a specific range for y

% To compare to analytic solutions and get feel for the average or variance
% of particle distribution over time... we need to simulate lots of
% particles.

nParticles = 1000;

% We need a 2D matrix to keep track of positions for EACH particle over
% time.
posMatrix = zeros(nParticles, nSteps);

% Run our full simulation, looping through all particles and all steps.
for j = 1:nParticles
    for i = 2:nSteps
        coinFlip = rand;
        if coinFlip>0.5
            posMatrix(j,i) = posMatrix(j,i-1) + 1; % hop to right
        else
            posMatrix(j,i) = posMatrix(j,i-1) - 1; % hop to left
        end
    end
end

% Plot all our particle trajectories over time
figure
plot(time,posMatrix')
xlabel('time step')
ylabel('position')
ylim([-nSteps/2, nSteps/2]); % force matlab to plot a specific range for y

% Note that the distribution of positions at each time point has a
% characteristic "width".  How could we quantify this width?
% We could.... calculate the displacement from starting position for each
% particle, and average these.
% How about the average of displacement^2, instead?

% Calculate the "mean squared displacement" = MSD
posSquare = posMatrix.^2; % square each element in our matrix of positions
MSD = mean(posSquare,1); % average down the first index (across all particles)

% Plot the MSD over time
figure
plot(time, MSD, '.', 'MarkerSize', 10)
xlabel('time step')
ylabel('mean squared displacement')

% MSD (units of length squared) is proportion to t. Pretty cool.




























