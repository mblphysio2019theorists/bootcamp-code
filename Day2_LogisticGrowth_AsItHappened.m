% Day 2 - Logistic Growth 
% In this exercise, we will solve numerically a differential equation for
% logistic growth. In part II, we'll write a simple stochastic simulation
% of logistic growth.

%% Part I. Solve (a.k.a. do a deterministic simulation of) logistic growth

% Background: Our bacterial cultures don't actually grow unbounded to the
% size of the earth.

% Remember from our board work:
% our logistic equation for growth
% dN/dt = r*N - r*N.^2/Nstar

clear variables
close all

% Define our biological parameters
r = 1/37; % growth rate in 1/min
N0 = 1; % initial condition: number of cells at time point 0
Nstar = 1e5; % carrying capacity = saturation, steady-state population size

% Define our simulation parameters
dt = 0.3; % time step in minutes. 
totalSteps = 3000; % total number of steps we'll take in our simulation

% Now, we'll simulate our logistic growth!

% Define a matlab array that will hold N(t)
N = zeros(totalSteps,1);

N(1) = N0; % set our first N (time 0) to be the starting population we decided on

% Use a "for loop" to calculate N at all time points.

for i = 2:totalSteps
    N(i) = N(i-1) + r*N(i-1)*dt - r*dt*N(i-1)^2/Nstar;
end

% Last but not least: we'll plot N as a function of time to see the results
% of our simulation.

timesSim = (0:dt:(totalSteps-1)*dt);

% Plot population size vs time
figure
plot(timesSim, N)
xlabel('Time (mins)')
ylabel('Number of cells')
title('Exponential population growth over time')


%% Part II. Stochastic "Gillespie" simulation of logistic growth

% Keep our biological and simulation parameters from above.

% Define a matlab array that will hold N(t)
stochN = zeros(totalSteps, 1);

stochN(1) = N0; % set our first row (time 0) to be the starting population we decided on

% Now, the fun part!  We want to turn our rate equation into code, and we
% want to incorporate this idea of a STOCHASTIC simulation.

% Use a for loop to calculate N at all time points
for i = 2:totalSteps
    pGrow = r*stochN(i-1) / (r*stochN(i-1) + r*stochN(i-1)^2/Nstar);
    %       poss of grow   /   poss of grow + poss of shrink
    pShrink = r*stochN(i-1)^2/Nstar / (r*stochN(i-1) + r*stochN(i-1)^2/Nstar);
    
    coin = rand; % trick: matlab's random number generator gives us a number between 0 and 1
    if coin < pGrow % e.g., if pGrow 60%, the coin will be less than 0.6, 60% of the time
        % then, grow!
        stochN(i) = stochN(i-1) + r*stochN(i-1)*dt; 
        % N now   =  old N      +  growth rate * timestep
        
    else % if not, shrink
        stochN(i) = stochN(i-1) - r*stochN(i-1)^2/Nstar*dt;
        
    end % end if statement
    
end % end for loops

% Last but not least: let's plot stochN and N over time to see our results

% Define an array of times
times = (0:dt:(totalSteps-1)*dt);

% Plot!
figure
hold on
plot(times, N) % deterministic prediction from Part I
plot(times, stochN) % stochastic simulation
hold off
xlabel('Time (mins)')
ylabel('Number of cells')
title('Simulations of logistic growth')
legend('Deterministic simulation', 'Stochastic simulation')




















