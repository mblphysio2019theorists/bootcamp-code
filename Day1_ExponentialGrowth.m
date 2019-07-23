% % Day 1 - Exponential Growth % %
% In this exercise, we will simulate exponential growth of a bacterial colony,
% segment images to quantify colony growth, and use chi-squared minimisation to 
% find the growth rate that best describes the data, given the assumption of exponential growth

%% Part I. Simulate the exponential growth of a colony
% What equation describes exponential growth?  How can we solve dN/dt = kN?
% a.k.a. Intro to solving differential equations numerically
% 
% Background: If cells double every 1/r min, we can describe the growth of a
% population of cells (N) as a function of time (t) using

% dN/dt = r * N(t) 

% where r has units of 1/min.
% This is a *differential equation* describing the *rate* of cell division.
% If we want to simulate and watch the evolution of the population over
% time, we need to solve this differential equation for N(t).
% But, integration is hard! Euler's method lets us take (1) a rate equation
% like this plus (2) specific initial conditions to find a *numerical
% solution* for population growth... i.e. to simulate/predict the dynamics we want.

clear variables
close all

% Define our biological parameters.
r = 1/30;            % growth rate in 1/min
N0 = 1;              % initial condition: number of cells at time point 0

% Define our simulation parameters.
dt = 5; % time step in minutes. Seems reasonable giving doubling time of around 30 min.
totalSteps = 50; % total number of steps we'll take in our simulation (could define as totalTime / dt)
% totalTime = totalSteps * dt;

% Next, we want to turn our rate equation into code. Useful trick: figure
% out the first time point manually.
N1 = N0 + r * N0 * dt; % after dt min, our population will have grown by:
% (previous pop size) * (growth rate) * length of time we run the clock
% before counting the cells. Take a minute to make sure you internalize.

% Now, we'll simulate all the other time points!

% Define a matlab array that will hold N(t).
N = zeros(totalSteps,1); % this is a tall vector, with 1 column. Check it out in your workspace.
N(1) = N0; % set our first row (time 0) to be the starting population we decided on.

% Use a "for loop" to calculate N at all time points. Review: for loops.
for i = 2:totalSteps
    N(i) = N(i-1) + r * N(i-1) * dt;
end

% Last but not least: we'll plot N(t) to see the results of our simulation.
% Define an array of times, so that we can plot N vs t.
timesSim = (0:dt:(totalSteps-1)*dt);
% Plot population size vs time
figure
plot(timesSim, N)
xlabel('Time (mins)')
ylabel('Number of cells')
title('Population growth over time')

% Bonus check: let's plot our results on a logscale to check that our
% growth is truly exponential.
figure
semilogy(timesSim,N)
xlabel('Time (mins)')
ylabel('Number of cells')
title('Log of population growth over time')

%% Part II. Theory meets experiment: Measure bacterial growth rate
% How can we extract growth rate from time-lapse images of bacteria?
% a.k.a. Review of image analysis

% First, let's import an image and figure out how to segment our bacteria.

% Select directory containing images
folder_name = uigetdir; % opens a GUI to let you navigate to your chosen folder
cd(folder_name); % makes that folder matlab's current working folder

image  = imread('colony_growth_25.tif'); % important random example image
imshow(image)        % Display the image

% Use imtool(image) to find suitable threshold. Second icon from left in imtool tool bar.
threshold = 100; % define threshold
imThresh = image > threshold; % do the thresholding
imshow(imThresh); % check our segmentation

% Side note: what is ImThresh, exactly?
a = [2,4,1,6,4,3];
a > 3; % check out 'a' to see example of binary mask

% Measure the area of the colony (a proxy for population size)
area = sum(sum(imThresh));

% Now, the real deal! Let's take this measurement at all time points.

% The "structure array" called fileList saves info about our images
fileList = dir([folder_name, '/colony_growth*.tif']); %check out fileList. See how we can use it to call images.

% Define arrays for our times and our area measurements.
area = zeros(1,length(fileList));    % This is a preallocated array of zeros

% We'll do the same for the time at which each image was taken
dt = 5;    % the interval between image acquisitions (in minutes)
times = [];       % this is another way to preallocate an array

% The big moment: we'll use a for loop to open, segment, and measure cell area for each of our images.
for i = 1:length(fileList)
    image = imread(['ColonyGrowthData/',fileList(i).name]);
    imThresh = image > threshold;
    area(i) = sum(sum(imThresh));    % sum over the columns and the rows to find the total area
    times(i) = dt*(i-1);    %this gives us the time at which each area is measured
    
    % add after, to see whether we're doing okay
    figure(3)
    imshow(imThresh)
    %drawnow
end 

% Finally, let's plot our measurent of area v time.
plot(times, area, '.')       % '.' causes the data to be plotted as points rather than a continuous line
xlabel('Time (mins)')
ylabel('Colony Area (pixels^2)')
title('Measured Colony Growth')

% Want to see a log plot, to check if it's exponential, like we did in Part I?
semilogy(times, area, '.')       % '.' causes the data to be plotted as points rather than a continuous line
xlabel('Time (mins)')
ylabel('Colony Area (pixels^2)')
title('Measured Colony Growth (log scale)')

% Note the noise in our early segmentations. We could use regionprops and
% throw out everything except the biggest continuous region... but we won't
% get into this today.

%% Part III. Chi-square minimization to find the bacterial division time
% Taking it to the finish line: how do we find a best fit for E coli's
% division time (how long it takes to make a new cell) based on our data?
% a.k.a. Use chi-squared minimization to find the growth rate that best describes our data, given the assumption of exponential growth

% Give background on board to the chi-square statistic. 

% First, let's look at our data again - on a log plot.
semilogy(times, area, '.')       % '.' causes the data to be plotted as points rather than a continuous line
xlabel('Time (mins)')
ylabel('Colony Area (pixels^2)')
title('Measured Colony Growth (log scale)')

% Let's fit the exponential part of the colony's growth -- not the first ten time points.
area(10:end); % here's how we take a subset of the array "area" (remove the semi-colon to print)

% For one time point, calculate the error between data and prediction
% define a possible range of growth rates (r's) to test
rRange = 0.01:0.00001:0.03;    % start : increment : end, centered around r = 0.02 (30 min doubling time)

% Evaluate the fit for one sample r, to help us develop our logic
rTest = 0.02;

dataForFit = area(10:end);     % Define the range of the data we will use, for simplicity's sake
initialCondition = dataForFit(1); 
timesForFit = dt*(0:length(dataForFit)-1);     % Time vector for fitting (dt is 5 mins - defined earlier)

% Now we're going to employ the chi^2 equation to calculate the difference 
% between our data and the prediction based on a growth rate of 0.02. 
theoryTest = initialCondition * exp(rTest* timesForFit);
chiSqTest = sum((theoryTest-dataForFit).^2); 

% Visualize how well our sample fit does.
semilogy(timesForFit,dataForFit,'k.')       % our data
hold on
semilogy(timesForFit,initialCondition*exp(rTest*timesForFit), 'r')          % best fit line
hold off
xlabel('time starting from time point 10 (mins)')
ylabel('colony area (pixels^2)')
legend('data', 'guess')

% Repeat the above, now looping through all values of  in rRange to compute the  statistic for each one.
for i = 1:length(rRange)
    theory = initialCondition * exp(rRange(i)*timesForFit);
    chiSq(i) = sum((theory - dataForFit).^2);
end

% Plot the chi-square for each value of r
semilogy(rRange, chiSq, 'o')
xlabel('Growth rate, r (mins^-1)')
ylabel('\chi^2 statistic')
title('\chi^2 for a range of growth rates')

% Zoom in and use data cursor to estimate the best fit of r.
rFit = 0.0186; % could also use find(chiSq == min(chiSq)). rRange(865) = ?

% This corresponds to a doubling time of
doublingTime = log(2)/rFit % COOOOOL. We did it. How should we explain this log2/r thing? N(t) = N*e^rt

% Now we can do our most accurate simulation of colony growth yet, using
% the growth rate we measured from data.
bestFit = initialCondition*exp(rFit*timesForFit);    

% Lastly, let's plot the fit of our data.
semilogy(times, area,'k.')       % our data
hold on

timeVectFit = timesForFit + 10 * dt;    % because we only data from timepoint 10 (50 mins), 
                                      % we need to start our time vector at 50 mins for plotting (frame 10 * 5 min intervals)
semilogy(timeVectFit,bestFit, 'r')    % best fit line
hold off
xlabel('Time (mins)')
ylabel('Colony Area (pixels^2)')
legend('Data', 'Fit')

% Fun bonus: why is rate of growth at the end of the experiment higher?
% Example of theory-inspired biological hypotheses. Why it's so great to
% have simple calculations and key theories in your biologist toolbox!

% Recap: return to goals at top of code. Give students a minute to
% annotate and absorb.



