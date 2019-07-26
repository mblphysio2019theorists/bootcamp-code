% Plot entropy maximization to understand equilibrium

% Goal: plot TOTAL entropy of the system as a function of L, where L is the
% number of particles on the right half of a two-sectioned lattice with
% omega sites

% Define our parameters
omega = 1e9; % number of 1 nm^3 cubes in an E coli
LTot = 1e4; % total number of particles, which out to be 10uM

% Define L-space
L = (1:100:LTot); % from 1 particle on right side to LTot on right side

SL = (LTot-L).*log(omega) - (LTot-L).*log(LTot-L) + (LTot-L); % define vector 
% of left side's entropy, for a range of possible L's

SR = L.*log(omega) - L.*log(L) + L; % define vector of right side entropy

STot = SL + SR; %total entropy

figure
plot(L,STot, 'g') % total entropy, in green
hold on
plot(L, SL, 'b') % left side entropy, in blue
plot(L,SR, 'k') % right side entropy, in black
hold off
xlabel('L, number of particles on right side')
ylabel('S, entropy')
legend('STot', 'SL', 'SR')





