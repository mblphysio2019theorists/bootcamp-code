% Plot probability of open channel as a function of the energy of the open
% state (where the channel as two possible microstates: open and closed)

% Define a range of E_open
E_open = (-20:0.1:10);

% Define E_closed
E_closed = -10; % in k_B*T units

figure
plot(E_open, exp(-E_open) ./ (exp(-E_open) + exp(-E_closed)) )
xlabel('E open')
ylabel('probability of open state')