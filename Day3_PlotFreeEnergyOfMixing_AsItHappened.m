% Day 3 - plot of deltaG = Flory-Huggins free energy of mixing

% Define our parameters
phi = (0 : 0.01 : 1 ); % fraction of particles that are 'A'
chi = (-2 : 1 : 3); % flory parameter = energetic cost of being next to a 
% particle that you don't like

N = 10; % degree of polymerization

% Create an array for deltaG
deltaG = zeros(1, length(phi));

% Free energy of mixing calculation

for k = 1:length(chi)
    for j = 1:length(phi)
        deltaG(j,k) = 1/N*phi(j)*log(phi(j)) + (1-phi(j))*log(1-phi(j))...
            + chi(k)*phi(j)*(1-phi(j));
    end
    
    plot(phi, deltaG)
    legendCell(k) = strcat('chi = ', string(num2cell( chi(k) )));
    
end

xlabel('phi, the fraction of particles are type A')
ylabel('deltaG, free energy of mixing')
legend(legendCell)








