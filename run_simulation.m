%-------------------------------------------------------------------------%
% This is an example simulation of the dissipative XYZ Heisenberg model in
% a square (checkerboard) lattice. The type of dissipation is spontaneous
% local dissipation. The heterodyne unravelling is used for the quantum 
% trajectories.
%
% In this example we calculate the following observables:
%   
%   m2x        Second moment in the x-direction (Classical + Quantum)
%   m2y        Second moment in the y-direction (Classical + Quantum)
%   m2z        Second moment in the z-direction (Classical + Quantum)
%
% NOTE: one can extract other observables from the simulation as well. 
%       This is left for another example.
%
% Each of the above observables are calculate for a number of quantum 
% trajectories. For a certain set of parameter values all independent 
% trajectories are stored in a cell array for each respective observable. 
% These cell arrays are then stored in the corresponding cell array element
% of the specific parameter set.
%
% These cell structures are then saved to a .mat file. This file can later 
% be loaded to be analyzed and e.g. combined with other (independent) 
% simulations on the same parameter set.
%
% More information on the method and parameters can be found in the arXiv
% article at arXiv:2209.13377v3. 
%-------------------------------------------------------------------------%

parpool([4 64])                                                             % Start a parallel pool

% System parameters:
gamma = 1;                                                                  % Dissipation strength (local emission)
Jxeff = 0.9*gamma;                                                          % Coupling strength in the x-direction
Jzeff = gamma;                                                              % Coupling strength in the z-direction
dt = 1e-3;                                                                  % Time step in the Euler time evolution
N = 4;                                                                      % 2D-lattice dimension (N x N)-lattice

effi = 1;                                                                   % Measurement efficiency
Jy_list = (1.01:0.01:1.15)*gamma;                                           % List of coupling strength in the y-direction
timearray = 0:0.5:10;                                                       % List of times at which result is saved

ntraj = 16;                                                                 % Number of simulated trajectories


% Storage cell arrays for each coupling strength Jy:
lenJy = length(Jy_list);            
m2xc = cell(1, lenJy);                                                      % Cell array to store the trajectories
m2yc = cell(1, lenJy);                                                      % Cell array to store the trajectories
m2zc = cell(1, lenJy);                                                      % Cell array to store the trajectories

rng('shuffle');                                                             % Shuffles on local worker only
saveseeds = zeros(1,lenJy);                                                 % Array to save seeds for RNG


% The simulation:
for i = 1:length(Jy_list)

    Jyeff = Jy_list(i)*gamma;                                               % Coupling strength in the y-direction

    % Temporary storage cell arrays for each coupling strength (Jy) value
    r1 = cell(1, ntraj);                                                    % Each element is a quantum trajectory
    r2 = cell(1, ntraj);                                                    % Each element is a quantum trajectory
    r3 = cell(1, ntraj);                                                    % Each element is a quantum trajectory
    
    
    tic

    % Shuffle rng on each parallel worker and save the seeds.
    seed_offset = randi(floor(intmax/10));
    saveseeds(i) = seed_offset;                                             % Save the seed for this coupling strength

    parfor j = 1:ntraj
        % Unique rng seed for each quantum trajectory
	    rng(j + seed_offset);                                               % Each rng seed is determined by the trajectory number.

        try
            % Initialise the system in a certain state:
            initial_state = random_initial_state(N, 'random');          

            % The quantum trajectory evolution:
            [m2x, m2y, m2z]=heterodyne_trajectory_no_2ndorder(Jxeff,...     % This evolution does not include the second order 
                Jyeff, Jzeff,gamma, timearray, dt, initial_state, N, effi); % Wiener noise term. To include the second order
                                                                            % noise terms use "heterodyne_trajectory_2ndorder(...)".
            % Store the results of each trajectory:
            r1{j} = m2x;
            r2{j} = m2y;
            r3{j} = m2z;
        catch
            r1{j} = NaN(1, length(timearray));                              % Trajectory result set to NaN
            r2{j} = NaN(1, length(timearray));                              % Trajectory result set to NaN
            r3{j} = NaN(1, length(timearray));                              % Trajectory result set to NaN
        end
    end
    toc

    % Save all trajectories for a given coupling strength in the
    % y-direction:
    m2xc{i} = r1;
    m2yc{i} = r2;
    m2zc{i} = r3;

    save('run_1.mat')
end
