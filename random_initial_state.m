function magn = random_initial_state(Nlength, type)
%-------------------------------------------------------------------------%
%   Prepares an initial state with specific values for the magnetization in
%   the x- and y-direction (no z-magnetization) and where all cumulants are 
%   equal to zero.
%   
%Parameters:
%   Nlength         Lattice dimension
%   type            Type of initial state
%                   'random': calculate a random state on the Bloch sphere
%                             for each spin in the system.
%                   else: calculates a random state on the Bloch sphere and
%                         duplicates this one for each spin in the system.
%-------------------------------------------------------------------------%

    Nspins = Nlength^2;                                                     % Number of spins

    if type == 'random'
        sign = (-1).^(rand(Nspins,2) < 0.5);                                % Random sign of the x and y magnetization 

        mx = rand(Nspins,1);                                                % Random x-magnetization
        my = sqrt(1 - mx.^2);                                               % Random y-magnetization

        magn = [mx my].*sign;                                               % Apply the random sign to obtain a random 
                                                                            % initial state on the Bloch sphere.
    else 
       mx = rand*(-1)^(rand < 0.5);                                         % Random x-magnetization
       my = sqrt(1 - mx^2)*(-1)^(rand<0.5);                                 % Random y-magnetization
        
       magn = [mx*ones(Nspins,1) my*ones(Nspins,1)];                        % Apply this Bloch sphere state to each spin
    end
end

