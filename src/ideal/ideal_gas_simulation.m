clear;
close all;

% Simulation parameters

P = 101325; % Pressure of particles / Pa
T = 293; % Thermodynamic temperature of particles / K
m = 4.65e-26; % Mass of particles / kg
N = 1500; % Number of particles to consider

R = 5.87e-10; % Radius of particles / m
t = 1.15e-10; % Total time to observe for / s
n_iter = 200; % Resolution of the simulation

%-------------------------------------------

% Setup

k_B = 1.38e-23; % Boltzmann's constant
v_init = sqrt(3*k_B*T/m); % Initial (uniform) speed of every particle (using 3/2 k_B T = 1/2 m v^2) / m s^-1
V = N*k_B*T/P; % Volume of container / m^3
side = nthroot(V,3); % Side length of the cubic box we're considering / m

dt = t/n_iter; % Time per step / s

s = side*rand(3,N); % Random x,y,z positions for each particle / m
        
angles = -pi + 2*pi*rand(2,N); % Two random angles per particle (between -pi and pi) / rad
unit_vect = [cos(angles(1,:)).*sin(angles(2,:)); ...
                sin(angles(1,:)).*sin(angles(2,:)); ...
                cos(angles(2,:))]; % A random uniformly distributed 3-dimensional unit vector per particle
v = v_init*unit_vect; % Random velocities with uniform specified magnitude / m s^-1



% Data storage

wallcollisions = 0; % Collision counters
particlecollisions = 0;
firsttimeparticlecollisions = 0;

s_history = zeros(3,N,n_iter+1); % 3D matrices to hold the entire history of particle positions and velocities
v_history = zeros(3,N,n_iter+1);

s_history(:,:,1) = s; % Fill in the initial values
v_history(:,:,1) = v;

animation = struct('cdata',[],'colormap',[]); % Movie array

% Main loop

for iter = 1:n_iter
    
    for i = 1:N % For each particle
        
        % Check for collision with wall
        
        for dim = 1:3 % Repeat for x,y,z
            
            if s(dim,i) <= R || side - s(dim,i) <= R % Check distance from each wall
                
                v(dim,i) = -v(dim,i); % Reverse velocity
                wallcollisions = wallcollisions + 1; % Increment collision counter
                
            end
            
        end
        
        % Check for collision with another particle
        
        for j = i+1:N % Check against all particles with a higher index (to avoid duplicate collision checks)

            if (s(1,j)-s(1,i))^2 + ...
                    (s(2,j)-s(2,i))^2 + ...
                    (s(3,j)-s(3,i))^2 <= 4*R^2 % Detect collision
                
                % Unit vector between the two particles' centres
                difference_vector = (s(:,j)-s(:,i));
                difference_magnitude = sqrt((s(1,j)-s(1,i))^2 + ...
                                            (s(2,j)-s(2,i))^2 + ...
                                            (s(3,j)-s(3,i))^2);
                direction = difference_vector/difference_magnitude;
                
                % Component velocities in the collision direction
                
                v_i = dot(v(:,i),direction)*direction;
                v_j = dot(v(:,j),direction)*direction;
                
                delta_v = v_j - v_i;
                
                % Swap these velocities in this direction
                
                v(:,i) = v(:,i) + delta_v;
                v(:,j) = v(:,j) - delta_v;
                
                particlecollisions = particlecollisions + 1; % Increment collision counter
                
                if iter == 1
                    
                    firsttimeparticlecollisions = firsttimeparticlecollisions + 1;
                    
                end
                
            end
        end
    end
    
    % Update all positions
    
    s = s + v*dt;
    
    % Add to history
    
    s_history(:,:,iter+1) = s;
    v_history(:,:,iter+1) = v;
    
end

speeds = sqrt((v_history(1,:,:)).^2+(v_history(2,:,:)).^2+(v_history(3,:,:)).^2); % Matrix of history of speeds

% The matrix 'speeds' can then be used to plot histograms etc. of the speed
% distribution at different times.