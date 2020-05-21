function impulse_response = monte_carlo(tx, rx, particle, scenario)
% This algorithm carries out the Monte Carlo simulation of an UWOC
% scenario. The simulation is carried out in a ray-basis and it uses a FIFO
% queue to iterate the dynamic tree structure derived from the scattering.

% Preparation of figure
if scenario.plot
    figure('Color','white');
    scatter(0,Inf);
    hold on;
    xlabel('Time (nanoseconds)', 'interpreter', 'latex', 'fontsize', 16);
    ylabel('Power (dBm)', 'interpreter', 'latex', 'fontsize', 16);
end

% Result will be provided as a Buffer
impulse_response = Buffer();

% Buffer initialization
buff = Buffer();

% We generate one direct ray and N random rays and store them
% Rx-Tx direction vector
direction = rx.position - tx.position;
direction = direction/sqrt(direction*direction');
ray = Ray(tx.position, direction, 1, 1, 0, 0);
% We multiply by the output power. The angle corresponds to the difference
% between tx.orientation and the current ray position.
ray.power = ray.power*outputPower(tx.type,tx.params, ...
    ray.getTheta(tx.orientation));

% We store it on the buffer
buff.push(ray);

% Now we randomly generate N random rays
for i = 1:scenario.N
    ray = Ray(tx.position, tx.orientation, 1, 0, 0, 0);
    ray = ray.generateRandomDirection(tx.type,tx.params);
    ray = ray.generateBasis();
    ray.power = ray.power/scenario.N;
    buff.push(ray);
end

iteration_counter = 1;

% Last buffer length
last_buff_len = length(buff.storage);
buffer_empty_rate = 0;

% We start counting time
tic;

% While there exist rays to process...
while(~isempty(buff.storage))
    
    if (~rem(iteration_counter, scenario.info_period))
        fprintf('Rays in buffer: %d\n', length(buff.storage));
        
        elapsed_time = toc;
        aux = (last_buff_len - length(buff.storage))...
            /elapsed_time;
        
        % We update the buffer length for the next report
        last_buff_len = length(buff.storage);
        
        if (aux > 0)
            buffer_empty_rate = 0.75*buffer_empty_rate + 0.25*aux;
            fprintf('Time to end: %d minutes\n',...
                fix(length(buff.storage)/buffer_empty_rate/60));
            fprintf('\n');
        end
        tic;
        
    end
    
    currentRay = buff.pop();
    
    %     fprintf('Ray status (travelling)\n');
    %     fprintf('Hop: %d\n',currentRay.hop);
    %     fprintf('Distance: %1.2f\n',currentRay.distance);
    %     fprintf('Power: %f dBm\n',10*log10(currentRay.power*1000));
    %     fprintf('------------\n\n');
    
    if (currentRay.hop == scenario.max_hops)
        iteration_counter = iteration_counter + 1;
        continue;
    elseif (currentRay.power < scenario.power_threshold)
        iteration_counter = iteration_counter + 1;
        continue;
    end
    
    % We calculate the following impact, and we act depending on the
    % ray nature (direct or random).
    delta = exprnd(1/scenario.beta);
    impact_deltas = currentRay.calculateImpact(scenario);
    
    % Plane index determines if it has impacted on the
    % x(1),y(2) or z(3) plane
    [~, plane_index] = min(impact_deltas);
    % This checks if the ray has impacted on the boundaries or previously
    if (delta == min([delta impact_deltas]))
        % It has impacted on a particle. We generate 1 direct ray and M
        % random rays following the particle's scattering phase function
        currentRay = currentRay.updateRay(delta);
        
        % We generate the direct ray by forcing the new direction
        nextRay = currentRay;
        nextRay.direction = rx.position - nextRay.position;
        nextRay.direction = nextRay.direction/norm(nextRay.direction);
        
        % We update the ray basis
        nextRay = nextRay.generateBasis();
        
        % We update its weigth
        nextRay.power = nextRay.power*...
            outputPower(particle.type,particle.params);
        
        % is direct, we define it
        nextRay.is_direct = 1;
        
        buff.push(nextRay);
        
        % Now the M random rays
        for i = 1:scenario.M
            nextRay = currentRay;
            nextRay = nextRay.generateRandomDirection(particle.type,particle.params);
            nextRay = nextRay.generateBasis();
            nextRay.power = nextRay.power/scenario.M;
            nextRay.is_direct = 0;
            buff.push(nextRay);
        end
        
    else
        % It has impacted on a boundary
        if (currentRay.is_direct)&&(plane_index == 3)
            
            % We update its parameters
            % Distance is the cumulative sum plus the final distance.
            % Power must be updated too
            currentRay = currentRay.updateRay(norm(rx.position - currentRay.position));
            
            % Absorption final power update
            currentRay.power = currentRay.power*...
                exp(-scenario.alfa*currentRay.distance);
            
            % It must be taken into account that this plot is tentative. In
            % essence, it is assuming an isotropic receiver or unitary
            % radius. Poynting vector integration module.
            if scenario.plot
                scatter(currentRay.distance/3e8*scenario.n*1e9, ...
                    10*log10(currentRay.power*1000), 5, 'k');
                drawnow;
            end
            
            impulse_response.push(currentRay);
            
        end
        
    end
    
    % We update the iteration counter
    iteration_counter = iteration_counter + 1;
end