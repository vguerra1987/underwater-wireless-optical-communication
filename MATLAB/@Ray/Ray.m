classdef Ray
   % This class is used for easy-handling of the ray-based simulation
    
    properties
        power = 0;
        position = [0,0,0];
        direction = [0,0,0];
        
        % Coordinate basis to properly pass from the direction-defined
        % coordinate system to the scenario's cartesian system.
        basis = randn(3,3);
        
        % Is direct flag
        is_direct = 0;
        
        % number of hops
        hop = 0;
        
        % travelled distance
        distance = 0;
    end
    
    methods
        % Constructor
        function obj = Ray(pos, dir, pow, is_direct, hop, distance)
            obj.power = pow;
            obj.direction = dir;
            obj.position = pos;
            obj = obj.generateBasis();
            obj.is_direct = is_direct;
            obj.hop = hop;
            obj.distance = distance;
        end
        
        % Generate basis
        % This method uses Gram-Schmidt orthogonalization
        function obj = generateBasis(obj)
            % z-component
            uz = obj.direction;
            uz = uz/sqrt(uz*uz');
            
            % x-component
            aux = rand(1,3);
            aux = aux/sqrt(aux*aux');            
            ux = aux - dot(aux,uz)*uz;
            ux = ux/sqrt(ux*ux');
            
            % y-component
            aux = rand(1,3);
            aux = aux/sqrt(aux*aux');            
            uy = aux - dot(aux,ux)*ux - dot(aux,uz)*uz;
            uy = uy/sqrt(uy*uy');
            
            % basis composition
            obj.basis = [ux;uy;uz]';
        end
        
        % This function generates a random direction respect to the current
        % propagation direction.
        function obj = generateRandomDirection(obj, type, params)
            % randomDir is in the ray's reference system
            randomDir = randomDirection(type,params);
            % We pass it to the simulation reference (cartesian)
            obj.direction = randomDir*obj.basis';
        end
        
        % This function updates travelled distance, position, power and
        % hops
        function obj = updateRay(obj, delta)
            obj.position = obj.position + ...
                delta*obj.direction;
            obj.distance = obj.distance + delta;
            if (delta > 0.1)
                obj.power = obj.power*1/delta^2;
            else
                obj.power = obj.power*0.5; % Hemispherical approximation
            end
            obj.hop = obj.hop + 1;
        end
        
        % Calculate impact 
        function impacts = calculateImpact(obj, scenario)
            % The scenario is formed by a cube (2x3 planes)
            
            %plane_interest_coordinate = position + lambda*v
            
            % xmin impact
            aux = scenario.limits(1);
            xmin = (aux - obj.position(1))/obj.direction(1);
            
            % xmax impact
            aux = scenario.limits(1+3);
            xmax = (aux - obj.position(1))/obj.direction(1);
            
            % ymin impact
            aux = scenario.limits(2);
            ymin = (aux - obj.position(2))/obj.direction(2);
            
            % ymax impact
            aux = scenario.limits(2+3);
            ymax = (aux - obj.position(2))/obj.direction(2);
            
            % zmin impact
            aux = scenario.limits(3);
            zmin = (aux - obj.position(3))/obj.direction(3);
            
            % zmax impact
            aux = scenario.limits(3+3);
            zmax = (aux - obj.position(3))/obj.direction(3);
            
            impacts = [min([xmin(xmin>0) xmax(xmax>0)]), ...
                       min([ymin(ymin>0) ymax(ymax>0)]), ...
                       min([zmin(zmin>0) zmax(zmax>0)])];
        end
        
        % Gets the theta elevation respecto to a given direction
        function theta = getTheta(obj, direction)
            % We express direction in terms of the current ray basis
            direction = direction*obj.basis;
            % We extract the theta angle
            theta = acos(direction(3));
        end
        
    end
    
    
end