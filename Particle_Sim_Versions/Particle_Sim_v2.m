clc; clear; close all

%% Initialise particles' state
nPart = 5; % no. of particles
initX = (40 - -40).*rand(1,nPart) + -40; % X-coordinate
initY = (40 - -40).*rand(1,nPart) + -40; % Y-coordinate
size = 10; % particle diameter
angle = (360 - 0).*rand(1,nPart) + 0; % angle w.r.t x-axis
itrSkip = int32(ceil(10 .* rand(1,nPart))); % iteration skip for having different speeds
discreteSpeed = 0.5 .* ones(1,nPart); % fixed discretized speed in units/iteration
box_xmax = 100; % max box width
box_ymax = 100; % max box breadth
time = 1000; % equal to iterations
isCollision = zeros(nPart,nPart); % collision detection matrix

%% Simulate particles
X = zeros(time,length(initX));
Y = zeros(time,length(initY));
X(1,:) = initX;
Y(1,:) = initY;

velX = discreteSpeed.*cosd(angle);
velY = discreteSpeed.*sind(angle);

xmax = box_xmax/2 - size/2;
xmin = -box_xmax/2 + size/2;
ymax = box_ymax/2 - size/2;
ymin = -box_ymax/2 + size/2;

p=plot(X(1,:),Y(1,:),'o','MarkerFaceColor','blue', 'MarkerSize', size*2.5);
axis equal;
axis([-box_xmax/2 box_xmax/2 -box_ymax/2 box_ymax/2]);
xticks(-box_xmax/2:10:box_xmax/2);
yticks(-box_ymax/2:10:box_ymax/2);
grid on;
drawnow;

% calculate particles' coordinates for all times
for i = 2:time
    % move particle
    X(i,:) = X(i-1,:) + (velX .* (mod(i,itrSkip) == 0));
    Y(i,:) = Y(i-1,:) + (velY .* (mod(i,itrSkip) == 0));

    % detect collision
    for m = 1:nPart-1
        for n = m+1:nPart
            dist = sqrt((abs(X(i,m) - X(i,n)))^2 + (abs(Y(i,m) - Y(i,n)))^2);
            if dist < size
                isCollision(m,n) = 1; isCollision(n,m) = 1;
                % disp("Collision detected!");
            else
                isCollision(m,n) = 0; isCollision(n,m) = 0;
            end
        end
    end

    % disp(X(i,:));
    % disp(Y(i,:));
    % disp(dist);
    % disp(isCollision);

    % collision with particle [TRY USING COMPLEX VECTORS]
    for m = 1:nPart-1
        for n = m+1:nPart
            if isCollision(m,n) == 1
                % disp("exchanging momentum")
                slope = (Y(i,m) - Y(i,n))/(X(i,m) - X(i,n));
                slopeDeg = atand(slope);
                % disp(slopeDeg)
                % disp(velX([m n]));
                % disp(velY([m n]));
                
                if slopeDeg >= 0
                    % resolve velX and velY along line-of-collision coordinates
                    velParallel = velX([m n]).*cosd(slopeDeg) + velY([m n]).*sind(slopeDeg);
                    velPerpendicular = velY([m n]).*cosd(slopeDeg) - velX([m n]).*sind(slopeDeg);
                    
                    % exchange momentum (parallel velocities swap if same mass)
                    velParallel([1 2]) = velParallel([2 1]);
                    
                    % resolve velocities in velX and velY directions
                    velX([m n]) = velParallel.*cosd(slopeDeg) - velPerpendicular.*sind(slopeDeg);
                    velY([m n]) = velParallel.*sind(slopeDeg) + velPerpendicular.*cosd(slopeDeg);
                    % disp(velX([m n]));
                    % disp(velY([m n]));
                else
                    slopeDeg = abs(slopeDeg);
                    % resolve velX and velY along line-of-collision coordinates
                    velParallel = velY([m n]).*sind(slopeDeg) - velX([m n]).*cosd(slopeDeg);
                    velPerpendicular = velY([m n]).*cosd(slopeDeg) + velX([m n]).*sind(slopeDeg);
                    
                    % exchange momentum (parallel velocities swap if same mass)
                    velParallel([1 2]) = velParallel([2 1]);
                    
                    % resolve velocities in velX and velY directions
                    velX([m n]) = velPerpendicular.*sind(slopeDeg) - velParallel.*cosd(slopeDeg);
                    velY([m n]) = velParallel.*sind(slopeDeg) + velPerpendicular.*cosd(slopeDeg);
                    % disp(velX([m n]));
                    % disp(velY([m n]));
                end
            end
            X(i,:) = X(i-1,:) + (velX .* (mod(i,itrSkip) == 0));
            Y(i,:) = Y(i-1,:) + (velY .* (mod(i,itrSkip) == 0));
        end
    end

    % collision with wall
    velX = velX.*(((X(i,:) > xmax) + (X(i,:) < xmin) - 0.5)/(-0.5));
    velY = velY.*(((Y(i,:) > ymax) + (Y(i,:) < ymin) - 0.5)/(-0.5));
    X(i,:) = X(i,:) - 2*((X(i,:) > xmax).*rem(X(i,:), xmax));
    Y(i,:) = Y(i,:) - 2*((Y(i,:) > ymax).*rem(Y(i,:), ymax));
    X(i,:) = X(i,:) + 2*((X(i,:) < xmin).*abs(rem(X(i,:), xmin)));
    Y(i,:) = Y(i,:) + 2*((Y(i,:) < ymin).*abs(rem(Y(i,:), ymin)));

    % % animate particles
    % w = waitforbuttonpress;
    % if w == 1
    %     p.XData = X(i,:);
    %     p.YData = Y(i,:);
    %     drawnow;
    % end
end

% animate particles
for i = 2:time
    w = waitforbuttonpress;
    if w == 1
        p.XData = X(i,:);
        p.YData = Y(i,:);
        drawnow;
    end
    % p.XData = X(i,:);
    % p.YData = Y(i,:);
    % drawnow;
end