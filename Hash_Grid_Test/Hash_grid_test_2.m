clc; clear; close all

%% Initialise particles' state
nPart = 1000; % no. of particles
size = 2; % particle diameter
box_xmax = 100; % max box width
box_ymax = 100; % max box breadth
time = 1000; % equal to iterations

initX = 2*(box_xmax/2-size*3)*rand(1,nPart) + -(box_xmax/2-size*3); % initial X-coordinates
initY = 2*(box_ymax/2-size*3)*rand(1,nPart) + -(box_ymax/2-size*3); % initial Y-coordinates
initAngle = 360*rand(1,nPart) + 0; % initial angle w.r.t x-axis
initSpeed = 0.5*ones(1,nPart); % initial speed in units/iteration
isCollision = zeros(nPart,nPart); % collision detection matrix
cellSize_x = box_xmax/floor(box_xmax/size); % Hash grid cell size
cellSize_y = box_ymax/floor(box_ymax/size);
cellNos_x = floor(box_xmax/size); % No. of cells
cellNos_y = floor(box_ymax/size);
collisionCtr = 0;

%% Simulate particles
X = zeros(time,length(initX));
Y = zeros(time,length(initY));
X(1,:) = initX;
Y(1,:) = initY;

velX = zeros(time,length(initX));
velY = zeros(time,length(initY));
velX(1,:) = initSpeed.*cosd(initAngle);
velY(1,:) = initSpeed.*sind(initAngle);
velX0 = velX(1,:);
velY0 = velY(1,:);

xmax = box_xmax/2 - size/2;
xmin = -box_xmax/2 + size/2;
ymax = box_ymax/2 - size/2;
ymin = -box_ymax/2 + size/2;

% calculate particles' coordinates for all times
for i = 2:time
    % move particle
    X(i,:) = X(i-1,:) + velX0;
    Y(i,:) = Y(i-1,:) + velY0;
    
    % collision with wall
    velX0 = velX0.*(((X(i,:) > xmax) + (X(i,:) < xmin) - 0.5)/(-0.5));
    velY0 = velY0.*(((Y(i,:) > ymax) + (Y(i,:) < ymin) - 0.5)/(-0.5));
    X(i,:) = X(i,:) - 2*((X(i,:) > xmax).*rem(X(i,:), xmax));
    Y(i,:) = Y(i,:) - 2*((Y(i,:) > ymax).*rem(Y(i,:), ymax));
    X(i,:) = X(i,:) + 2*((X(i,:) < xmin).*abs(rem(X(i,:), xmin)));
    Y(i,:) = Y(i,:) + 2*((Y(i,:) < ymin).*abs(rem(Y(i,:), ymin)));

    % Reset grid
    lookup = zeros(cellNos_x,cellNos_y);

    % Assign particles to grid
    for p = 1:nPart
        cellX = ceil((box_xmax/2 + X(i,p))/cellSize_x); 
        cellY = ceil((box_ymax/2 + Y(i,p))/cellSize_y);
        lookup(cellX,cellY) = p;
    end


    % Check neighbouring cells for each particle and register collision
    for m = 1:nPart
        cellX = ceil((box_xmax/2 + X(i,m))/cellSize_x); 
        cellY = ceil((box_ymax/2 + Y(i,m))/cellSize_y);

        for x = (cellX-1):(cellX+1)
            for y = (cellY-1):(cellY+1)
                if x >= 1 && x <= cellNos_x && y >= 1 && y <= cellNos_y
                    if lookup(x,y) ~= 0
                        n = lookup(x, y);
                        collisionCtr = collisionCtr + 1;
                        if (sqrt((X(i,m) - X(i,n))^2 + (Y(i,m) - Y(i,n))^2) < size) && m~=n
                            slopeVect = [(X(i,m) - X(i,n)) (Y(i,m) - Y(i,n))];
                            slopeVect = slopeVect/norm(slopeVect); % unit vector
                            velVect_m = [velX0(m) velY0(m)];
                            velVect_n = [velX0(n) velY0(n)];

                            % resolve velX and velY along line-of-collision coordinates
                            velParallel = [dot(velVect_m,slopeVect)*slopeVect; dot(velVect_n,slopeVect)*slopeVect];
                            velPerpendicular = [velVect_m; velVect_n] - velParallel;

                            % exchange momentum (parallel velocities swap if same mass)
                            velParallel([1 2],:) = velParallel([2 1],:);

                            % check if velocities separate particles
                            if dot(slopeVect, velParallel(1,:) - velParallel(2,:)) < 0
                                velParallel([1 2],:) = velParallel([2 1],:);
                            end

                            % resolve velocities in velX and velY directions
                            velX0([m n]) = [velParallel(1,1)+velPerpendicular(1,1), velParallel(2,1)+velPerpendicular(2,1)];
                            velY0([m n]) = [velParallel(1,2)+velPerpendicular(1,2), velParallel(2,2)+velPerpendicular(2,2)];
                        end
                    end
                end
            end
        end
    end

    % % collision with particle
    % for m = 1:nPart-1
    %     for n = m+1:nPart
    %         collisionCtr = collisionCtr + 1;
    %         if sqrt((X(i,m) - X(i,n))^2 + (Y(i,m) - Y(i,n))^2) < size
    %             slopeVect = [(X(i,m) - X(i,n)) (Y(i,m) - Y(i,n))];
    %             slopeVect = slopeVect/norm(slopeVect); % unit vector
    %             velVect_m = [velX0(m) velY0(m)];
    %             velVect_n = [velX0(n) velY0(n)];
    % 
    %             % resolve velX and velY along line-of-collision coordinates
    %             velParallel = [dot(velVect_m,slopeVect)*slopeVect; dot(velVect_n,slopeVect)*slopeVect];
    %             velPerpendicular = [velVect_m; velVect_n] - velParallel;
    % 
    %             % exchange momentum (parallel velocities swap if same mass)
    %             velParallel([1 2],:) = velParallel([2 1],:);
    % 
    %             % check if velocities separate particles
    %             if dot(slopeVect, velParallel(1,:) - velParallel(2,:)) < 0
    %                 velParallel([1 2],:) = velParallel([2 1],:);
    %             end
    % 
    %             % resolve velocities in velX and velY directions
    %             velX0([m n]) = [velParallel(1,1)+velPerpendicular(1,1), velParallel(2,1)+velPerpendicular(2,1)];
    %             velY0([m n]) = [velParallel(1,2)+velPerpendicular(1,2), velParallel(2,2)+velPerpendicular(2,2)];
    %         end
    %     end
    % end
    
    % update velocity array
    velX(i,:) = velX0;
    velY(i,:) = velY0;
end

%% Evaluate thermo-physical quantities
speed = zeros(time,length(initX));
for i = 1:time
    speed(i,:) = sqrt(velX(i,:).^2 + velY(i,:).^2);
end


%% Animate particles and make video

% v = VideoWriter('gas');
% v.Quality = 100;
% open(v);

p=plot(X(1,:),Y(1,:),'o','MarkerFaceColor','blue', 'MarkerSize', size*2.5);
axis equal;
axis([-box_xmax/2 box_xmax/2 -box_ymax/2 box_ymax/2]);
xticks(-box_xmax/2:10:box_xmax/2);
yticks(-box_ymax/2:10:box_ymax/2);
grid on;
drawnow;
for i = 2:time
    % w = waitforbuttonpress;
    % if w == 1
    %     p.XData = X(i,:);
    %     p.YData = Y(i,:);
    %     drawnow;
    % end

    p.XData = X(i,:);
    p.YData = Y(i,:);
    drawnow;

    % frame = getframe(gcf);
    % writeVideo(v, frame);
end
% close(v);

%% Plot thermo-physical quantities
