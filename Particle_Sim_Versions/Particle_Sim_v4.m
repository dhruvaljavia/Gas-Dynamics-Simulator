clc; clear; close all

%% Initialise particles' state
nPart = 20; % no. of particles
size = 5; % particle diameter
box_xmax = 100; % max box width
box_ymax = 100; % max box breadth
time = 5000; % equal to iterations

initX = randi(int32([-(box_xmax/2-size*3) (box_xmax/2-size*3)]/size), 1, nPart)*size; % X-coordinate
initY = randi(int32([-(box_ymax/2-size*3) (box_ymax/2-size*3)]/size), 1, nPart)*size; % Y-coordinate
angle = (360 - 0).*rand(1,nPart) + 0; % angle w.r.t x-axis
discreteSpeed = 0.5 .* ones(1,nPart); % fixed discretized speed in units/iteration
isCollision = zeros(nPart,nPart); % collision detection matrix

%% Simulate particles
X = zeros(time,length(initX));
Y = zeros(time,length(initY));
X(1,:) = initX;
Y(1,:) = initY;

velX = zeros(time,length(initX));
velY = zeros(time,length(initY));
velX(1,:) = discreteSpeed.*cosd(angle);
velY(1,:) = discreteSpeed.*sind(angle);

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
    X(i,:) = X(i-1,:) + velX(i-1,:);
    Y(i,:) = Y(i-1,:) + velY(i-1,:);

    % detect collision
    for m = 1:nPart-1
        for n = m+1:nPart
            dist = sqrt((abs(X(i,m) - X(i,n)))^2 + (abs(Y(i,m) - Y(i,n)))^2);
            if dist < size
                isCollision(m,n) = 1;
            else
                isCollision(m,n) = 0;
            end
        end
    end

     % collision with wall
    velX(i,:) = velX(i-1,:).*(((X(i,:) > xmax) + (X(i,:) < xmin) - 0.5)/(-0.5));
    velY(i,:) = velY(i-1,:).*(((Y(i,:) > ymax) + (Y(i,:) < ymin) - 0.5)/(-0.5));
    X(i,:) = X(i,:) - 2*((X(i,:) > xmax).*rem(X(i,:), xmax));
    Y(i,:) = Y(i,:) - 2*((Y(i,:) > ymax).*rem(Y(i,:), ymax));
    X(i,:) = X(i,:) + 2*((X(i,:) < xmin).*abs(rem(X(i,:), xmin)));
    Y(i,:) = Y(i,:) + 2*((Y(i,:) < ymin).*abs(rem(Y(i,:), ymin)));

    % collision with particle
    for m = 1:nPart-1
        for n = m+1:nPart
            if isCollision(m,n) == 1
                slopeVect = [(X(i,m) - X(i,n)) (Y(i,m) - Y(i,n))];
                slopeVect = slopeVect/norm(slopeVect); % unit vector
                velVect_m = [velX(i-1,m) velY(i-1,m)];
                velVect_n = [velX(i-1,n) velY(i-1,n)];
                
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
                velX(i,[m n]) = [velParallel(1,1)+velPerpendicular(1,1), velParallel(2,1)+velPerpendicular(2,1)];
                velY(i,[m n]) = [velParallel(1,2)+velPerpendicular(1,2), velParallel(2,2)+velPerpendicular(2,2)];
            end
        end
    end

    % animate particles
    % w = waitforbuttonpress;
    % if w == 1
    %     p.XData = X(i,:);
    %     p.YData = Y(i,:);
    %     drawnow;
    % end
end

%% Evaluate thermo-physical quantities



%% Animate particles and make video

% v = VideoWriter('gas');
% v.Quality = 100;
% open(v);

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