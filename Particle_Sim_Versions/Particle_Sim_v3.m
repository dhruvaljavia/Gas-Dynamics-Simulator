clc; clear; close all

%% Initialise particles' state
nPart = 7; % no. of particles
initX = (40 - -40).*rand(1,nPart) + -40; % X-coordinate
initY = (40 - -40).*rand(1,nPart) + -40; % Y-coordinate
size = 10; % particle diameter
angle = (360 - 0).*rand(1,nPart) + 0; % angle w.r.t x-axis
discreteSpeed = 0.5 .* ones(1,nPart); % fixed discretized speed in units/iteration
box_xmax = 100; % max box width
box_ymax = 100; % max box breadth
time = 10000; % equal to iterations
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
    X(i,:) = X(i-1,:) + velX;
    Y(i,:) = Y(i-1,:) + velY;

    % detect collision
    for m = 1:nPart-1
        for n = m+1:nPart
            dist = sqrt((abs(X(i,m) - X(i,n)))^2 + (abs(Y(i,m) - Y(i,n)))^2);
            if dist < size
                isCollision(m,n) = 1; isCollision(n,m) = 1;
            else
                isCollision(m,n) = 0; isCollision(n,m) = 0;
            end
        end
    end

    % collision with particle
    for m = 1:nPart-1
        for n = m+1:nPart
            if isCollision(m,n) == 1
                slopeVect = [(X(i,m) - X(i,n)) (Y(i,m) - Y(i,n))];
                slopeVect = slopeVect/norm(slopeVect); % unit vector
                velVect_m = [velX(m) velY(m)];
                velVect_n = [velX(n) velY(n)];
                
                % resolve velX and velY along line-of-collision coordinates
                velParallel = [dot(velVect_m,slopeVect)*slopeVect; dot(velVect_n,slopeVect)*slopeVect];
                velPerpendicular = [velVect_m; velVect_n] - velParallel;

                % exchange momentum (parallel velocities swap if same mass)
                velParallel([1 2],:) = velParallel([2 1],:);

                % resolve velocities in velX and velY directions
                velX([m n]) = [velParallel(1,1)+velPerpendicular(1,1), velParallel(2,1)+velPerpendicular(2,1)];
                velY([m n]) = [velParallel(1,2)+velPerpendicular(1,2), velParallel(2,2)+velPerpendicular(2,2)];
            end
        end
    end

    % collision with wall
    velX = velX.*(((X(i,:) > xmax) + (X(i,:) < xmin) - 0.5)/(-0.5));
    velY = velY.*(((Y(i,:) > ymax) + (Y(i,:) < ymin) - 0.5)/(-0.5));
    X(i,:) = X(i,:) - 2*((X(i,:) > xmax).*rem(X(i,:), xmax));
    Y(i,:) = Y(i,:) - 2*((Y(i,:) > ymax).*rem(Y(i,:), ymax));
    X(i,:) = X(i,:) + 2*((X(i,:) < xmin).*abs(rem(X(i,:), xmin)));
    Y(i,:) = Y(i,:) + 2*((Y(i,:) < ymin).*abs(rem(Y(i,:), ymin)));

    % animate particles
    % w = waitforbuttonpress;
    % if w == 1
    %     p.XData = X(i,:);
    %     p.YData = Y(i,:);
    %     drawnow;
    % end
end

% animate particles and make video
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
    frame = getframe(gcf);
    % writeVideo(v, frame);
end
% close(v);