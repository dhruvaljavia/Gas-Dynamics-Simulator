clc; clear; close all

%% Initialise particles' state
nPart = 20; % no. of particles
initX = (40 - -40).*rand(1,nPart) + -40; % X-coordinate
initY = (40 - -40).*rand(1,nPart) + -40; % Y-coordinate
size = 3; % particle diameter
angle = (360 - 0).*rand(1,nPart) + 0; % angle w.r.t x-axis
speed = (3 - 0.5).*rand(1,nPart) + 0.5; % speed in units/iteration
box_xmax = 100; % max box width
box_ymax = 100; % max box breadth
time = 1000; % equal to iterations

%% Simulate particles
X = zeros(time,length(initX));
Y = zeros(time,length(initY));
X(1,:) = initX;
Y(1,:) = initY;

velX = speed.*cosd(angle);
velY = speed.*sind(angle);

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

    % collision with wall
    velX = velX.*(((X(i,:) > xmax) + (X(i,:) < xmin) - 0.5)/(-0.5));
    velY = velY.*(((Y(i,:) > ymax) + (Y(i,:) < ymin) - 0.5)/(-0.5));
    X(i,:) = X(i,:) - 2*((X(i,:) > xmax).*rem(X(i,:), xmax));
    Y(i,:) = Y(i,:) - 2*((Y(i,:) > ymax).*rem(Y(i,:), ymax));
    X(i,:) = X(i,:) + 2*((X(i,:) < xmin).*abs(rem(X(i,:), xmin)));
    Y(i,:) = Y(i,:) + 2*((Y(i,:) < ymin).*abs(rem(Y(i,:), ymin)));
end

% animate particles
for i = 2:time
    p.XData = X(i,:);
    p.YData = Y(i,:);
    drawnow;
end