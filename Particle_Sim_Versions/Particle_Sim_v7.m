clc; clear; close all

%% Choose type of simulation and simulation settings
twoParticleMix = 1;
innerPartitions = 1;
heatTransfer = 0;
reactiveFlow = 0;
playAnimationOnly = 1;
makeVideo = 0;

%% Enter simulation parameters 
% Particle A
massPart_A = 2; % mass of particles 
nPart_A = 50; % no. of particles
size_A = 2; % particle diameter
speed_A = 0.5; % initial speed in units/iteration

% Particle B
massPart_B = 5;
nPart_B = 10;
size_B = 5;
speed_B = 0.5;

box_xmax = 100; % max box width
box_ymax = 100; % max box breadth
time = 1000; % equal to iterations

rectangle('Position',[-box_xmax/2,-box_ymax/2,box_xmax,box_ymax]);
axis equal;
axis([-box_xmax/2 box_xmax/2 -box_ymax/2 box_ymax/2]);
xticks(-box_xmax/2:10:box_xmax/2);
yticks(-box_ymax/2:10:box_ymax/2);
grid on;
title('Select initial region of A (and B) particles');
rect_A = drawrectangle(Color="b");
if twoParticleMix == 1
    rect_B = drawrectangle(Color="r");
end
pause;
initLoc_A = rect_A.Position;
if twoParticleMix ==1
    initLoc_B = rect_B.Position;
else
    initLoc_B = [0 0 0 0];
end

title('Select end points of walls');
ctr=1;
while 1
    wall = drawline(Color="k");
    wallPts(1:2,1:2,ctr) = wall.Position;
    b = waitforbuttonpress;
    c = get(gcf,'CurrentCharacter');
    if c == 'x'
        break;
    end
    ctr=ctr+1;
end
close;


%% Initialise particles' state

if twoParticleMix == 0
    nPart_B = 0;
end

size = [size_A size_B];
massPart = [massPart_A massPart_B];
initX_A = initLoc_A(3)*rand(1,nPart_A) + initLoc_A(1); % initial X-coordinates
initY_A = initLoc_A(4)*rand(1,nPart_A) + initLoc_A(2); % initial Y-coordinates
initX_B = initLoc_B(3)*rand(1,nPart_B) + initLoc_B(1);
initY_B = initLoc_B(4)*rand(1,nPart_B) + initLoc_B(2);
initX = [initX_A initX_B];
initY = [initY_A initY_B];
initAngle = 360*rand(1,nPart_A+nPart_B) + 0; % initial angle w.r.t x-axis
initSpeed = [speed_A*ones(1,nPart_A) speed_B*ones(1,nPart_B)];
isCollision = zeros(nPart_A+nPart_B,nPart_A+nPart_B); % collision detection matrix

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

xmax = box_xmax/2 - [size_A/2*ones(1,nPart_A) size_B/2*ones(1,nPart_B)];
xmin = -xmax;
ymax = box_ymax/2 - [size_A/2*ones(1,nPart_A) size_B/2*ones(1,nPart_B)];
ymin = -ymax;

% calculate particles' coordinates for all times
for i = 2:time
    % move particle
    X(i,:) = X(i-1,:) + velX0;
    Y(i,:) = Y(i-1,:) + velY0;

    % detect inter-particle collision
    for m = 1:(nPart_A+nPart_B)-1
        for n = m+1:(nPart_A+nPart_B)
            dist = sqrt((abs(X(i,m) - X(i,n)))^2 + (abs(Y(i,m) - Y(i,n)))^2);
            if dist < size((m>nPart_A)+1)/2 + size((n>nPart_A)+1)/2
                isCollision(m,n) = 1;
            else
                isCollision(m,n) = 0;
            end
        end
    end

    % collision with particle
    for m = 1:(nPart_A+nPart_B)-1
        for n = m+1:(nPart_A+nPart_B)
            if isCollision(m,n) == 1
                slopeVect = [(X(i,m) - X(i,n)) (Y(i,m) - Y(i,n))];
                slopeVect = slopeVect/norm(slopeVect); % unit vector
                velVect_m = [velX0(m) velY0(m)];
                velVect_n = [velX0(n) velY0(n)];
                
                % resolve velX and velY along line-of-collision coordinates
                velParallel = [dot(velVect_m,slopeVect)*slopeVect; dot(velVect_n,slopeVect)*slopeVect];
                velParallelMag = [dot(velVect_m,slopeVect); dot(velVect_n,slopeVect)];
                velPerpendicular = [velVect_m; velVect_n] - velParallel;

                % exchange momentum
                m1 = massPart((m>nPart_A)+1);
                m2 = massPart((n>nPart_A)+1);
                oldVelParallel = velParallel([1 2],:);
                v1f = ((m1-m2)*velParallelMag(1) + 2*m2*velParallelMag(2))/(m1+m2);
                v2f = ((m2-m1)*velParallelMag(2) + 2*m1*velParallelMag(1))/(m1+m2);
                velParallel = [v1f*slopeVect; v2f*slopeVect];

                % check if velocities separate particles
                if dot(slopeVect, velParallel(1,:) - velParallel(2,:)) < 0
                    velParallel([1 2],:) = oldVelParallel;
                end

                % resolve velocities in velX and velY directions
                velX0([m n]) = [velParallel(1,1)+velPerpendicular(1,1), velParallel(2,1)+velPerpendicular(2,1)];
                velY0([m n]) = [velParallel(1,2)+velPerpendicular(1,2), velParallel(2,2)+velPerpendicular(2,2)];
            end
        end
    end

    % collision with box wall
    velX0 = velX0.*(((X(i,:) > xmax) + (X(i,:) < xmin) - 0.5)/(-0.5));
    velY0 = velY0.*(((Y(i,:) > ymax) + (Y(i,:) < ymin) - 0.5)/(-0.5));
    X(i,:) = X(i,:) - 2*((X(i,:) > xmax).*rem(X(i,:), xmax));
    Y(i,:) = Y(i,:) - 2*((Y(i,:) > ymax).*rem(Y(i,:), ymax));
    X(i,:) = X(i,:) + 2*((X(i,:) < xmin).*abs(rem(X(i,:), xmin)));
    Y(i,:) = Y(i,:) + 2*((Y(i,:) < ymin).*abs(rem(Y(i,:), ymin)));

    % collision with inside walls
    
    
    % update velocity array
    velX(i,:) = velX0;
    velY(i,:) = velY0;
end

%% Evaluate thermo-physical quantities
speedArray = zeros(time,length(initX));
for i = 1:time
    speedArray(i,:) = sqrt(velX(i,:).^2 + velY(i,:).^2);
end

isBottomWallCollision = zeros(time,length(initX));
for i = 2:time
    isBottomWallCollision(i,:) = velY(i-1,:)./velY(i,:) == -1 & velY(i-1,:) < 0;
end
pressure = zeros(time,1);
for i = 2:time
    pressure(i) = (massPart(1)*dot(velY(i,1:nPart_A),isBottomWallCollision(i,1:nPart_A))*2)/box_ymax; % (m*dv/dt)/L where dt=1 for 1 itr
    if twoParticleMix == 1
        pressure(i) = pressure(i) + (massPart(2)*dot(velY(i,nPart_A+1:end),isBottomWallCollision(i,nPart_A+1:end))*2)/box_ymax;
    end
end

%% Plot thermo-physical quantities
% figure();
% plot(pressure)
% mean(pressure)

%% Animate particles

if playAnimationOnly == 1
    figure();
    p_A = plot(X(1,1:nPart_A),Y(1,1:nPart_A),'o','MarkerFaceColor','blue', 'MarkerSize', size_A*2.5);
    hold on;
    p_B = plot(X(1,nPart_A+1:end),Y(1,nPart_A+1:end),'o','MarkerFaceColor','red', 'MarkerSize', size_B*2.5);
    axis equal;
    axis([-box_xmax/2 box_xmax/2 -box_ymax/2 box_ymax/2]);
    xticks(-box_xmax/2:10:box_xmax/2);
    yticks(-box_ymax/2:10:box_ymax/2);
    grid on;
    drawnow;
    
    title("Press any key to start");
    pause;

    for i = 2:time
        % w = waitforbuttonpress;
        % if w == 1
        %     p_A.XData = X(i,1:nPart_A);
        %     p_A.YData = Y(i,1:nPart_A);
        %     p_B.XData = X(i,nPart_A+1:end);
        %     p_B.YData = Y(i,nPart_A+1:end);
        %     drawnow;
        % end
    
        p_A.XData = X(i,1:nPart_A);
        p_A.YData = Y(i,1:nPart_A);
        if twoParticleMix == 1
            p_B.XData = X(i,nPart_A+1:end);
            p_B.YData = Y(i,nPart_A+1:end);
        end
        drawnow;
    end
end

%% Make video

if makeVideo == 1
    figure();
    p_A = plot(X(1,1:nPart_A),Y(1,1:nPart_A),'o','MarkerFaceColor','blue', 'MarkerSize', size_A*2.5);
    hold on;
    p_B = plot(X(1,nPart_A+1:end),Y(1,nPart_A+1:end),'o','MarkerFaceColor','red', 'MarkerSize', size_B*2.5);
    axis equal;
    axis([-box_xmax/2 box_xmax/2 -box_ymax/2 box_ymax/2]);
    xticks(-box_xmax/2:10:box_xmax/2);
    yticks(-box_ymax/2:10:box_ymax/2);
    grid on;
    drawnow;
    
    v = VideoWriter('gas');
    v.Quality = 100;
    open(v);
    
    for i = 2:time
        p_A.XData = X(i,1:nPart_A);
        p_A.YData = Y(i,1:nPart_A);
        if twoParticleMix == 1
            p_B.XData = X(i,nPart_A+1:end);
            p_B.YData = Y(i,nPart_A+1:end);
        end
        drawnow;
        
        frame = getframe(gcf);
        writeVideo(v , frame);
    end
    
    close(v);
end