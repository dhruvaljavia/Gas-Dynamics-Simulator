clc; clear; close all

%% Choose type of simulation and simulation settings
twoParticleMix = 1;
innerPartitions = 1;
heatTransfer = 0;
reactiveFlow = 0;
playAnimationOnly = 1;
makeVideo = 0;

%% Input simulation parameter values - by typing
% Particle A
massPart_A = 5; % mass of particles 
nPart_A = 50; % no. of particles
size_A = 5; % particle diameter
speed_A = 0.5; % initial speed in units/iteration

% Particle B
massPart_B = 3;
nPart_B = 70;
size_B = 3;
speed_B = 0.5;

% Environment and post-processing settings
box_xmax = 300; % max box width
box_ymax = 50; % max box breadth
time = 1000; % equal to iterations
movingAverageWindow = 100; % moving average window size for calculation of physical parameters

%% Input simulation parameter values - interactively

% Initial particle region
fig = figure();
fig.WindowState = 'maximized';
rectangle('Position',[-box_xmax/2,-box_ymax/2,box_xmax,box_ymax]);
axis equal;
axis([-box_xmax/2 - 10, box_xmax/2 + 10, -box_ymax/2 - 10, box_ymax/2 + 10]);
xticks(-box_xmax/2:10:box_xmax/2);
yticks(-box_ymax/2:10:box_ymax/2);
grid on;
title('Select initial region of A (and B) particles');
rect_A = drawrectangle(Color="blue");
if twoParticleMix == 1
    rect_B = drawrectangle(Color="red");
end
pause;
initLoc_A = rect_A.Position;
if twoParticleMix ==1
    initLoc_B = rect_B.Position;
else
    initLoc_B = [0 0 0 0];
end

% Wall locations
if innerPartitions == 1
    title('Select end points of walls');
    ctr=1;
    while 1
        wall = drawline(Color="black");
        wallPts(1:2,1:2,ctr) = wall.Position;
        b = waitforbuttonpress;
        c = get(gcf,'CurrentCharacter');
        if c == 'x'
            break;
        end
        ctr=ctr+1;
    end
    noOfWalls = size(wallPts,3);
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
if innerPartitions == 1
    isCollisionWallEnd = zeros(nPart_A+nPart_B,2,noOfWalls); % wall end collision detection matrix
    isCollisionWallEdge = zeros(noOfWalls,nPart_A+nPart_B); % wall edge collision detection matrix
    isWallJitter = zeros(noOfWalls,nPart_A+nPart_B); % wall jitter detection matrix
end

% Making the wall box
if innerPartitions == 1
    wallNormalVect = zeros(noOfWalls,2); % unit wall normal vector
    for w = 1:noOfWalls
        wallNormalVect(w,:) = [-(wallPts(2,2,w)-wallPts(1,2,w)), wallPts(2,1,w)-wallPts(1,1,w)];
        wallNormalVect(w,:) = wallNormalVect(w,:)/norm(wallNormalVect(w,:));
    end
    addVect_A = (size(1)/2)*wallNormalVect;
    wallBox_A = zeros(4,2,noOfWalls);
    if twoParticleMix == 1
        addVect_B = (size(2)/2)*wallNormalVect;
        wallBox_B = zeros(4,2,noOfWalls);
    end
    for w = 1:noOfWalls
        wallBox_A(:,:,w) = [wallPts(1,:,w)+addVect_A(w,:); wallPts(1,:,w)-addVect_A(w,:); wallPts(2,:,w)-addVect_A(w,:); wallPts(2,:,w)+addVect_A(w,:)];
        if twoParticleMix == 1
            wallBox_B(:,:,w) = [wallPts(1,:,w)+addVect_B(w,:); wallPts(1,:,w)-addVect_B(w,:); wallPts(2,:,w)-addVect_B(w,:); wallPts(2,:,w)+addVect_B(w,:)];
        end
    end
end

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

    % detect inter-particle and inside wall collision and change velocities
    for m = 1:(nPart_A+nPart_B)-1
        for n = m+1:(nPart_A+nPart_B)
            distPart = sqrt((X(i,m) - X(i,n))^2 + (Y(i,m) - Y(i,n))^2);
            isCollisionPart = distPart < size((m>nPart_A)+1)/2 + size((n>nPart_A)+1)/2;
            if isCollisionPart == 1
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
        
        if innerPartitions == 1
            for w = 1:noOfWalls
                distWallEnd1 = sqrt((wallPts(1,1,w) - X(i,m))^2 + (wallPts(1,2,w) - Y(i,m))^2);
                distWallEnd2 = sqrt((wallPts(2,1,w) - X(i,m))^2 + (wallPts(2,2,w) - Y(i,m))^2);
                isCollisionWallEnd(m,:,w) = [distWallEnd1 distWallEnd2] < (size((m>nPart_A)+1)/2)*ones(1,2);
            end
        end
    end
    if innerPartitions == 1
        for w = 1:noOfWalls
            distWallEnd1 = sqrt((wallPts(1,1,w) - X(i,end))^2 + (wallPts(1,2,w) - Y(i,end))^2);
            distWallEnd2 = sqrt((wallPts(2,1,w) - X(i,end))^2 + (wallPts(2,2,w) - Y(i,end))^2);
            isCollisionWallEnd(end,:,w) = [distWallEnd1 distWallEnd2] < (size(2)/2)*ones(1,2);
            
            in = inpolygon(X(i,1:nPart_A),Y(i,1:nPart_A),wallBox_A(:,1,w),wallBox_A(:,2,w));
            isWallJitter(w,1:nPart_A) = isCollisionWallEdge(w,1:nPart_A) & in;
            isCollisionWallEdge(w,1:nPart_A) = in;
            if twoParticleMix == 1
                in = inpolygon(X(i,nPart_A+1:end),Y(i,nPart_A+1:end),wallBox_B(:,1,w),wallBox_B(:,2,w));
                isWallJitter(w,nPart_A+1:end) = isCollisionWallEdge(w,nPart_A+1:end) & in;
                isCollisionWallEdge(w,nPart_A+1:end) = in;
            end
            
            % remove overlapping wall detection for end and edge regions
            isCollisionWallEnd(:,1,w) = (~(isCollisionWallEnd(:,1,w) & (isCollisionWallEdge(w,:))') & isCollisionWallEnd(:,1,w));
            isCollisionWallEnd(:,2,w) = (~(isCollisionWallEnd(:,2,w) & (isCollisionWallEdge(w,:))') & isCollisionWallEnd(:,2,w));
        end
    end

    % collision with box wall
    velX0 = velX0.*((((X(i,:) > xmax) | (X(i,:) < xmin)) - 0.5)/(-0.5));
    velY0 = velY0.*((((Y(i,:) > ymax) | (Y(i,:) < ymin)) - 0.5)/(-0.5));
    X(i,:) = X(i,:) - 2*((X(i,:) > xmax).*(X(i,:) - xmax));
    Y(i,:) = Y(i,:) - 2*((Y(i,:) > ymax).*(Y(i,:) - ymax));
    X(i,:) = X(i,:) + 2*((X(i,:) < xmin).*abs(xmin - X(i,:)));
    Y(i,:) = Y(i,:) + 2*((Y(i,:) < ymin).*abs(ymin - Y(i,:)));

    % collision with inside walls
    if innerPartitions == 1
        for m = 1:nPart_A+nPart_B
            for w = 1:noOfWalls
                if isCollisionWallEnd(m,1,w) == 1
                    slopeVect = [(X(i,m) - wallPts(1,1,w)) (Y(i,m) - wallPts(1,2,w))];
                    slopeVect = slopeVect/norm(slopeVect); % unit vector
                    velVect = [velX0(m) velY0(m)];
                    
                    % resolve velX and velY along line-of-collision coordinates
                    velParallel = dot(velVect,slopeVect)*slopeVect;
                    velParallelMag = dot(velVect,slopeVect);
                    velPerpendicular = velVect - velParallel;

                    % reverse parallel velocity
                    velParallel = -velParallel;

                    % check if new velocity moves particle away from wall
                    if dot(velParallel,slopeVect) < 0
                        velParallel = -velParallel;
                    end

                    % resolve velocities in velX and velY directions
                    velX0(m) = velParallel(1)+velPerpendicular(1);
                    velY0(m) = velParallel(2)+velPerpendicular(2);
                end
                if isCollisionWallEnd(m,2,w) == 1
                    slopeVect = [(X(i,m) - wallPts(2,1,w)) (Y(i,m) - wallPts(2,2,w))];
                    slopeVect = slopeVect/norm(slopeVect); % unit vector
                    velVect = [velX0(m) velY0(m)];
                    
                    % resolve velX and velY along line-of-collision coordinates
                    velParallel = dot(velVect,slopeVect)*slopeVect;
                    velParallelMag = dot(velVect,slopeVect);
                    velPerpendicular = velVect - velParallel;

                    % reverse parallel velocity
                    velParallel = -velParallel;

                    % check if new velocity moves particle away from wall
                    if dot(velParallel,slopeVect) < 0
                        velParallel = -velParallel;
                    end

                    % resolve velocities in velX and velY directions
                    velX0(m) = velParallel(1)+velPerpendicular(1);
                    velY0(m) = velParallel(2)+velPerpendicular(2);
                end

                if isCollisionWallEdge(w,m) == 1
                    velVect = [velX0(m) velY0(m)];
                    wallVect = [wallPts(2,1,w)-wallPts(1,1,w) wallPts(2,2,w)-wallPts(1,2,w)];
                    wallVect = wallVect/norm(wallVect); % unit vector
                    
                    % reflect the velocity from wall
                    velReflected = 2*dot(velVect,wallVect)*wallVect - velVect;

                    % check if jitter is occuring
                    if isWallJitter(w,m) == 1
                        velReflected = -velReflected;
                    end

                    % resolve velocities in velX and velY directions
                    velX0(m) = velReflected(1);
                    velY0(m) = velReflected(2);
                end
            end
        end
    end
    
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
figure();
plot(pressure)
hold on
plot(movmean(pressure,movingAverageWindow));
title("Pressure exerted on the bottom box wall")
xlabel("Time");
ylabel("Pressure");
legend("Localised instantaneous pressure","Moving wall-averaged pressure");
grid on;

%% Animate particles and make video

if playAnimationOnly == 1 || makeVideo == 1
    
    fig = figure();
    fig.WindowState = 'maximized';
    drawnow; % Force MATLAB to update the figure size
    
    % Set axis properties first
    axis equal;
    axis([-box_xmax/2 - 10, box_xmax/2 + 10, -box_ymax/2 - 10, box_ymax/2 + 10]);
    xticks(-box_xmax/2:10:box_xmax/2);
    yticks(-box_ymax/2:10:box_ymax/2);
    grid on;
    hold on;
    
    % Plot walls (if any)
    if innerPartitions == 1
        for w = 1:noOfWalls
            line(wallPts(:,1,w), wallPts(:,2,w), 'Color', 'black','LineWidth',1.5);
        end
    end
    rectangle('Position',[-box_xmax/2,-box_ymax/2,box_xmax,box_ymax]);
    
    % SCALING CALCULATION
    ax = gca;
    drawnow; % Force the figure to render to get accurate measurements
    
    % Get axis position in pixels (more reliable than points)
    set(ax, 'Units', 'pixels');
    axPos = get(ax, 'Position');
    set(ax, 'Units', 'normalized'); % Reset to default
    
    % Get axis data limits
    xLimits = xlim(ax);
    yLimits = ylim(ax);
    
    % Calculate pixels per data unit
    pixelsPerDataUnit_X = axPos(3) / diff(xLimits); % pixels per x-unit
    pixelsPerDataUnit_Y = axPos(4) / diff(yLimits); % pixels per y-unit
    
    % Use the smaller scaling to maintain circular particles in axis equal
    pixelsPerDataUnit = min(pixelsPerDataUnit_X, pixelsPerDataUnit_Y);
    
    % Convert pixels to points (72 points = 1 inch, typical screen = 96 DPI)
    screenDPI = get(0, 'ScreenPixelsPerInch'); % Get screen DPI
    pointsPerPixel = 72 / screenDPI;
    
    % Calculate marker sizes in points^2 (area units for scatter)
    diameter_A_pixels = size_A * pixelsPerDataUnit;
    diameter_B_pixels = size_B * pixelsPerDataUnit;
    
    diameter_A_points = diameter_A_pixels * pointsPerPixel;
    diameter_B_points = diameter_B_pixels * pointsPerPixel;
    
    % Convert diameter to area (scatter uses area, not diameter)
    markerSize_A = (diameter_A_points)^2;
    markerSize_B = (diameter_B_points)^2;

    % Create scatter plots (not regular plot)
    s_A = scatter(X(1,1:nPart_A), Y(1,1:nPart_A), markerSize_A, 'blue', 'filled', 'MarkerEdgeColor', 'blue');
    
    if twoParticleMix == 1
        s_B = scatter(X(1,nPart_A+1:end), Y(1,nPart_A+1:end), markerSize_B, 'red', 'filled', 'MarkerEdgeColor', 'red');
    end
    
    title("Press any key to start");
    pause;
    
    if makeVideo == 1
        v = VideoWriter('gas');
        v.Quality = 100;
        open(v);
    end

    for i = 2:time
        % w = waitforbuttonpress;
        % if w == 1
        %     % Update scatter positions (not markersize)
        %     set(s_A, 'XData', X(i,1:nPart_A), 'YData', Y(i,1:nPart_A));
        %     if twoParticleMix == 1
        %         set(s_B, 'XData', X(i,nPart_A+1:end), 'YData', Y(i,nPart_A+1:end));
        %     end
        %     drawnow;
        % end
    
        % Update scatter positions (not markersize)
        set(s_A, 'XData', X(i,1:nPart_A), 'YData', Y(i,1:nPart_A));
        if twoParticleMix == 1
            set(s_B, 'XData', X(i,nPart_A+1:end), 'YData', Y(i,nPart_A+1:end));
        end
        drawnow;
        
        if makeVideo == 1
            frame = getframe(gcf);
            writeVideo(v , frame);
        end
    end
    
    if makeVideo == 1
        close(v);
    end
end