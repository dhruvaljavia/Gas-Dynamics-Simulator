clc; clear; close all

%% Choose type of simulation and simulation settings
innerPartitions = 1;
gravity = 0;
playAnimationOnly = 1;
makeVideo = 0;
probeRegion = 0;

%% Input simulation parameter values - by typing
% Particle
massPart = 1; % mass of particles [amu OR gm/mol]
size = 1; % particle diameter [angstrom]
inlet_velocity = 0.5; % inlet fluid flow velocity [angstrom/femtosecond]

% Environment and post-processing settings
partExcess = 20; % Excess particle outside box at inlet
partSpace = 1.2; % particle spacing factor (space = factor*size)
g = 9.8; % acceleration due to gravity (delta_V/timeStep) [m/s2]
box_xmax = 400; % max box width [angstrom]
box_ymax = 200; % max box breadth [angstrom]
timeTotal = 1000; % [femtoseconds]
timeStep = 1; % dicretized time [femtoseconds]
simSpeedUp = 1; % simulation speed up [>= 1, integer]
movingAverageWindow = 100; % moving average window size for calculation of physical parameters

%% Input simulation parameter values - interactively

% Wall locations
fig = figure();
fig.WindowState = 'maximized';
rectangle('Position',[-box_xmax/2,-box_ymax/2,box_xmax,box_ymax],'LineWidth',1.5);
axis equal;
axis([-box_xmax/2 - 10, box_xmax/2 + 10, -box_ymax/2 - 10, box_ymax/2 + 10]);
xticks(-box_xmax/2:10:box_xmax/2);
yticks(-box_ymax/2:10:box_ymax/2);
grid on;
if innerPartitions == 1
    title('Select end points of walls');
    ctr=1;
    while 1
        wall = drawline(Color="black");
        wallPts(1:2,1:2,ctr) = wall.Position;
        pause;
        c = get(gcf,'CurrentCharacter');
        if c == 'x'
            break;
        end
        ctr=ctr+1;
    end
    noOfWalls = ctr;
end

% Probe region
if probeRegion == 1
    title('Select probe region for concentration and internal energy')
    rect = drawrectangle(Color="green");
    pause;
    probeRect = rect.Position;
    probeRectCoord_X = [probeRect(1) probeRect(1) (probeRect(1)+probeRect(3)) (probeRect(1)+probeRect(3))];
    probeRectCoord_Y = [probeRect(2) (probeRect(2)+probeRect(4)) (probeRect(2)+probeRect(4)) probeRect(2)];
end

close;

%% Initialise particles' state

itr = floor(timeTotal/timeStep);
nX = floor((box_xmax/size)/partSpace) + partExcess;
nY = floor((box_ymax/size)/partSpace);
nPart = nX*nY;
xloc = linspace(-box_xmax/2 + size/2 - (nX-1)*partSpace*size,-box_xmax/2 + size/2,nX);
yloc = linspace(-box_ymax/2 + size/2,box_ymax/2 - size/2,nY);
dY = yloc(2) - yloc(1);
ctrY = 1;
minX = min(xloc)-partSpace*size;
[initX, initY] = meshgrid(xloc,yloc);
initX = reshape(initX,1,[]);
initY = reshape(initY,1,[]);
initAngle = zeros(1,nPart); % initial angle w.r.t x-axis
initSpeed = inlet_velocity*ones(1,nPart); % initial speed in units/iteration
if innerPartitions == 1
    isCollisionWallEnd = zeros(nPart,2,noOfWalls); % wall end collision detection matrix
    isCollisionWallEdge = zeros(noOfWalls,nPart); % wall edge collision detection matrix
    isWallJitter = zeros(noOfWalls,nPart); % wall jitter detection matrix
end
cellSize_x = box_xmax/floor(box_xmax/size); % Hash grid cell size
cellSize_y = box_ymax/floor(box_ymax/size);
cellNos_x = floor(box_xmax/size); % No. of cells
cellNos_y = floor(box_ymax/size);
cellSearchRange = 1;

% Making the wall box
if innerPartitions == 1
    wallNormalVect = zeros(noOfWalls,2); % unit wall normal vector
    for w = 1:noOfWalls
        wallNormalVect(w,:) = [-(wallPts(2,2,w)-wallPts(1,2,w)), wallPts(2,1,w)-wallPts(1,1,w)];
        wallNormalVect(w,:) = wallNormalVect(w,:)/norm(wallNormalVect(w,:));
    end
    addVect = (size(1)/2)*wallNormalVect;
    wallBox = zeros(4,2,noOfWalls);
    for w = 1:noOfWalls
        wallBox(:,:,w) = [wallPts(1,:,w)+addVect(w,:); wallPts(1,:,w)-addVect(w,:); wallPts(2,:,w)-addVect(w,:); wallPts(2,:,w)+addVect(w,:)];
    end
end

%% Simulate particles
X = zeros(itr,length(initX));
Y = zeros(itr,length(initY));
X(1,:) = initX;
Y(1,:) = initY;

velX = zeros(itr,length(initX));
velY = zeros(itr,length(initY));
velX(1,:) = initSpeed.*cosd(initAngle);
velY(1,:) = initSpeed.*sind(initAngle);
velX0 = velX(1,:);
velY0 = velY(1,:);

xmax = box_xmax/2 - size/2;
xmin = -xmax;
ymax = box_ymax/2 - size/2;
ymin = -ymax;

% calculate particles' coordinates for all times
for i = 2:itr
    % move particle and particle append flag (if windTunnel is on)
    X(i,:) = X(i-1,:) + velX0*timeStep;
    Y(i,:) = Y(i-1,:) + velY0*timeStep;
    minX = minX + inlet_velocity*timeStep;

    % Indices for particles inside the box
    partI = find(X(i,:)>=xmin & X(i,:)<=xmax);

    % collision with box wall
    velY0 = velY0.*((((Y(i,:) > ymax) | (Y(i,:) < ymin)) - 0.5)/(-0.5));
    Y(i,:) = Y(i,:) - 2*((Y(i,:) > ymax).*(Y(i,:) - ymax));
    Y(i,:) = Y(i,:) + 2*((Y(i,:) < ymin).*(ymin - Y(i,:)));

    % check inside wall edge collision
    if innerPartitions == 1
        for w = 1:noOfWalls
            in = inpolygon(X(i,partI),Y(i,partI),wallBox(:,1,w),wallBox(:,2,w));
            isWallJitter(w,partI) = isCollisionWallEdge(w,partI) & in;
            isCollisionWallEdge(w,partI) = in;
        end
    end

    % reset grid
    lookup = zeros(cellNos_x,cellNos_y);
    lookupCtr = ones(cellNos_x,cellNos_y);

    % assign particles to grid
    for p = partI
        cellX = ceil((box_xmax/2 + X(i,p))/cellSize_x); 
        cellY = ceil((box_ymax/2 + Y(i,p))/cellSize_y);
        lookup(cellX,cellY,lookupCtr(cellX,cellY)) = p;
        lookupCtr(cellX,cellY) = lookupCtr(cellX,cellY) + 1;
    end

    % check neighbouring cells for each particle and register collision
    for m = partI
        if innerPartitions == 1
            for w = 1:noOfWalls
                distWallEnd1 = sqrt((wallPts(1,1,w) - X(i,m))^2 + (wallPts(1,2,w) - Y(i,m))^2);
                distWallEnd2 = sqrt((wallPts(2,1,w) - X(i,m))^2 + (wallPts(2,2,w) - Y(i,m))^2);
                isCollisionWallEnd(m,:,w) = [distWallEnd1 distWallEnd2] < (size/2)*ones(1,2);

                % remove overlapping wall detection for end and edge regions
                isCollisionWallEnd(m,1,w) = (~(isCollisionWallEnd(m,1,w) & isCollisionWallEdge(w,m)) & isCollisionWallEnd(m,1,w));
                isCollisionWallEnd(m,2,w) = (~(isCollisionWallEnd(m,2,w) & isCollisionWallEdge(w,m)) & isCollisionWallEnd(m,2,w));
                
                if isCollisionWallEnd(m,1,w) == 1
                    slopeVect = [(X(i,m) - wallPts(1,1,w)) (Y(i,m) - wallPts(1,2,w))];
                    slopeVect = slopeVect/norm(slopeVect); % unit vector
                    velVect = [velX0(m) velY0(m)];
                    
                    % resolve velX and velY along line-of-collision coordinates
                    velParallel = dot(velVect,slopeVect)*slopeVect;
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

                    % reflect X and Y coordinates
                    d = size/2 - sqrt((X(i,m) - wallPts(1,1,w))^2 + (Y(i,m) - wallPts(1,2,w))^2);
                    X(i,m) = X(i,m) + 2*slopeVect(1)*d;
                    Y(i,m) = Y(i,m) + 2*slopeVect(2)*d;
                end
                if isCollisionWallEnd(m,2,w) == 1
                    slopeVect = [(X(i,m) - wallPts(2,1,w)) (Y(i,m) - wallPts(2,2,w))];
                    slopeVect = slopeVect/norm(slopeVect); % unit vector
                    velVect = [velX0(m) velY0(m)];
                    
                    % resolve velX and velY along line-of-collision coordinates
                    velParallel = dot(velVect,slopeVect)*slopeVect;
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

                    % reflect X and Y coordinates
                    d = size/2 - sqrt((X(i,m) - wallPts(2,1,w))^2 + (Y(i,m) - wallPts(2,2,w))^2);
                    X(i,m) = X(i,m) + 2*slopeVect(1)*d;
                    Y(i,m) = Y(i,m) + 2*slopeVect(2)*d;
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

                    % reflect X and Y coordinates
                    wallEndToPartVect = [(X(i,m) - wallPts(1,1,w)) (Y(i,m) - wallPts(1,2,w))];
                    wallEdgeToPartVect = wallEndToPartVect - dot(wallEndToPartVect,wallVect)*wallVect;
                    wallEdgeToPart = norm(wallEdgeToPartVect);
                    normalVect = wallEdgeToPartVect/wallEdgeToPart;
                    d = size/2 - wallEdgeToPart;
                    X(i,m) = X(i,m) + 2*normalVect(1)*d;
                    Y(i,m) = Y(i,m) + 2*normalVect(2)*d;
                end
            end
        end

        cellX = ceil((box_xmax/2 + X(i,m))/cellSize_x); 
        cellY = ceil((box_ymax/2 + Y(i,m))/cellSize_y);

        for x = (cellX-cellSearchRange):(cellX+cellSearchRange)
            for y = (cellY-cellSearchRange):(cellY+cellSearchRange)
                if x >= 1 && x <= cellNos_x && y >= 1 && y <= cellNos_y
                    if lookup(x,y,1) ~= 0
                        for n = lookup(x,y,lookupCtr(x,y)-1)
                            if m~=n
                                distPart = sqrt((X(i,m) - X(i,n))^2 + (Y(i,m) - Y(i,n))^2);
                                isCollisionPart = distPart < size;
                                if (isCollisionPart == 1)
                                    slopeVect = [(X(i,m) - X(i,n)) (Y(i,m) - Y(i,n))];
                                    slopeVect = slopeVect/norm(slopeVect); % unit vector
                                    velVect_m = [velX0(m) velY0(m)];
                                    velVect_n = [velX0(n) velY0(n)];
                    
                                    % resolve velX and velY along line-of-collision coordinates
                                    velParallel = [dot(velVect_m,slopeVect)*slopeVect; dot(velVect_n,slopeVect)*slopeVect];
                                    velPerpendicular = [velVect_m; velVect_n] - velParallel;
                    
                                    % exchange momentum
                                    oldVelParallel = velParallel([1 2],:);
                                    velParallel([1 2],:) = velParallel([2 1],:);
                    
                                    % check if velocities separate/join particles
                                    if (isCollisionPart == 1) && (dot(slopeVect, velParallel(1,:) - velParallel(2,:)) < 0)
                                        velParallel([1 2],:) = oldVelParallel;
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
        end
    end

    % add -Y velocity increase due to gravity
    if gravity == 1
        velY0(partI) = velY0(partI) - (g*(10^10)/(10^30))*timeStep;
    end

    % Detect exited particles and transfer them at inlet
    partO = find(X(i,:)>xmax);
    for p = partO
        if ctrY > nY
            minX = minX - partSpace*size;
            ctrY = 1;
        end
        X(i,p) = minX;
        Y(i,p) = ymin + (ctrY-1)*dY;
        velX0(p) = inlet_velocity;
        velY0(p) = 0;
        ctrY = ctrY + 1;
    end

    % update velocity array
    velX(i,:) = velX0;
    velY(i,:) = velY0;

    disp("Simulation Completion Status: " + i/itr*100 + "%");
end

%% Evaluate thermo-physical quantities

% calculate particle speed, pressure on top wall and region concentration/internal energy
speedArray = zeros(itr,length(initX));
if probeRegion == 1
    concentration = zeros(itr);
    energy = zeros(itr);
end
for i = 1:itr
    % speed
    speedArray(i,:) = sqrt(velX(i,:).^2 + velY(i,:).^2);

    if probeRegion == 1
        % region concentration/internal energy
        partI = find(X(i,:)>=xmin & X(i,:)<=xmax);
        in = inpolygon(X(i,partI),Y(i,partI),probeRectCoord_X,probeRectCoord_Y);
        concentration(i) = sum(in)/(probeRect(3)*probeRect(4)*10^-20); % particles/m2
        energy(i) = sum(0.5 * massPart * (speedArray(i,partI).*in).^2 * 0.001/6.023/10^23 * 10^10); % J
    end
end

%% Plot thermo-physical quantities
if probeRegion == 1
    figure()
    plot(concentration);
    title("Particle concentration")
    xlabel("Time (in fs)");
    ylabel("Concentration (in particles/m2)");
    grid on;
end

if probeRegion == 1
    figure()
    plot(energy);
    title("Internal energy")
    xlabel("Time (in fs)");
    ylabel("Energy (in J)");
    grid on;
end

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
    rectangle('Position',[-box_xmax/2,-box_ymax/2,box_xmax,box_ymax],'LineWidth',1.5);
    
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
    diameter_pixels = size * pixelsPerDataUnit;
    
    diameter_points = diameter_pixels * pointsPerPixel;
    
    % Convert diameter to area (scatter uses area, not diameter)
    markerSize = (diameter_points)^2;

    % Create scatter plots (not regular plot)
    partI = find(X(1,:)>=xmin & X(1,:)<=xmax);
    s = scatter(X(1,partI), Y(1,partI), markerSize, 'blue', 'filled', 'MarkerEdgeColor', 'blue');
    
    title("Press any key to start");
    pause;
    title("Animation started");

    if makeVideo == 1
        v = VideoWriter('gas');
        v.Quality = 100;
        open(v);
    end

    for i = 2:simSpeedUp:itr
        partI = find(X(i,:)>=xmin & X(i,:)<=xmax);
        
        % w = waitforbuttonpress;
        % if w == 1
        %     % Update scatter positions (not markersize)
        %     set(s, 'XData', X(i,partI), 'YData', Y(i,partI));
        %     drawnow;
        % end
    
        % Update scatter positions (not markersize)
        set(s, 'XData', X(i,partI), 'YData', Y(i,partI));
        drawnow;

        if makeVideo == 1
            frame = getframe(gcf);
            writeVideo(v , frame);
        end
    end
    
    title("Animation ended");

    if makeVideo == 1
        close(v);
    end
end