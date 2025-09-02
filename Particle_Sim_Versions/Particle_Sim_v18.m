clc; clear; close all

%% Choose type of simulation and simulation settings
twoParticleMix = 1;
innerPartitions = 0;
heatTransferBottom = 0;
reactiveGas = 1;
windTunnel = 0;
gravity = 0;
playAnimationOnly = 1;
makeVideo = 0;
selfReact_A = 1;
selfReact_B = 1;
probeRegion = 0;

%% Input simulation parameter values - by typing
% Particle A
massPart_A = 2; % mass of particles [amu OR gm/mol]
nPart_A = 50; % no. of particles
size_A = 2; % particle diameter [angstrom]
speed_A = 0.5; % initial speed [angstrom/femtosecond]

% Particle B
massPart_B = 1; % [amu OR gm/mol]
nPart_B = 100;
size_B = 1; % [angstrom]
speed_B = 0.5; % [angstrom/femtosecond]

% Environment and post-processing settings
g = 9.8; % acceleration due to gravity (delta_V/timeStep) [m/s2]
heatIn = 0.000005; % in units of energy/time [W]
criticalEnergy = 0; % Minimum total K.E. of two atoms along line-of-centre for reaction to happen [kcal/mol]
box_xmax = 100; % max box width [angstrom]
box_ymax = 100; % max box breadth [angstrom]
timeTotal = 1000; % [femtoseconds]
timeStep = 1; % dicretized time [femtoseconds]
simSpeedUp = 1; % simulation speed up [>= 1, integer]
movingAverageWindow = 100; % moving average window size for calculation of physical parameters

%% Input simulation parameter values - interactively

% Initial particle region
fig = figure();
fig.WindowState = 'maximized';
rectangle('Position',[-box_xmax/2,-box_ymax/2,box_xmax,box_ymax],'LineWidth',1.5);
axis equal;
axis([-box_xmax/2 - 10, box_xmax/2 + 10, -box_ymax/2 - 10, box_ymax/2 + 10]);
xticks(-box_xmax/2:10:box_xmax/2);
yticks(-box_ymax/2:10:box_ymax/2);
grid on;
title('Select initial region of blue (and red) particles');
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

if twoParticleMix == 0
    nPart_B = 0;
end

itr = floor(timeTotal/timeStep);
size = [size_A size_B];
massPart = [massPart_A massPart_B];
massPartArray = [massPart_A*ones(1,nPart_A) massPart_B*ones(1,nPart_B)];
initX_A = initLoc_A(3)*rand(1,nPart_A) + initLoc_A(1); % initial X-coordinates
initY_A = initLoc_A(4)*rand(1,nPart_A) + initLoc_A(2); % initial Y-coordinates
initX_B = initLoc_B(3)*rand(1,nPart_B) + initLoc_B(1);
initY_B = initLoc_B(4)*rand(1,nPart_B) + initLoc_B(2);
initX = [initX_A initX_B];
initY = [initY_A initY_B];
initAngle = 360*rand(1,nPart_A+nPart_B) + 0; % initial angle w.r.t x-axis
initSpeed = [speed_A*ones(1,nPart_A) speed_B*ones(1,nPart_B)];
isTopWallCollision = zeros(itr,length(initX));
if innerPartitions == 1
    isCollisionWallEnd = zeros(nPart_A+nPart_B,2,noOfWalls); % wall end collision detection matrix
    isCollisionWallEdge = zeros(noOfWalls,nPart_A+nPart_B); % wall edge collision detection matrix
end
isMolecule = zeros(1,nPart_A+nPart_B);
molCtr = 0; % molecule counter
speciesHistory = zeros(itr,5); % B, R, BB, RR, BR
speciesHistory(:,1) = nPart_A;
speciesHistory(:,2) = nPart_B;
cellSize_x = box_xmax/floor(box_xmax/max(size)); % Hash grid cell size
cellSize_y = box_ymax/floor(box_ymax/max(size));
cellNos_x = floor(box_xmax/max(size)); % No. of cells
cellNos_y = floor(box_ymax/max(size));
cellSearchRange = reactiveGas + 1;

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

xmax = box_xmax/2 - [size_A/2*ones(1,nPart_A) size_B/2*ones(1,nPart_B)];
xmin = -xmax;
ymax = box_ymax/2 - [size_A/2*ones(1,nPart_A) size_B/2*ones(1,nPart_B)];
ymin = -ymax;

% calculate particles' coordinates for all times
for i = 2:itr
    % move particle
    X(i,:) = X(i-1,:) + velX0*timeStep;
    Y(i,:) = Y(i-1,:) + velY0*timeStep;

    % collision with box wall
    velX0 = velX0.*((((X(i,:) > xmax) | (X(i,:) < xmin)) - 0.5)/(-0.5));
    velY0 = velY0.*((((Y(i,:) > ymax) | (Y(i,:) < ymin)) - 0.5)/(-0.5));
    X(i,:) = X(i,:) - 2*((X(i,:) > xmax).*(X(i,:) - xmax));
    Y(i,:) = Y(i,:) - 2*((Y(i,:) > ymax).*(Y(i,:) - ymax));
    X(i,:) = X(i,:) + 2*((X(i,:) < xmin).*(xmin - X(i,:)));
    Y(i,:) = Y(i,:) + 2*((Y(i,:) < ymin).*(ymin - Y(i,:)));
    isTopWallCollision(i,:) = velY(i-1,:)./velY0 == -1 & velY0 < 0;

     % add heat energy to particles collided with bottom wall
    if heatTransferBottom == 1
        isBottomWallCollision = velY(i-1,:)./velY0 == -1 & velY(i-1,:) < 0;
        velY0 = isBottomWallCollision .* sqrt(velY0.^2 + 2*(heatIn*timeStep*10^-15)*(1000*(6.023*10^23)*(10^-10))./massPartArray) + ~isBottomWallCollision.*velY0;
    end

    % check inside wall edge collision
    if innerPartitions == 1
        for w = 1:noOfWalls
            in = inpolygon(X(i,1:nPart_A),Y(i,1:nPart_A),wallBox_A(:,1,w),wallBox_A(:,2,w));
            isCollisionWallEdge(w,1:nPart_A) = in;
            if twoParticleMix == 1
                in = inpolygon(X(i,nPart_A+1:end),Y(i,nPart_A+1:end),wallBox_B(:,1,w),wallBox_B(:,2,w));
                isCollisionWallEdge(w,nPart_A+1:end) = in;
            end
        end
    end

    % reset grid
    lookup = zeros(cellNos_x,cellNos_y);
    lookupCtr = ones(cellNos_x,cellNos_y);

    % assign particles to grid
    for p = 1:(nPart_A+nPart_B)
        cellX = ceil((box_xmax/2 + X(i,p))/cellSize_x); 
        cellY = ceil((box_ymax/2 + Y(i,p))/cellSize_y);
        lookup(cellX,cellY,lookupCtr(cellX,cellY)) = p;
        lookupCtr(cellX,cellY) = lookupCtr(cellX,cellY) + 1;
    end

    % check neighbouring cells for each particle and register collision
    for m = 1:(nPart_A+nPart_B)
        if innerPartitions == 1
            for w = 1:noOfWalls
                distWallEnd1 = sqrt((wallPts(1,1,w) - X(i,m))^2 + (wallPts(1,2,w) - Y(i,m))^2);
                distWallEnd2 = sqrt((wallPts(2,1,w) - X(i,m))^2 + (wallPts(2,2,w) - Y(i,m))^2);
                isCollisionWallEnd(m,:,w) = [distWallEnd1 distWallEnd2] < (size((m>nPart_A)+1)/2)*ones(1,2);

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
                    d = size((m>nPart_A)+1)/2 - sqrt((X(i,m) - wallPts(1,1,w))^2 + (Y(i,m) - wallPts(1,2,w))^2);
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
                    d = size((m>nPart_A)+1)/2 - sqrt((X(i,m) - wallPts(2,1,w))^2 + (Y(i,m) - wallPts(2,2,w))^2);
                    X(i,m) = X(i,m) + 2*slopeVect(1)*d;
                    Y(i,m) = Y(i,m) + 2*slopeVect(2)*d;
                end

                if isCollisionWallEdge(w,m) == 1
                    velVect = [velX0(m) velY0(m)];
                    wallVect = [wallPts(2,1,w)-wallPts(1,1,w) wallPts(2,2,w)-wallPts(1,2,w)];
                    wallVect = wallVect/norm(wallVect); % unit vector
                    
                    % reflect the velocity from wall
                    velReflected = 2*dot(velVect,wallVect)*wallVect - velVect;

                    % resolve velocities in velX and velY directions
                    velX0(m) = velReflected(1);
                    velY0(m) = velReflected(2);
                    
                    % reflect X and Y coordinates
                    wallEndToPartVect = [(X(i,m) - wallPts(1,1,w)) (Y(i,m) - wallPts(1,2,w))];
                    wallEdgeToPartVect = wallEndToPartVect - dot(wallEndToPartVect,wallVect)*wallVect;
                    wallEdgeToPart = norm(wallEdgeToPartVect);
                    normalVect = wallEdgeToPartVect/wallEdgeToPart;
                    d = size((m>nPart_A)+1)/2 - wallEdgeToPart;
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
                                isCollisionPart = distPart < size((m>nPart_A)+1)/2 + size((n>nPart_A)+1)/2;
                                if (isCollisionPart == 1) || (isCollisionPart == 0 && reactiveGas == 1 && isMolecule(m) == isMolecule(n) && isMolecule(m) ~= 0)
                                    slopeVect = [(X(i,m) - X(i,n)) (Y(i,m) - Y(i,n))];
                                    slopeVect = slopeVect/norm(slopeVect); % unit vector
                                    velVect_m = [velX0(m) velY0(m)];
                                    velVect_n = [velX0(n) velY0(n)];
                                    m1 = massPart((m>nPart_A)+1);
                                    m2 = massPart((n>nPart_A)+1);
                    
                                    % resolve velX and velY along line-of-collision coordinates
                                    velParallel = [dot(velVect_m,slopeVect)*slopeVect; dot(velVect_n,slopeVect)*slopeVect];
                                    velParallelMag = [dot(velVect_m,slopeVect); dot(velVect_n,slopeVect)];
                                    velPerpendicular = [velVect_m; velVect_n] - velParallel;
                    
                                    % check if reaction occured
                                    if reactiveGas == 1 && isMolecule(m) == 0 && isMolecule(n) == 0
                                        if 0.5*(m1*(norm(velParallel(1,:)))^2 + m2*(norm(velParallel(2,:)))^2) >= (criticalEnergy * (4184*1000*(10^-10)))
                                            molCtr = molCtr + 1;
                                            isMolecule(m) = molCtr; isMolecule(n) = molCtr;
                                            if (selfReact_A == 0 && m <= nPart_A && n <= nPart_A) || (selfReact_B == 0 && m > nPart_A && n > nPart_A)
                                                molCtr = molCtr - 1;
                                                isMolecule(m) = 0; isMolecule(n) = 0;
                                            end
                                            
                                            % Update species count
                                            if isMolecule(m) ~= 0 && isMolecule(n) ~= 0
                                                if m <= nPart_A && n <= nPart_A
                                                    speciesHistory(i:end,3) = speciesHistory(i:end,3) + 1;
                                                    speciesHistory(i:end,1) = speciesHistory(i:end,1) - 2;
                                                elseif m > nPart_A && n > nPart_A
                                                    speciesHistory(i:end,4) = speciesHistory(i:end,4) + 1;
                                                    speciesHistory(i:end,2) = speciesHistory(i:end,2) - 2;
                                                else
                                                    speciesHistory(i:end,5) = speciesHistory(i:end,5) + 1;
                                                    speciesHistory(i:end,1) = speciesHistory(i:end,1) - 1;
                                                    speciesHistory(i:end,2) = speciesHistory(i:end,2) - 1;
                                                end
                                            end
                                        end
                                    end
                    
                                    % exchange momentum
                                    oldVelParallel = velParallel([1 2],:);
                                    v1f = ((m1-m2)*velParallelMag(1) + 2*m2*velParallelMag(2))/(m1+m2);
                                    v2f = ((m2-m1)*velParallelMag(2) + 2*m1*velParallelMag(1))/(m1+m2);
                                    velParallel = [v1f*slopeVect; v2f*slopeVect];
                    
                                    % check if velocities separate/join particles
                                    if (isCollisionPart == 1) && (dot(slopeVect, velParallel(1,:) - velParallel(2,:)) < 0)
                                        velParallel([1 2],:) = oldVelParallel;
                                    elseif (isCollisionPart == 0 && reactiveGas == 1 && isMolecule(m) == isMolecule(n) && isMolecule(m) ~= 0) && (dot(slopeVect, velParallel(1,:) - velParallel(2,:)) > 0)
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
        velY0 = velY0 - (g*(10^10)/(10^30))*timeStep;
    end

    % update velocity array
    velX(i,:) = velX0;
    velY(i,:) = velY0;

    disp("Simulation Completion Status: " + i/itr*100 + "%");
end

%% Evaluate thermo-physical quantities

% calculate particle speed, pressure on top wall and region particle concentration/internal energy
speedArray = zeros(itr,length(initX));
pressure = zeros(itr,1);
if probeRegion == 1
    concentration = zeros(itr,2);
    energy = zeros(itr,2);
end
for i = 1:itr
    % speed
    speedArray(i,:) = sqrt(velX(i,:).^2 + velY(i,:).^2);

    % pressure on top wall
    if i>1
        pressure(i) = (massPart(1)*dot(velY(i-1,1:nPart_A),isTopWallCollision(i,1:nPart_A))*2)/timeStep/box_ymax * (0.001/(6.023*10^23)*10^30); % (m*dv/dt)/L [N/m]
        if twoParticleMix == 1
            pressure(i) = pressure(i) + (massPart(2)*dot(velY(i-1,nPart_A+1:end),isTopWallCollision(i,nPart_A+1:end))*2)/timeStep/box_ymax * (0.001/(6.023*10^23)*10^30);
        end
    end
    
    if probeRegion == 1
        % region concentration/internal energy
        in_A = inpolygon(X(i,1:nPart_A),Y(i,1:nPart_A),probeRectCoord_X,probeRectCoord_Y);
        concentration(i,1) = sum(in_A)/(probeRect(3)*probeRect(4)*10^-20); % particles/m2
        energy(i,1) = sum(0.5 * massPart_A * (speedArray(i,1:nPart_A).*in_A).^2 * 0.001/6.023/10^23 * 10^10); % J
        if twoParticleMix == 1
            in_B = inpolygon(X(i,nPart_A+1:end),Y(i,nPart_A+1:end),probeRectCoord_X,probeRectCoord_Y);
            concentration(i,2) = sum(in_B)/(probeRect(3)*probeRect(4)*10^-20); % particles/m2
            energy(i,2) = sum(0.5 * massPart_B * (speedArray(i,nPart_A+1:end).*in_B).^2 * 0.001/6.023/10^23 * 10^10); % J
        end
    end
end

%% Plot thermo-physical quantities
% figure();
% plot(pressure)
% hold on
% plot(movmean(pressure,movingAverageWindow));
% title("Pressure exerted on the top box wall")
% xlabel("Time (in fs)");
% ylabel("Pressure (in N/m)");
% legend("Localised instantaneous pressure","Moving wall-averaged pressure");
% grid on;

% figure()
% h = histogram(speedArray(1,:));
% title("Particle speed distribution")
% xlabel("Speed (in Å/fs)");
% ylabel("No. of particles");
% pause;
% for i = 2:itr
%     pause(0.05);
%     h.Data = speedArray(i,:);
% end

figure()
plot(speciesHistory);
legend("Blue monoatomic gas","Red monoatomic gas", "Blue diatomic gas", "Red diatomic gas", "Blue-red diatomic gas");
title("Chemical species time history");
xlabel("Time (in fs)");
ylabel("No. of species");
grid on;


if probeRegion == 1
    figure()
    plot(concentration(:,1));
    hold on;
    plot(concentration(:,2));
    legend("Blue particle","Red particle");
    title("Particle concentration")
    xlabel("Time (in fs)");
    ylabel("Concentration (in particles/m2)");
    grid on;
end

if probeRegion == 1
    figure()
    plot(energy(:,1));
    hold on;
    plot(energy(:,2));
    legend("Blue particle","Red particle");
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
    if heatTransferBottom == 1
        line([-box_xmax/2 box_xmax/2], [-box_ymax/2 -box_ymax/2], 'Color', 'yellow','LineWidth',1.5);
    end
    
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
    title("Animation started");

    if makeVideo == 1
        v = VideoWriter('gas');
        v.Quality = 100;
        open(v);
    end

    for i = 2:simSpeedUp:itr
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
    
    title("Animation ended");

    if makeVideo == 1
        close(v);
    end
end