clc;
clear;

numExecutions = 8; % Number of executions
numMaps = 7; % Number of maps generated each time
minDistance = 50; % Minimum distance

result = zeros(10, 50);
for executionIndex = 1:numExecutions
    rng('shuffle');
    numObstacles = 10;
    goalCount = 5;
    
    data = struct();
   
    [startEndPoints, omappings, data] = generate_problems(numMaps, numObstacles, goalCount, minDistance, data);

    for mapIndex = 1:numMaps
        omap = omappings{mapIndex};
        mapStruct = data.maps(mapIndex);
        startPose = mapStruct.start;
        goalPoints = mapStruct.goals;
        
        try
            [dscore,k] = calculate_dscore(startPose, goalPoints, omap);
            result(1, numMaps*(executionIndex-1)+mapIndex) = dscore;
            result(6, numMaps*(executionIndex-1)+mapIndex) = k;
            
        catch
           result(6, numMaps*(executionIndex-1)+mapIndex) = 1;
        end
        
        try
            [d,k] = calculateDQuantities(startPose, goalPoints, omap);
            result(2, numMaps*(executionIndex-1)+mapIndex) = d;
            result(7, numMaps*(executionIndex-1)+mapIndex) = k;
            
        catch
           result(7, numMaps*(executionIndex-1)+mapIndex) = 1;

        end
        
        try
            [dd,k] = cal1(startPose, goalPoints, omap);
            result(3, numMaps*(executionIndex-1)+mapIndex) = dd;
            result(8, numMaps*(executionIndex-1)+mapIndex) = k;
            
        catch
            result(8, numMaps*(executionIndex-1)+mapIndex) = 1;

        end
        
        try
            [ddd,k] = cal2(startPose, goalPoints, omap);
            result(4, numMaps*(executionIndex-1)+mapIndex) = ddd;
            result(9, numMaps*(executionIndex-1)+mapIndex) = k;
           
        catch
            result(9, numMaps*(executionIndex-1)+mapIndex) = 1;

        end
        
        try
            [dddd,k] = cal3(startPose, goalPoints, omap);
            result(5, numMaps*(executionIndex-1)+mapIndex) = dddd;
            result(10, numMaps*(executionIndex-1)+mapIndex) = k;
            
        catch
           result(10, numMaps*(executionIndex-1)+mapIndex) = 1;

        end
    end
    fprintf('第%d次了\n', numMaps*(executionIndex-1)+mapIndex);
    overallAvgPathLength = mean(avgPathLengths);
    fprintf('Overall Average Path Length: %.2f meters\n', overallAvgPathLength);
end

function [dquantity,k] = calculate_dscore(startPose, goalPoints, omap)
    % Consider unknown spaces to be unoccupied
   
    firstPoint = goalPoints(1, 1:3);

    % Calculate distances between each point and the first point
    distances = sqrt(sum((goalPoints(:, 1:3) - firstPoint).^2, 2));
    
    % Set the distance of the first point to infinity to exclude itself
    distances(1) = inf;
    
    % Find the index of the closest point
    [~, closestIndex] = min(distances);
    % Output the closest point to the first point
    decoyPose3 = goalPoints(closestIndex, :);
    

    % Define UAV state space
    ss = ExampleHelperUAVStateSpace("MaxRollAngle", pi/6, ...
                                    "AirSpeed", 6, ...
                                    "FlightPathAngleLimit", [-0.1 0.1], ...
                                    "Bounds", [-20 520; -20 520; 10 100; -pi pi]);

    % Create validator for the map
    sv = validatorOccupancyMap3D(ss, "Map", omap);
    sv.ValidationDistance = 0.1;

    % Create planners for the paths
    planner1 = plannerRRTStar(ss, sv);
    planner1.MaxConnectionDistance = 50;
    planner1.GoalBias = 0.10;
    planner1.MaxIterations = 400;

    planner2 = plannerRRTStar(ss, sv);
    planner2.MaxConnectionDistance = 50;
    planner2.GoalBias = 0.10;
    planner2.MaxIterations = 400;

    planner3 = plannerRRTStar(ss, sv);
    planner3.MaxConnectionDistance = 50;
    planner3.GoalBias = 0.10;
    planner3.MaxIterations = 400;

    % Calculate the new point using the given formula
    a = norm(startPose(1:3) - goalPoints(1, 1:3));
    b = norm(startPose(1:3) - decoyPose3(1:3));
    c = norm(decoyPose3(1:3) - goalPoints(1, 1:3));

    d = (c + a - b) / 2;
    directionVec = (goalPoints(1, 1:3) - decoyPose3(1:3)) / c;
    newPoint = goalPoints(1, 1:3) - d * directionVec;
    newPoint = [newPoint pi/2];

    % Plan paths
    rng(50, "twister");
    [pthObj1, solnInfo1] = plan(planner1, startPose, goalPoints(1, :));

    rng(50, "twister");
    [pthObj2, solnInfo2] = plan(planner2, startPose, newPoint);
    
    rng(50, "twister");
    [pthObj3, solnInfo3] = plan(planner3, newPoint, goalPoints(1, :));

    if (solnInfo1.IsPathFound && solnInfo2.IsPathFound && solnInfo3.IsPathFound)
        % Calculate path lengths
        pathLen1 = pathLength(pthObj1);
        pathLen2 = pathLength(pthObj2);
        pathLen3 = pathLength(pthObj3);
        
        % Compute overlap lengths
        overlapThreshold = 1; % Define a distance threshold
        overlapLength = 0;

        % Calculate overlap between path1 and path2
        [xENU1, yENU1, zENU1] = exampleHelperSimulateUAV(pthObj1.States, ss.AirSpeed, pathLen1 / ss.AirSpeed);
        [xENU2, yENU2, zENU2] = exampleHelperSimulateUAV(pthObj2.States, ss.AirSpeed, pathLen2 / ss.AirSpeed);

        for i = 1:length(xENU1) - 1
            segmentStart1 = [xENU1(i), yENU1(i), zENU1(i)];
            segmentEnd1 = [xENU1(i+1), yENU1(i+1), zENU1(i+1)];

            for j = 1:length(xENU2) - 1
                segmentStart2 = [xENU2(j), yENU2(j), zENU2(j)];
                segmentEnd2 = [xENU2(j+1), yENU2(j+1), zENU2(j+1)];

                if minDistBetweenSegments(segmentStart1, segmentEnd1, segmentStart2, segmentEnd2) < overlapThreshold
                    overlapLength = overlapLength + norm(segmentEnd2 - segmentStart2);
                    break;
                end
            end
        end

        % Calculate overlap between path1 and path3
        [xENU3, yENU3, zENU3] = exampleHelperSimulateUAV(pthObj3.States, ss.AirSpeed, pathLen3 / ss.AirSpeed);

        for i = 1:length(xENU1) - 1
            segmentStart1 = [xENU1(i), yENU1(i), zENU1(i)];
            segmentEnd1 = [xENU1(i+1), yENU1(i+1), zENU1(i+1)];

            for j = 1:length(xENU3) - 1
                segmentStart2 = [xENU3(j), yENU3(j), zENU3(j)];
                segmentEnd2 = [xENU3(j+1), yENU3(j+1), zENU3(j+1)];

                if minDistBetweenSegments(segmentStart1, segmentEnd1, segmentStart2, segmentEnd2) < overlapThreshold
                    overlapLength = overlapLength + norm(segmentEnd2 - segmentStart2);
                    break;
                end
            end
        end
 
        % Compute dscore
        dquantity = max((pathLen1 - overlapLength) / pathLen1,0);
        %dcost = (pathLen2 + pathLen3 - pathLen1) / pathLen1;
        %dscore = dquantity / dcost;
        fprintf('dissimulation dquantity is: %.2f\n', dquantity);
        maxProbGoals = run_rrt_second([], goalPoints, planner2, pthObj2);
        n = length(goalPoints) - 3;
        k=mode(maxProbGoals(1:n));
    fprintf(' Max Probability Goal Index = %d\n', k);
   
    end   
end
function [startEndPoints, omappings, data] = generate_problems(numMaps, numObstacles, goalCount, minDistance, data)
    omappings = cell(numMaps, 1);
    startEndPoints = zeros(goalCount + 1, 4, numMaps);

    mapLength = 500;
    mapWidth = 500;

    for mapIndex = 1:numMaps
        % Create 3D occupancy map
        omap3D = occupancyMap3D;
        
        % Generate obstacles
        obstacleNumber = 1;
        while obstacleNumber <= numObstacles
            width = randi([1 50],1);
            length = randi([1 50],1);
            height = randi([1 150],1);
            xPosition = randi([0 mapWidth-width],1);
            yPosition = randi([0 mapLength-length],1);

            [xObstacle,yObstacle,zObstacle] = meshgrid(xPosition:xPosition+width,yPosition:yPosition+length,0:height);
            xyzObstacles = [xObstacle(:) yObstacle(:) zObstacle(:)];

            checkIntersection = false;
            for i = 1:size(xyzObstacles,1)
                if checkOccupancy(omap3D,xyzObstacles(i,:)) == 1
                    checkIntersection = true;
                    break
                end
            end
            if checkIntersection
                continue
            end

            setOccupancy(omap3D,xyzObstacles,1)
            
            obstacleNumber = obstacleNumber + 1;
        end

        % Set ground as an obstacle
        [xGround,yGround,zGround] = meshgrid(0:mapWidth,0:mapLength,0);
        xyzGround = [xGround(:) yGround(:) zGround(:)];
        setOccupancy(omap3D,xyzGround,1)

        % Set unknown spaces to be unoccupied
        omap3D.FreeThreshold = omap3D.OccupiedThreshold;

        % Generate start point
        while true
            startX = randi([1, mapWidth]);
            startY = randi([1, mapLength]);
            startZ = randi([1, 50]);
            startTheta = rand() * 2 * pi;

            % Ensure start point is in free space
            if ~checkOccupancy(omap3D, [startX, startY, startZ])
                break;
            end
        end
        
        startEndPoints(1,:,mapIndex) = [startX, startY, startZ, startTheta];

        mapStruct = struct();
        mapStruct.start = [startX, startY, startZ, startTheta];
        mapStruct.goals = [];
        
        % Generate the goal points
        for goal = 1:goalCount
            while true
                endX = randi([1, mapWidth]);
                endY = randi([1, mapLength]);
                endZ = randi([1, 50]);
                endTheta = rand() * 2 * pi;

                % Calculate distance between start and end points
                distanceToStart = sqrt((endX - startX)^2 + (endY - startY)^2 + (endZ - startZ)^2);

                % Calculate distances to all previously placed goals
                distancesToGoals = arrayfun(@(g) sqrt((endX - startEndPoints(g+1, 1, mapIndex))^2 + ...
                                                      (endY - startEndPoints(g+1, 2, mapIndex))^2 + ...
                                                      (endZ - startEndPoints(g+1, 3, mapIndex))^2), ...
                                            1:goal-1);

                % Check if the distance is at least minDistance and the point is in free space
                if distanceToStart >= minDistance && all(distancesToGoals >= minDistance) && ~checkOccupancy(omap3D, [endX, endY, endZ])
                    break;
                end
            end

            startEndPoints(goal + 1, :, mapIndex) = [endX, endY, endZ, endTheta];
            mapStruct.goals = [mapStruct.goals; [endX, endY, endZ, endTheta]];
            
        end

        omappings{mapIndex} = omap3D;
        data.maps(mapIndex) = mapStruct;
    end
end
function [dquantity1,k] = calculateDQuantities(startPose, goalPoints, omap)
    firstPoint = goalPoints(1, 1:3);

    % Calculate distances between each point and the first point
    distances = sqrt(sum((goalPoints(:, 1:3) - firstPoint).^2, 2));
    
    % Set the distance of the first point to infinity to exclude itself
    distances(1) = inf;
    
    % Find the index of the closest point
    [~, closestIndex] = min(distances);
    % Output the closest point to the first point
    newPoint = goalPoints(closestIndex, :);

    % Define UAV state space
    ss = ExampleHelperUAVStateSpace("MaxRollAngle", pi/6, ...
                                "AirSpeed", 6, ...
                                "FlightPathAngleLimit", [-0.1 0.1], ...
                                "Bounds", [-20 520; -20 520; 10 100; -pi pi]);
% Create validator for the map
sv = validatorOccupancyMap3D(ss, "Map", omap);
sv.ValidationDistance = 0.1;

% Create planner for the first path (startPose to goalPose)
planner1 = plannerRRT(ss, sv);
planner1.MaxConnectionDistance = 50;
planner1.GoalBias = 0.10;
planner1.MaxIterations = 400;


% Plan path from startPose to goalPose



% Create planner for the second path (startPose to decoyPose3 to goalPose)
planner2 = plannerRRT(ss, sv);
planner2.MaxConnectionDistance = 50;
planner2.GoalBias = 0.10;
planner2.MaxIterations = 400;

% Create planner for the third path (decoyPose3 to goalPose)
planner3 = plannerRRT(ss, sv);
planner3.MaxConnectionDistance = 50;
planner3.GoalBias = 0.10;
planner3.MaxIterations = 400;

% Plan paths
    rng(50, "twister");
    [pthObj1, solnInfo1] = plan(planner1, startPose, goalPoints(1, :));

    rng(50, "twister");
    [pthObj2, solnInfo2] = plan(planner2, startPose, newPoint);

    rng(50, "twister");
    [pthObj3, solnInfo3] = plan(planner3, newPoint, goalPoints(1, :));

if (solnInfo1.IsPathFound && solnInfo2.IsPathFound && solnInfo3.IsPathFound)
    % Calculate path lengths
    pathLen1 = pathLength(pthObj1);
    pathLen2 = pathLength(pthObj2);
    pathLen3 = pathLength(pthObj3);

    % Simulate UAV trajectory based on fixed-wing guidance model
    timeToReachGoal1 = pathLen1 / ss.AirSpeed;
    waypoints1 = pthObj1.States;
    [xENU1, yENU1, zENU1] = exampleHelperSimulateUAV(waypoints1, ss.AirSpeed, timeToReachGoal1);

    % Simulate path from startPose to decoyPose3
    waypoints2 = pthObj2.States;
    timeToReachGoal2 = pathLen2 / ss.AirSpeed;
    [xENU2, yENU2, zENU2] = exampleHelperSimulateUAV(waypoints2, ss.AirSpeed, timeToReachGoal2);

    % Simulate path from decoyPose3 to goalPose
    waypoints3 = pthObj3.States;
    timeToReachGoal3 = pathLen3 / ss.AirSpeed;
    [xENU3, yENU3, zENU3] = exampleHelperSimulateUAV(waypoints3, ss.AirSpeed, timeToReachGoal3);

    % Calculate overlap length between paths
    overlapThreshold = 1;
    overlapLength = 0;

    % Calculate overlap between path1 and path2
    for i = 1:length(xENU1) - 1
        segmentStart1 = [xENU1(i), yENU1(i), zENU1(i)];
        segmentEnd1 = [xENU1(i+1), yENU1(i+1), zENU1(i+1)];

        for j = 1:length(xENU2) - 1
            segmentStart2 = [xENU2(j), yENU2(j), zENU2(j)];
            segmentEnd2 = [xENU2(j+1), yENU2(j+1), zENU2(j+1)];

            segmentLength1 = norm(segmentEnd2 - segmentStart2);

            if minDistBetweenSegments(segmentStart1, segmentEnd1, segmentStart2, segmentEnd2) < overlapThreshold
                overlapLength = overlapLength + segmentLength1;
                break;
            end
        end
    end

    % Calculate overlap between path1 and path3
    for i = 1:length(xENU1) - 1
        segmentStart1 = [xENU1(i), yENU1(i), zENU1(i)];
        segmentEnd1 = [xENU1(i+1), yENU1(i+1), zENU1(i+1)];

        for j = 1:length(xENU3) - 1
            segmentStart2 = [xENU3(j), yENU3(j), zENU3(j)];
            segmentEnd2 = [xENU3(j+1), yENU3(j+1), zENU3(j+1)];

            segmentLength1 = norm(segmentEnd2 - segmentStart2);

            if minDistBetweenSegments(segmentStart1, segmentEnd1, segmentStart2, segmentEnd2) < overlapThreshold
                overlapLength = overlapLength + segmentLength1;
                break;
            end
        end
    end

    % Calculate dquantity
    dquantity1 = max((pathLen1 - overlapLength) / pathLen1,0);
  fprintf('simulation dquantity is: %.2f\n', dquantity1);
  maxProbGoals = run_rrt_second([], goalPoints, planner2, pthObj2);
        n = length(goalPoints) - 3;
        k=mode(maxProbGoals(1:n));
    fprintf(' Max Probability Goal Index = %d\n', k);

    
  
end
end
function [dquantity2,k] = cal1(startPose, goalPoints, omap)

    newPoint = mean(goalPoints);

    % Define UAV state space
    ss = ExampleHelperUAVStateSpace("MaxRollAngle", pi/6, ...
                                "AirSpeed", 6, ...
                                "FlightPathAngleLimit", [-0.1 0.1], ...
                                "Bounds", [-20 520; -20 520; 10 100; -pi pi]);

% Create validator for the map
sv = validatorOccupancyMap3D(ss, "Map", omap);
sv.ValidationDistance = 0.1;

% Create planner for the first path (startPose to goalPose)
planner1 = plannerRRT(ss, sv);
planner1.MaxConnectionDistance = 50;
planner1.GoalBias = 0.10;
planner1.MaxIterations = 400;



% Create planner for the second path (startPose to decoyPose3 to goalPose)
planner2 = plannerRRT(ss, sv);
planner2.MaxConnectionDistance = 50;
planner2.GoalBias = 0.10;
planner2.MaxIterations = 400;

% Create planner for the third path (decoyPose3 to goalPose)
planner3 = plannerRRT(ss, sv);
planner3.MaxConnectionDistance = 50;
planner3.GoalBias = 0.10;
planner3.MaxIterations = 400;

% Plan paths
    rng(50, "twister");
    [pthObj1, solnInfo1] = plan(planner1, startPose, goalPoints(1, :));

    rng(50, "twister");
    [pthObj2, solnInfo2] = plan(planner2, startPose, newPoint);

    rng(50, "twister");
    [pthObj3, solnInfo3] = plan(planner3, newPoint, goalPoints(1, :));

if (solnInfo1.IsPathFound && solnInfo2.IsPathFound && solnInfo3.IsPathFound)
    % Calculate path lengths
    pathLen1 = pathLength(pthObj1);
    pathLen2 = pathLength(pthObj2);
    pathLen3 = pathLength(pthObj3);

    % Simulate UAV trajectory based on fixed-wing guidance model
    timeToReachGoal1 = pathLen1 / ss.AirSpeed;
    waypoints1 = pthObj1.States;
    [xENU1, yENU1, zENU1] = exampleHelperSimulateUAV(waypoints1, ss.AirSpeed, timeToReachGoal1);

    % Simulate path from startPose to decoyPose3
    waypoints2 = pthObj2.States;
    timeToReachGoal2 = pathLen2 / ss.AirSpeed;
    [xENU2, yENU2, zENU2] = exampleHelperSimulateUAV(waypoints2, ss.AirSpeed, timeToReachGoal2);

    % Simulate path from decoyPose3 to goalPose
    waypoints3 = pthObj3.States;
    timeToReachGoal3 = pathLen3 / ss.AirSpeed;
    [xENU3, yENU3, zENU3] = exampleHelperSimulateUAV(waypoints3, ss.AirSpeed, timeToReachGoal3);

    % Calculate overlap length between paths
    overlapThreshold = 1;
    overlapLength = 0;

    % Calculate overlap between path1 and path2
    for i = 1:length(xENU1) - 1
        segmentStart1 = [xENU1(i), yENU1(i), zENU1(i)];
        segmentEnd1 = [xENU1(i+1), yENU1(i+1), zENU1(i+1)];

        for j = 1:length(xENU2) - 1
            segmentStart2 = [xENU2(j), yENU2(j), zENU2(j)];
            segmentEnd2 = [xENU2(j+1), yENU2(j+1), zENU2(j+1)];

            segmentLength1 = norm(segmentEnd2 - segmentStart2);

            if minDistBetweenSegments(segmentStart1, segmentEnd1, segmentStart2, segmentEnd2) < overlapThreshold
                overlapLength = overlapLength + segmentLength1;
                break;
            end
        end
    end

    % Calculate overlap between path1 and path3
    for i = 1:length(xENU1) - 1
        segmentStart1 = [xENU1(i), yENU1(i), zENU1(i)];
        segmentEnd1 = [xENU1(i+1), yENU1(i+1), zENU1(i+1)];

        for j = 1:length(xENU3) - 1
            segmentStart2 = [xENU3(j), yENU3(j), zENU3(j)];
            segmentEnd2 = [xENU3(j+1), yENU3(j+1), zENU3(j+1)];

            segmentLength1 = norm(segmentEnd2 - segmentStart2);

            if minDistBetweenSegments(segmentStart1, segmentEnd1, segmentStart2, segmentEnd2) < overlapThreshold
                overlapLength = overlapLength + segmentLength1;
                break;
            end
        end
    end
    % Calculate dquantity
    dquantity2 = max((pathLen1 - overlapLength) / pathLen1,0);
  fprintf('centroid1 dquantity is: %.2f\n', dquantity2);
  maxProbGoals = run_rrt_second([], goalPoints, planner2, pthObj2);
        n = length(goalPoints) - 3;
        k=mode(maxProbGoals(1:n));
    fprintf(' Max Probability Goal Index = %d\n', k);
end
end
function [dquantity3,k] = cal2(startPose, goalPoints, omap)

    newPoint = mean(goalPoints(2:length(goalPoints),:));

    % Define UAV state space
    ss = ExampleHelperUAVStateSpace("MaxRollAngle", pi/6, ...
                                "AirSpeed", 6, ...
                                "FlightPathAngleLimit", [-0.1 0.1], ...
                                "Bounds", [-20 520; -20 520; 10 100; -pi pi]);

% Create validator for the map
sv = validatorOccupancyMap3D(ss, "Map", omap);
sv.ValidationDistance = 0.1;

% Create planner for the first path (startPose to goalPose)
planner1 = plannerRRT(ss, sv);
planner1.MaxConnectionDistance = 50;
planner1.GoalBias = 0.10;
planner1.MaxIterations = 400;



% Create planner for the second path (startPose to decoyPose3 to goalPose)
planner2 = plannerRRT(ss, sv);
planner2.MaxConnectionDistance = 50;
planner2.GoalBias = 0.10;
planner2.MaxIterations = 400;

% Create planner for the third path (decoyPose3 to goalPose)
planner3 = plannerRRT(ss, sv);
planner3.MaxConnectionDistance = 50;
planner3.GoalBias = 0.10;
planner3.MaxIterations = 400;

% Plan paths
    rng(50, "twister");
    [pthObj1, solnInfo1] = plan(planner1, startPose, goalPoints(1, :));

    rng(50, "twister");
    [pthObj2, solnInfo2] = plan(planner2, startPose, newPoint);

    rng(50, "twister");
    [pthObj3, solnInfo3] = plan(planner3, newPoint, goalPoints(1, :));

if (solnInfo1.IsPathFound && solnInfo2.IsPathFound && solnInfo3.IsPathFound)
    % Calculate path lengths
    pathLen1 = pathLength(pthObj1);
    pathLen2 = pathLength(pthObj2);
    pathLen3 = pathLength(pthObj3);

    % Simulate UAV trajectory based on fixed-wing guidance model
    timeToReachGoal1 = pathLen1 / ss.AirSpeed;
    waypoints1 = pthObj1.States;
    [xENU1, yENU1, zENU1] = exampleHelperSimulateUAV(waypoints1, ss.AirSpeed, timeToReachGoal1);

    % Simulate path from startPose to decoyPose3
    waypoints2 = pthObj2.States;
    timeToReachGoal2 = pathLen2 / ss.AirSpeed;
    [xENU2, yENU2, zENU2] = exampleHelperSimulateUAV(waypoints2, ss.AirSpeed, timeToReachGoal2);

    % Simulate path from decoyPose3 to goalPose
    waypoints3 = pthObj3.States;
    timeToReachGoal3 = pathLen3 / ss.AirSpeed;
    [xENU3, yENU3, zENU3] = exampleHelperSimulateUAV(waypoints3, ss.AirSpeed, timeToReachGoal3);

    % Calculate overlap length between paths
    overlapThreshold = 1;
    overlapLength = 0;

    % Calculate overlap between path1 and path2
    for i = 1:length(xENU1) - 1
        segmentStart1 = [xENU1(i), yENU1(i), zENU1(i)];
        segmentEnd1 = [xENU1(i+1), yENU1(i+1), zENU1(i+1)];

        for j = 1:length(xENU2) - 1
            segmentStart2 = [xENU2(j), yENU2(j), zENU2(j)];
            segmentEnd2 = [xENU2(j+1), yENU2(j+1), zENU2(j+1)];

            segmentLength1 = norm(segmentEnd2 - segmentStart2);

            if minDistBetweenSegments(segmentStart1, segmentEnd1, segmentStart2, segmentEnd2) < overlapThreshold
                overlapLength = overlapLength + segmentLength1;
                break;
            end
        end
    end

    % Calculate overlap between path1 and path3
    for i = 1:length(xENU1) - 1
        segmentStart1 = [xENU1(i), yENU1(i), zENU1(i)];
        segmentEnd1 = [xENU1(i+1), yENU1(i+1), zENU1(i+1)];

        for j = 1:length(xENU3) - 1
            segmentStart2 = [xENU3(j), yENU3(j), zENU3(j)];
            segmentEnd2 = [xENU3(j+1), yENU3(j+1), zENU3(j+1)];

            segmentLength1 = norm(segmentEnd2 - segmentStart2);

            if minDistBetweenSegments(segmentStart1, segmentEnd1, segmentStart2, segmentEnd2) < overlapThreshold
                overlapLength = overlapLength + segmentLength1;
                break;
            end
        end
    end
    % Calculate dquantity
    dquantity3 = max((pathLen1 - overlapLength) / pathLen1,0);
    fprintf('centroid2 dquantity is: %.2f\n', dquantity3);
    maxProbGoals = run_rrt_second([], goalPoints, planner2, pthObj2);
    n = length(goalPoints) - 3;
    k=mode(maxProbGoals(1:n));
    fprintf(' Max Probability Goal Index = %d\n', k);
 
    
end
end
function [dquantity4,k] = cal3(startPose, goalPoints, omap)
    firstPoint = goalPoints(1, 1:3);

    % Calculate distances between each point and the first point
    distances = sqrt(sum((goalPoints(:, 1:3) - firstPoint).^2, 2));
    
    % Set the distance of the first point to infinity to exclude itself
    distances(1) = inf;
    
    % Find the index of the closest point
    [~, closestIndex] = min(distances);
    % Output the closest point to the first point
    decoyPose3 = goalPoints(closestIndex, :);
    newPoint = (goalPoints(1,:)+decoyPose3)/2;

    % Define UAV state space
    ss = ExampleHelperUAVStateSpace("MaxRollAngle", pi/6, ...
                                "AirSpeed", 6, ...
                                "FlightPathAngleLimit", [-0.1 0.1], ...
                                "Bounds", [-20 520; -20 520; 10 100; -pi pi]);

% Create validator for the map
sv = validatorOccupancyMap3D(ss, "Map", omap);
sv.ValidationDistance = 0.1;

% Create planner for the first path (startPose to goalPose)
planner1 = plannerRRT(ss, sv);
planner1.MaxConnectionDistance = 50;
planner1.GoalBias = 0.10;
planner1.MaxIterations = 400;



% Create planner for the second path (startPose to decoyPose3 to goalPose)
planner2 = plannerRRT(ss, sv);
planner2.MaxConnectionDistance = 50;
planner2.GoalBias = 0.10;
planner2.MaxIterations = 400;

% Create planner for the third path (decoyPose3 to goalPose)
planner3 = plannerRRT(ss, sv);
planner3.MaxConnectionDistance = 50;
planner3.GoalBias = 0.10;
planner3.MaxIterations = 400;

% Plan paths
    rng(50, "twister");
    [pthObj1, solnInfo1] = plan(planner1, startPose, goalPoints(1, :));

    rng(50, "twister");
    [pthObj2, solnInfo2] = plan(planner2, startPose, newPoint);

    rng(50, "twister");
    [pthObj3, solnInfo3] = plan(planner3, newPoint, goalPoints(1, :));

if (solnInfo1.IsPathFound && solnInfo2.IsPathFound && solnInfo3.IsPathFound)
    % Calculate path lengths
    pathLen1 = pathLength(pthObj1);
    pathLen2 = pathLength(pthObj2);
    pathLen3 = pathLength(pthObj3);

    % Simulate UAV trajectory based on fixed-wing guidance model
    timeToReachGoal1 = pathLen1 / ss.AirSpeed;
    waypoints1 = pthObj1.States;
    [xENU1, yENU1, zENU1] = exampleHelperSimulateUAV(waypoints1, ss.AirSpeed, timeToReachGoal1);

    % Simulate path from startPose to decoyPose3
    waypoints2 = pthObj2.States;
    timeToReachGoal2 = pathLen2 / ss.AirSpeed;
    [xENU2, yENU2, zENU2] = exampleHelperSimulateUAV(waypoints2, ss.AirSpeed, timeToReachGoal2);

    % Simulate path from decoyPose3 to goalPose
    waypoints3 = pthObj3.States;
    timeToReachGoal3 = pathLen3 / ss.AirSpeed;
    [xENU3, yENU3, zENU3] = exampleHelperSimulateUAV(waypoints3, ss.AirSpeed, timeToReachGoal3);

    % Calculate overlap length between paths
    overlapThreshold = 1;
    overlapLength = 0;

    % Calculate overlap between path1 and path2
    for i = 1:length(xENU1) - 1
        segmentStart1 = [xENU1(i), yENU1(i), zENU1(i)];
        segmentEnd1 = [xENU1(i+1), yENU1(i+1), zENU1(i+1)];

        for j = 1:length(xENU2) - 1
            segmentStart2 = [xENU2(j), yENU2(j), zENU2(j)];
            segmentEnd2 = [xENU2(j+1), yENU2(j+1), zENU2(j+1)];

            segmentLength1 = norm(segmentEnd2 - segmentStart2);

            if minDistBetweenSegments(segmentStart1, segmentEnd1, segmentStart2, segmentEnd2) < overlapThreshold
                overlapLength = overlapLength + segmentLength1;
                break;
            end
        end
    end

    % Calculate overlap between path1 and path3
    for i = 1:length(xENU1) - 1
        segmentStart1 = [xENU1(i), yENU1(i), zENU1(i)];
        segmentEnd1 = [xENU1(i+1), yENU1(i+1), zENU1(i+1)];

        for j = 1:length(xENU3) - 1
            segmentStart2 = [xENU3(j), yENU3(j), zENU3(j)];
            segmentEnd2 = [xENU3(j+1), yENU3(j+1), zENU3(j+1)];

            segmentLength1 = norm(segmentEnd2 - segmentStart2);

            if minDistBetweenSegments(segmentStart1, segmentEnd1, segmentStart2, segmentEnd2) < overlapThreshold
                overlapLength = overlapLength + segmentLength1;
                break;
            end
        end
    end
    % Calculate dquantity
    dquantity4 = max((pathLen1 - overlapLength) / pathLen1,0);
  fprintf('centroid3 dquantity is: %.2f\n', dquantity4);
 maxProbGoals = run_rrt_second([], goalPoints, planner2, pthObj2);
        n = length(goalPoints) - 3;
        k=mode(maxProbGoals(1:n));
    fprintf(' Max Probability Goal Index = %d\n', k);
end
end
function d = minDistBetweenSegments(A, B, C, D)
        AB = B - A;
        CD = D - C;
        AC = C - A;

        AB_squared = dot(AB, AB);
        CD_squared = dot(CD, CD);

        AB_dot_AC = dot(AB, AC);
        AB_dot_CD = dot(AB, CD);
        CD_dot_AC = dot(CD, AC);

        denom = AB_squared * CD_squared - AB_dot_CD * AB_dot_CD;
        t = (AB_dot_AC * CD_squared - CD_dot_AC * AB_dot_CD) / denom;
        u = (CD_dot_AC * AB_squared - AB_dot_AC * AB_dot_CD) / denom;

        t = max(0, min(1, t));
        u = max(0, min(1, u));

        closestPointOnAB = A + t * AB;
        closestPointOnCD = C + u * CD;
        d = norm(closestPointOnAB - closestPointOnCD);
    end
function maxProbGoals = run_rrt_second(~, goals, planner2, path1)
    numPointsToSelect = 280;
    numGoals = size(goals, 1);

    % Initialize cell arrays to store metrics
    pathLengths = cell(numPointsToSelect, numGoals);
    angles = cell(numPointsToSelect, numGoals);
    orientationOffsets = cell(numPointsToSelect, numGoals);
    sld_values = cell(numPointsToSelect, numGoals);
    angularDifference = cell(numPointsToSelect, numGoals);

        % Use the original path1 directly
        numStates = size(path1.States, 1);
        selectedIndices = round(linspace(1, numStates, numPointsToSelect));
        selectedStates = path1.States(selectedIndices, 1:4);

        % Initialize arrays to store the max probability goals and values
        maxProbGoals = zeros(numPointsToSelect, 1);
        maxProbValues = zeros(numPointsToSelect, 1);

        % Running RRT from selected indices and calculating probabilities
        for i = 1:numPointsToSelect
            % Arrays to store weights
            weights_point = zeros(1, numGoals);
            weights_angle = zeros(1, numGoals);
            weights_orientation = zeros(1, numGoals);
            weights_angdiff = zeros(1, numGoals);
            
            startpose = [selectedStates(i,:), 0];
            start = startpose(1:4);

            for j = 1:numGoals
                goal = goals(j, :);
                [path, ~] = plan(planner2, start, goal);

                % Path length
                pathLengths{i, j} = sum(sqrt(sum(diff(path.States(:, 1:3)).^2, 2)));

                % SLD
                sld_values{i, j} = norm(selectedStates(i, 1:3) - goals(j, 1:3));

                % Angle
                goalvector = goal - start;
                referenceDirection = goal(1:3);
                angle = atan2(norm(cross(referenceDirection, goalvector(1:3))), dot(referenceDirection, goalvector(1:3)));
                angles{i, j} = rad2deg(angle);

                % Orientation
                agentOrientation = selectedStates(i, 4);
                goalOrientation = goals(j, 4);
                orientationOffset = rad2deg(goalOrientation - agentOrientation);
                orientationOffset = mod(orientationOffset, 360);
                if orientationOffset == 0
                    orientationOffset = 0.1;
                end
                if orientationOffset > 180
                    orientationOffset = 360 - orientationOffset;
                end
                orientationOffsets{i, j} = orientationOffset;

                % Angular difference
                agentOrientation = rad2deg(selectedStates(i, 4));
                goalVector = goal(1:3) - selectedStates(i, 1:3);
                angle_difference = acosd(dot(goalVector, [cosd(agentOrientation), sind(agentOrientation), 0]) / norm(goalVector));
                if angle_difference > 180
                    angle_difference = 360 - angle_difference;
                end
                angularDifference{i, j} = angle_difference;

                % Compute weights
                weights_point(j) = 1 / pathLengths{i, j};
                weights_angle(j) = 1 / angles{i, j};
                weights_orientation(j) = 1 / orientationOffsets{i, j};
                weights_angdiff(j) = 1 / angularDifference{i, j};
            end

            % Combine weights
            weight_combined = weights_point .* weights_angle .* weights_orientation .* weights_angdiff;
            normalized_weights = weight_combined / sum(weight_combined);
            probabilities = normalized_weights / sum(normalized_weights);

            % Find the maximum probability goal
            [maxProb, maxIndex] = max(probabilities);
            maxProbGoals(i) = maxIndex;
            maxProbValues(i) = maxProb;
        end
    end
