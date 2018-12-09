clc
clear

% 1. Import Blade Points (from Scan Points/)
%        Should be STL, n x 3 x 3 matrix (use fcnSTLREAD if needed)
% 2. Manually rotate the scan about 3-axis
%        Rotation axis should end up in global Z, spanwise in global Y
% 3. Translate it so the axis of rotation is at [0 0 0]
% 4. Set the root and tip Y values to crop to the blade
% 5. Set number of spanwise stations
% 6. Compare to known (r/R, c/R, beta) geometry

% matPOINTS = fcnSTLREAD('CAD Geom/master_airscrew.stl');
blade_name = 'MAE 11x7'

%% Import Blade
load(['Scan Points/', blade_name, '.mat'])
matVLST = unique(cat(1, matPOINTS(:,:,1), matPOINTS(:,:,2), matPOINTS(:,:,3)),'rows');

% scatter3(matVLST(:,1), matVLST(:,2), matVLST(:,3), 2,'xk')
% xlabel('X-Dir','Fontsize',15)
% ylabel('Y-Dir','Fontsize',15)
% zlabel('Z-Dir','Fontsize',15)
% axis equal
% grid minor
% box on

%% Rotate Blade so hub is parallel with x-z plane
% Rotate about y-axis
u_x = 0;
u_y = 1;
u_z = 0;
theta = deg2rad(11.3);

R = [ cos(theta) + u_x.^2.*(1 - cos(theta)) u_x.*u_y.*(1 - cos(theta)) - u_z.*sin(theta) u_x.*u_z.*(1 - cos(theta)) + u_y.*sin(theta); ...
    u_y.*u_x.*(1 - cos(theta)) + u_z.*sin(theta) cos(theta) + u_y.^2.*(1 - cos(theta)) u_y.*u_z.*(1 - cos(theta)) - u_x.*sin(theta); ...
    u_z.*u_x.*(1 - cos(theta)) - u_y.*sin(theta) u_z.*u_y.*(1 - cos(theta)) + u_x.*sin(theta) cos(theta) + u_z.^2.*(1 - cos(theta))];

for i = 1:size(matVLST)
    matVLST(i,:) = R*matVLST(i,:)';
end

% Rotate about x-axis
u_x = 1;
u_y = 0;
u_z = 0;
theta = deg2rad(2.5);

R = [ cos(theta) + u_x.^2.*(1 - cos(theta)) u_x.*u_y.*(1 - cos(theta)) - u_z.*sin(theta) u_x.*u_z.*(1 - cos(theta)) + u_y.*sin(theta); ...
    u_y.*u_x.*(1 - cos(theta)) + u_z.*sin(theta) cos(theta) + u_y.^2.*(1 - cos(theta)) u_y.*u_z.*(1 - cos(theta)) - u_x.*sin(theta); ...
    u_z.*u_x.*(1 - cos(theta)) - u_y.*sin(theta) u_z.*u_y.*(1 - cos(theta)) + u_x.*sin(theta) cos(theta) + u_z.^2.*(1 - cos(theta))];

for i = 1:size(matVLST)
    matVLST(i,:) = R*matVLST(i,:)';
end

% Rotate about z-axis
u_x = 0;
u_y = 0;
u_z = 1;
theta = deg2rad(180);

R = [ cos(theta) + u_x.^2.*(1 - cos(theta)) u_x.*u_y.*(1 - cos(theta)) - u_z.*sin(theta) u_x.*u_z.*(1 - cos(theta)) + u_y.*sin(theta); ...
    u_y.*u_x.*(1 - cos(theta)) + u_z.*sin(theta) cos(theta) + u_y.^2.*(1 - cos(theta)) u_y.*u_z.*(1 - cos(theta)) - u_x.*sin(theta); ...
    u_z.*u_x.*(1 - cos(theta)) - u_y.*sin(theta) u_z.*u_y.*(1 - cos(theta)) + u_x.*sin(theta) cos(theta) + u_z.^2.*(1 - cos(theta))];

for i = 1:size(matVLST)
    matVLST(i,:) = R*matVLST(i,:)';
end

% Translate
matVLST = matVLST + [0.2992 0.1879 -1.435];

%% Crop to Single Blade
root = 0.03; % Manually set
tip = 0.275; % Manually set
Radius = max(matVLST(:,2)); % True radius (sort of) before chopping off the tip
matVLST(matVLST(:,2) < root,:) = [];
matVLST(matVLST(:,2) > tip,:) = [];

%% Stations
hFig1 = figure(1);
clf(1);
scatter3(matVLST(:,1), matVLST(:,2), matVLST(:,3), 2,'xk')
xlabel('X-Dir','Fontsize',15)
ylabel('Y-Dir','Fontsize',15)
zlabel('Z-Dir','Fontsize',15)
axis equal
grid minor
box on

num_stations = 20;
station_width = 0.005;

station_location = linspace(root, tip, num_stations);
hold on
for i = 1:num_stations
    idx = station_location(i) - station_width <= matVLST(:,2) & matVLST(:,2) <= station_location(i) + station_width;
    station(i).points = matVLST(idx,:);
    station(i).points(:,2) = station_location(i);
    scatter3(station(i).points(:,1), station(i).points(:,2), station(i).points(:,3), 2, 'xr');
    [fitresult, ~] = createFit(station(i).points(:,1), station(i).points(:,3));
    station(i).camberline = fitresult;
    station(i).leading_edge = min(station(i).points(:,1));
    station(i).trailing_edge = max(station(i).points(:,1));
end
hold off

%% Chordwise Stations

% This part needs to be fixed
for i = 1:num_stations
    station(i).chord_stations = [station(i).leading_edge, station(i).trailing_edge];
    station(i).chord_z = [station(i).camberline(station(i).leading_edge), station(i).camberline(station(i).trailing_edge)]';
    tmp = [station(i).leading_edge - station(i).trailing_edge, station(i).chord_z(1) - station(i).chord_z(end)];
    chord_actual(i) = norm(tmp);
    tmp = tmp./sqrt(sum(tmp.^2,2));
    station(i).pitch = acosd(dot([-1 0], tmp, 2));
    
    %     hFig20 = figure(20);
    %     clf(20);
    %     scatter(station(i).points(:,1), station(i).points(:,3),20,'xk')
    %     grid minor
    %     box on
    %     axis equal
    %
    %     hold on
    %     tmp2 = linspace(station(i).leading_edge, station(i).trailing_edge, 100);
    %     plot(tmp2, station(i).camberline(tmp2), ':^m')
    %     plot(station(i).chord_stations, station(i).chord_z, '-sr')
    %     hold off
    
    beta(i) = station(i).pitch;
end

%% Spanwise
hold on
for i = 1:num_stations - 1
    x = [station(i).leading_edge; station(i).trailing_edge; station(i+1).trailing_edge; station(i+1).leading_edge;];
    y = [repmat(station_location(i), 2, 1); repmat(station_location(i+1),2,1)];
    z = [station(i).camberline(x(1:2)); station(i+1).camberline(x(3:4))];
    patch(x,y,z,'b','LineWidth',2,'EdgeAlpha',0.6,'FaceAlpha',0.2);
end
hold off

%% Compare
try
    A = dlmread(['Known Geometry/', blade_name '.txt'],'',1,0);
    
    hFig2 = figure(2);
    clf(2);
    xlabel('r/R','FontSize',15);
    yyaxis left
    ylabel('c/R','FontSize',15);
    hold on
    plot(A(:,1), A(:,2), '-ok')
    plot(station_location./Radius, chord_actual./Radius, '-b^')
    
    yyaxis right
    ylabel('Beta (Degrees)','FontSize',15);
    plot(A(:,1), A(:,3), '--ok')
    plot(station_location./Radius, beta, '--rs')
    
    grid minor
    box on
    
    legend('Illinois r/R','My r/R (Actual Chord)','Illinois Beta','My Beta','Location','SouthWest')
catch
    disp('Could not find known geometry to compare to.')
end



