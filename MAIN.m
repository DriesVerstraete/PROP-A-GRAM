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
    
    tmp1 = repmat(station(i).points, 1, 1, size(station(i).points,1));
    tmp2 = repmat(reshape(station(i).points', 1, 3, size(station(i).points,1)), size(station(i).points,1),1,1);
    lens = sqrt(sum((tmp2 - tmp1).^2, 2));
    [idx1,~] = find(lens >= max(max(lens)));
    tmp1 = station(i).points(idx1,:);
    [~,idx] = sort(tmp1(:,1),1,'ascend');
    tmp1 = tmp1(idx,:);
    
    station(i).leading_edge = tmp1(1,:);
    station(i).trailing_edge = tmp1(2,:);
    
%     hFig20 = figure(20);
%     clf(20);
%     scatter(station(i).points(:,1), station(i).points(:,3),20,'xk')
%     hold on
%     scatter(tmp1(:,1), tmp1(:,3), 'or')
%     hold off
%     grid minor
%     box on
%     axis equal
    
    %%
%     LE = tmp1(1,:);
%     TE = tmp1(2,:);
%     chord = norm(LE - TE);
%         
%     p = [0.0288    0.3368    0.1159    0.2447    0.2544    0.0302   -0.9551    0.0010    0.0108    4.3094   11.1424    0.3994];
%     truncation = 0.2;
%     eps = acosd(dot([-1 0], (LE(1:2:3) - TE(1:2:3))./chord));
%     z = [p, truncation, eps];
%     
%     nvars = length(z);
%     lb = [0.005 0.25    0.03    -1  0.25    0.03 -1  -0.01    0.008  -5 10 0.32 0 -10];
%     ub = [0.04  0.4     0.15     1   0.4     0.10  1   0.01    0.012  5 20  0.40 0.2 70];
%    
%     options = optimoptions('ga', 'Display', 'iter', 'InitialPopulation',z, 'UseParallel', true, 'MaxGenerations', 1000, 'StallGenLimit', 50, 'MutationFcn','mutationadaptfeasible', 'CreationFcn', 'gacreationlinearfeasible');
%     [x,fval,exitflag,output,population,scores] = ga({@fcnOBJECTIVE, station(i).points(:,1:2:3), LE, chord},nvars,[],[],[],[],lb,ub,[],[],options);
%     
%     [~, foil] = fcnOBJECTIVE(x, station(i).points(:,1:2:3), LE, chord);
%     plot(foil(:,1), foil(:,2),'--om')
 
end
hold off

%% Chordwise Stations

% This part needs to be fixed
for i = 1:num_stations
    chord_actual(i) = norm(station(i).leading_edge - station(i).trailing_edge);
    tmp = (station(i).leading_edge - station(i).trailing_edge)./chord_actual(i);
    station(i).pitch = acosd(dot([-1 0 0], tmp, 2));
    
%     hFig20 = figure(20);
%     clf(20);
%     scatter(station(i).points(:,1), station(i).points(:,3),20,'xk')
%     grid minor
%     box on
%     axis equal
        
    beta(i) = station(i).pitch;
end

%% Spanwise
% hold on
% for i = 1:num_stations - 1
%     x = [station(i).leading_edge; station(i).trailing_edge; station(i+1).trailing_edge; station(i+1).leading_edge;];
%     y = [repmat(station_location(i), 2, 1); repmat(station_location(i+1),2,1)];
%     z = [station(i).camberline(x(1:2)); station(i+1).camberline(x(3:4))];
%     patch(x,y,z,'b','LineWidth',2,'EdgeAlpha',0.6,'FaceAlpha',0.2);
% end
% hold off

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




