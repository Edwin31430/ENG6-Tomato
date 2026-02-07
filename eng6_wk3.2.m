clear; clc; close all;

%% ================================================================
%  1. LOAD DATA
% ================================================================
T = readtable('data_NXTgen_merged_titles_clean.xlsx.xlsx');
P = readtable('Overview_of_sensors.xlsx.xlsx');

%% ================================================================
%  2. CLEAN METADATA
% ================================================================
P.S_n = string(P.S_n);
P.S_n = upper(strtrim(P.S_n));

validIDs = [
    "0013A20041E3D7B5"
    "0013A20041E2FBE9"
    "0013A20041E273BA"
    "0013A20041E273D5"
    "0013A20041E3D7F4"
    "04A39"
    "04B0F"
    "04ACD"
    "04A15"
    "007B9"
    "00817"
    "00818"
    "0081D"
];

P = P(ismember(P.S_n, validIDs), :);

P.RowNumber        = str2double(string(P.RowNumber));
P.Height_m_        = str2double(string(P.Height_m_));
P.DistanceFromPath = str2double(string(P.DistanceFromPath));

hasCoords = isfinite(P.RowNumber) & isfinite(P.Height_m_);
P = P(hasCoords,:);

fprintf("Metadata sensors with coordinates: %d\n", height(P));

%% ================================================================
%  3. FIND TEMPERATURE COLUMNS
% ================================================================
names = T.Properties.VariableNames;
tempMask = contains(names,"Temperature__C_") | contains(names,"AirTemperature__C_");
tempNames = names(tempMask);

%% ================================================================
%  4. BUILD SENSOR MAP
% ================================================================
map = table();
map.SensorID   = string.empty(0,1);
map.ColName    = string.empty(0,1);
map.X          = zeros(0,1);
map.Y          = zeros(0,1);
map.Z          = zeros(0,1);

for i = 1:height(P)
    id = P.S_n(i);

    if strlength(id)==5
        hit = tempNames(contains(tempNames,lower(id)));
    else
        suffix = lower(extractAfter(id,strlength(id)-3));
        hit = tempNames(contains(tempNames,"sensor_"+suffix));
    end

    if isempty(hit); continue; end

    map = [map; {id,string(hit(1)),P.DistanceFromPath(i),P.RowNumber(i),P.Height_m_(i)}];
end

map.Properties.VariableNames = ["SensorID","ColName","X","Y","Z"];
disp(map)

%% ================================================================
%  5. EXTRACT TEMPERATURE MATRIX
% ================================================================
Tmat = T{:,map.ColName};

%% ================================================================
%  6. TIME INDEX
% ================================================================
t0 = 1000;
Trow = Tmat(t0,:);

%% ================================================================
%  7. BUILD CONSISTENT VECTORS
% ================================================================
Tv_all = Trow(:);
y_all  = map.Y(:);
z_all  = map.Z(:);

valid = isfinite(Tv_all) & isfinite(y_all) & isfinite(z_all);

Tv   = Tv_all(valid);
yPos = y_all(valid);
zPos = z_all(valid);

fprintf("Valid sensors: %d\n",numel(Tv));

% Merge duplicates
[YZuniq,~,ic] = unique([yPos zPos],'rows');
Tv   = accumarray(ic,Tv,[],@mean);
yPos = YZuniq(:,1);
zPos = YZuniq(:,2);

%% ================================================================
%  8. INTERPOLATION GRID
% ================================================================
yq = linspace(min(yPos),max(yPos),80);
zq = linspace(min(zPos),max(zPos),40);
[YI,ZI] = meshgrid(yq,zq);

%% ================================================================
%  9. INTERPOLATION METHODS (LINEAR, SPLINE, KRIGING)
% ================================================================

% --- Linear ---
F_lin = scatteredInterpolant(yPos,zPos,Tv,'linear','nearest');
T_lin = F_lin(YI,ZI);

% --- Spline (via griddata) ---
T_spl = griddata(yPos,zPos,Tv,YI,ZI,'cubic');

% --- Kriging (GPR) ---
gprMdl = fitrgp([yPos zPos],Tv,'KernelFunction','squaredexponential');
T_krig = reshape(predict(gprMdl,[YI(:) ZI(:)]),size(YI));


%% ================================================================
%  10. VISUALIZATION
% ================================================================

figure;

subplot(1,3,1)
contourf(YI,ZI,T_lin,20,'LineColor','none'); colorbar
hold on; scatter(yPos,zPos,60,Tv,'filled','k')
title('Linear'); xlabel('Row'); ylabel('Height')

subplot(1,3,2)
contourf(YI,ZI,T_spl,20,'LineColor','none'); colorbar
hold on; scatter(yPos,zPos,60,Tv,'filled','k')
title('Spline')

subplot(1,3,3)
contourf(YI,ZI,T_krig,20,'LineColor','none'); colorbar
hold on; scatter(yPos,zPos,60,Tv,'filled','k')
title('Kriging (GPR)')

sgtitle(sprintf('Spatial interpolation with height at t = %d',t0))

%% ================================================================
%  10.5 RMSE COMPARISON (LEAVE-ONE-OUT)
% ================================================================

nS = numel(Tv);

rmse_lin  = zeros(nS,1);
rmse_spl  = zeros(nS,1);
rmse_krig = zeros(nS,1);

for i = 1:nS
    idx = true(nS,1);
    idx(i) = false;

    y_sub = yPos(idx);
    z_sub = zPos(idx);
    T_sub = Tv(idx);

    % ---- Linear ----
    F = scatteredInterpolant(y_sub,z_sub,T_sub,'linear','nearest');
    rmse_lin(i) = (F(yPos(i),zPos(i)) - Tv(i)).^2;

    % ---- Spline ----
    T_pred_spl = griddata(y_sub,z_sub,T_sub,yPos(i),zPos(i),'cubic');
    rmse_spl(i) = (T_pred_spl - Tv(i)).^2;

    % ---- Kriging ----
    gpr = fitrgp([y_sub z_sub],T_sub,'KernelFunction','squaredexponential');
    T_pred_k = predict(gpr,[yPos(i) zPos(i)]);
    rmse_krig(i) = (T_pred_k - Tv(i)).^2;
end

RMSE_linear  = sqrt(mean(rmse_lin,'omitnan'));
RMSE_spline  = sqrt(mean(rmse_spl,'omitnan'));
RMSE_kriging = sqrt(mean(rmse_krig,'omitnan'));

fprintf('\nInterpolation RMSE:\n');
fprintf('Linear  : %.4f\n',RMSE_linear);
fprintf('Spline  : %.4f\n',RMSE_spline);
fprintf('Kriging : %.4f\n',RMSE_kriging);

figure;
bar([RMSE_linear RMSE_spline RMSE_kriging])
set(gca,'XTickLabel',{'Linear','Spline','Kriging'})
ylabel('RMSE (°C)')
title('Interpolation RMSE Comparison')
grid on

%% ================================================================
%  10.7 SPLINE-BASED SENSOR PLACEMENT (UNCERTAINTY MAP)
% ================================================================

% Build uncertainty proxy: distance to nearest sensor
distMap = zeros(size(YI));

for k = 1:numel(YI)
    d = sqrt((yPos - YI(k)).^2 + (zPos - ZI(k)).^2);
    distMap(k) = min(d);
end

% Best candidate = max distance
[~,idxBest] = max(distMap(:));
bestRow    = YI(idxBest);
bestHeight = ZI(idxBest);

fprintf('\nRecommended new sensor (Spline support):\n');
fprintf('Row = %.2f , Height = %.2f\n',bestRow,bestHeight);

%% ---- Plot spline + uncertainty + placement ----
figure;

contourf(YI,ZI,T_spl,20,'LineColor','none'); 
colorbar
hold on;

scatter(yPos,zPos,80,Tv,'filled','k')
plot(bestRow,bestHeight,'rp','MarkerSize',18,'MarkerFaceColor','r')

title('Spline interpolation with recommended sensor placement')
xlabel('Row')
ylabel('Height')
legend('Spline temperature','Sensors','Recommended new sensor','Location','best')

%% ================================================================
%  11. PCA ANALYSIS (MULTICOLLINEARITY / DOMINANT MODES)
% ================================================================

% Build full matrix: time × sensors
Tfull = Tmat;   % already Ntime x Nsensors

% Remove rows with NaN
goodRows = all(isfinite(Tfull),2);
Tfull = Tfull(goodRows,:);

% Normalize
Tnorm = (Tfull - mean(Tfull)) ./ std(Tfull);

% PCA
[coeff,score,latent,~,explained] = pca(Tnorm);

fprintf('\nPCA Variance Explained:\n');
disp(explained(1:5))

%% Plot variance
figure;
pareto(explained)
xlabel('Principal Component')
ylabel('Variance Explained (%)')
title('PCA Variance Explained')

%% Spatial loading plot (PC1, PC2)
figure;
subplot(1,2,1)
scatter(map.Y,map.Z,80,coeff(:,1),'filled'); colorbar
title('PC1 Spatial Loading')
xlabel('Row'); ylabel('Height')

subplot(1,2,2)
scatter(map.Y,map.Z,80,coeff(:,2),'filled'); colorbar
title('PC2 Spatial Loading')
xlabel('Row'); ylabel('Height')

sgtitle('Spatial structure from PCA')

%% ================================================================
%  12. SENSOR PLACEMENT OPTIMIZATION (LEAVE-ONE-OUT)
% ================================================================

nS = numel(Tv);
rmse_LOO = zeros(nS,1);

for i = 1:nS
    idx = true(nS,1);
    idx(i) = false;

    y_sub = yPos(idx);
    z_sub = zPos(idx);
    T_sub = Tv(idx);

    F = scatteredInterpolant(y_sub, z_sub, T_sub, 'linear', 'nearest');

    T_pred = F(yPos(i), zPos(i));
    rmse_LOO(i) = sqrt((T_pred - Tv(i)).^2);
end

%% ---- FORCE CONSISTENT SIZES ----
mapValid = map(valid,:);

SensorID = mapValid.SensorID(:);
RowPos   = yPos(:);
Height   = zPos(:);
TempVal  = Tv(:);
RMSEval  = rmse_LOO(:);

N = min([numel(SensorID), numel(RowPos), numel(Height), numel(TempVal), numel(RMSEval)]);

SensorID = SensorID(1:N);
RowPos   = RowPos(1:N);
Height   = Height(1:N);
TempVal  = TempVal(1:N);
RMSEval  = RMSEval(1:N);

%% ---- BUILD TABLE ----
placementTable = table( ...
    SensorID, ...
    RowPos, ...
    Height, ...
    TempVal, ...
    RMSEval, ...
    'VariableNames', {'SensorID','RowNumber','Height','Temperature','RMSE_if_removed'} );

disp('Sensor placement importance:')
disp(placementTable)

%% ---- PLOT ----
figure;
bar(RMSEval)
xticks(1:N)
xticklabels(SensorID)
xtickangle(45)
ylabel('RMSE when removed')
title('Sensor placement importance (leave-one-out)')




