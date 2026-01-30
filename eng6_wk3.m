%% ================================================================
%  1. LOAD DATA
% ================================================================
clear; clc; close all;

T = readtable('data_NXTgen_merged_titles_clean.xlsx.xlsx');
P = readtable('Overview_of_sensors.xlsx.xlsx');

%% ================================================================
%  2. CLEAN METADATA (coordinates + serial numbers)
% ================================================================
P.S_n = string(P.S_n);
P.S_n = strtrim(P.S_n);
P.S_n = replace(P.S_n, " ", "");
P.S_n = upper(P.S_n);
P.S_n(ismissing(P.S_n)) = "";
P.S_n(P.S_n == "NAN") = "";

% Keep only the sensors you actually use
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

% Coordinates
P.RowNumber        = str2double(strrep(strtrim(string(P.RowNumber)), ',', '.'));
P.Height_m_        = str2double(strrep(strtrim(string(P.Height_m_)), ',', '.'));
P.DistanceFromPath = str2double(strrep(strtrim(string(P.DistanceFromPath)), ',', '.'));

hasCoords = isfinite(P.RowNumber) & isfinite(P.Height_m_) & isfinite(P.DistanceFromPath);
P = P(hasCoords,:);

fprintf("Metadata sensors with coordinates: %d\n", height(P));

%% ================================================================
%  3. FIND TEMPERATURE COLUMNS IN MEASUREMENT DATA
% ================================================================
names = T.Properties.VariableNames;
tempMask = contains(names, "Temperature__C_") | contains(names, "AirTemperature__C_");
tempNames = names(tempMask);

%% ================================================================
%  4. BUILD SENSOR MAPPING TABLE
% ================================================================
map = table();
map.SensorID   = string.empty(0,1);
map.ColName    = string.empty(0,1);
map.X          = zeros(0,1);
map.Y          = zeros(0,1);
map.Z          = zeros(0,1);
map.SourceType = string.empty(0,1);

for i = 1:height(P)
    id = P.S_n(i);
    if id == ""
        continue
    end

    % ---------- Aranet (5‑char IDs) ----------
    if strlength(id) == 5
        hit = tempNames(contains(tempNames, lower(id)) & contains(tempNames,"Temperature__C_"));
        if isempty(hit)
            fprintf("No Aranet temperature column found for ID %s\n", id);
            continue
        end
        for h = hit'
            colName = string(h);
            newRow = {
                string(id), ...
                colName, ...
                P.DistanceFromPath(i), ...
                P.RowNumber(i), ...
                P.Height_m_(i), ...
                "Aranet" ...
            };
            map = [map; newRow];
        end
    end

    % ---------- Custom Zigbee (16‑char IDs) ----------
    if strlength(id) == 16
        % last 3 hex chars as suffix (e.g. 7B5, BE9, 3BA, 3D5, 7F4)
        suffix = lower(extractAfter(id, strlength(id)-3));

        % Prefer AirTemperature__C_, fallback to Temperature__C_
        hitAir = tempNames(contains(tempNames, "sensor_" + suffix) & ...
                           contains(tempNames, "AirTemperature__C_"));
        hitTemp = tempNames(contains(tempNames, "sensor_" + suffix) & ...
                            contains(tempNames, "Temperature__C_"));

        if ~isempty(hitAir)
            hit = hitAir(1);
        elseif ~isempty(hitTemp)
            hit = hitTemp(1);
        else
            fprintf("No custom temperature column found for ID %s (suffix %s)\n", id, suffix);
            continue
        end

        colName = string(hit);
        newRow = {
            string(id), ...
            colName, ...
            P.DistanceFromPath(i), ...
            P.RowNumber(i), ...
            P.Height_m_(i), ...
            "Custom" ...
        };
        map = [map; newRow];
    end
end

disp("Mapped sensors:");
disp(map)

if isempty(map)
    error("No sensors were successfully mapped. Check ID patterns and column names.");
end

%% ================================================================
%  5. EXTRACT TEMPERATURE MATRIX FOR MAPPED SENSORS
% ================================================================
Tmat = T{:, map.ColName};   % Ntime × Nsensors

%% ================================================================
%  6. CHOOSE TIME INDEX FOR ANALYSIS
% ================================================================
t0 = 1000;   % choose a time index with good coverage
if t0 > size(Tmat,1)
    error("t0 (%d) exceeds number of time steps (%d).", t0, size(Tmat,1));
end

Trow = Tmat(t0,:);

valid = isfinite(Trow);
yPos = map.Y(valid);          % use RowNumber as spatial axis
xPos = map.X(valid);          % DistanceFromPath (constant or nearly)
zPos = map.Z(valid);
Tv   = Trow(valid);

fprintf("Valid sensors at t=%d: %d\n", t0, numel(Tv));

if numel(Tv) < 3
    error("Not enough valid sensors at t=%d for interpolation.", t0);
end

%% ================================================================
%  7. 1D INTERPOLATION ALONG ROW DIRECTION (Y) WITH DUPLICATE MERGE
% ================================================================
% Sort by Y
[ySorted, idxSort] = sort(yPos);
Tsorted = Tv(idxSort);

% Merge duplicate Y positions by averaging temperatures
[uniqueY, ~, ic] = unique(ySorted);
Tmerged = accumarray(ic, Tsorted, [], @mean);

ySorted = uniqueY;
Tsorted = Tmerged;

yq = linspace(min(ySorted), max(ySorted), 100)';  % query positions

% Linear interpolation
T_lin = interp1(ySorted, Tsorted, yq, 'linear', 'extrap');

% Spline interpolation
T_spl = interp1(ySorted, Tsorted, yq, 'spline', 'extrap');

% Kriging / GPR in 1D (Y only)
gprMdl = fitrgp(ySorted, Tsorted, 'KernelFunction','squaredexponential');
T_krig = predict(gprMdl, yq);

%% ================================================================
%  8. BUILD PSEUDO‑2D GRID FOR VISUALIZATION
% ================================================================
% Extrude the 1D profile along a fake X‑axis to get a 2D map
xq = [0 1];                      % two columns, arbitrary width
[XI, YI] = meshgrid(xq, yq);

ZI_linear = repmat(T_lin, 1, numel(xq));
ZI_spline = repmat(T_spl, 1, numel(xq));
ZI_krig   = repmat(T_krig, 1, numel(xq));

%% ================================================================
%  9. VISUALIZATION — SPATIAL MAPS (PSEUDO‑2D)
% ================================================================
figure;
subplot(1,3,1)
contourf(XI, YI, ZI_linear, 20, 'LineColor','none'); colorbar
title('Linear interpolation (along row)')
hold on; scatter(zeros(size(ySorted)), ySorted, 60, Tsorted, 'filled', 'k')
ylabel('RowNumber')

subplot(1,3,2)
contourf(XI, YI, ZI_spline, 20, 'LineColor','none'); colorbar
title('Spline interpolation (along row)')
hold on; scatter(zeros(size(ySorted)), ySorted, 60, Tsorted, 'filled', 'k')

subplot(1,3,3)
contourf(XI, YI, ZI_krig, 20, 'LineColor','none'); colorbar
title('Kriging (GPR) (along row)')
hold on; scatter(zeros(size(ySorted)), ySorted, 60, Tsorted, 'filled', 'k')

sgtitle(sprintf('Temperature interpolation along row direction at t = %d', t0))

%% ================================================================
%  10. COMPARISON AGAINST ACTUAL SENSOR VALUES
% ================================================================
% Interpolate back at sensor Y positions
T_lin_at_sensors  = interp1(yq, T_lin,  ySorted, 'linear');
T_spl_at_sensors  = interp1(yq, T_spl,  ySorted, 'linear');
T_krig_at_sensors = interp1(yq, T_krig, ySorted, 'linear');

err_lin  = T_lin_at_sensors  - Tsorted;
err_spl  = T_spl_at_sensors  - Tsorted;
err_krig = T_krig_at_sensors - Tsorted;

figure;
subplot(1,3,1)
scatter(Tsorted, T_lin_at_sensors); hold on; plot(xlim,xlim,'k--')
title(sprintf('Linear (RMSE=%.3f)', sqrt(mean(err_lin.^2))))
xlabel('Measured'); ylabel('Interpolated')

subplot(1,3,2)
scatter(Tsorted, T_spl_at_sensors); hold on; plot(xlim,xlim,'k--')
title(sprintf('Spline (RMSE=%.3f)', sqrt(mean(err_spl.^2))))
xlabel('Measured'); ylabel('Interpolated')

subplot(1,3,3)
scatter(Tsorted, T_krig_at_sensors); hold on; plot(xlim,xlim,'k--')
title(sprintf('Kriging (RMSE=%.3f)', sqrt(mean(err_krig.^2))))
xlabel('Measured'); ylabel('Interpolated')

sgtitle('Interpolation accuracy at sensor locations (row direction)')

%% ================================================================
%  11. ROBUSTNESS TEST — DIFFERENT SENSOR CONFIGURATIONS
% ================================================================
% Use unique Y positions BEFORE selecting subsets
[uniqueY_all, ~, ic_all] = unique(ySorted);
Tunique_all = accumarray(ic_all, Tsorted, [], @mean);

% Now we have unique sensors
nUnique = numel(uniqueY_all);

configs = {
    1:nUnique, ...
    randperm(nUnique, max(3, round(nUnique*0.7))), ...
    randperm(nUnique, max(2, round(nUnique*0.4))) ...
};

figure;
for c = 1:numel(configs)
    idx = configs{c};
    y_c = uniqueY_all(idx);
    T_c = Tunique_all(idx);

    yq_c = linspace(min(uniqueY_all), max(uniqueY_all), 100)';
    T_lin_c = interp1(y_c, T_c, yq_c, 'linear', 'extrap');

    ZI_lin_c = repmat(T_lin_c, 1, numel(xq));

    subplot(1, numel(configs), c)
    contourf(XI, YI, ZI_lin_c, 20, 'LineColor','none'); colorbar
    hold on; scatter(zeros(size(y_c)), y_c, 60, T_c, 'filled', 'k')
    title(sprintf('Config: %d sensors', numel(idx)))
    ylabel('RowNumber')
end

sgtitle('Robustness: different sensor configurations (row direction)')

%% ================================================================
%  ROBUSTNESS: RMSE FOR ALL CONFIGURATIONS (4, 3, 2 sensors)
% ================================================================

% Unique sensors (after merging duplicates)
[uniqueY, ~, ic] = unique(ySorted);
Tunique = accumarray(ic, Tsorted, [], @mean);

nUnique = numel(uniqueY);
configSizes = [4, 3, 2];   % all possible sizes
RMSE_results = zeros(numel(configSizes), 1);

for k = 1:numel(configSizes)
    n = configSizes(k);

    % Random subset of size n
    idx = randperm(nUnique, n);
    y_c = uniqueY(idx);
    T_c = Tunique(idx);

    % Interpolate using linear method
    yq_c = linspace(min(uniqueY), max(uniqueY), 100)';
    T_lin_c = interp1(y_c, T_c, yq_c, 'linear', 'extrap');

    % Evaluate RMSE at the chosen sensor positions
    T_interp_at_sensors = interp1(yq_c, T_lin_c, y_c, 'linear');
    RMSE_results(k) = sqrt(mean((T_interp_at_sensors - T_c).^2));
end

% Display results
robustnessTable = table(configSizes', RMSE_results, ...
    'VariableNames', {'SensorCount', 'RMSE'});
disp(robustnessTable)

