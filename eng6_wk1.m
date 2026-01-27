clear; clc; close all;

%% 1. Load data
T = readtable('greenhouse_data.xlsx.xlsx');

%% 2. Convert time column
if ~isdatetime(T.Time)
    T.Time = datetime(T.Time,'InputFormat','dd-MM-yyyy HH:mm:ss');
end
time = T.Time;

%% 3. Sampling time
Ts = 5;  % 5 minutes (numeric for iddata)

%% 4. List of signals (output + inputs)
var_names = { ...
    'Meetbox_Temperature_Temperature___C__Averages', ...                
    'x1_Kas_Verwarmingstemperatuur_Setpoint_Temperature___C__Averages', ... 
    'x1_Kas_VentilatieLuweZijde_VentLee_Percentage____Averages', ...        
    'x1_Kas_VentilatieWindzijde_VentWind_Percentage____Averages', ...       
    'x1_Kas_Doek_EnergyCurtain_Percentage____Averages', ...                 
    'Gr_1_Meteo_OutsideTemperature_Temperature___C__Averages', ...         
    'WeatherDataTomatoWorld_GlobalRadiation_W_m2__Averages', ...            
    'Gr_1_Meteo_OutsideWindSpeed_WindSpeed_m_s__Averages', ...              
    'Gr_1_Meteo_OutsideRH_Humidity____Averages' ...                          
    };

%% 5. Convert all variables to numeric
data_numeric = zeros(height(T), length(var_names));

for k = 1:length(var_names)
    col = T.(var_names{k});
    if iscell(col)
        data_numeric(:,k) = str2double(col); % convert cell → numeric
    else
        data_numeric(:,k) = col;            % already numeric
    end
end

% Assign output and input
Y = data_numeric(:,1);        % inside temperature
U = data_numeric(:,2:end);    % inputs (8 signals)

%% 6. Remove rows with NaNs
idx = all(~isnan([Y U]),2);
Y = Y(idx);
U = U(idx,:);

%% 7. Preprocess: detrend and normalize numeric signals
Y = detrend(Y);
U = detrend(U);

% Normalize to zero mean, unit variance
Y = (Y - mean(Y)) / std(Y);
U = (U - mean(U)) ./ std(U);

%% 8. Create iddata object
data = iddata(Y, U, Ts);
data.OutputName = {'InsideTemperature'};
data.InputName  = {'Heating','VentLee','VentWind','Curtain','OutsideTemp','OutsideRad','OutsideWind','OutsideRH'};
data.TimeUnit = 'minutes';

%% 9. Plot for sanity check
figure;
plot(data);
title('Greenhouse Identification Data');

%% 10. Excitation check
disp('Standard deviation of inputs:')
disp(std(U));

%{
%% 11. Dead-time estimate example (Heating -> Temp)
[xc,lags] = xcorr(Y-mean(Y), U(:,1)-mean(U(:,1)), 40, 'coeff');
[~,i] = max(xc);
deadtime_samples = lags(i);
deadtime_minutes = deadtime_samples * Ts; % using numeric Ts
disp(['Estimated dead time (Heating -> Temp): ', num2str(deadtime_minutes),' minutes']);
%}

% 12. Dead time estimate for the rest of the variables
maxLag = 40;  % up to 200 minutes
nu = size(U,2);
nk = zeros(1,nu);

figure;
for i = 1:nu
    [xc,lags] = xcorr(Y, U(:,i), maxLag, 'coeff');
    
    subplot(4,2,i)
    stem(lags, xc,'filled')
    title(['xcorr: ', data.InputName{i}])
    xlabel('Lag (samples)')
    ylabel('Correlation')
    
    % Use positive lags only
    posIdx = lags >= 0;
    [~,idx] = max(xc(posIdx));
    nk(i) = lags(posIdx);
    nk(i) = nk(i)(idx);
end

disp('Estimated dead-times (samples):')
disp(array2table(nk,'VariableNames',data.InputName))

disp('Estimated dead-times (minutes):')
disp(array2table(nk*Ts,'VariableNames',data.InputName))

%% 13. Split data into train/validation/test
N = length(Y);
data_train = data(1:round(0.6*N));
data_val   = data(round(0.6*N)+1:round(0.8*N));
data_test  = data(round(0.8*N)+1:end);

disp('Week 1 data setup complete. Training, validation, and test sets created.');

%% 14. Implementing ARX
na = 2;
nb = repmat(2,1,nu);

sys_arx = arx(data_train, [na nb nk]);

figure;
compare(data_val, sys_arx);
title('ARX – Validation Data');
