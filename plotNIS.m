function plotNIS(csv_dir)
files = dir(strcat(csv_dir, '/*.csv'));
varNames = {'std_a', 'std_y', 'nq_lidar', 'nq_radar', 'pass_lidar', ...
            'pass_radar'};
varTypes = {'double','double','double','double','logical','logical'};
results = table('Size', [length(files), 6], 'VariableNames', varNames,...
                'VariableTypes', varTypes);
for i  = 1:length(files)
    file = files(i);
    filename = strcat(file.folder,'\', file.name);
    data = readData(filename);
    lidarNIS = table2array(data(data.sensorType=='LASER', {'nis'}));
    radarNIS = table2array(data(data.sensorType=='RADAR', {'nis'}));
    params = sscanf(file.name, 'a%f_%f.csv');
    makePlots(lidarNIS, radarNIS, params);
    
    % Test hyptothesis using method described on p20:
    % http://www.robots.ox.ac.uk/~ian/Teaching/Estimation/LectureNotes2.pdf
     
    % Determine if pass at p=0.05
    numLidar = length(lidarNIS);
    nqLidar = numLidar * mean(lidarNIS) / 2.0;
    passLidar = chi2inv(0.025, numLidar) <= nqLidar & ...
                nqLidar <= chi2inv(0.975, numLidar);
    numRadar = length(radarNIS);
    nqRadar = numRadar * mean(radarNIS) / 3.0;
    passRadar = chi2inv(0.025, numRadar) <= nqRadar & ...
                nqRadar <= chi2inv(0.975, numRadar);
    results(i,:) = {params(1), params(2), nqLidar, nqRadar, passLidar, ...
                    passRadar};
end
disp(results)
end
function makePlots(lidarNIS, radarNIS, params)
figure;
subplot(2,1,1);
histogram(lidarNIS, 'Normalization', 'probability')
x = 0:0.2:15;
y = chi2pdf(x, 2);
hold on;
plot(x,y);
title(['Lidar, std\_a: ', ...
       num2str(params(1)), ...
       ' std\_yawdd: ', ...
       num2str(params(2))])
hold off;
subplot(2,1,2);
histogram(radarNIS, 'Normalization', 'probability')
x = 0:0.2:15;
y = chi2pdf(x, 3);
hold on;
plot(x,y);
title(['Radar, std\_a: ', ...
       num2str(params(1)), ...
       ' std\_yawdd: ', ...
       num2str(params(2))])
end

function data = readData(filename)
formatSpec = '%C%f';
data = readtable(filename, 'Format', formatSpec);
data.Properties.VariableNames = {'sensorType' 'nis'};
end