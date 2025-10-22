%% get files from data dir and assign them to a struct
rawDir = fullfile(pwd, 'data_raw');
L = dir(rawDir);
L = L(~[L.isdir]);
L = L(endsWith({L.name},{'.csv','.CSV','.xlsx','.XLSX','.xls','.XLS'}));

Takes = struct();
for i = 1:numel(L)
    fpath = fullfile(L(i).folder, L(i).name);
    [~, base, ext] = fileparts(fpath);
    field = matlab.lang.makeValidName(base);

    %data starts from cell A8 in every table
    if any(strcmpi(ext,{'.xlsx','.xls'}))
        T = readtable(fpath, 'Range','A8', 'VariableNamingRule','preserve');
        % handle NaNs
        T = standardizeMissing(T, {'','NaN','nan'});
        for v = 1:width(T)
            if ~isnumeric(T.(v))
                T.(v) = str2double(string(T.(v)));
            end
        end
    else
        opts = detectImportOptions(fpath, 'VariableNamingRule','preserve');
        opts.VariableNamesLine = 8;          % headers row
        opts.DataLines         = [9 Inf];    % data rows
        opts = setvaropts(opts, opts.VariableNames, 'TreatAsMissing', {'','NaN','nan'});
        opts = setvartype(opts, opts.VariableNames, 'double');  % force numeric
        T    = readtable(fpath, opts);
    end

    Takes.(field) = T;
end

names = fieldnames(Takes);

%% basic 3D plotting

scale = 0.001;   % to make mm into meters

% helper func
toNum = @(v) (isnumeric(v) .* v) + (~isnumeric(v)) .* str2double(string(v));

figure('Color','k'); hold on
for i = 1:numel(names)
    T = Takes.(names{i});
    vars = string(T.Properties.VariableNames);

    % Marker1 starts triples (x,y,x values to loop trhough)
    startCol = find(vars=="Marker1",1); if isempty(startCol), startCol = 4; end
    nTriplets = floor((width(T) - startCol + 1)/3);

    for k = 0:nTriplets-1
        Xv = toNum(T{:, startCol + 3*k     });
        Yv = toNum(T{:, startCol + 3*k + 1 });
        Zv = toNum(T{:, startCol + 3*k + 2 });

        % scaling data to meters
        X = scale .* Xv;
        Y = scale .* Yv;
        Z = scale .* Zv;

        % swapping the axis to make Z vertical
        X_= X;
        X = Z;
        Z = X_;

        % filling missing data with fillmissing functionm
        X = fillmissing(X, 'nearest');
        Y = fillmissing(Y, "nearest");
        Z = fillmissing(Z, "nearest");

        % smoothing using a moving average filter with window sizes
        windowSize = 50;
        X = smoothdata(X, 'movmean', windowSize);
        Y = smoothdata(Y, 'movmean', windowSize);
        Z = smoothdata(Z, 'movmean', windowSize);


        % filtering out the remaining incomplete paths

        min_length = 4.5; % minimum distance travelled have to be 4.5m
        floor_height = 1.5; % to filter out noice on the floor

        % keep only real flights witht the above filters
        keepThisMarker = (max(Z) > floor_height) && ((max(Y) - min(Y)) > min_length);
        if ~keepThisMarker, continue; end

        good = isfinite(Y) & isfinite(X) & isfinite(Z);
        
        if any(good)
            plot3(X(good), Y(good), Z(good), '-', 'LineWidth', 1.8);
        end
    end
end
view(3); box on; grid on
xlabel('X(m)'); ylabel('Y(m)'); zlabel('Z(m)');
title('Filtered and noice removed results')
hold off


%% speed analysis

fs_default   = 120;  
min_length   = 4.5;       
floor_height = 1.5;        
scale        = 0.001;  % same constrains as before  

results = table();          % will collect per-take, per-marker measurements

for i = 1:numel(names)
    T = Takes.(names{i});
    vars = string(T.Properties.VariableNames);
    
    fs = fs_default;

    % --- marker triplets start (as in your code) ---
    startCol  = find(vars=="Marker1",1); if isempty(startCol), startCol = 4; end
    nTriplets = floor((width(T) - startCol + 1)/3);

    for k = 0:nTriplets-1
        % raw columns
        Xv = T{:, startCol + 3*k     };
        Yv = T{:, startCol + 3*k + 1 };
        Zv = T{:, startCol + 3*k + 2 };
        if ~isnumeric(Xv), Xv = str2double(string(Xv)); end
        if ~isnumeric(Yv), Yv = str2double(string(Yv)); end
        if ~isnumeric(Zv), Zv = str2double(string(Zv)); end

        % scaling and swapping axises
        X = scale.*Zv;  Y = scale.*Yv;  Z = scale.*Xv;

        % simple fill + smooth
        X = fillmissing(X,'nearest');  Y = fillmissing(Y,'nearest');  Z = fillmissing(Z,'nearest');
        X = smoothdata(X,'movmean',50); Y = smoothdata(Y,'movmean',50); Z = smoothdata(Z,'movmean',50);

        % keep only valid flights (same rules)
        keepThisMarker = (max(Z) > floor_height) && ((max(Y) - min(Y)) > min_length);
        if ~keepThisMarker, continue; end

        % valid segment mask
        good = isfinite(X) & isfinite(Y) & isfinite(Z);

        % horozontal distance
        % calculating the arc length on XY axis.
        Xg = X(good); Yg = Y(good);
        if numel(Xg) < 2, continue; end

        dXY = hypot(diff(Xg), diff(Yg)); % per-step horizontal distance
        dist_horiz = sum(dXY); % total horizontal distance (m)
        duration = length(Y)/fs_default; % seconds
        if duration <= 0, continue; end

        avg_speed  = dist_horiz / duration;

        % store a row
        results = [results; table( ...
            string(names{i}), k+1, duration, dist_horiz, avg_speed, ...
            'VariableNames', {'take','markerTriplet','duration_s','distance_m','avg_speed_mps'})]; %#ok<AGROW>
    end
end

disp('--- Speed summary (horizontal XZ) ---');
disp(results);

%saving the table for future analysis
if ~exist('results','dir'), mkdir('results'); end
writetable(results, fullfile('results','speed_summary.csv'));


%% removing arm swings
fs_default   = 120;         
min_length   = 4.5;         
floor_height = 1.5;       
scale        = 0.001;      
windowSize   = 50;          % same metrics

results = table(); % avg speed only

for i = 1:numel(names)
    T = Takes.(names{i});
    vars = string(T.Properties.VariableNames);

    startCol  = find(vars=="Marker1",1); if isempty(startCol), startCol = 4; end
    nTriplets = floor((width(T) - startCol + 1)/3);

    for k = 0:nTriplets-1
        % raw columns -> numeric
        Xv = T{:, startCol + 3*k     };
        Yv = T{:, startCol + 3*k + 1 };
        Zv = T{:, startCol + 3*k + 2 };
        if ~isnumeric(Xv), Xv = str2double(string(Xv)); end
        if ~isnumeric(Yv), Yv = str2double(string(Yv)); end
        if ~isnumeric(Zv), Zv = str2double(string(Zv)); end

        % scale + swap 
        X = scale.*Zv;  Y = scale.*Yv;  Z = scale.*Xv;

        % simple fill as before
        X = fillmissing(X,'nearest');  Y = fillmissing(Y,'nearest');  Z = fillmissing(Z,'nearest');
        X = smoothdata(X,'movmean',windowSize);
        Y = smoothdata(Y,'movmean',windowSize);
        Z = smoothdata(Z,'movmean',windowSize);

        keepThisMarker = (max(Z) > floor_height) && ((max(Y) - min(Y)) > min_length);
        if ~keepThisMarker, continue; end

        good = isfinite(X) & isfinite(Y) & isfinite(Z);
        if ~any(good), continue; end

        % remove arm-swing via simple speed gate
        % horizontal speed on the plane you plot (X–Y), uniform sampling at fs_default
        vx = [0; diff(X)] * fs_default;
        vy = [0; diff(Y)] * fs_default;
        vH = sqrt(vx.^2 + vy.^2);

        % only starting to plot when v > 20% v max (max speed)
        vmax   = max(vH(good),[],'omitnan');
        if ~isfinite(vmax) || vmax == 0, continue; end
        enterFrac = 0.20;   % start when >20% vmax
        exitFrac  = 0.12;   % end when >12% vmax
        sIdx = find(vH >= enterFrac*vmax, 1, 'first');
        eIdx = find(vH >= exitFrac *vmax,  1, 'last');
        if isempty(sIdx) || isempty(eIdx) || eIdx <= sIdx, continue; end

        % restrict to that in-flight window and to finite samples
        segMask = false(size(good)); segMask(sIdx:eIdx) = true;
        keep    = good & segMask;

        if nnz(keep) < 2, continue; end

        Xg = X(keep);  Yg = Y(keep);
        % horizontal distance 
        dist_m    = sum(hypot(diff(Xg), diff(Yg)));
        duration_s= nnz(keep) / fs_default;
        if duration_s <= 0, continue; end

        avg_speed = dist_m / duration_s;

        results = [results; table( ...
            string(names{i}), k+1, avg_speed, ...
            'VariableNames', {'take','markerTriplet','avg_speed_mps'})]; %#ok<AGROW>
    end
end

disp('--- Average speed (visible in-flight segments only) ---');
disp(results);

% optional save again
if ~exist('results','dir'), mkdir('results'); end
writetable(results, fullfile('results','avg_speed_visible.csv'));


%% plotting wothout the arm swings
fs_default   = 120;    
scale        = 0.001;  
windowSize   = 50;      
min_length   = 4.5;    
floor_height = 1.5;    

% Speed gate
enterFrac = 0.20;       % start when horizontal speed > 20% of vmax
exitFrac  = 0.12;       % end when speed > 12% of vmax

results = table();      % will hold: take, markers, duration_s and avg_speed_mps

toNum = @(v) (isnumeric(v) .* v) + (~isnumeric(v)) .* str2double(string(v));

figure('Color','k'); hold on
for i = 1:numel(names)
    T = Takes.(names{i});
    vars = string(T.Properties.VariableNames);

    % Time make it 120 fps if time col is not present
    tIdx = find(contains(lower(vars),'time'), 1);
    if ~isempty(tIdx)
        t = T.(vars(tIdx));
    else
        t = (0:height(T)-1)'/fs_default;
    end

    % Marker triplets (assume Marker1 starts)
    startCol  = find(vars=="Marker1",1); if isempty(startCol), startCol = 4; end
    nTriplets = floor((width(T) - startCol + 1)/3);

    for k = 0:nTriplets-1
        
        Xv = toNum(T{:, startCol + 3*k     });
        Yv = toNum(T{:, startCol + 3*k + 1 });
        Zv = toNum(T{:, startCol + 3*k + 2 });

        % --- scale + swap 
        X = scale.*Zv;   Y = scale.*Yv;   Z = scale.*Xv;

        X = fillmissing(X,'nearest');   Y = fillmissing(Y,'nearest');   Z = fillmissing(Z,'nearest');
        X = smoothdata(X,'movmean',windowSize);
        Y = smoothdata(Y,'movmean',windowSize);
        Z = smoothdata(Z,'movmean',windowSize);

        % --- keep only "real" flights (your two filters)
        keepThisMarker = (max(Z) > floor_height) && ((max(Y) - min(Y)) > min_length);
        if ~keepThisMarker, continue; end

        good = isfinite(X) & isfinite(Y) & isfinite(Z);
        if ~any(good), continue; end

        %remove arm swing via a simple speed gate
        vx = [0; diff(X)] * fs_default;
        vy = [0; diff(Y)] * fs_default;
        vH = sqrt(vx.^2 + vy.^2);

        vmax   = max(vH(good),[],'omitnan');
        if ~isfinite(vmax) || vmax==0, continue; end

        sIdx = find(vH >= enterFrac*vmax, 1, 'first');
        eIdx = find(vH >= exitFrac *vmax, 1, 'last');
        if isempty(sIdx) || isempty(eIdx) || eIdx <= sIdx, continue; end

        segMask = false(size(good)); segMask(sIdx:eIdx) = true;
        keep    = good & segMask;
        if nnz(keep) < 2, continue; end

        

        % visible in-flight paths
        Xg = X(keep);  Yg = Y(keep);  Zg = Z(keep);  tg = t(keep);

        % plotting the curve
        plot3(Xg, Yg, Zg, '-', 'LineWidth', 1.8);

        % add markers + label (use take name + marker index as the name)
        plot3(Xg(1),  Yg(1),  Zg(1),  'o', 'MarkerSize',4, 'MarkerFaceColor','w', 'MarkerEdgeColor','w');
        plot3(Xg(end),Yg(end),Zg(end),'s', 'MarkerSize',5, 'MarkerFaceColor','w', 'MarkerEdgeColor','w');
        lbl = sprintf('%s | m%d', names{i}, k+1);
        text(Xg(end), Yg(end), Zg(end), ['  ' lbl], 'Color','w', 'FontSize',8, 'Interpreter','none');

        % metrics: average speed + duration
        distance_m = sum(hypot(diff(Xg), diff(Yg)));   % horizontal on X–Y (your chart)
        duration_s = numel(Xg) / fs_default;              % uniform sampling
        if duration_s > 0
            avg_speed_mps = distance_m / duration_s;
            results = [results; table(string(names{i}), k+1, duration_s, distance_m, avg_speed_mps, ...
                'VariableNames', {'take','markerTriplet','duration_s','distance_m','avg_speed_mps'})]; %#ok<AGROW>
        end
    end
end

view(3); box on; grid on
ax = gca; ax.GridColor='w'; ax.XColor='w'; ax.YColor='w'; ax.ZColor='w';
xlabel('X (m)'); ylabel('Y(m)'); zlabel('Z(m)');
title('Flight Paths without arm swings')
hold off

disp(results);

% save for later analysis
if ~exist('results','dir'), mkdir('results'); end
writetable(results, fullfile('results','avg_speed_and_duration_visible.csv'));

%% visualizing and analizing the resutls with result table

figDir = fullfile('results','figs');
if ~exist(figDir,'dir'), mkdir(figDir); end

%% Figure 1: Average speed leaderboard
[vals, ord] = sort(results.avg_speed_mps, 'descend');
labs = cellstr(results.take(ord));
bar(vals); box on; grid on
ylabel('Average speed (m/s)'); title('Average speed by flight');
% marker names on the bottom in a 45 degree angle
set(gca,'XTick',1:numel(labs),'XTickLabel',labs,'XTickLabelRotation',45);


%% Figure 2: Distance leaderboard
[valsD, ordD] = sort(results.distance_m, 'descend');
labsD = cellstr(results.take(ordD));
bar(valsD); box on; grid on
ylabel('Distance (m)'); title('Distance by flight (in-flight segment)');
set(gca,'XTick',1:numel(labsD),'XTickLabel',labsD,'XTickLabelRotation',45);


%% figure 3 : distance leaderboard 

[valsT, ordT] = sort(results.duration_s, 'descend');
labsT = cellstr(results.take(ordT));
bar(valsT); box on; grid on
ylabel('Flight time (s)'); title('Flight time by flight');
set(gca,'XTick',1:numel(labsT),'XTickLabel',labsT,'XTickLabelRotation',45);

%% Speed vs Distance scatterplot

scatter(results.avg_speed_mps, results.distance_m, 60, 'filled'); grid on; box on
xlabel('Average speed (m/s)'); ylabel('Distance (m)'); title('Distance vs Average speed (in-flight segment)');

% gettign the range and avg speed for flights
range_speed = max(results.avg_speed_mps) - min(results.avg_speed_mps);
range_distance = max(results.distance_m) - min(results.distance_m);

dx = 0.01*range_speed; dy = 0.01*range_distance;
for i=1:height(results)
    text(results.avg_speed_mps(i)+dx, results.distance_m(i)+dy, results.take(i), 'FontSize',8, 'Interpreter','none');
end

%% Figure 5: Distance vs Flight time

scatter(results.duration_s, results.distance_m, 60, 'filled'); grid on; box on
xlabel('Flight time (s)'); ylabel('Distance (m)'); title('Distance vs Flight time (in-flight segment)');

range_distance = max(results.distance_m) - min(results.distance_m);
range_duration = max(results.duration_s) - min(results.duration_s);

dx = 0.01*range_duration; dy = 0.01*range_distance;
for i=1:height(results)
    text(results.duration_s(i)+dx, results.distance_m(i)+dy, results.take(i), 'FontSize',8, 'Interpreter','none');
end


%% logging fastest and most travelled flight
[~, iFast] = max(results.avg_speed_mps);
[~, iFar ] = max(results.distance_m);
fprintf('\nFastest flight: %s (%.2f m/s, %.2f m, %.2f s)\n', ...
    results.take(iFast), results.avg_speed_mps(iFast), ...
    results.distance_m(iFast), results.duration_s(iFast));
fprintf('Farthest flight: %s (%.2f m, %.2f m/s, %.2f s)\n\n', ...
    results.take(iFar), results.distance_m(iFar), ...
    results.avg_speed_mps(iFar), results.duration_s(iFar));


%% Histograms: avg speed, distance, duration
figDir = fullfile('results','figs');
if ~exist(figDir,'dir'), mkdir(figDir); end

% derive distance if its not already in the table
if ~ismember('distance_m', results.Properties.VariableNames)
    results.distance_m = results.avg_speed_mps .* results.duration_s;
end

% a simple bin count that works for small N
nbins = max(3, min(10, floor(sqrt(height(results)))));

% 1 Average speed histogram
fhs1 = figure;
histogram(results.avg_speed_mps, nbins);
grid on; box on
xlabel('Average speed (m/s)');
ylabel('Count');
title('Distribution of average speeds');
exportgraphics(fhs1, fullfile(figDir,'hist_avg_speed.png'), 'Resolution',300);

% 2) Distance histogram
fhs2 = figure;
histogram(results.distance_m, nbins);
grid on; box on
xlabel('Distance (m)');
ylabel('Count');
title('Distribution of distances');
% saving the graph
exportgraphics(fhs2, fullfile(figDir,'hist_distance.png'), 'Resolution',300);

% 3) Flight time histogram
fhs3 = figure;
histogram(results.duration_s, nbins);
grid on; box on
xlabel('Flight time (s)');
ylabel('Count');
title('Distribution of flight times');
exportgraphics(fhs3, fullfile(figDir,'hist_duration.png'), 'Resolution',300);

% quick console stats
fprintf('\n-- Summary stats --\n');
fprintf('Avg speed:  mean %.2f m/s,  sd %.2f\n', mean(results.avg_speed_mps,'omitnan'), std(results.avg_speed_mps,'omitnan'));
fprintf('Distance:   mean %.2f m,    sd %.2f\n', mean(results.distance_m,'omitnan'),   std(results.distance_m,'omitnan'));
fprintf('Duration:   mean %.2f s,    sd %.2f\n\n', mean(results.duration_s,'omitnan'),  std(results.duration_s,'omitnan'));



%% === Landing positions (ground projection) ===
% Uses: Takes, names, results (with .take, .markerTriplet, .avg_speed_mps...)
% Re-derives each visible in-flight segment to get the true landing point.
fs_default   = 120;
scale        = 0.001;
windowSize   = 50;
min_length   = 4.5;
floor_height = 1.5;
enterFrac    = 0.20;   % start when horiz speed > 20% vmax
exitFrac     = 0.12;   % end when > 12% vmax

toNum = @(v) (isnumeric(v).*v) + (~isnumeric(v)).*str2double(string(v));

% ensure cols exist in results
if ~ismember('x_end', results.Properties.VariableNames)
    results.x_end = nan(height(results),1);
end
if ~ismember('y_end', results.Properties.VariableNames)
    results.y_end = nan(height(results),1);
end

for r = 1:height(results)
    takeField = char(results.take(r));
    if ~isfield(Takes, takeField), continue; end

    T = Takes.(takeField);
    k = results.markerTriplet(r) - 1; 

    % time (same method)
    vars = string(T.Properties.VariableNames);
    tIdx = find(contains(lower(vars),'time'), 1);
    if ~isempty(tIdx), t = T.(vars(tIdx));
    else,              t = (0:height(T)-1)'/fs_default;
    end

    % raw triplet to numeric
    startCol = find(vars=="Marker1",1); if isempty(startCol), startCol = 4; end
    Xv = toNum(T{:, startCol + 3*k     });
    Yv = toNum(T{:, startCol + 3*k + 1 });
    Zv = toNum(T{:, startCol + 3*k + 2 });

    % scale + swap (your convention: XY = ground, Z = up)
    X = scale.*Zv;   Y = scale.*Yv;   Z = scale.*Xv;

    % same simple clean-up
    X = fillmissing(X,'nearest');  Y = fillmissing(Y,'nearest');  Z = fillmissing(Z,'nearest');
    X = smoothdata(X,'movmean',windowSize);
    Y = smoothdata(Y,'movmean',windowSize);
    Z = smoothdata(Z,'movmean',windowSize);

    % simple valid flights filters (same as your charts)
    keepMarker = (max(Z) > floor_height) && ((max(Y) - min(Y)) > min_length);
    if ~keepMarker, continue; end

    good = isfinite(X) & isfinite(Y) & isfinite(Z);
    if ~any(good), continue; end

    % speed gate (remove arm swing)
    vx = [0; diff(X)] * fs_default;
    vy = [0; diff(Y)] * fs_default;
    vH = sqrt(vx.^2 + vy.^2);

    vmax = max(vH(good),[],'omitnan');
    if ~isfinite(vmax) || vmax==0, continue; end

    sIdx = find(vH >= enterFrac*vmax, 1, 'first');
    eIdx = find(vH >= exitFrac *vmax, 1, 'last');
    if isempty(sIdx) || isempty(eIdx) || eIdx <= sIdx, continue; end

    segMask = false(size(good)); segMask(sIdx:eIdx) = true;
    keepSeg = good & segMask;
    if nnz(keepSeg) < 2, continue; end

    % landing location = last point of the visible in-flight segment
    Xg = X(keepSeg);  Yg = Y(keepSeg);
    results.x_end(r) = Xg(end);
    results.y_end(r) = Yg(end);
end

% plot the landing plots
figure;
scatter(results.x_end, results.y_end, 60, 'filled'); grid on; box on
axis equal
xlabel('Ground X (m)'); ylabel('Ground Y (m)');
title('Landing positions (ground projection)');

% label points with take names 
rngX = (max(results.x_end) - min(results.x_end)); if rngX==0, rngX=1; end
rngY = (max(results.y_end) - min(results.y_end)); if rngY==0, rngY=1; end
dx = 0.01*rngX; dy = 0.01*rngY;
for i = 1:height(results)
    if isfinite(results.x_end(i)) && isfinite(results.y_end(i))
        text(results.x_end(i)+dx, results.y_end(i)+dy, string(results.take(i)), 'Interpreter','none', 'FontSize',8);
    end
end

% optional save
figDir = fullfile('results','figs');
if ~exist(figDir,'dir'), mkdir(figDir); end
exportgraphics(gcf, fullfile(figDir,'landing_positions.png'), 'Resolution',300);
