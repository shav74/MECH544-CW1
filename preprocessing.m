%% section1
clc; clear;

fname = 'C:\Users\aruns\Documents\MATLAB\Pablo\data_raw\Take 2025-10-06 10.01.47 AM.csv';  % or .xlsx
[~,~,ext] = fileparts(fname);
if any(strcmpi(ext,{'.xlsx','.xls'}))
    T = readtable(fname, 'Range','A8', 'VariableNamingRule','preserve');
else
    opts = detectImportOptions(fname, 'VariableNamingRule','preserve');
    opts.VariableNamesLine = 8;
    opts.DataLines = [9 Inf];
    T = readtable(fname, opts);
end


%% section 2
% === minimal 3D with white lines + axis labels ===
%T = take1;              % <-- your table variable
scale = 0.001;              % 0.001 if mm, 0.01 if cm, 1 if m

vars = string(T.Properties.VariableNames);
startCol = find(vars=="Marker1",1); if isempty(startCol), startCol = 4; end
nTriplets = floor((width(T) - startCol + 1)/3);

figure('Color','k'); hold on                 % black background (optional)
for k = 0:nTriplets-1
    % safely convert all types to numeric
    Xr = T{:, startCol + 3*k     }; X = scale * str2double(string(Xr));
    Yr = T{:, startCol + 3*k + 1 }; Y = scale * str2double(string(Yr));
    Zr = T{:, startCol + 3*k + 2 }; Z = scale * str2double(string(Zr));

    good = isfinite(X) & isfinite(Y) & isfinite(Z);
    if any(good)
        plot3(X(good), Y(good), Z(good), '-', 'Color','w', 'LineWidth',1.2);
    end
end
view(3); box on; grid on
xlabel('X'); ylabel('Y'); zlabel('Z');
hold off

%% section 3

clc;clear;

%% Build a struct of tables: Takes.<file_basename>
rawDir = fullfile(pwd, 'data_raw');
L = dir(rawDir);
L = L(~[L.isdir]);                                    % files only
L = L(endsWith({L.name},{'.csv','.CSV','.xlsx','.XLSX','.xls','.XLS'}));

Takes = struct();
for i = 1:numel(L)
    fpath = fullfile(L(i).folder, L(i).name);
    [~, base, ext] = fileparts(fpath);
    field = matlab.lang.makeValidName(base);          % safe struct field

    if any(strcmpi(ext,{'.xlsx','.xls'}))
        T = readtable(fpath, 'Range','A8', 'VariableNamingRule','preserve');
    else
        opts = detectImportOptions(fpath, 'VariableNamingRule','preserve');
        opts.VariableNamesLine = 8;       % row 8 = headers
        opts.DataLines         = [9 Inf]; % data from row 9
        T = readtable(fpath, opts);
    end

    Takes.(field) = T;
end

names = fieldnames(Takes);  % ← quick check of what you loaded

% Set scale if needed (e.g., 0.001 for mm → m; 0.01 for cm → m)
scale = 1;

figure('Color','k'); hold on
for i = 1:numel(names)
    T = Takes.(names{i});

    vars = string(T.Properties.VariableNames);
    % assume columns are Marker1, Marker2, ... in X,Y,Z triplets
    startCol = find(vars=="Marker1",1); if isempty(startCol), startCol = 4; end
    nTriplets = floor((width(T) - startCol + 1)/3);

    for k = 0:nTriplets-1
        toNum = @(v) str2double(string(v));

            Xv = toNum(T{:, startCol + 3*k     });
            Yv = toNum(T{:, startCol + 3*k + 1 });
            Zv = toNum(T{:, startCol + 3*k + 2 });
            
            X = scale .* Xv;
            Y = scale .* Yv;
            Z = scale .* Zv;
            
            good = isfinite(X) & isfinite(Y) & isfinite(Z);
            plot3(X(good), Y(good), Z(good), '-', 'Color','w', 'LineWidth', 1.2);

        end
end

view(3); box on; grid on;
xlabel('X'); ylabel('Y'); zlabel('Z');
title('all data points in 3D')
hold off


