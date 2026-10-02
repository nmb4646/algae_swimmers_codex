clear; clc; close all;

um_per_px = 0.48; 
fps = 15;            
dt = 1 / fps;         
summaryResults = []; % Initialize array for velocities + D

for i = 1:50
    % This searches for any file named "1", "1.xlsx", "1.csv", etc.
    searchPattern = dir(fullfile(pwd, [num2str(i) '*']));
    
    if isempty(searchPattern)
        fprintf('Loop %d: No file starting with "%d" found in %s\n', i, i, pwd);
        continue;
    end
    
    % Use the first match found
    actualFilename = searchPattern(1).name;
    
    try
        % Read the file
        data = readtable(actualFilename, 'FileType', 'spreadsheet');
        
        % Calculate coordinates and velocity
        x_um = data.X * um_per_px;
        y_um = data.Y * um_per_px;
        dist = sqrt(diff(x_um).^2 + diff(y_um).^2);
        avg_vel = mean(dist) / dt;
        
        % ===== ADD: Diffusion coefficient =====
        N = length(x_um);
        maxLag = floor(N/4);
        MSD = zeros(maxLag,1);
        tau = (1:maxLag)' * dt;
        
        for lag = 1:maxLag
            dx = x_um(1+lag:end) - x_um(1:end-lag);
            dy = y_um(1+lag:end) - y_um(1:end-lag);
            MSD(lag) = mean(dx.^2 + dy.^2);
        end
        
        % Linear fit (early region)
        fitRange = 1:round(maxLag/4);
        p = polyfit(tau(fitRange), MSD(fitRange), 1);
        slope = p(1);
        D = slope / 4; % 2D diffusion
        % ===== END ADD =====
        
        % Store for final table
        summaryResults = [summaryResults; i, avg_vel, D];
        
        % Create and show the figure
        figure('Name', ['File ' actualFilename], 'NumberTitle', 'off');
        plot(x_um, y_um, '-ob', 'MarkerFaceColor', 'b');
        axis equal; grid on;
        xlabel('X (\mum)'); ylabel('Y (\mum)');
        title(sprintf('Trajectory %s | Vel: %.2f \\mum/s | D: %.4f \\mum^2/s', ...
            actualFilename, avg_vel, D));
        
        fprintf('✅ Processed: %s | Vel = %.2f um/s | D = %.4f um^2/s\n', ...
            actualFilename, avg_vel, D);
        
    catch ME
        fprintf('❌ Error reading %s: %s\n', actualFilename, ME.message);
    end
end

% Create the summary table at the end
if ~isempty(summaryResults)
    FinalTable = table(summaryResults(:,1), summaryResults(:,2), summaryResults(:,3), ...
        'VariableNames', {'File_ID', 'Avg_Velocity_um_s', 'Diffusion_um2_s'});
    disp(FinalTable);
else
    disp('No data was processed. Check if your column headers are exactly "X" and "Y".');
end