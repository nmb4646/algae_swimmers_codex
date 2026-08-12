clear;close;clc;
%%% Modeling based on the analytical diffusion coefficients of 
% Ebbens et al. (2010) and
% Marine et al. (2013) derived from Langevin equations
% for motors with rotational diffusion (microalgae) and constant rotational velocity (picoalgae)
% -- main assumption: at long timescales, pico and microalgae behave diffusively. 
% -- this may be false for microalgae, which, if trapped, may be subdiffusive.
% -- for picoalgae, since the looping motion has no preferred drift orientation, there
% should be no ballistic bias at long timescales and the assumption of long-time diffusive transport seems rational.
%%%

avg_type = "geo";


%%% Picoalgae diffusion coefficients

dir_w = './data/algae_motion/buffer/'; Nw = 36;
dir_m = './data/algae_motion/mucus/';  Nm = 25;
dt = 1/100; pix = .32e-6;

filtertime = 100; % Number of timesteps duration below which we disregard sample

D_pw = []; % Diffusion coefficients for picoalgae in water
D_pm = []; % Diffusion coefficients for picoalgae in mucus

for n=1:Nw
    filename = dir_w + sprintf("%i.csv",n);
    warning('off','all');
    T = readtable(filename);
    warning('on','all');
    % Extract X and Y columns (make sure names match exactly)

    x = T.X; x = x-x(1);x = x*pix;
    y = T.Y; y = y-y(1);y = y*pix;
    
    if length(x)<filtertime
        continue
    end


    % Combine into a 2xN matrix
    A = [x, y];   % transpose to make row vectors
    %smoothing
    go = 3;gw=11;
    x = sgolayfilt(x,go,gw); y = sgolayfilt(y,go,gw);
    dx = gradient(x)/dt;dy = gradient(y)/dt;
    speed = sqrt(dot([dx,dy],[dx,dy],2));
    theta = atan2(dy, dx);

    % unwrap to avoid jumps at ±pi
    theta_unwrapped = unwrap(theta);
    
    % total turning angle
    total_turn = theta_unwrapped(end) - theta_unwrapped(1);
    
    % mean rotational velocity  
    omega = total_turn/(length(x)*dt);
  
    % mean speed
    v = mean(speed);


    D_pw = [D_pw; D2(1e-6,v,omega)];
end

for n=1:Nm
    filename = dir_m + sprintf("%i.csv",n);
    warning('off','all');
    T = readtable(filename);
    warning('on','all');
    % Extract X and Y columns (make sure names match exactly)
    maxT = length(T.X); 
    x = T.X; x = x-x(1);x = x*pix;
    y = T.Y; y = y-y(1);y = y*pix;

    if length(x)<filtertime
        continue
    end
    
    % Combine into a 2xN matrix
    A = [x, y];   % transpose to make row vectors
    %smoothing
    go = 3;gw=11;
    x = sgolayfilt(x,go,gw); y = sgolayfilt(y,go,gw);
    dx = gradient(x)/dt;dy = gradient(y)/dt;
    speed = sqrt(dot([dx,dy],[dx,dy],2));
    theta = atan2(dy, dx);

    % unwrap to avoid jumps at ±pi
    theta_unwrapped = unwrap(theta);
    
    % total turning angle
    total_turn = theta_unwrapped(end) - theta_unwrapped(1);
    
    % mean rotational velocity  
    omega = total_turn/(length(x)*dt);
  
    % mean speed
    v = mean(speed);


    D_pm = [D_pm; D2(1e-6,v,omega)];
end


[p,h,stats] = ranksum(D_pw,D_pm);
swarmchart(0*D_pw,D_pw,"DisplayName","pico, water");hold on;swarmchart(1+0*D_pm,D_pm,"DisplayName","pico, mucus")

if avg_type == "median"
disp("Median water diffusion coefficient, picoalgae: " + median(D_pw))
disp("Median mucus diffusion coefficient, picoalgae: " + median(D_pm))
elseif avg_type == "geo"
disp("Mean water diffusion coefficient, picoalgae: " + 10^mean(log10(D_pw))+" m^2/s"); 
disp("Mean mucus diffusion coefficient, picoalgae: " + 10^mean(log10(D_pm))+" m^2/s");
end
 
set(gca,'YScale','log')

%disp(mean(log10(D_pw)));disp(mean(log10(D_pm)))



%%% Microalgae diffusion coefficients

dir_w = './data/algae_motion/buffer_micro/'; Nw = 20;
dir_m = './data/algae_motion/mucus_micro/'; Nm = 20;
dt = 1/16.5; pix = .24e-6;

filtertime = 100; % Number of timesteps duration below which we disregard sample

D_mw = []; % Diffusion coefficients for picoalgae in water
D_mm = []; % Diffusion coefficients for picoalgae in mucus

for n=1:Nw
    filename = dir_w + sprintf("Microalgae in water (%i).csv",n);
    warning('off','all');
    T = readtable(filename);
    warning('on','all');
    % Extract X and Y columns (make sure names match exactly)

    x = T.X; x = x-x(1);x = x*pix;
    y = T.Y; y = y-y(1);y = y*pix;
    
    if length(x)<filtertime
        continue
    end


    % Combine into a 2xN matrix
    A = [x, y];   % transpose to make row vectors
    %smoothing
    go = 3;gw=11;
    x = sgolayfilt(x,go,gw); y = sgolayfilt(y,go,gw);
    dx = gradient(x)/dt;dy = gradient(y)/dt;
    speed = sqrt(dot([dx,dy],[dx,dy],2));
    theta = atan2(dy, dx);

    % unwrap to avoid jumps at ±pi
    theta_unwrapped = unwrap(theta);
    
    % total turning angle
    total_turn = theta_unwrapped(end) - theta_unwrapped(1);
    
    % mean rotational velocity  
    omega = total_turn/(length(x)*dt);
  
    % mean speed
    v = mean(speed);


    D_mw = [D_mw; D1(10e-6,v)];
end

for n=1:Nm
    filename = dir_m + sprintf("Microalgae in mucus (%i).csv",n);
    warning('off','all');
    T = readtable(filename);
    warning('on','all');
    % Extract X and Y columns (make sure names match exactly)
    maxT = length(T.X); 
    x = T.X; x = x-x(1);x = x*pix;
    y = T.Y; y = y-y(1);y = y*pix;

    if length(x)<filtertime
        continue
    end
    
    % Combine into a 2xN matrix
    A = [x, y];   % transpose to make row vectors
    %smoothing
    go = 3;gw=11;
    x = sgolayfilt(x,go,gw); y = sgolayfilt(y,go,gw);
    dx = gradient(x)/dt;dy = gradient(y)/dt;
    speed = sqrt(dot([dx,dy],[dx,dy],2));
    theta = atan2(dy, dx);

    % unwrap to avoid jumps at ±pi
    theta_unwrapped = unwrap(theta);
    
    % total turning angle
    total_turn = theta_unwrapped(end) - theta_unwrapped(1);
    
    % mean rotational velocity  
    omega = total_turn/(length(x)*dt);
  
    % mean speed
    v = mean(speed);


    D_mm = [D_mm; D1(10e-6,v)];
end


[p,h,stats] = ranksum(D_mw,D_mm);
swarmchart(0*D_mw,D_mw,"DisplayName","micro, water");hold on;swarmchart(1+0*D_mm,D_mm,"DisplayName","micro, mucus")

if avg_type == "median"
disp("Median water diffusion coefficient, microalgae: " + median(D_mw))
disp("Median mucus diffusion coefficient, microalgae: " + median(D_mm))
elseif avg_type == "geo"
disp("Mean water diffusion coefficient, microalgae: " + 10^mean(log10(D_mw))+" m^2/s")
disp("Mean mucus diffusion coefficient, microalgae: " + 10^mean(log10(D_mm))+" m^2/s")
end

legend;
set(gca,'YScale','log')
title("Diffusion Coefficients, swarm visualization")

% Bar chart
figure;


if avg_type == "median"
b = bar(["Picoalgae robot","Microalgae robot"],[median(D_pw)/median(D_pm),median(D_mw)/median(D_mm)]);
elseif avg_type == "geo"
b = bar(["Picoalgae robot","Microalgae robot"],[(10^mean(log10(D_pw)))/(10^mean(log10(D_pm))),(10^mean(log10(D_mw)))/(10^mean(log10(D_mm)))]);
end
b.FaceColor = 'flat';
b.CData(2,:) = [.7,.3,.1];

ylabel("Dw/Dm")
set(gca,"YScale","log")
title("Diffusion coefficients in different media")

set(gca,"FontSize",16)
