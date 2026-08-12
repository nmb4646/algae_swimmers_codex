% File name
clear; clc; close all;

blue = [.06,.43,.74];

dir_w = './data/algae_motion/buffer_micro/';
dir_m = './data/algae_motion/mucus_micro/';


n = 19; dt = 1/16.5; pix = .24e-6;


filename = dir_w + sprintf("Microalgae in water (%i).csv",n);
warning('off','all');
T = readtable(filename);
warning('on','all');
% Extract X and Y columns (make sure names match exactly)
maxT = length(T.X); disp("Timesteps: " + maxT)
x = T.X; x = x-x(1);x = x*pix;%x = [x;x+x(end)];x = [x;x+x(end)];x = [x;x+x(end)];
y = T.Y; y = y-y(1);y = y*pix; %y = [y;y+y(end)];y = [y;y+y(end)];y = [y;y+y(end)];

%[x,y] = repeat_traj_rot4(x,y);[x,y] = repeat_traj_rot4(x,y);


% Combine into a 2xN matrix
A = [x, y];   % transpose to make row vectors
%smoothing
go = 3;gw=11;
x = sgolayfilt(x,go,gw); y = sgolayfilt(y,go,gw);
dx = gradient(x)/dt;dy = gradient(y)/dt;
ddx = gradient(dx)/dt;ddy = gradient(dy)/dt;
speed = sqrt(dot([dx,dy],[dx,dy],2));
theta = atan2(dy, dx);


% unwrap to avoid jumps at ±pi
theta_unwrapped = unwrap(theta);

% total turning angle
total_turn = theta_unwrapped(end) - theta_unwrapped(1);

% winding number
wn = total_turn / (2*pi);

% mean rotational velocity

omega = total_turn/(length(x)*dt);

% mean speed

v = mean(speed);



% calculate diffusion from analytical expression
kb = 1.38e-23;
Temp = 293;
mu = .001;
a = 10e-6;

Do = kb*Temp/(6*pi*mu*a);
Dr = kb*Temp/(8*pi*mu*(a^3));

Dl = Do + (v^2)/(4*Dr); % rotationally diffusive motor 
%Dl = Do + (v^2)/(2*Dr*(1+ (omega/Dr)^2)); % circle-swimming motor (assume constant speed, spherical shape)

disp("Analytical value: " + Dl)


time = (1:20*length(x))*dt;

[D, alpha, msd, tau]=compute_msd_log(x,y,dt,maxT*.5,1);
%scatter(tau,msd)


plot(time,Do*time + (time.^2)*(v^2),"--","DisplayName","Short-time ballistic motion","LineWidth",2);
xline(maxT*dt)

%xlim([0,3]);ylim([0 3e-9])
set(gca,"XScale","log")
set(gca,"YScale","log")