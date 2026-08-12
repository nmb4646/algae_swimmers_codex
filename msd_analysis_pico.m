% File name
clear; clc; close all;

blue = [.06,.43,.74];

dir_w = './data/algae_motion/buffer/';
dir_m = './data/algae_motion/mucus/';

use_mucus=true;



dt = 1/100; pix = .32e-6;
tiledlayout(4,1)

if use_mucus
    ns = [2,3,7,10];
else
    ns= [2,10,15,34];
end

for n=  ns
    nexttile
    if use_mucus
    filename = dir_m + sprintf("%i.csv",n);
    else
        filename = dir_w + sprintf("%i.csv",n);
    end
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
    a = 1e-6;
    
    Do = kb*Temp/(6*pi*mu*a);
    Dr = kb*Temp/(8*pi*mu*(a^3));
    
    %Dl = Do + (v^2)/(4*Dr); % rotationally diffusive motor 
    Dl = Do + (v^2)/(2*Dr*(1+ (omega/Dr)^2)); % circle-swimming motor (assume constant speed, spherical shape)
    
    disp("Analytical value: " + Dl)
    
    
    time = (1:20*length(x))*dt;
    
    [D, alpha, msd, tau]=compute_msd_log(x,y,dt,maxT,1);

    %Reference lines
    scatter(tau,msd); hold on;
    tau_ref = tau(round(1));
    msd_ref = msd(round(1));
    %plot(tau(1:5),msd(1)*(tau(1:5)/tau(1)).^2);
    plot(tau, msd_ref*(tau/tau_ref).^1, 'k:',  'LineWidth', 1);
    plot(tau, msd_ref*(tau/tau_ref).^2, 'k--', 'LineWidth', 1);

    %plot(time,Do*time + (time.^2)*(v^2),"--","DisplayName","Short-time ballistic motion","LineWidth",2);
    %xline(maxT*dt)
    %xline(maxT*dt*.2,'--',Color='r',LineWidth=1.5)
    msd_an = 4*Do*tau ...
        + (2*v^2*Dr*tau)/(Dr^2 + omega^2) ...
        + (2*v^2*(omega^2 - Dr^2))/((Dr^2 + omega^2)^2) ...
        + (2*v^2*exp(-Dr*tau))/((Dr^2 + omega^2)^2) .* ...
          ((Dr^2 - omega^2)*cos(omega*tau) - 2*omega*Dr*sin(omega*tau)); plot(tau,msd_an,Color=[.3,.9,1],LineWidth=2)
    legend('off')
    
    xlim([0,maxT*dt]);
    ylim([min(msd) 2*max([max(msd_an),max(msd)])])
    set(gca,"XScale","log")
    set(gca,"YScale","log")

    ylabel('MSD(\tau)',FontSize=14)

end
xlabel('\tau',FontSize=18)
set(gcf,"Position",[700,100,1000,1000])


if use_mucus
    nexttile(1);title('Mean squared displacement, picoalgae in mucus',FontSize=17)
else
nexttile(1);title('Mean squared displacement, picoalgae in water',FontSize=17)
end