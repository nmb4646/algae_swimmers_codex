function[D_out] = D1(a,v)
% Effective diffusivity of a swimmer with translational motion at speed v and rotional
% diffusivity
kb = 1.38e-23;
T = 293;
mu = .001;

Do = kb*T/(6*pi*mu*a);
Dr = kb*T/(8*pi*mu*(a^3));

D_out = Do + (v^2)/(4*Dr);
end

