kb = 1.38e-23;
T = 293;
mu = .001;
a = 1e-6;
%v = 50e-6;

Do = kb*T/(6*pi*mu*a);
Dr = kb*T/(8*pi*mu*(a^3));

Dl = Do + (v^2)/(4*Dr); % Rotationally diffusive motor
Dl = Do + (v^2)/(2*Dr*(1+ (omega/Dr)^2)) % Rotating motor with angular and translational diffusion