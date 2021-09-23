%% conceptCheck9.m

psi0 = 40 * pi / 180;
theta0 = 30 * pi / 180;
phi0 = 80 * pi / 180;

x0 = [psi0; theta0; phi0];

%%

clc

T = 42;
dt = 0.1;
t = transpose(0:dt:T);

omega_t = [sin(0.1*t) 0.01*ones(length(t),1) cos(0.1*t)] * 20 * pi / 180;

xt = zeros(length(t), size(x0,1));

xt(1,:) = x0';

for ii = 2:length(t)
    
    fprintf("iter: %i\t", ii);
    
    theta_i_minus_1 = xt(ii-1,2);
    phi_i_minus_1 = xt(ii-1,3);
    
    fprintf("theta(t)=%f\t", theta_i_minus_1)
    fprintf("phi(t)=%f\t", phi_i_minus_1)
    
    eulDot = eulerDeriv(theta_i_minus_1, phi_i_minus_1,...
        transpose(omega_t(ii-1,:)));
    
    xt(ii,:) = xt(ii-1,:) + transpose(eulDot) * dt;
    
    fprintf("\n");
    
end

%%

figure
subplot(2,1,1)
plot(t,omega)
grid on
subplot(2,1,2)
plot(t,xt)
grid on

%% 


function eulDot = eulerDeriv(theta, phi, omega)


eulDot = 1/cos(theta) * ...
    [0 sin(phi) cos(phi);
    0 cos(phi)*cos(theta) -sin(phi)*cos(theta);
    cos(theta) sin(phi)*sin(theta) cos(phi)*sin(theta)] * ...
    omega;



end

% function xt = myIntegrator(derivFunction, x0, t)
% 
% xt = zeros(length(t), size(x0,1));
% 
% end

