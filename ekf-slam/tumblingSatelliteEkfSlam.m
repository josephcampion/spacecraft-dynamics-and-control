%% tumblingSatelliteEkfSlam.m


%%


Jx = 10;
Jy = 50;
Jz = 50;
inertia = [Jx; Jy; Jz];

T = 30;
dt = 1e-3;

%%

% h = figure;


%% Simulate dynamics

t = transpose(0:dt:T);
n = size(t,1);

quat_t = zeros(n,4);
quat_t(1,1) = 1.0;

omega_t = zeros(n,3);

tau_t = [ones(n,1) cos(t) -2 * sin(t)];

for ii = 1:(n-1)
    
    omega_t0 = transpose(omega_t(ii,:));
    omega_t1 = updateDynamics(inertia, omega_t0, tau_t(ii,:), dt);
    
    omega_t(ii+1,:) = transpose(omega_t1);
    
    quat_t0 = transpose(quat_t(ii,:));
    quat_t1 = updateAttitude(quat_t0, omega_t0, dt);
    % re-normalize:
    quat_t1 = quat_t1 / norm(quat_t1);
    quat_t(ii+1,:) = transpose(quat_t1);
   
end


% Plot results
figure
subplot(2,1,1)
plot(t,omega_t)
grid on
subplot(2,1,2)
plot(t,quat_t)
grid on

%% Simulate stars

% N = 10;
% 
% phi_stars = zeros(N,1);
% theta_stars = zeros(N,1);
% 
% for ii = 1:N
%     
%     % azimuthal angles
%     phi_stars(ii) = 2 * pi * rand();
%     
%     % polar angles
%     theta_stars(ii) = pi * rand();
%     
% end
% 
% figure
% plotStars(phi_stars, theta_stars)
% grid on
% axis equal


%%



%% Helper functions

function plotStars(phi, theta)

n = length(phi);


if (length(theta) ~= n)
    
    fprintf("Azimuthal angle vector and polar angle vector need to have the same length.\n");
    
    return;
end

for ii = 1:n
    
    sx = sin(theta(ii)) * cos(phi(ii));
    sy = sin(theta(ii)) * sin(phi(ii));
    sz = cos(theta(ii));
    
    plot3(sx, sy, sz, 'b*')
    hold on
    
end


end

%% Satellite Dynamics

function omega_t1 = updateDynamics(inertia, omega_t0, torque_t0, dt)

Jx = inertia(1);
Jy = inertia(2);
Jz = inertia(3);

wx_t0 = omega_t0(1);
wy_t0 = omega_t0(2);
wz_t0 = omega_t0(3);

tau_x_t0 = torque_t0(1);
tau_y_t0 = torque_t0(2);
tau_z_t0 = torque_t0(3);

wx_t1 = wx_t0 + ((Jy - Jz) * wy_t0 * wz_t0 + tau_x_t0) / Jx  * dt;
wy_t1 = wy_t0 + ((Jz - Jx) * wz_t0 * wx_t0 + tau_y_t0) / Jy * dt;
wz_t1 = wz_t0 + ((Jx - Jy) * wx_t0 * wy_t0 + tau_z_t0) / Jz * dt;

omega_t1 = [wx_t1; wy_t1; wz_t1];

end

%% Update attitude (quaternion)

function quat_t1 = updateAttitude(quat_t0, omega_t0, dt)

A = [0 -transpose(omega_t0); omega_t0 skewMat(omega_t0)];

quat_t1 = (eye(4) + A * dt) * quat_t0;

end

%% Skew matrix

function skew_mat = skewMat(omega)

skew_mat = [0          -omega(1)    omega(2);
            omega(2)    0          -omega(3);
            omega(2)    omega(3)    0];

end

