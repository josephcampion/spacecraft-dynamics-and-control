% Euler Angle Addition and Subtraction

clc

%% Problem 1

% 3-2-1 rotation matrix (DCM) of Euler angles below:

theta1 =      10 /180*pi;
theta2 =    20 /180*pi;
theta3 =      30 /180*pi;

DCM321 = dirCosMat(theta1, theta2, theta3, 3, 2, 1)

% 3-1-3 rotation matrix (DCM) of Euler angles below:

% psi2 =      38.4812 /180*pi;
% theta2 =    -9.84655 /180*pi;
% phi2 =      17.4952 /180*pi;

psi2 =      40.6423 /180*pi;
theta2 =    35.5313 /180*pi;
phi2 =      -36.0524 /180*pi;

% 
DCM313 = dirCosMat(psi2, theta2, phi2, 3, 1, 3)


%% Problem 2

DCM_BN = dirCosMat(theta1, theta2, theta3, 3, 2, 1)

DCM_RN = dirCosMat(-5*pi/180, 5*pi/180, 5*pi/180, 3, 2, 1) 

DCM_BR = DCM_BN * transpose(DCM_RN)


%% %%%%%%%%%%%%%%%%%%%%%% Helper functions: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function DCM = dirCosMat(theta1, theta2, theta3, n1, n2, n3)

Ri = dcmHelper(theta1, n1);
Rj = dcmHelper(theta2, n2);
Rk = dcmHelper(theta3, n3);

DCM = Rk * Rj * Ri;

end

function Rx = myRotX(q)

Rx = [1 0 0;
    0 cos(q) sin(q);
    0 -sin(q) cos(q)];

end

function Ry = myRotY(q)

Ry = [cos(q) 0 -sin(q);
    0 1 0;
    sin(q) 0 cos(q)];

end

function Rz = myRotZ(q)

Rz = [cos(q) sin(q) 0;
    -sin(q) cos(q) 0;
    0 0 1];

end

function R = dcmHelper(q, n)

R = eye(3);

if n == 1
    R = myRotX(q);
elseif n == 2
    R = myRotY(q);
elseif n == 3
    R = myRotZ(q);
else
    fprintf("n = 1, 2, or 3");
end

end

