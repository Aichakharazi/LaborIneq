function ww = equationSystem(x,alp, VrsV,u_mV, LYV, LMV)
ww = [ % Unpack parameters
 

    % Given values
x(1) - (1 - u_mV);
% System of equations
x(3) - (x(2))^(alp/(alp-1));
%x(4) - VrsV / u_mV;
%x(5) - VrsV /x(4);

x(3) - (x(4))^(-alp);
x(2) - (x(4))^(1-alp);
x(2) - (LYV+LMV)/x(1);
];
