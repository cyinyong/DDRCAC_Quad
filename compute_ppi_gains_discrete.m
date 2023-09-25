% Description:
%   This function computes the P-PI controller gains of an integrator plant
%   based on the specified Wn and Damping
% ----------------------------- %
% Input: 
%   sys_tf: Integrator plant in transfer function
%   wn1:     Wn of pos loop
%   wn2:     Wn of vel loop
%   zeta2:   Damping of vel loop

function [Kp1, Kp2, Ki2, r_d] = compute_ppi_gains_discrete(sys_tf, wn1, wn2, zeta, Ts)

b = sys_tf.num{:}(end);

% Compute required poles
if (zeta <= 1)
    j1 = complex(-zeta*wn2, wn2*sqrt(1-zeta^2));
    j2 = complex(-zeta*wn2, -wn2*sqrt(1-zeta^2));
else
    j1 = -zeta*wn2+wn2*sqrt(zeta^2-1);
    j2 = -zeta*wn2-wn2*sqrt(zeta^2-1);
end

r = [-wn1 j1 j2];
r_d = exp(r.*Ts);

p = poly(r_d);

syms kp ki kd z s;

cc = (b*kp*(ki + kd*s))/(s^3 + b*ki*kp + b*ki*s + b*kd*s^2 + b*kd*kp*s);
cc_d = subs(cc,s,1/Ts*log(z));
[~,D] = numden(cc_d);
d1 = subs(D,z,r_d(1));
d2 = subs(D,z,r_d(2));
d3 = subs(D,z,r_d(3));

eqns = [d1 == 0, d2 == 0, d3 == 0];

% eqns = [(-3 - kd*b*Ts/2 + (ki*b+kp*kd)*(Ts/2)^2 + 3*ki*kp*(Ts/2)^3)/(1 + kd*b*Ts/2 + (ki*b+kp*kd)*(Ts/2)^2 + ki*kp*(Ts/2)^3), ...
%         (3 - kd*b*Ts/2 - (ki*b+kp*kd)*(Ts/2)^2 + 3*ki*kp*(Ts/2)^3)/(1 + kd*b*Ts/2 + (ki*b+kp*kd)*(Ts/2)^2 + ki*kp*(Ts/2)^3), ...
%         (-1 + kd*b*Ts/2 + (ki*b+kp*kd)*(Ts/2)^2 + 3*ki*kp*(Ts/2)^3)/(1 + kd*b*Ts/2 + (ki*b+kp*kd)*(Ts/2)^2 + ki*kp*(Ts/2)^3)]; 

S = solve(eqns,[kp ki kd]);
Kp = (double(S.kp));
Kd = (double(S.kd));
Ki = (double(S.ki));

if (~isreal(Kp) || ~isreal(Ki) || ~isreal(Kd))
    disp("ppi unable to meet poles")
end

Kp1 = abs(Kp(1));
Kp2 = abs(Kd(1));
Ki2 = abs(Ki(1));

