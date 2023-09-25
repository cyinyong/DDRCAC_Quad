% Description:
%   This function computes the P-PI controller gains of an integrator plant
%   based on the specified Wn and Damping
% ----------------------------- %
% Input: 
%   sys_tf: Integrator plant in transfer function
%   wn1:     Wn of pos loop
%   wn2:     Wn of vel loop
%   zeta2:   Damping of vel loop

function [K, r] = compute_fsb_gains_v2(sys_A, sys_B, r)

% Compute required poles
% if (zeta <= 1)
%     j1 = complex(-zeta*wn1, wn1*sqrt(1-zeta^2));
%     j2 = complex(-zeta*wn1, -wn1*sqrt(1-zeta^2));
% else
%     j1 = -zeta*wn1+wn1*sqrt(zeta^2-1);
%     j2 = -zeta*wn1-wn1*sqrt(zeta^2-1);
% end
% r = [j1 j2 -wn2];

% if (zeta <= 1)
%     j1 = complex(-zeta*wn2, wn2*sqrt(1-zeta^2));
%     j2 = complex(-zeta*wn2, -wn2*sqrt(1-zeta^2));
% else
%     j1 = -zeta*wn2+wn2*sqrt(zeta^2-1);
%     j2 = -zeta*wn2-wn2*sqrt(zeta^2-1);
% end
% r = [-wn1 j1 j2];


K = place(sys_A,sys_B,r);

