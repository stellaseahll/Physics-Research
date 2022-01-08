function [h,cw] = getBipartition(rho)
%SN 06.02.18: Get bipartition with reduced density matrices of hot mode and
%cold*work modes for 3-qubit system

h = zeros(2);
h(1,1) = rho(1,1) + rho(2,2) + rho(3,3) + rho(4,4);
h(2,2) = rho(5,5) + rho(6,6) + rho(7,7) + rho(8,8);
h(1,2) = rho(1,5) + rho(2,6) + rho(3,7) + rho(4,8);
h(2,1) = rho(5,1) + rho(6,2) + rho(7,3) + rho(8,4);

cw = rho(1:4,1:4) + rho(5:8,5:8);
% cw = zeros(4);
% cw(1,2) = rho(1,2) + rho(5,6);
% cw(1,3) = rho(1,3) + rho(5,7);
% cw(1,4) = rho(1,4) + rho(5,8);
% cw(2,3) = rho(2,3) + rho(6,7);
% cw(2,4) = rho(2,4) + rho(6,8);
% cw(3,4) = rho(3,4) + rho(7,8);
% 
% cw(1,1) = rho(1,1) + rho(5,5);
% cw(2,2) = rho(2,2) + rho(6,6);
% cw(3,3) = rho(3,3) + rho(7,7);
% cw(4,4) = rho(4,4) + rho(8,8);
