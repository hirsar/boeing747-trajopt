% ------------------------------------------------------
%  Derivation of TPBVP for B747 Trajectory Optimization
%  -----------------------------------------------------
%  by Hiren Saravana, March 2020
%  -----------------------------------------------------

syms X Z U W th Q rho CD S CL m FT Iyy MAC p0 p1 p2 p3 p4 p5 p6 g

a = atan(W/U);
D = 0.5*rho*(U^2+W^2)*CD*S;
L = 0.5*rho*(U^2+W^2)*CL*S;
dX = U*cos(th)+W*sin(th);
dZ = -U*sin(th)+W*cos(th);
dU = -W*Q-g*sin(th)-cos(a)*D/m+sin(a)*L/m+FT/m;
dW = U*Q+g*cos(th)-sin(a)*D/m-cos(a)*L/m;
dth = Q;
dQ = MAC/Iyy;

H = p0*(FT^2+MAC^2)+p1*dX+p2*dZ+p3*dU+p4*dW+p5*dth+p6*dQ;
dHdFT = diff(H, FT);
dHdMAC = diff(H, MAC);
dp1 = -diff(H,X);
dp2 = -diff(H,Z);
dp3 = -diff(H,U);
dp4 = -diff(H,W);
dp5 = -diff(H,th);
dp6 = -diff(H,Q);
