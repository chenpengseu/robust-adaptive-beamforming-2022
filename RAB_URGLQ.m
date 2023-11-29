%--------------------RAB-------------------%
clc;
clear all;
close all;

%--------------------Parameters-------------------%
vecH = @(MAT) MAT(:).';
M = 10;     % the number of sensors
L = 80;     % the number of snapshots
thetai1 = -30;      %actual direction of interference 1
thetai2 = 40;       %actual direction of interference 2
thetas1 = 10;       %actual direction of desired signal
theta = deg2rad([thetas1,thetai1,thetai2]);
SNR = 20;
INR = 20;
p1 = 10^(SNR/10);
p2 = 10^(INR/10); 
I = 20; %the number of discretizations
d = 0.5; %the spacing of sensors
steerVec = @(angTmp) exp(1j*2*pi*[0:1:M-1].'*d*sin(vecH(angTmp))); % steering vector
steerVec2 = @(angTmp) exp(1j*angTmp);

%--------------------DOA error--------------------%
As = steerVec(theta(1));
Ai1 = steerVec(theta(2));
Ai2 = steerVec(theta(3));
thetaedi1 = (thetai1 + 8*rand(1)-4);
thetaedi2 = (thetai2 + 8*rand(1) -4);
thetaeds1 = (thetas1 + 8*rand(1) -4);
thetaed = deg2rad([thetaeds1,thetaedi1,thetaedi2]); %% estimated directions  
Aed2 = 1*exp(1j*2*pi*[0:1:M-1].'*d*sin(thetaed(2)));
Aed3 = 1*exp(1j*2*pi*[0:1:M-1].'*d*sin(thetaed(3)));
Aed1 = 1*exp(1j*2*pi*[0:1:M-1].'*d*sin(thetaed(1)));
%--------------------DOA error--------------------%

%--------------------Optimal--------------------%
signal = sqrt(p1/2) * (randn(L, 1) + 1j* randn(L, 1)); 
Xs = As * signal.'; 
inter = sqrt(p2/2) * (randn(L, 2) + 1j* randn(L, 2)); 
Xi = Ai1 * inter(:,1).' + Ai2 * inter(:,2).';  
Xn = (randn(L, M) + 1j* randn(L, M)).';
Xx=Xs+Xi+Xn; 
Rxx = Xx*Xx'/L;
Rii = (Xi+Xn)*(Xi+Xn)'/L;
wopt = inv(Rii)*As/(As'*inv(Rii)*As);
%--------------------Optimal--------------------%

%--------------------Proposed algorithm--------------------%
[un,pn_c] = eig(Rxx);
[pn,indeln] = sort(diag(pn_c),'ascend');
un = un(:,indeln);
Rsi = 100*trace(Rxx)*(Aed1*Aed1') + pn(1)*eye(M);
[p,u] = eig(Rsi);
[u_sort,index] = sort(diag(u),'descend');
p_sort = p(:,index);
B = (eye(M)-p_sort(:,1)*p_sort(:,1)')/pn(1);
Rxx2 = B'*Rxx*B;
Rin1 = Rxx2 + pn(1)*(eye(M) - B'*B);
f0 = steerVec(deg2rad(thetaedi1-8*sqrt(15)/5))*steerVec(deg2rad(thetaedi1-8*sqrt(15)/5))'/(steerVec(deg2rad(thetaedi1-8*sqrt(15)/5))'*inv(Rin1)*steerVec(deg2rad(thetaedi1-8*sqrt(15)/5)));
f1 = steerVec(deg2rad(thetaedi1))*steerVec(deg2rad(thetaedi1))'/(steerVec(deg2rad(thetaedi1))'*inv(Rin1)*steerVec(deg2rad(thetaedi1)));
f2 = steerVec(deg2rad(thetaedi1+8*sqrt(15)/5))*steerVec(deg2rad(thetaedi1+8*sqrt(15)/5))'/(steerVec(deg2rad(thetaedi1+8*sqrt(15)/5))'*inv(Rin1)*steerVec(deg2rad(thetaedi1+8*sqrt(15)/5)));
Rr1 = 16/2*(5/9*f0+8/9*f1+5/9*f2);
f0 = steerVec(deg2rad(thetaedi2-8*sqrt(15)/5))*steerVec(deg2rad(thetaedi2-8*sqrt(15)/5))'/(steerVec(deg2rad(thetaedi2-8*sqrt(15)/5))'*inv(Rin1)*steerVec(deg2rad(thetaedi2-8*sqrt(15)/5)));
f1 = steerVec(deg2rad(thetaedi2))*steerVec(deg2rad(thetaedi2))'/(steerVec(deg2rad(thetaedi2))'*inv(Rin1)*steerVec(deg2rad(thetaedi2)));
f2 = steerVec(deg2rad(thetaedi2+8*sqrt(15)/5))*steerVec(deg2rad(thetaedi2+8*sqrt(15)/5))'/(steerVec(deg2rad(thetaedi2+8*sqrt(15)/5))'*inv(Rin1)*steerVec(deg2rad(thetaedi2+8*sqrt(15)/5)));
Rr2 = 16/2*(5/9*f0+8/9*f1+5/9*f2);
R_pro = (Rr1 + Rr2) + pn(1)*eye(M);
As_pro = Aed1;
cvx_begin
    variable e(M,1) complex;
    zj = quad_form(e+As_pro,inv(R_pro));
    zj2 = quad_form(e+As_pro,R_pro);
    zj3 = quad_form(As_pro,R_pro);
    minimize(zj);
    subject to 
        (As_pro'*e)==0;
        (zj2) <= (zj3);
cvx_end
As_pro = As_pro + e;
w_prop = inv(R_pro)*As_pro/(As_pro'*inv(R_pro)*As_pro);
phi = deg2rad([-90:90]);
for i =1:181
    a1 = steerVec(phi(i));
    y15(i) = w_prop'*a1;
end
y15 = 20*log10(abs(y15)/max(abs(y15)));
%--------------------------Proposed algorithm--------------------------%

%-------------------------------SINR Calculate-----------------------------------%
Xso = zeros(L,1);
Xio = zeros(L,1);
Xno = zeros(L,1);
for i=1:L
    Xso(i) = w_prop'*Xs(:,i);
    Xio(i) = w_prop'*Xi(:,i);
    Xno(i) = w_prop'*Xn(:,i);
end
SINRprop = 10*log10(sum(abs(Xso).^2)/sum(abs(Xio).^2+abs(Xno).^2));
for i=1:L
    Xso(i) = wopt'*Xs(:,i);
    Xio(i) = wopt'*Xi(:,i);
    Xno(i) = wopt'*Xn(:,i);
end
SINRopt = 10*log10(sum(abs(Xso).^2)/sum(abs(Xio).^2+abs(Xno).^2));
