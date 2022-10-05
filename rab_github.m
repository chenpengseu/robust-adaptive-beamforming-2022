%--------------------RAB-------------------%
clc;
clear all;
close all;

%--------------------Parameters-------------------%
vecH = @(MAT) MAT(:).';
M = 10;     % the number of sensors
L = 80;     % the number of snapshots
SINRoptz = [];
SINRvolz = [];
SINRausz  = [];
SINRpropz = [];                                                                                                                                                           
SINRlinearz = [];
SINRsubz = [];
SINRmepsz = [];

%--------------------Parameters-------------------%
for counterz=1:300
    clc;
    close all;
    counterz = counterz + 1; %counter
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
%     As = steerVec(theta(1));
%     Ai1 = steerVec(theta(2));
%     Ai2 = steerVec(theta(3));
%     thetaedi1 = (-34+8*rand(1));
%     thetaedi2 = (36+8*rand(1));
%     thetaeds1 = (6 + 8*rand(1));
%     thetaed = deg2rad([thetaeds1,thetaedi1,thetaedi2]); %% estimated directions  
%     Aed2 = 1*exp(1j*2*pi*[0:1:M-1].'*d*sin(thetaed(2)));
%     Aed3 = 1*exp(1j*2*pi*[0:1:M-1].'*d*sin(thetaed(3)));
%     Aed1 = 1*exp(1j*2*pi*[0:1:M-1].'*d*sin(thetaed(1)));
    %--------------------DOA error--------------------%

    %--------------------SV mismatch--------------------%
    thetaedi1 = -30;
    thetaedi2 = 40;
    thetaeds1 = 10;
    thetaed = deg2rad([thetaeds1,thetaedi1,thetaedi2]);
    ammis = unifrnd(0,sqrt(0.03))*ones(1,M);
    phmis = deg2rad(unifrnd(0,360,1,M));
    mise = (ammis.* steerVec2(phmis)).'; %% mismatch due to sv random error
     As = steerVec(theta(1)) + mise;
    Ai1 = steerVec(theta(2)) + mise;
    Ai2 = steerVec(theta(3)) + mise;
    Aed2 = steerVec(thetaed(2));
    Aed3 = steerVec(thetaed(3));
    Aed1 = steerVec(thetaed(1));
    %--------------------SV mismatch--------------------%

    %--------------------Gain and phase perturbations--------------------%
%     thetaedi1 = -30;
%     thetaedi2 = 40;
%     thetaeds1 = 10;
%     thetaed = deg2rad([thetaeds1,thetaedi1,thetaedi2]);
%     randzy = round(rand*9)+1;
%     As = zeros(M,1);
%     Ai1 = zeros(M,1);
%     Ai2 = zeros(M,1);
%     for i=1:M
%         fudu = normrnd(1,0.05^2);
%         xiangwei = (normrnd(0,(0.025*pi)^2));
%         As(i) = (fudu*exp(1j*(2*pi*(i-1)*d*sin(theta(1)+xiangwei))));
%         Ai1(i) = (fudu*exp(1j*(2*pi*(i-1)*d*sin(theta(2)+xiangwei))));
%         Ai2(i) = (fudu*exp(1j*(2*pi*(i-1)*d*sin(theta(3)+xiangwei))));
%     end
% %     As = As.';
% %     Ai1 = Ai1.';
% %     Ai2 = Ai2.';
%     Aed2 = steerVec(thetaed(2));
%     Aed3 = steerVec(thetaed(3));
%     Aed1 = steerVec(thetaed(1));
    %--------------------Gain and phase perturbations--------------------%

    %--------------------Optimal--------------------%
    signal = sqrt(p1/2) * (randn(L, 1) + 1j* randn(L, 1)); 
    Xs = As * signal.'; 
    inter = sqrt(p2/2) * (randn(L, 2) + 1j* randn(L, 2)); 
    Xi = Ai1 * inter(:,1).' + Ai2 * inter(:,2).';  
    Xn = (randn(L, M) + 1j* randn(L, M)).';
    Xx=Xs+Xi+Xn; 
    Rxx = Xx*Xx'/L;
    Rii = (Xi+Xn)*(Xi+Xn)'/L;
    wopt = inv(Rii)*Aed1/(Aed1'*inv(Rii)*Aed1);
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

    %--------------------------AUS--------------------------%
    sita2 = 0.3;
    detatheta = fullfact(2*ones(1,10))-1;
    for i=1:1024
        for j = 1:10
            if(detatheta(i,j)==0)
                detatheta(i,j) = -1;
            end
        end
    end
    R = zeros(M,M,I*2);
    counter = 1;
    for i = (thetaedi1)-7.9:0.8:(thetaedi1)+7.9   
        for j = 1:1024
            a = steerVec(deg2rad(i))+sita2/sqrt(M)*detatheta(j);
            R(:,:,counter) = a*a'/(a'*inv(Rxx)*a)+ R(:,:,counter);
        end
        counter = counter + 1;
    end
    for i = (thetaedi2)-7.9:0.8:(thetaedi2)+7.9 
        for j = 1:1024
            a = steerVec(deg2rad(i))+sita2/sqrt(M)*detatheta(j);
            R(:,:,counter) = a*a'/(a'*inv(Rxx)*a)+ R(:,:,counter);
        end
        counter = counter + 1;
    end
    R_aus = zeros(M,M);
    for i=1:I
        R_aus = R(:,:,i)+R_aus;
    end
    R_aus = 1/2*R_aus + pn(1)*eye(M);
    
    sita4 = 0.3;
    Ri = zeros(M,M,I);
    counter2 = 1;
    for i = thetaeds1-7.9:0.8:thetaeds1+7.9 
        for j = 1:1024
            a = steerVec(deg2rad(i))+sita4/sqrt(M)*detatheta(j);
            Ri(:,:,counter2) = a*a'/(a'*inv(Rxx)*a)+ Ri(:,:,counter2);
        end
        counter2 = counter2 + 1;
    end
    Rs_aus = zeros(M,M);
    for i=1:I
        Rs_aus = Ri(:,:,i)+Rs_aus;
    end
    Rs_aus = 1/2*Rs_aus + pn(1)*eye(M);
    [V,D] = eig(Rs_aus);
    [D_sort,index1] = sort(diag(D),'descend');
    V_sort = V(:,index1);
    sum1 = 0;
    sum2 = 0;
    sum12 = 0;
    for i=1:M
        if(sum12<=0.9)
            sum1 = 0;
            sum2 = 0;
            for j=1:i
                sum1 = (abs(D_sort(j)))^2+sum1;
            end
            for j=i:M
                sum2 = (abs(D_sort(j)))^2+sum2;
            end
            sum12 = sum1/sum2;
            counter3 = i;
        else
            break;
        end
    end
    Vs = V_sort(1:counter3,:);
    [V2,D2] = eig(Rxx);
    [D2_sort,index2] = sort(diag(D2),'descend');
    V2_sort = V2(:,index2);
    Es = V2_sort(:,1:3).';
    PP = Es'*Es*Vs'*Vs;
    [V3,D3] = eig(PP);
    maxD3 = max(max(D3));
    [x1,y1] = find(D3==maxD3);
    As_aus = V3(:,y1);
    w_aus = inv(R_aus)*As_aus/(As_aus'*inv(R_aus)*As_aus);
    phi = deg2rad([-90:90]);
    for i =1:181
        a1 = steerVec(phi(i));
        y_asu(i) = w_aus'*a1;
    end
    y_asu = 20*log10(abs(y_asu)/max(abs(y_asu)));
    %--------------------------AUS--------------------------%
    %--------------------------VOLUME--------------------------%
    sita2 = 0.3;
    detatheta = fullfact(2*ones(1,10))-1;
    for i=1:1024
        for j = 1:10
            if(detatheta(i,j)==0)
                detatheta(i,j) = -1;
            end
        end
    end
    R = zeros(M,M,2*I);
    counter = 1;
    for i = (thetaedi1)-7.9:0.8:(thetaedi1)+7.9   
        for j = 1:1024
            a = steerVec(deg2rad(i))+sita2/sqrt(M)*detatheta(j);
            R(:,:,counter) = a*a'/(a'*inv(Rxx)*a)+ R(:,:,counter);
        end
        counter = counter + 1;
    end
    for i = (thetaedi2)-7.9:0.8:(thetaedi2)+7.9 
        for j = 1:1024
            a = steerVec(deg2rad(i))+sita2/sqrt(M)*detatheta(j);
            R(:,:,counter) = a*a'/(a'*inv(Rxx)*a)+ R(:,:,counter);
        end
        counter = counter + 1;
    end
    R_vol = zeros(M,M);
    for i=1:I
        R_vol = R(:,:,i)+R_vol;
    end
    R_vol = 1/2*R_vol + pn(1)*eye(M);
    R = inv(Rxx); 
    As_vol = Aed1;
    cvx_begin
        variable e(M,1) complex;
        zj = quad_form(e+As_vol,R);
        zj2 = quad_form(e+As_vol,R_vol);
        zj3 = quad_form(As_vol,R_vol);
        minimize(zj);
        subject to 
            (As_vol'*e)==0;
            (zj2) <= (zj3);
    cvx_end
    As_vol = As_vol + e;
    w_vol = inv(R_vol)*As_vol/(As_vol'*inv(R_vol)*As_vol);
    phi = deg2rad([-90:90]);
    for i =1:181
        a1 = steerVec(phi(i));
        y_vol(i) = w_vol'*a1;
    end
    y_vol = 20*log10(abs(y_vol)/max(abs(y_vol)));
    %--------------------------VOLUME--------------------------%
    %--------------------------LINEAR--------------------------%
    Rinlinear = zeros(M,M);
    for i=thetaedi1 - 7.9 : 0.8 : thetaedi1 + 7.9
        Rinlinear = (steerVec(deg2rad(i))*steerVec(deg2rad(i))')/(steerVec(deg2rad(i))'*inv(Rxx)*steerVec(deg2rad(i))) + Rinlinear;
    end
    for i=thetaedi2 - 7.9 : 0.8 : thetaedi2 + 7.9
        Rinlinear = (steerVec(deg2rad(i))*steerVec(deg2rad(i))')/(steerVec(deg2rad(i))'*inv(Rxx)*steerVec(deg2rad(i))) + Rinlinear;
    end
     R = inv(Rxx); 
    As_linear = Aed1;
    cvx_begin
        variable e(M,1) complex;
        zj = quad_form(e+As_linear,R);
        zj2 = quad_form(e+As_linear,Rinlinear);
        zj3 = quad_form(As_linear,Rinlinear);
        minimize(zj);
        subject to 
            (As_linear'*e)==0;
            (zj2) <= (zj3);
    cvx_end
    As_linear = As_linear + e;
    wlinear = inv(Rinlinear)*As_linear/(As_linear'*inv(Rinlinear)*As_linear);
    phi = deg2rad([-90:90]);
    for i =1:181
        a1 = steerVec(phi(i));
        ylinear(i) = wlinear'*a1;
    end
    ylinear = 20*log10(abs(ylinear)/max(abs(ylinear)));
    %--------------------------LINEAR--------------------------%
    %--------------------------SUB--------------------------%
    Rinsub1 = zeros(M,M);
    Rinsub2 = zeros(M,M);
    for i=thetaedi1 - 7.9 : 0.8 : thetaedi1 + 7.9
        Rinsub1 = (steerVec(deg2rad(i))*steerVec(deg2rad(i))')/(steerVec(deg2rad(i))'*inv(Rxx)*steerVec(deg2rad(i))) + Rinsub1;
    end
    [V_sub,D_sub] = eig(Rinsub1);
    [D_subsort,index_sub1] = sort(diag(D_sub),'descend');
    V_sort = V_sub(:,index_sub1);
    v_sum = 0;
    for i =1 :M
        v_sum = v_sum + D_subsort(i);
        if((v_sum/sum(D_subsort))>0.9)
            counter = i;
            break;
        end
    end
    Vsub1 = V_sort(:,1:counter);
    [V2_sub,D2_sub] = eig(Rxx);
    [D2sub_sort,index_sub2] = sort(diag(D2_sub),'descend');
    V2sub_sort = V2_sub(:,index_sub2);
    Es_sub = V2sub_sort(:,1:3);
    PP_sub1 = Vsub1*Vsub1'*Es_sub*Es_sub';
    [V3_sub,D3_sub] = eig(PP_sub1);
    maxD3_sub = max(max(D3_sub));
    [x1_sub,y1_sub] = find(D3_sub==maxD3_sub);
    inter_1 = V3_sub(:,y1_sub);

    for i=thetaedi2 - 7.9 : 0.8 : thetaedi2 + 7.9
        Rinsub2 = (steerVec(deg2rad(i))*steerVec(deg2rad(i))')/(steerVec(deg2rad(i))'*inv(Rxx)*steerVec(deg2rad(i))) + Rinsub2;
    end
    [V_sub2,D_sub2] = eig(Rinsub2);
    [D_subsort2,index_sub2] = sort(diag(D_sub2),'descend');
    V_subsort2 = V_sub2(:,index_sub2);
    v_sum2 = 0;
    for i =1 :M
        v_sum2 = v_sum2 + D_subsort2(i);
        if((v_sum2/sum(D_subsort2))>0.9)
            counter2 = i;
            break;
        end
    end
    Vssub2 = V_subsort2(:,1:counter2);
    PP_sub2 = Vssub2*Vssub2'*Es_sub*Es_sub';
    [V4_sub,D4_sub] = eig(PP_sub2);
    maxD4_sub = max(max(D4_sub));
    [x2_sub,y2_sub] = find(D4_sub==maxD4_sub);
    inter_2 = V4_sub(:,y2_sub);
    powersub_1 = 1/(inter_1'*inv(Rxx)*inter_1);
    powersub_2 = 1/(inter_2'*inv(Rxx)*inter_2);
    R_sub = powersub_1 * inter_1*inter_1' + powersub_2 * inter_2 * inter_2'+ pn(1)*eye(M);
    a_sub = sqrt(M)*un(:,M);
    w_sub = inv(R_sub)*a_sub/(a_sub'*inv(R_sub)*a_sub);
    phi = deg2rad([-90:90]);
    for i =1:181
        a1 = steerVec(phi(i));
        y_sub(i) = w_sub'*a1;
    end
    y_sub = 20*log10(abs(y_sub)/max(abs(y_sub)));
    %--------------------------SUB--------------------------%
  
    %--------------------------MEPS--------------------------%
    u1 = zeros(M,1);
    u1(1) = 1;
    p_meps = zeros(M,M);
    R_meps_d = zeros(M,M);
    for i = thetaedi1 - 7.9 : 0.8 : thetaedi1 + 7.9
        p_meps = p_meps + 0.8*steerVec(deg2rad(i))*steerVec(deg2rad(i))'/((1/(u1.'*inv(Rxx)*u1))*(steerVec(deg2rad(i))'*inv(Rxx)*u1*u1'*inv(Rxx)'*steerVec(deg2rad(i))));
    end
    for i = thetaedi2 - 7.9 : 0.8 : thetaedi2 + 7.9
        p_meps = p_meps + 0.8*steerVec(deg2rad(i))*steerVec(deg2rad(i))'/((1/(u1.'*inv(Rxx)*u1))*(steerVec(deg2rad(i))'*inv(Rxx)*u1*u1'*inv(Rxx)'*steerVec(deg2rad(i))));
    end

    for i = thetaeds1 - 7.9 : 0.8 : thetaeds1 + 7.9
        R_meps_d = R_meps_d + 2*steerVec(deg2rad(i))*steerVec(deg2rad(i))'/((1/(u1.'*inv(Rxx)*u1))*(steerVec(deg2rad(i))'*inv(Rxx)*u1*u1'*inv(Rxx)'*steerVec(deg2rad(i))));
    end
    a_meps = R_meps_d*Aed1;
    p_meps = p_meps + pn(1)*eye(M);
    w_meps = inv(p_meps)*a_meps/(a_meps'*inv(p_meps)*a_meps);
    phi = deg2rad([-90:90]);
    for i =1:181
        a1 = steerVec(phi(i));
        y_meps(i) = w_meps'*a1;
    end
    y_meps = 20*log10(abs(y_meps)/max(abs(y_meps)));
    %--------------------------MEPS--------------------------%

    %-------------------------------SINR Calculate-----------------------------------%
    Xso = zeros(L,1);
    Xio = zeros(L,1);
    Xno = zeros(L,1);
    for i=1:L
        Xso(i) = w_aus'*Xs(:,i);
        Xio(i) = w_aus'*Xi(:,i);
        Xno(i) = w_aus'*Xn(:,i);
    end
    SINRaussvp = 10*log10(sum(abs(Xso).^2)/sum(abs(Xio).^2+abs(Xno).^2));
    for i=1:L
        Xso(i) = wopt'*Xs(:,i);
        Xio(i) = wopt'*Xi(:,i);
        Xno(i) = wopt'*Xn(:,i);
    end
    SINRopt = 10*log10(sum(abs(Xso).^2)/sum(abs(Xio).^2+abs(Xno).^2));
    for i=1:L
        Xso(i) = w_vol'*Xs(:,i);
        Xio(i) = w_vol'*Xi(:,i);
        Xno(i) = w_vol'*Xn(:,i);
    end
    SINRvol = 10*log10(sum(abs(Xso).^2)/sum(abs(Xio).^2+abs(Xno).^2));

    for i=1:L
        Xso(i) = w_prop'*Xs(:,i);
        Xio(i) = w_prop'*Xi(:,i);
        Xno(i) = w_prop'*Xn(:,i);
    end
    SINRprop = 10*log10(sum(abs(Xso).^2)/sum(abs(Xio).^2+abs(Xno).^2));
    for i=1:L
        Xso(i) = wlinear'*Xs(:,i);
        Xio(i) = wlinear'*Xi(:,i);
        Xno(i) = wlinear'*Xn(:,i);
    end
    SINRlinearqc = 10*log10(sum(abs(Xso).^2)/sum(abs(Xio).^2+abs(Xno).^2));
    for i=1:L
        Xso(i) = w_sub'*Xs(:,i);
        Xio(i) = w_sub'*Xi(:,i);
        Xno(i) = w_sub'*Xn(:,i);
    end
    SINRsub = 10*log10(sum(abs(Xso).^2)/sum(abs(Xio).^2+abs(Xno).^2));
     for i=1:L
        Xso(i) = w_meps'*Xs(:,i);
        Xio(i) = w_meps'*Xi(:,i);
        Xno(i) = w_meps'*Xn(:,i);
    end
    SINRmeps = 10*log10(sum(abs(Xso).^2)/sum(abs(Xio).^2+abs(Xno).^2));
    SINRoptz = [SINRopt,SINRoptz];
    SINRvolz = [SINRvol,SINRvolz];
    SINRausz = [SINRaussvp,SINRausz];
    SINRpropz = [SINRprop,SINRpropz];
    SINRlinearz = [SINRlinearqc,SINRlinearz];
    SINRsubz = [SINRsub,SINRsubz];
    SINRmepsz = [SINRmeps,SINRmepsz];
end
s1 = zeros(8,1);
s1(1) = sum(SINRoptz)/300;
s1(2) = sum(SINRpropz)/300;
s1(3) = sum(SINRausz)/300;
s1(4) = sum(SINRvolz)/300;
s1(5) = sum(SINRlinearz)/300;
s1(6) = sum(SINRsubz)/300;
s1(7) = sum(SINRmepsz)/300;

