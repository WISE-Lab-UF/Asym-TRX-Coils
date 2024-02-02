clc;
clear all;
close all;

%% Inputs
k    = 0.1;
Lrx  = 200e-9;   % can be modified based on RX Area
R1   = 0.5;
R2   = R1;
% M    = sqrt(k*k*L1*L2);
fopt = 20e6;
zp   = 50; 

%% finding Lreq
Lreq = sqrt((R1+zp)*(R2+zp))/(2*pi*fopt*sqrt(1-k^2));
Ltx  = (Lreq^2)/Lrx;

%% RX coil parameters
n_rx = 1:10;
w_rx = (1e-3)*(0.1:.01:2);
ri_rx= (1e-3)*(0.1:.01:2);
seg  = 8;
ur   = 0.999994;
u    = 4*pi*1e-7;
Coct = [1.07, 2.29, 0.00, 0.19];
err  = 1e-10;
m    = 1;
mm   = 1;

for i= 1:length(n_rx)
    for j= 1:length(w_rx)
        for k= 1: length(ri_rx)
            delr(m) = 2*w_rx(j);   
               D(m) = w_rx(j)+2*((ri_rx(k)+n_rx(i)*delr(m))*cos(pi/seg));
             Din(m) = 2*ri_rx(k)*cos(pi/seg)-w_rx(j);
               b(m) = (D(m)-Din(m))./(D(m)+Din(m));
            davg(m) = 0.5*(D(m)+Din(m));
               L(m) = Coct(1,1)*u*ur*n_rx(i)^2*0.5*davg(m)*(log(Coct(1,2)/b(m))+(Coct(1,3))*b(m)+(Coct(1,4))*(b(m))^2);
             tol(m) = abs((L(m))-Lrx);
             if tol(m) < err && Din(m) > 3*w_rx(j) && D(m) < 20e-3
                LL(mm,1:6)=[n_rx(i) (1e3)*w_rx(j) (1e3)*ri_rx(k) (1e9)*L(m) (1e3)*(D(m))^1 (1e6)*(D(m))^2];
                mm=mm+1;
             end
            m=m+1;
        end
    end
end 
Final = sortrows(LL,6);             

%% RX coil parameters
% Similar as above replacing Lrx by Ltx
