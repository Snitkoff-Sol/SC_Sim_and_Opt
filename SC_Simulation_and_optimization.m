clc
clear all
close all
%---BADS simulation of surface confiend species----
rng default
Num=1;
%-------------Call experimental data-------------------------------
Exp=load(''); % three columns potential|current|time
[C,ia,ic]=unique(Exp(:,3)-Exp(1,3)); % In case there are non unique points
t_e=linspace(0,C(end),2^17); 
Jex=interp1(C,Exp(ia,2),t_e'); % To ensure constant sampling rate
f=0.119; % Applied sinewave frequency in Hz
%----------------------------------------------------------------------




%--------Setting up the Fourier transform---------------------------------------------
Jex_fft = fft(Jex); 
Y=zeros(length(t_e),5);
%Sampling frequency
fs=length(t_e)./t_e(end);
%Window band
band=0.02;
% Length of the input signal
L = length(t_e);
% Define the frequency axis
x_f = (0:L-1)*fs/L';
% Define the width of the Gaussian mask
std=10;
alpha=(L-1)/(2*std);
% Define the Gaussian mask
g=gausswin(L,alpha);
% Define the frequency band to extract
Nh=7; %Number of harmonics
H=3:Nh;
f_low=H.*f-band/2;  % Lower frequency boundary
f_high = H*f+band/2;  % Upper frequency boundary
%------------------------------------------------------------------


%--------Extract harmonics-----------------------------------------
for i=1:length(H)
mask = (x_f >= f_low(i)) & (x_f <= f_high(i)); % Rectangular 
Conv_mask=conv(mask,g,'same')'; % Convolute with gaussian
Norm_Mask=Conv_mask./max(Conv_mask);
% Apply the mask to the FFT of the input signal
Jex_fft_masked = Jex_fft.*Norm_Mask;
% Perform IFFT on the masked FFT to obtain the extracted frequency band
Jex_ifft = ifft((Jex_fft_masked));
%YHarC=bandpass(Jfc,wind,df,'ImpulseResponse','iir','Steepness',0.7);
Y(:,i)=2*abs(hilbert(real(Jex_ifft)));
end
%-----------------------------------------------------------------




%------Set up parameters for BADS-----------------------------------
% Please refer to https://github.com/acerbilab/bads for details on BADS

% P=[gamma1,k1,E1]
P0=[3.16E-8,10,0.8];
LB=[1.0E-9,0.1,0.5];
UB=[2.756E-7,100,0.85];
PLB=[1E-8,1,0.52];
PUB=[2.7E-7,50,0.82];
funfun=@(P)Onesite(P,Y,t_e,Num);
%-----------Run BADS--------------------------------------------------
[P,Fval]=bads(funfun,P0,LB,UB,PLB,PUB,[]);
%----------------------------------------------------------------------




%----------Simulation and objective function------
function Output=Onesite(P,Y,t_e,Num)

%----Values for the samples----------------
Rcell=0.715;      %Ohm
CDL=0.2863;       %F/cm2


%----Constant Parameters--------
F=96485; %Faraday's constant
R=8.31; %Gas constant 
T=353.15; %Temperature in K
b=R.*T./F; %Thermal voltage in V
%-----Experimental Paramteres----------------
Amp=130/1000; %Amplitude in V
v=0.476/1000; %Scan rate V/s
f=0.119; %Frequency in Hz
%f=0.358; %Frequency in Hz
Ru=Rcell(Num); %Uncompensated resistance Ohm
Cdl=CDL(Num);% Double layer F/cm2
Ei=0.55; %Initial potential in V
A=1; %Area of electrode in cm2


%----------Parameters to evaluate-------------
gamma1=P(1); % Total number of active sites in mol/cm2 
k1=P(2); % Rate constant
E1=P(3); % Redox potential
%----Non-Dimensionalizing--------
T=t_e.*f; %Time
Xsii=Ei./b;%Initial potential
K1=k1./f; %Rate constant
E1=E1./b; % Redox potential
Am=Amp/b; % Amplitude
Ni=v./(f.*b); % Scan rate
E=Xsii+Ni.*T+Am.*sin(2*pi*T); % Potential
RC=Cdl*Ru*f*A; 
TCA=linspace(0,2./f,2^12).*f; % Time for CA
J_Norm1=f*gamma1*F;
Cdl1=Cdl*f*b/(J_Norm1);
%-------Solving the set of ODEs------------------------
opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
Am1=0; % Zero amplitude for CA
[~,y1]=ode15s(@(t,y)One_site(t,y,K1,E1,Ni,Am1,Xsii,Cdl1,RC),TCA,[1,Xsii],opts);
% Solution for odes
[~,y]=ode15s(@(t,y)One_site(t,y,K1,E1,Ni,Am,Xsii,Cdl1,RC),T,[y1(end,1),y1(end,2)],opts);

% ---------Total current-----
J(:)=((E(:)-y(:,2)).*b./Ru).*1000;


%-------------Loss function evaluation---------------
 



%---Setup before harmonuc extraction-------------
%FFT the current
J_fft = fft(J'); 

%Preallocate memory

%Sampling frequency
fs=length(t_e)./t_e(end);
%Window band
band=0.02;
% Length of the input signal
L = length(t_e);
% Define the frequency axis
x_f = (0:L-1)*fs/L';
% Define the width of the Gaussian mask
std=10;
alpha=(L-1)/(2*std);
% Define the Gaussian mask
g=gausswin(L,alpha);
% Max harmonic to evaluate
Nh=7;
% Harmonics to evaluate 3 to 7
H=3:Nh;
f_low=H.*f-band/2;  % Lower frequency boundary
f_high = H*f+band/2;  % Upper frequency boundary
%--------Extract harmonics------------
% start the loss function at 0
Loss=0;
for i=1:length(H)

% Create the Gaussian mask
mask = (x_f >= f_low(i)) & (x_f <= f_high(i));
Conv_mask=conv(mask,g,'same')';
Norm_Mask=Conv_mask./max(Conv_mask);
% Apply the mask to the FFT of the input signal

J_fft_masked = J_fft.*Norm_Mask;
% Perform IFFT on the masked FFT to obtain the extracted frequency band
J_ifft = ifft((J_fft_masked));
%YHarC=bandpass(Jfc,wind,df,'ImpulseResponse','iir','Steepness',0.7);
I_sim=2*abs(hilbert(real(J_ifft)));
Loss=Loss+(sum((I_sim(:)-Y(:,i)).^2)./sum(Y(:,i).^2));
end

Output=sqrt(Loss);
end


%------Simulation for the SC species -------------
function dydx = One_site(t,y,K1,E1,Ni,Am,Xsii,Cdl1,RC)

E=Xsii+Ni.*t+Am.*sin(2*pi*t);
xsi1=(y(2)-E1)./(2);
r1=K1*((1-y(1)).*exp(-xsi1)-(y(1)).*exp(xsi1));
rE=(E-y(2))/(RC)+(r1./Cdl1);
dydx=[r1;rE];
end