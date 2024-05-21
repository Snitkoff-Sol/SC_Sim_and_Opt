clear y y1
close all


%Code for ploting and comparing the output from the data optimization
%scheme

Rcell=[]; % Resistance Ohm 
CDL=[]; % Double layer capacitance F/cm2
P=[]; % Best fit parameters


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
Ru=Rcell; %Uncompensated resistance Ohm
Cdl=CDL;% Double layer F/cm2
Ei=0.55; %Initial potential in V
A=1; %Area of electrode in cm2

%---------------Call Experimental Data---------------
Exp=load(''); % three columns potential|current|time
[C,ia,ic]=unique(Exp(:,3)-Exp(1,3)); % In case there are non unique points
t_e=linspace(0,C(end),2^17); 
Jex=interp1(C,Exp(ia,2),t_e'); % To ensure constant sampling rate


%----------Kinetic Parameters-------------
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


%-------Solving the odes--------------------
opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
[t1,y1]=ode15s(@(t,y)One_site(t,y,K1,E1,Ni,0,Xsii,Cdl1,RC),TCA,[0,Xsii],opts);


opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
[t,y]=ode15s(@(t,y)One_site(t,y,K1,E1,Ni,Am,Xsii,Cdl1,RC),T,[y1(end,1),y1(end,2)],opts);

% Cell current
J(:)=((E(:)-y(:,2)).*b./Ru).*1000;



tdc=Ei+t_e'.*v;


Nh=7; %Number of harmonics 1st to 8th



%---Setup before loop-------------
%FFT the current
J_fft = fft(J'); 
Jex_fft = fft(Jex); 
%Preallocate memory
Harmonic=zeros(length(t),Nh);
Harmonic_ex=zeros(length(t_e),Nh);
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

%--------Extract harmonics------------
for i=1:Nh
H=i+1;
% % Define the frequency band to extract
f_low = H*f-band/2;  % Lower frequency boundary
f_high = H*f+band/2;  % Upper frequency boundary
% Create the Gaussian mask
mask = (x_f >= f_low) & (x_f <= f_high);
Conv_mask=conv(mask,g,'same')';
Norm_Mask=Conv_mask./max(Conv_mask);
% Apply the mask to the FFT of the input signal

J_fft_masked = J_fft.*Norm_Mask;
Jex_fft_masked = Jex_fft.*Norm_Mask;
% Perform IFFT on the masked FFT to obtain the extracted frequency band
J_ifft = ifft((J_fft_masked));
Jex_ifft = ifft((Jex_fft_masked));
%YHarC=bandpass(Jfc,wind,df,'ImpulseResponse','iir','Steepness',0.7);
Harmonic(:,i)=2*abs(hilbert(real(J_ifft)));
Harmonic_ex(:,i)=2*abs(hilbert(real(Jex_ifft)));
end



for j=1:6
    figure(j)
    k=num2str(j+1);
    plot(tdc,Harmonic(:,j),tdc,Harmonic_ex(:,j)),title([k '^{th} Harmonic']), legend('Simulation','Experiment'), xlabel('Cell Voltage/V'), ylabel('i/mA')
end



%------Simulation for the SC species -------------
function dydx = One_site(t,y,K1,E1,Ni,Am,Xsii,Cdl1,RC)

E=Xsii+Ni.*t+Am.*sin(2*pi*t);
xsi1=(y(2)-E1)./(2);
r1=K1*((1-y(1)).*exp(-xsi1)-(y(1)).*exp(xsi1));
rE=(E-y(2))/(RC)+(r1./Cdl1);
dydx=[r1;rE];
end