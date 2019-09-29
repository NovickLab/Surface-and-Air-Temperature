function [Tprofile,zz,stab]=Textrap_ForSharing(doy,hr,Ta,H,ustar,pres,zm,zd,zo,U,Taero,maxz,ustarthresh,isDS,isGS,isnight,isday,h,sitenum);


cp = 1006; %specific heat capacity of dry air, j/kg/K
rho_w = 1000; %density of water, kg/m^3
k=0.4; %von Karman constant, unitless
pres=98;
rho_a = 101.3*10^3/287.04./(Ta+273.15); %density of dry air kg/m^3
S = 2508./(Ta+237.3).^2.*exp(17.3*Ta./(Ta+237)); %slope of saturate water vapor function, kPa/KSTst
Lv = 2500-2.360*Ta; %latent heat of vaporization, kJ/kg
gamma = cp.*pres./(0.622.*Lv*1000); %psychrometeric constant, kPa/K

Ta_K=Ta+273.15;  %degress kelvin
Cp=1005; %specific heat capacity of dry air
OL=-rho_a.*(ustar+eps).^3./(.4*9.81*((H+eps)./(Cp*Ta_K))); %Obuhkov lenght
stab=(zm)./OL; %atmospheric stability parameter

psiH=6*log(1+stab); %diabatic correction for heat, stable conditions
psiM=psiH;     %diabatic correction for momentum, stable conditions

%when stab<0, conditions are unstable, and the psiH and psiM require
%modification
psiH(stab<0)=-2*log((1+(1-16*stab(stab<0)).^0.5)/2);  %diabatic correction for heat, unstable conditions
psiM(stab<0)=0.6*psiH(stab<0); %diabatic correction for momentum, unstable conditions


Ta(Ta<-20)=NaN; Ta(Ta>50)=NaN; %exclude unreasonable values
Taero(Taero<-20)=NaN; Taero(Taero>50)=NaN;

%set the points of extrapolation (representing elevation in meters)
zz=[1:1:20 25:5:55 60:10:170];

Tprofile=ones(length(Taero),length(zz))*NaN; %set up the profile vector
for i=1:length(zz)
Tprofile(:,i)=[Taero-H./k./rho_a./cp./ustar.*(log((zz(i))./zo)+psiH)];
end





