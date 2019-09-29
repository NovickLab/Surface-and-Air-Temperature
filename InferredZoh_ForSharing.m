function [Zoh,Zom]=InferredZoh(H,ustar,pres,zm,zd,U,Ta,Taero,isday);

cp = 1006; %specific heat capacity of dry air, j/kg/K
rho_w = 1000; %density of water, kg/m^3
k=0.4; %von Karman constant, unitless
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



%Invert Equation 4 to find Zoh
lnterm=(Ta-Taero).*(0.4*rho_a*cp.*ustar)./(-H)-psiH;

kk=find(lnterm>0&lnterm<50&ustar>0.01&(stab)<1&stab>-2&H>20&isday==1); %apply some data filters...
    %avoid very high ln(zo/z), avoid low H, avoid highly stable or unstable conditions (-2<stab<1)
    %limit analysis to daytime
  
    %find the zoh
 Zoh=nanmean((zm-zd))./(exp(nanmean(lnterm(kk))));

    
%Invert Equation 3 to find Zom
lnterm_m=U./ustar.*0.4-psiM;
kk=find(lnterm_m>0&lnterm_m<50&ustar>0.01&(stab)<1&stab>-2&H>20&isday==1); %apply same data filters...
Zom=nanmean(zm-zd)./(exp(nanmean(lnterm_m(kk))));





