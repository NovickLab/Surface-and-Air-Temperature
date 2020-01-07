clear all

%%%%%%SECTION 1 - IMPORT DATA, make day/night indicators, indicators for
%%%%%%dormant and growing season


%import data; 
MyData=importdata('Duke_OF2005to2008.csv');
zm = 2.7; %measurement height, m
h=0.5; % canopy height
ele=153; %elevation, m
zd=0.6*h; %zero plane displacement, m

%extract key variables from the dataset
ExtractSiteData_ForSharing




%make a vector 'isDS' that equals 1 during the dormant season, 0 otherwise
isDS=ones(size(doy)); isDS(doy>100&doy<300)=0;

%make a vector 'isGS' that equals 1 during the growing season, 0 otherwise
isGS=zeros(size(doy)); isGS(doy>150&doy<270)=1;

%make variables to identify daytime and nighttime conditions
isnight=zeros(size(doy)); isnight(hr<7|hr>20)=1;
 isday=zeros(size(isnight));  isday(hr>9&hr<17)=1;
% 



 
 %%%%%%SECTION 2 - FIND THE AERODYNAMIC TEMPERATURE 
%this approach relies on the assumption that, when H is
%small (|H|<30 W/m2), the Ta and Taero are the same 
%The output of the loop is the ensemble mean ratio of Ta/Tsurf for every
%hour of the day, separately for the dormant and growing season
%this ratio is X_T in Eq. 5 of Novick and Katul


Taero=Ts*NaN; %set up an empty Taero vectro)
for j=1:24 %cycle through 24 hours in the day
gg=find(abs(H)<30&abs(Ta)<60&abs(Ts)<60&hr==(j-1)&isDS==1&ustar>0); %find dormant season data collected when H is small, Ta and Ts make sense, 
%and |Ts| is greater than 2 (to avoid instabilities when Ts is near zero)
%the variable 'gg' is an indexing variable that stores the locations of the
%suitable data
XTDS(j)=nanmean(Ta(gg)+273.15)./nanmean(Ts(gg)+273.15); %save the ratio X_T for each hour of the day, during the growing season
                                                        %note conversion to
                                                        %degrees K
sXTDS(j)=nanmean(Ta(gg)+273.15)./nanmean(Ts(gg)+273.15)./sqrt(length(gg)); %find the standard error of the mean for X_T

hh=find(hr==(j-1)&isDS==1); %find all dormant season data for the jth hour
Taero(hh)=XTDS(j).*(Ts(hh)+273.15)-273.15; %find the Taero for all time periods as the product of X_T and Tsurf
                                            %calcuate in deg K, convert
                                            %back to deg C

%now, repeat for the growing season
gg=find(abs(H)<30&abs(Ta)<60&abs(Ts)<60&hr==(j-1)&isGS==1&abs(Ts)>2&ustar>0); 
XTGS(j)=nanmean(Ta(gg)+273.15)./nanmean(Ts(gg)+273.15);
sXTGS(j)=nanmean(Ta(gg)+273.15)./nanmean(Ts(gg)+273.15)./sqrt(length(gg)); 
hh=find(hr==(j-1)&isGS==1); 
Taero(hh)=XTGS(j).*(Ts(hh)+273.15)-273.15; 
%note that Taero is still NaN during the shoulder seasons
end


%This ensures that Taero is bounded by Ta and Tsurf
Taero(Ts>Ta&Taero<Ta)=Ta(Ts>Ta&Taero<Ta); %if Ts>Ta, and Taero<Ta, Taero=Ta
Taero(Ts<Ta&Taero<Ts)=Ts(Ts<Ta&Taero<Ts); %If Ts<Ta, and Taero<Ts, Taero=Ts
Taero(Ts<Ta&Taero>Ta)=Ta(Ts<Ta&Taero>Ta); %if Ts<Ta, and Taero>Ta, Taero=Ta
Taero(Ta<Ts&Taero>Ts)=Ts(Ta<Ts&Taero>Ts); %if Ta<Ts, and Taero>Ts, Taero=Ts
Taero(Taero<-50)=NaN; %convert -9999s to NaNs for better visualization    
    
   
%find the mean Ts, Ta, and Taero for each hour of the day, for visualization purposes
for j=1:24
    hh=find(hr==(j-1)&isDS==1&abs(Ta)<50&abs(Ts)<60&abs(Taero)<50);
hTsDS(j)=nanmean(Ts(hh));
hTaDS(j)=nanmean(Ta(hh));
hTaeroDS(j)=nanmean(Taero(hh));
hh=find(hr==(j-1)&isGS==1&abs(Ta)<50&abs(Ts)<60&abs(Taero)<50);
hTsGS(j)=nanmean(Ts(hh));
hTaGS(j)=nanmean(Ta(hh));
hTaeroGS(j)=nanmean(Taero(hh));
end

%let's have a look
figure(3)
clf
subplot(1,3,1)
plot(XTDS,'b-')
hold on
plot(XTGS,'g-')
xlabel('Time of Day')
ylabel('X_T')
legend('dormant season','growing season')
subplot(1,3,2)
plot(hTsGS,'b-')
hold on
plot(hTaGS,'r-')
plot(hTaeroGS,'k-')
xlabel('Time of day')
ylabel('T (deg C)')
legend('Tsurf','Ta','Taero')
title('Growing season')
subplot(1,3,3)
plot(hTsDS,'b-')
hold on
plot(hTaDS,'r-')
plot(hTaeroDS,'k-')
xlabel('Time of day')
ylabel('T (deg C)')
legend('Tsurf','Ta','Taero')
title('Dormant season')


 



%%%%%%SECTION 3 - FIND THE ROUGHNLESS LENGTH FOR HEAT (Zo,h) 

%now that we have a guess of Taero, we can estimate the roughness length for heat (zo,h) by inverting the
%diabatic profile equation (eq 4 in Novick & Katul)

%we are going to do that by ensemble averaging within bins delineated by u*
%in Novick and Katul, this is done inside a bootstrap that incorporates
%uncertaintly in Taero...we're going to leave out the bootstrapping in this
%code, for clarity

sortedustar=sort(ustar(ustar>0.01&ustar<1.5));  %sort ustar, after first excluding -999s, zeros Nans, and very large values)
nn=length(sortedustar); %find the length of the sorted ustar

%now, we are going to make the edges of our u* bins
ustarbins=[0.01]; %the minimum value of our lowest u* bin
for j=1:15
    ustarbins(j+1)=sortedustar(ceil(j/15*nn)); %this codes divides the data into quantiles (15 of them)
end


%this loop cycles through the 15 u* bins, calling the subroutine to infer
%zo within each loop; Do it separately for dormant and growing season
%for informational purposes, we will also infer the momentum roughness
%lenght (zom)
%assume the zero plane displacement is 0.6 times canopy height
zoh=ones(size(ustar))*NaN; clear zohGS_binned zohDS_binned zomGS_binned zomDS_binned
for j=1:(length(ustarbins)-1)
    
    %growing season first
       gg=find(ustar>ustarbins(j)&ustar<=ustarbins(j+1)&isGS==1); %find data in the jth bin, growing season first
    [zohGS_binned(j),zomGS_binned(j)]=InferredZoh_ForSharing(H(gg),ustar(gg),Pa(gg),zm,0.6*h,WS(gg),Ta(gg),Taero(gg),isday(gg)); %find zom and zoh in each bin
    ustarplotGS(j)=nanmean(ustar(gg)); %find the mean ustar in each bin, for plotting
   
    %Make a continuous Zo,h time series
    zoh(gg)=zohGS_binned(j); 
    
    %dormant season
     gg=find(ustar>ustarbins(j)&ustar<=ustarbins(j+1)&isDS==1); %find data in the jth bin, growing season first
    [zohDS_binned(j),zomDS_binned(j)]=InferredZoh_ForSharing(H(gg),ustar(gg),Pa(gg),zm,0.6*h,WS(gg),Ta(gg),Taero(gg),isday(gg)); %find zom and zoh in each bin
    ustarplotDS(j)=nanmean(ustar(gg)); %find the mean ustar in each bin, for plotting
    %add to the continuous Zo,h time series
    zoh(gg)=zohDS_binned(j); 
      
end


%let's take a look
%the zom will likely increase as a function of u*
%the zoh should be less than, or approximately equal to, the zom
%in short ecosytems, this may not be the case, as the zoh and zom are very
%small, and small measurement errors can affect the ratio of the two
figure(4)
clf
plot(ustarplotGS,zohGS_binned./h,'k-','LineWidth',2)
hold on
plot(ustarplotGS,zomGS_binned./h,'k--','LineWidth',2)
plot(ustarplotDS,zohDS_binned./h,'k-','color',[.7 .7 .7],'LineWidth',2)
plot(ustarplotDS,zomDS_binned./h,'k-','color',[.7 .7 .7],'LineWidth',2)
legend('zo_h','zo_m')
xlabel('u^* (m/s)')
ylabel('zo_h or zo_m')
title('Growing Season')  
legend('zo_h/h, GS','zo_m/h, GS','zo_h/h, DS','zo_m/h, DS')
xlabel('u^* (m/s)')
ylabel('zo_h or zo_m')
title('Dormant Season')






%%%%%%SECTION 4 - EXTRAPOLATE THE TEMPERATURE OVER THE FIRST 150m of the
%%%%%%surface layer...this approach conceptually collapses the vertical
%%%%%%structure of the canopy, permitting a more apples-to-apples comarison

%%%%This is done in a subroutine..the main engine is Eq. 4 of Novick &
%%%%Katul
clear  Textrap
 
  %outputs of the subroutine are extrapolated temperature (deg C), the associated
 %elevation (m), and the atmospheric stability parameter
[Textrap,zz,stab]=Textrap_ForSharing(doy,hr,Ta,H,ustar,Pa,zm,zd,zoh./h,WS,Taero,0,0.0,isDS,isGS,isnight,isday,1);
   
 %find time periods when all the temperature metrics are good, and
 %-5<stab<1
Goodtas=zeros(size(Ta));
Goodtas(abs(Ta)<60&abs(Taero)<60&abs(Ts)<60&abs(Textrap(:,35))<60)=1;


%Note that, in the Novick & Katul paper, additional data were screened to
%ensure that temeprature metrics from the pine and hardwood site were also
%available; to recreate the profiles for the grass field shown in Novick & Katul Figure 7,
%you can directly load the QAQC variable called 'Goodtasall', which was
%used in the manuscript, by uncommenting the following line of code

%Goodtas=load('Goodtasall.txt');

%Plot the profiles
figure(5)
clf
subplot(2,2,1)
gg=find(isDS==1&Goodtas==1&isnight==1);
plot(nanmean(Textrap(gg,:)),zz,'r-','LineWidth',3)
hold on
plot(nanmean(Textrap(gg,:))-nanstd(Textrap(gg,:))./sqrt(length(gg)),zz,'r-','LineWidth',1)
plot(nanmean(Textrap(gg,:))+nanstd(Textrap(gg,:))./sqrt(length(gg)),zz,'r-','LineWidth',1)
xlabel('T_{extrap}, deg C')
ylabel('z (m)')
title('Dormant season, night')


subplot(2,2,2)
gg=find(isDS==1&Goodtas==1&isday==1);
plot(nanmean(Textrap(gg,:)),zz,'r-','LineWidth',3)
hold on
plot(nanmean(Textrap(gg,:))-nanstd(Textrap(gg,:))./sqrt(length(gg)),zz,'r-','LineWidth',1)
plot(nanmean(Textrap(gg,:))+nanstd(Textrap(gg,:))./sqrt(length(gg)),zz,'r-','LineWidth',1)
xlabel('T_{extrap}, deg C')
ylabel('z (m)')
title('Dormant season, day')

subplot(2,2,3)
gg=find(isGS==1&Goodtas==1&isnight==1);
plot(nanmean(Textrap(gg,:)),zz,'r-','LineWidth',3)
hold on
plot(nanmean(Textrap(gg,:))-nanstd(Textrap(gg,:))./sqrt(length(gg)),zz,'r-','LineWidth',1)
plot(nanmean(Textrap(gg,:))+nanstd(Textrap(gg,:))./sqrt(length(gg)),zz,'r-','LineWidth',1)
xlabel('T_{extrap}, deg C')
ylabel('z (m)')
title('Growing season, night')

subplot(2,2,4)
gg=find(isGS==1&Goodtas==1&isday==1);
plot(nanmean(Textrap(gg,:)),zz,'r-','LineWidth',3)
hold on
plot(nanmean(Textrap(gg,:))-nanstd(Textrap(gg,:))./sqrt(length(gg)),zz,'r-','LineWidth',1)
plot(nanmean(Textrap(gg,:))+nanstd(Textrap(gg,:))./sqrt(length(gg)),zz,'r-','LineWidth',1)
xlabel('T_{extrap}, deg C')
ylabel('z (m)')
title('Growing season, Day')







