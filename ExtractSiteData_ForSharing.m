
Tstampstr=num2str(MyData.data(:,1));

year=str2num(Tstampstr(:,1:4)); 
month=str2num(Tstampstr(:,5:6)); 
day=str2num(Tstampstr(:,7:8)); 
hr=str2num(Tstampstr(:,9:10));   %start of averaging period is 0 minutes

[yrout,doy]=makeDOY(year,month,day); doy=doy';
cdoy=year+doy./366; 


Ta=MyData.data(:,6);


WS=MyData.data(:,12);



ustar=MyData.data(:,11);



H=MyData.data(:,5);


H(WS<-100)=-9999;



Le=MyData.data(:,4);

Pa=MyData.data(:,13);


Rlout=MyData.data(:,17);


Rsout=MyData.data(:,15);

Rsin=MyData.data(:,14);


%calculate surface T

[Ts]=calculateTs(cdoy,hr,Rsin,Rsout,Rlout);



Ts(abs(Ts)>50)=NaN;


%use barometric formula to convert to potential temperature
%P=101.325*exp(-0.00012*z);
Tapot=real(Ta.*(101.3./(101.325*exp(-0.00012*(ele+2.7)))).^(2/7));




Ta=Tapot; 
