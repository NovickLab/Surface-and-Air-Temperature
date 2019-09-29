function [Ts]=calculateTs(cdoy,hr,Rsin,Rsout,Rlout)

ucdoy=unique(cdoy);

for i=1:length(ucdoy)
    gg=find(cdoy==ucdoy(i)&hr>10&hr<15);
   RSrati=(Rsout(gg)./Rsin(gg));
albedo(cdoy==ucdoy(i))=nanmean(RSrati(RSrati>0&RSrati<0.3));
end
albedo=albedo';
es=-0.16*albedo+ 0.99;


Ts=(Rlout./es./(5.67*10^-8)).^(0.25)-273.15;
Ts=real(Ts);
Ts(isnan(Ts)==1)=-9999;
%Ts(isreal(Ts)==0)=-9999;