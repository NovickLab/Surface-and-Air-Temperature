function [yearout,doyout]=makeDOY(year,month,day)


uyear=unique(year);
doyout=[];
for j=1:length(uyear)
   yj=year(year==uyear(j));
   mj=month(year==uyear(j));
   dj=day(year==uyear(j));
   

if floor(year(1)/4)==year(1)/4
    leap=1;
else
    leap=0;
end

if leap==1
mdays=[0 31 60 91 121 152 192 213 244 274 305 335 366];
else
mdays=[0 31 59 90 120 151 181 212 243 273 304 334 365];
end

size(mj)
for i=1:12
    cc=find(mj==i);
    doy(cc)=dj(cc)+mdays(i);
end

size(doy)
doyout(year==uyear(j))=doy;
yearout(year==uyear(j))=uyear(j);
clear doy
end