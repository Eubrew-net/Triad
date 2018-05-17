% Daily model for midday
% O3= Ai*Ii + B*(t-to) +C(t-to)^2
% to=solar noon
%

sum_day=[];
for dia=1:310
r_date=datenum(2015,0,0)+dia;
noon_m=solar_noon(r_date, -16.499);
r_date=r_date+noon_m/24/60;
diag=[];
for j=1:3

  if dia>140 & dia<157 & j==3
      
  cdata=1;
  else 
    jd=find(diaj(summary{j}(:,1))==dia);
    idx=zeros(size(jd,1),3);
    idx(:,j)=1;
    diag_=[summary{j}(jd,6),idx,summary{j}(jd,1)-r_date,(summary{j}(jd,1)-r_date).^2,(summary{j}(jd,1)-r_date).^3];
    diag=[diag;diag_];
  end
end
[B,BINT,R,RINT,STATS] = regress(diag(:,1),diag(:,2:end));
sum_day=[sum_day;[dia,B',STATS]];

end
sum_day(sum_day==0)=nan;
sum_day(abs(sum_day)>1000)=nan;

r1=mean(sum_day(:,2:4),2);
plot(sum_day(:,1),100*matdiv(matadd(sum_day(:,2:4),-r1),r1))
