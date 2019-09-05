cloud=[];aod15=[];
for ano=2016:2019
    j=ano-2015;
    s1_=strrep('/Users/aredondas/CODE/rbcce.aemet.es/iberonesia/RBCC_E/2019/Triad/BSRN/cloudScreening.txt','2019',num2str(ano))
    if exist(s1_)
        s=load(s1_);
        cloud=[cloud;s];
    else
        disp('cloud');
        disp('ano');
    end
    
    s1_=strrep('/Users/aredondas/CODE/rbcce.aemet.es/iberonesia/RBCC_E/2019/Triad/BSRN/190101_191231_Izana.lev15','19',num2str(ano-2000))
    
    if exist(s1_)
        s1=fileread(s1_);
        if ano~=2016  % quito cabeceras
          for i=1:5 % 5 cabeceras
            [a,s1]=strtok(s1,char(10));disp(a)
          end
        end
        aod15=[aod15,s1];
    else
        disp('aod error');
        disp('ano');
    end
    
end

filewrite('aod_2016_2019.lev15',aod15);
f1 = fopen('cloudScreening.txt','w');
fprintf(f1,'%%date dayj clear_AM clear_PM\n');
fprintf(f1,'%f %03d %d %d\n',cloud')
fclose(f1)

 [aod,aod_m,aeronet]=read_aeronet('aod_2016_2019.lev15');
 aod_m(:,1)=fix(aod_m(:,1))
 aod_cloud=scan_join(aod_m,cloud);
 
 figure
 %errorbar(aod_cloud(:,1),aod_cloud(:,2),aod_cloud(:,3),'-')
 plot(aod_cloud(:,1),aod_cloud(:,2),'-')
  hold on
 gscatter(aod_cloud(:,1),aod_cloud(:,2),{aod_cloud(:,6),aod_cloud(:,7)},'','o')
 datetick
 grid
 title('AOD 340, colors: 0=clear 1=cloudy (AM,PM) ')
 