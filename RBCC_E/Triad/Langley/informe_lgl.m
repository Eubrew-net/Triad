% Buscar y sustituir este directorio por el apropiado
%root_path='/Users/aredondas/CODE/rbcce.aemet.es/iberonesia/RBCC_E/2019/Triad/'

function info_log=informe_lgl(n)

persistent in;
persistent brw;
persistent ano;
persistent root_path;

root_path='/Users/aredondas/CODE/rbcce.aemet.es/iberonesia/RBCC_E/Triad/Langley'
brw=[157,183,185];
tic
for ano=2016:2019
    pf=strrep('/Users/aredondas/CODE/rbcce.aemet.es/iberonesia/RBCC_E/2019/Triad/Langley/load_data.m','2019',num2str(ano));
    options_pub.outputDir=strrep('/Users/aredondas/CODE/rbcce.aemet.es/iberonesia/RBCC_E/2019/Triad/latex','2019',num2str(ano))
    options_pub.showCode=true;
    if exist(pf,'file')
        fd=strrep('/Users/aredondas/CODE/rbcce.aemet.es/iberonesia/RBCC_E/2019/Triad/Langley','2019',num2str(ano));
        wd=cd(fd)
        publish(pf,options_pub);
        disp(['OK ',' ',pf]);
    end
    for in=n(1):n(2)
         try
           pf=strrep(strrep('/Users/aredondas/CODE/rbcce.aemet.es/iberonesia/RBCC_E/2019/Triad/Langley/ozone/Langley157.m','157',num2str(brw(in))),'2019',num2str(ano));
            options_pub.outputDir=strrep('/Users/aredondas/CODE/rbcce.aemet.es/iberonesia/RBCC_E/2019/Triad/latex','2019',num2str(ano))
            options_pub.showCode=true;
            if exist(pf,'file')
                fd=strrep('/Users/aredondas/CODE/rbcce.aemet.es/iberonesia/RBCC_E/2019/Triad/Langley/ozone/','2019',num2str(ano));
                wd=cd(fd)
                publish(pf,options_pub);
                disp(['OK ',' ',pf]);
            end
        catch
            close all;
            disp(['ERROR',' '])
            w_d=cd('/Users/aredondas/CODE/rbcce.aemet.es/iberonesia/RBCC_E/Triad/Langley')
        end
        toc
    end
    toc
end




