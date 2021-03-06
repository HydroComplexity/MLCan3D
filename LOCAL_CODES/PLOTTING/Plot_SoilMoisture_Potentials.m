

vl_mn = min( min(min(volliq_store)) );
vl_mx = max( max(max(volliq_store)) );
rpp_mn = min( min(min(rppMPA_store)) );
rpp_mx = max( max(max(rppMPA_store)) );
smp_mn = min( min(min(smpMPA_store)) );
smp_mx = max( max(max(smpMPA_store)) );
figure(fignum); clf
    subplot(4,1,1)
        pcolor(timevect, -zns, volliq_store) 
        colormap(flipdim(jet,1))
        caxis([vl_mn vl_mx])
        shading interp
        colorbar
        ylabel('z [m]')
        title('SOIL MOISTURE [m^3 m^{-3}]')
    subplot(4,1,2)
        pcolor(timevect, -zns, rppMPA_store)
        caxis([rpp_mn rpp_mx])
        shading interp
        colorbar
        ylabel('z [m]')
        title('ROOT PRESSURE POTENTIAL [MPa]')
    subplot(4,1,3)
        pcolor(timevect, -zns, smpMPA_store)
        caxis([smp_mn smp_mx])
        shading interp
        colorbar
        ylabel('z [m]')
        title('SOIL MOISTURE POTENTIAL [MPa]')
    subplot(4,1,4)
        pcolor(timevect, -zns, rppMPA_store-smpMPA_store)
        %caxis([smp_mn smp_mx])
        shading interp
        colorbar
        ylabel('z [m]')
        title('(ROOT - SOIL) MOISTURE POTENTIAL [MPa]')

