%  Plot the srElevation structure

if srGeometry.tf_latlon
    subplot(121)
    pcolor(srElevation.longitude, srElevation.latitude, srElevation.data' )
    shading interp
    axis image
    hold on
    contour(srElevation.longitude, srElevation.latitude, srElevation.data' ,'k')
    plot(srStation.longitude, srStation.latitude,'ob', ...
        'MarkerSize',10, ...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','r')
    
    subplot(122)
    contour(srElevation.longitude,srElevation.latitude,srElevation.data' );axis image;hold on;
    plot(srStation.longitude,srStation.latitude,'ob','MarkerSize',10,...
        'MarkerEdgeColor','k','MarkerFaceColor','r');
    
elseif ~srGeometry.tf_latlon
    subplot(121)
    pcolor(srElevation.easting, srElevation.northing, srElevation.data' )
    shading interp
    axis image
    hold on
    contour(srElevation.easting, srElevation.northing, srElevation.data' ,'k')
    plot(srStation.longitude, srStation.latitude,'ob', ...
        'MarkerSize',10, ...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','r')
    
    subplot(122)
    contour(srElevation.easting,srElevation.northing,srElevation.data' );axis image;hold on;
    plot(srStation.easting,srStation.northing,'ob','MarkerSize',10,...
        'MarkerEdgeColor','k','MarkerFaceColor','r');
end