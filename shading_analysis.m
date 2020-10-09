function [shadow_score] = shading_analysis(b,lat,pix_length,scale,force_run,output_path_save)
% SHADING_ANALYSIS() analyzes the shading of rooftops over the course of a
% full year by directly raytracing from the Sun's position on the two
% solstices and an equinox.
% Output is saved into "output/"+b.name+"_shadows_"+length+".mat", and is
% not re-generated while this file exists.
%
% Returns:
% --------
%   shadow_score : double (2D map, size of b.Z)
%       Relative amount of time that a region spends in SHADE.
%       Lower scores (i.e. more sunlight) are better.

% Shadow study - Full Year
    dates = [81, 172, 264, 355];     % Mar 21, Jun 21, Sep 21, Dec 21 (Mar 21, Sep 21 are identical)
    hours = 1:24;               % Hourly checks
    
    % Generate building
    name = b.name;
    output_shadows_name           = output_path_save+b.name+"_shadows_"+pix_length+".mat";
    
    % Shadow study!
    shadows = ones([size(b.Z),numel(dates)*numel(hours)],'logical');
%     tic

    if ~isfile(output_shadows_name) || force_run==true
        disp("shading_analysis() : Analyzing shadows throughout the year for "+name+".");
        for i = 1:numel(dates)
            day = dates(i);
            dec = 23.45*sind(360/365*(284+day));    % (approximate) Solar declination; ranging from 23.45 to -23.45
            H_angle_rise = -acosd( -tand(lat) * tand(dec) );    % Hour angle at sunrise
            hour_rise = H_angle_rise / 60 / 0.25 + 12.00;       % Hour at sunrise

            for j = 1:numel(hours)
                hour_AST = hours(j);
                if abs(hour_AST - 12.00) < (abs(hour_rise - 12) - 1.5)
                    % Perform shadow analysis
                    % Get solar angles
                    H_angle = 0.25*(hour_AST - 12.00)*60;   % Hour angle; i.e. 0.25deg per minute since Noon
                    dec = 23.45*sind(360/365*(284+day+hour_AST/24));    % Solar declination; ranging from 23.45 to -23.45
                    alt_noon = 90-lat+dec;                  % Altitude of sun in sky AT SOLAR NOON
                    alt = asind( cosd(lat)*cosd(dec).*cosd(H_angle) + sind(lat)*sind(dec) ); % Altitude
                    azi = asind( sind(H_angle).*cosd(dec)/cosd(alt) );                       % Azimuth (deg CW from South)
                    % Perform shading!
%                     disp("Shading "+b.name+" on Day "+day+", at Hour "+hour_AST+"...");
                    if day==264
                        % Mar 21 and Sep 21 are identical; recycle!
                        shadows(:,:,(i-1)*numel(hours)+j) = shadows(:,:,(i-3)*numel(hours)+j);
                    else
                        shadows(:,:,(i-1)*numel(hours)+j) = shading_raytrace(b,alt,azi,scale);  % 1 where shadow exists; 0 where illuminated.
                    end
                else
                    % Skip
                end
            end
        end
        save(output_shadows_name,'shadows');
        disp("shading_analysis() : "+name+" shading analysis complete!");
    end
    load(output_shadows_name,'shadows');
    shadow_score = sum(shadows,3) / size(shadows,3);
%     toc
end