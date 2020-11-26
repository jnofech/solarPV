function generate_tables(name,clrsky,folder_output)
    filename = name + '_info.mat';                 % (!) Joseph edit
    global num alpha_p  
%     clrsky=readtable('Clear-sky-index.xlsx');
%     clrsky=table2array(clrsky);
%         clrsky = [0;clrsky];        % (!) Joseph edit. I assume MATLAB skips the first line of the '.xlsx' file by treating it as a header.
    % matfiles=dir(fullfile('*.mat')); %% access to mat files
%         folder_output = 'output/transfer/';                 % (!) Joseph edit
        matfiles=dir(fullfile(folder_output,filename));       % (!) Joseph edit
    l_mat=length(matfiles); %% length of matfiles
%     names=strings(l_mat,1);
    roofs_num=0;
    
    for j=1:l_mat %% loop for buildings
%         fid=fopen(matfiles(j).name); %% opening the mat files in a loop
        fid=fopen(filename); %% opening the mat files in a loop
%         h=matfiles(j).name; %% extracting names of buidlings
%         h=h(1:end-9); %% purification of the names
        
%         names(j,1)=h; %% names matrix
    %     q=load(matfiles(j).name); %% loading each building table
            name;                                        % (!) Joseph edit
            q=load(strcat(folder_output,filename)); %% loading each building table      % (!) Joseph edit
        roof_data=q.table;
        imax=height(roof_data);
        v_array=roof_data.polygons; %% extracting vertices arrays
        type=roof_data.type; %% extracting type of each blob
        slope=roof_data.slope; %% extracting slope of each blob
        pix_size=q.pix2m; %% use pixel size to convert coordinates to the meter
        num_panels=0; %% defining a variable for each building panels number
        for i=1:imax
            v=cell2mat(v_array(i)); %% vertices
            v(isnan(v))=0;
            if length(v)<3
                v(2,1)=0;
                v(2,2)=0;
            end
            if type(i)=="flat"

                l_v=length(v);
                for k=1:l_v
                    x_v(k)=v(k,1)*pix_size;
                    y_v(k)=-1*v(k,2)*pix_size;
                end
                x_v(l_v+1)=v(1,1)*pix_size;
                y_v(l_v+1)=-1*v(1,2)*pix_size;
    %             i
    %             figure
    %             plot(x_v,y_v,'LineWidth',2) % polygon 
    %             axis equal
    %             hold on            
                [num, numbers_panel, alpha_p, alpha_s, az_solar, G, G_cs, Generation, omega, max_Gen]=Panels_number(x_v, y_v,clrsky);
                num;
                num_panels=num_panels+num;    
                numbers_panel;
                roofs_num=roofs_num+1;            
%                 T(roofs_num,1)=names(j,1);
                T(roofs_num,1)=name;
                T(roofs_num,2)=i;
                T(roofs_num,3)=num;
                T(roofs_num,4)=max_Gen;            
                for k=1:l_v
                    x_v(:)=0;
                    y_v(:)=0;
                end
            end       
        end   
        num_panels;
        num_panels=0;

    end
    names_var(1,1)="Building";names_var(1,2)="Roof";names_var(1,3)="Number";names_var(1,4)="Generation(kWh)";
    T1=array2table(T,'VariableNames',names_var);
    % writetable(T1,'Genration-number.xlsx')
        writetable(T1,folder_output+name+'_Genration-number.xlsx')
end

function [num, numbers_panel, alpha_p, alpha_s, az_solar, G, G_cs, Generation, omega, max_Gen]=Panels_number(x_v, y_v,clrsky)    
    x_min=min(x_v)+1.0;
    x_max=max(x_v);
    y_min=min(y_v)+1.0;
    y_max=max(y_v);
    dist_max=sqrt((x_min-x_max)^2+(y_min-y_max)^2);
    delta_x=0.992;
    delta_y=2.078;
    eff=0.204;
    P=420;
    i_r_variable=[0:0.1:1];
    beta=[0:10:90];    
    l_beta=length(beta);
    l_i_r_variable=length(i_r_variable);
    i_r=zeros(l_beta, l_i_r_variable);
    for i=1:l_beta
        beta_tilt=beta(i);        
        [alpha_p, G_cs, alpha_s, az_solar, omega]=Power (beta_tilt, clrsky); 
        l_alpha_p=length(alpha_p);
        for j=1:l_i_r_variable
            i_r_constant=0.6;           
            x_q = x_min-dist_max:delta_x:x_min+dist_max;          
            i_r_constant=i_r_constant+i_r_variable(j);
            y_q(1)=y_min-dist_max;
            k=1;
            while (y_q(k)<y_min+dist_max)
                k=k+1;
                if rem(k,2)==0
                    y_q(k)=y_q(k-1)+delta_y*cosd(beta(i));
                else                    
                    %y_q(k)=y_q(k-1)+delta_y+i_r_constant-delta_y*cosd(beta(i));
                    y_q(k)=y_q(k-1)+i_r_constant;
                end
            end
            y_q(end)=[];
            %SEA_new=atand(delta_y*sind(beta(i))/(delta_y+i_r_constant-delta_y*cosd(beta(i))));
            SEA_new=atand(delta_y*sind(beta(i))/i_r_constant);            
                for kk=1:length(alpha_p)
                    if alpha_p(kk)>0 & alpha_p(kk)<SEA_new & alpha_s(kk)>0
                        shadow_length(kk)=delta_y.*sind(beta(i)).*sind(SEA_new-alpha_p(kk))./(sind(beta(i)+alpha_p(kk)).*sind(SEA_new));
                    else
                        shadow_length(kk)=0.0;
                    end
                end
                shadow_length=max(shadow_length,0);
                shadow_length=min(shadow_length,delta_y);
                       
            theta=0;
            R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
            [X,Y]=meshgrid(x_q,y_q);
            XY=[X(:) Y(:)];
            rotXY=XY*R';
            Xqr = reshape(rotXY(:,1), size(X,1), []);
            Yqr = reshape(rotXY(:,2), size(Y,1), []);
            [in,on] = inpolygon(Xqr,Yqr,x_v,y_v);
%             hold on
            in = inpolygon(Xqr,Yqr,x_v,y_v);
            in(on)=[0];
            in(:,end+1)=0;
            in(end+1,:)=0;
            in_mat=in;            
            n=0;
            m=0;
            n_x=length(x_q);
            n_y=length(y_q); 
            %figure
            for ii=1:2:n_y
                for jj=1:n_x
                    if in_mat(ii,jj)==1 & in_mat(ii+1,jj)==1 & in_mat(ii,jj+1)==1 & in_mat(ii+1,jj+1)==1
                        x=[Xqr(ii,jj) Xqr(ii,jj+1) Xqr(ii+1,jj+1) Xqr(ii+1,jj) Xqr(ii,jj)];
                        y=[Yqr(ii,jj) Yqr(ii,jj+1) Yqr(ii+1,jj+1) Yqr(ii+1,jj) Yqr(ii,jj)];
                        n=n+1;
                    else
                        n=n;
                    end
                end
            end
            for b=1:n_x
                x_q(b)=0;
            end
            for b=1:n_y
                y_q(b)=0;
            end
            numbers_panel(i,j)=n;
            num=n;            
            G=min(G_cs*delta_y*delta_x*num*eff/1000,num*P/1000);
            %G_firstrow=min(G_cs*delta_y*delta_x*m*eff/1000,num_first*P/1000);
            Generation(i,j)=sum(G.*(delta_y-shadow_length)/(delta_y*60));
            Gen_per_number(i,j)=Generation(i,j)./num;
            %h=sum(G)
            
        end
    end
    numbers_panel;
    num=max(numbers_panel(:));
    Generation;
    max_Gen=max(Generation(:));
    Gen_per_number;
    max_Gen_per_number=max(Gen_per_number(:));
    %%---------------------------------------------------------------------
    %%---------------------------------------------------------------------

end

function [alpha_p, G_cs, alpha_s, az_solar, omega]=Power (beta_tilt, clrsky)
    SM=7; LM=103; phi=53.625447;
    G_sc=1366.1; %% Solar Constant
    beta=beta_tilt; %% tilt angle
    az=0; %% azimuth angle
    DOY=1:1:366; %% Days of Year
    DOY_m=1:1/1440:(367-(1/1440)); %% Days of Year (minutely)
    EH=0:1/60:(24-(1/60)); %% Exact hour (standard time) for a day
    EH_new=EH;
    for i=1:365
        EH_new=[EH_new,EH]; 
    end
    B=360*(DOY_m-81)/365; %% A coefficient used in the calculations
    EoT=9.87*sind(2*B)-7.53*cosd(B)-1.5*sind(B); %% Equation of Time
    %%ST=EH_new+((SM-LM).*4+EoT)/60; %% Solar Time *****Excel???
    if DOY_m>=72 & DOY_m<=310 %% Solar Time based on Excel data
        ST=EH_new+(4*6+EoT)/60;
    else
        ST=EH_new+(4*7+EoT)/60;
    end
    %%decl=23.45*sind(B); %% Declination angle 
    decl=-23.45*cosd(0.986*(DOY_m+10.5)); %% Declination angle based on excel
    omega=(ST-12)*15; %% Hour Angle
    zenith=acosd(cosd(phi).*cosd(omega).*cosd(decl)+sind(phi).*sind(decl));
    alpha_s=90-zenith;
    %teta=acosd(sind(decl).*sind(phi).*cosd(beta)-sind(decl).*cosd(phi).*sind(beta).*cosd(az) ...
    %    +cosd(decl).*cosd(phi).*cosd(beta).*cosd(omega)+cosd(decl).*sind(phi).*sind(beta).*cosd(az).*cosd(omega) ...
    %    +cosd(decl).*sind(beta).*sind(az).*sind(omega)); %% Calculating the Incident Angle
    teta=acosd(sind(decl).*sind(phi).*cosd(beta)-sind(decl).*cosd(phi).*sind(beta).*cosd(az) ...
        +cosd(decl).*cosd(phi).*cosd(beta).*cosd(omega)+cosd(decl).*sind(phi).*sind(beta).*cosd(az).*cosd(omega) ...
        +cosd(decl).*sind(beta).*sind(az).*sind(omega)); %% Calculating the Incident Angle
    az_solar=sign(omega).*abs(acosd((cosd(zenith).*sind(phi)-sind(decl))./(sind(zenith).*cosd(phi)))); %% solar azimuth angle
    alpha_p=atand(tand(90-zenith)./cosd(az_solar-az)); %% profile angle
    rb=cosd(teta)./cosd(zenith); %% The ratio of beam radiation on the tilted surface to that on a horizontal surface
    for i=1:length(rb)
        if(rb(i)<0)
            rb(i)=0.0;
        else
            rb(i)=rb(i);
        end
    end
    %G_ex=G_sc*(1.000110+0.034221*cosd(B)+0.001280*sind(B)+ ...
    %    0.000719*cosd(2*B)+0.000077*sind(2*B)); %% Extraterrestrial radiation incidient on the plane normal to the radiation
    G_ex=G_sc*(1.000110+0.034221*cosd(B)+0.001280*sind(B)- ...
        0.000719*cosd(2*B)+0.000077*sind(2*B)); %% (Excel)Extraterrestrial radiation incidient on the plane normal to the radiation
    GHI_cs=0.7.*G_ex.*cosd(zenith); %% clear-sky global horizontal irradiance
    for i=1:length(GHI_cs)
        if (GHI_cs(i)<0)
            GHI_cs(i)=0;
        else
            GHI_cs(i)=GHI_cs(i);
        end
    end
    Beam=rb.*GHI_cs; %% Direct or Beam Radiation on Tilted Surface (clear sky)
    G_diff=14.29+21.04*((pi/2)-(zenith*pi/180));
    for i=1:length(G_diff)
        if(G_diff(i)<0)
            G_diff(i)=0.0;
        else
            G_diff(i)=G_diff(i);
        end
    end
    R_diff=(3+cosd(2*beta))/4; %% coefficient for diffuse on the tilted surface
    Diffuse=G_diff.*R_diff; %% Diffuse Radiation on Tilted Surface (clear sky)
    Reflected=0.5*GHI_cs*(1-cosd(beta))/2; %% Reflected Radiation on Tilted Surface (clear sky)
    G_cs=(Beam+Diffuse+Reflected).*clrsky'; %% Global Irradiance on Tilted Surface     

end