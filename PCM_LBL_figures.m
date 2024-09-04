%% Isolate Tsurf for Warmest Gases

gas_interest = {'H4O','HNO','SO2','CH6','NH3'}; %5 warmest gases
fieldnames_array = fieldnames(data);
Index=false(size(fieldnames_array));
pressure_name = {'5e4','2e5'};

for k=1:length(gas_interest)
    for j=1:length(fieldnames_array)
        if contains(fieldnames_array{j},gas_interest{k})
            Index(j)=true; %identify which mat files are associated with gases of interest
        end
    end
end

goi_idx = find(Index);

temp_conc_array = cell(numel(conc_ppm_array),numel(gas_interest),3);

for l=1:length(goi_idx)

    ppm_cell = extractBetween(fieldnames_array{goi_idx(l)},"_","ppm",'Boundaries','Inclusive');
    ppm_str = char(ppm_cell);
    ppm_str = ppm_str(2:end);

    for j=1:length(conc_ppm_array)
        for m=1:length(gas_interest)
          

                clear col_num
                clear row_num
                clear z_num

                if strcmp(ppm_str,conc_ppm_array{j})
                    row_num = j;
                end

                if contains(fieldnames_array{goi_idx(l)},gas_interest{m})
                    col_num = m;
                end

                z_var = find(strcmp(fieldnames_array{goi_idx(l)}(end-2:end),pressure_name));

                if isempty(z_var)
                    z_num = 2;
                else
                    if z_var == 2
                        z_num = 3;
                    elseif z_var == 1
                        z_num = 1;
                    end
                end

                if exist('col_num')
                    if exist('row_num')
                        if exist('z_num')



                            if z_num == 2
                                temp_conc_array{row_num,col_num,2} = data.(fieldnames_array{goi_idx(l)}).Tsurf_final;
                            elseif z_num == 3
                                temp_conc_array{row_num,col_num,3} = data.(fieldnames_array{goi_idx(l)}).Tsurf_final;
                            elseif z_num == 1
                                temp_conc_array{row_num,col_num,1} = data.(fieldnames_array{goi_idx(l)}).Tsurf_final;
                            end


                        end
                        
                    end
                end

        end
    end
end

Tsurf_table1 = cell2table(temp_conc_array(:,:,1),'VariableNames',gas_interest,'RowNames',conc_ppm_array);
Tsurf_table2 = cell2table(temp_conc_array(:,:,2),'VariableNames',gas_interest,'RowNames',conc_ppm_array);
Tsurf_table3 = cell2table(temp_conc_array(:,:,3),'VariableNames',gas_interest,'RowNames',conc_ppm_array);

%% Plot Surface Temp vs Conc for all Gas Species

gas_name_array1 = {'H_2O_2','SO_2','NH_3','C_2H_4','HNO_3','O_3','N_2O','CO','CH_4','NO_2', 'HBr', 'OCS', 'H_2CO', 'HCN','H_2S','C_2H_6'};
gas_name_array2 = {'O_3','N_2O','CO','CH_4','NO_2', 'HBr', 'OCS', 'H_2CO', 'HCN','H_2S','C_2H_6'};
conc_ppm1 = [0.001, 0.005, 0.01, 0.025, 0.05, 0.075, 0.1, 0.5, 1, 5, 10, 25, 50, 75, 100, 150, 200, 250, 350, 500];
%conc_ppm2 = [0.1, 0.5, 1, 5, 10, 25, 50, 75, 100, 150, 200, 250, 350, 500];
%t_surf_idx = ~cellfun(@isempty,Tsurf_data_table{:,7});

figure()
hold on;

hf1 = subplot(1,3,1);

for k = 1:5
    hold on;
    plot(conc_ppm1,Tsurf_data_table{:,k},'LineWidth',2.5)
    set(gca,'XScale','log')
    legend(gas_name_array1{1:5},'location','northwest')
    set(gca,'FontSize',18,'linewidth',2.5)
    ylabel('Surface Temperature [K]')
    title('(a)', 'FontSize', 24);
    %set(gca,'XTick',([10^{-3} 10^{-2} 10^{-1} 10 10^1 10^2]))
    %xticks([10^{-3} 10^{-2} 10^{-1} 10 10^1 10^2])
    %xticklabels('10^{-3}','10^{-2}','10^{-1}','10^2','10^3')
    %xticklabels('manual')
    set(gca,'XLim',[0 1000])
    xticks([10^-3 10^-2 10^-1 10^0 10^1 10^2 10^3])
    %xticklabels('10^{-3}','10^{-2}','10^{-1}','10^{0}','10^{1}','10^{2}','10^{3}')

end

hf2 = subplot(1,3,2);

for j = 6:10
    hold on;
    plot(conc_ppm1,Tsurf_data_table{:,j},'LineWidth',2.5)
    set(gca,'XScale','log')
    legend(gas_name_array1{6:10})
     set(gca,'FontSize',18,'linewidth',2.5)
    xlabel('Gas Concentration [ppmv]')
    title('(b)', 'FontSize', 24);
    set(gca,'XLim',[0 1000])
    %xticklabels('manual')
    xticks([10^-3 10^-2 10^-1 10^0 10^1 10^2 10^3])
    %xticklabels('10^{-1}','10^{0}','10^{1}','10^{2}','10^{3}')

end

%title('Surface Temperature vs. Gas Concentration')
%subtitle('1 bar atmosphere of CO_2 and H_2O')

hf3 = subplot(1,3,3);

for f = 11:16
    hold on;
    plot(conc_ppm1,Tsurf_data_table{:,f},'LineWidth',2.5)
    set(gca,'XScale','log')
    legend(gas_name_array1{11:16})
    set(gca,'FontSize',18,'linewidth',2.5)
    title('(c)', 'FontSize', 24);
    set(gca,'XLim',[0 1000])
    %xticklabels('manual')
    xticks([10^-3 10^-2 10^-1 10^0 10^1 10^2 10^3])
    %xticklabels('10^{-1}','10^{0}','10^{1}','10^{2}','10^{3}')

end

%set(gca,'XScale','log')
%set(gca,'YLim',[210 330])
%set(gca,'FontSize',14)

%ylabel('Surface Temperature [K]')
%xlabel('Gas Concentration [ppm]')
linkaxes([hf1,hf2,hf3],'y')
%legend(gas_array)

%% Plot Surface Temp vs Conc for Gas Species of Interest

color_vec = ["#ca0020","#EC7745","#FFEF00","#92c5de","#0571b0"];
%conc_ppm = [0.1, 0.5, 1, 5, 10, 25, 50, 75, 100, 150, 200, 250, 350, 500];

figure()
hold on;
% 
% for h = 1:numel(gas_interest)
% 
%     plot(conc_ppm,Tsurf_table1.(h),"LineWidth",1.5,"LineStyle",":","Color",color_vec(h)) 
%     plot(conc_ppm,Tsurf_table2.(h),"LineWidth",1.5,"LineStyle","--","Color",color_vec(h)) 
%     plot(conc_ppm,Tsurf_table3.(h),"LineWidth",1.5,"LineStyle","-","Color",color_vec(h)) 
% 
% end


hf1=subplot(1,3,1);

for h = 1:numel(gas_interest)

       plot(conc_ppm1,Tsurf_table1.(h),"LineWidth",2.5,"LineStyle","-","Color",color_vec(h))   
       hold on;
       set(gca,'XScale','log')
       set(gca,'YLim',[210 330])
       set(gca,'FontSize',20,'linewidth',2.5)
       ylabel('Surface Temperature [K]')
       title('(a) 0.5 bar')
       xticks([10^-3 10^-2 10^-1 10^0 10^1 10^2 10^3])
       set(gca,'XLim',[0 1000])

end

hf2=subplot(1,3,2);

for g = 1:numel(gas_interest)

       plot(conc_ppm1,Tsurf_table2.(g),"LineWidth",2.5,"LineStyle","-","Color",color_vec(g)) 
       hold on;
        
       set(gca,'XScale','log')
       set(gca,'YLim',[210 330])
       set(gca,'FontSize',20,'linewidth',2.5)
       xlabel('Gas Concentration [ppmv]')
       title('(b) 1 bar')
       xticks([10^-3 10^-2 10^-1 10^0 10^1 10^2 10^3])
       set(gca,'XLim',[0 1000])
end

hf3=subplot(1,3,3);


for v = 1:numel(gas_interest)

       plot(conc_ppm1,Tsurf_table3.(v),"LineWidth",2.5,"LineStyle","-","Color",color_vec(v))  
       hold on;
       set(gca,'XScale','log')
       set(gca,'YLim',[210 330])
       set(gca,'FontSize',20,'linewidth',2.5)
       title('(c) 2 bar')
       xticks([10^-3 10^-2 10^-1 10^0 10^1 10^2 10^3])
       set(gca,'XLim',[0 1000])
end

%plot(conc_ppm,Tsurf_table2.(h),"LineWidth",1.5,"LineStyle","-.","Color",color_vec(h))
 %      plot(conc_ppm,Tsurf_table3.(h),"LineWidth",1.5,"LineStyle","--","Color",color_vec(h))

%xlabel('Gas Concentration [ppm]')
%ylabel('Surface Temperature [K]')
%title('CO_2-H_2O')
goi_name = {'H_2O_2','HNO_3','SO_2','C_2H_4','NH_3'};
legend(goi_name,'Location','best')

%set(gca,'XScale','log')
%set(gca,'YLim',[210 330])
%
linkaxes([hf1,hf2,hf3],'')




%% Plot Absorption Cross Section vs Wavenumber

nS = 2000;
nLay = 50;
nrows = 4;
ncols = 5;

sigma_dir = dir ('/Users/jasonjorge/Thesis/Mars_Code/PCM_LBL/example_run/saved_sigma_data/*lw.dat');

for k = 1:numel(sigma_dir)

    read_gas_sigma = readtable(['/Users/jasonjorge/Thesis/Mars_Code/PCM_LBL/example_run/saved_sigma_data/',sigma_dir(k).name]);

    if numel(read_gas_sigma) == 2000
        gas_sigma = read_gas_sigma;
    else
        gas_sigma = read_gas_sigma.(3);
    end

    %gas_sigma = read_gas_sigma.(3);
    %gas_sigma_final = reshape(gas_sigma,nS,[]);
    %gas_sigma_final = gas_sigma_final';
    %gas_sigma_surf = gas_sigma_final(:,150);
    gas_data_array.(sigma_dir(k).name(1:end-4)) = gas_sigma; 

end

nu_vec = linspace(1,2500,2000);
fieldname_array = fieldnames(gas_data_array);

fig = figure();

for n = 1:length(fieldnames(gas_data_array))

    ax1 = subplot(nrows,ncols,n);
       
    plot(nu_vec,gas_data_array.(fieldname_array{n}))
    set(ax1,'YScale','log')
    set(ax1,'YLim',[10e-26 10e-17])
    set(ax1,'xLim',[0 2000])
    title(fieldname_array{n}(7:9))
end

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Absorption Cross Section [cm^2/molecule]');
xlabel(han,'Wavenumber [cm^{-1}]');

set(han.YLabel,'Position',[-0.0784 0.5000 0])

x = gas_data_array.nu_lw;

plot(nu_vec,gas_data_array.sigma_H4O_lw)



%% Set-up Blackbody Eqns for OLR Comp.

PCM_H4O_normTIPS = load('10ppm-H4O-.mat');
PCM_H4O_adjTIPS = load('10ppm-H4O-adjTIPS.mat');

c = 3e10; %  cm/s 
h = 6.626e-27;
kB = 1.38e-16;

nu_cm = (0:0.1:2500); %cm^-1 ❤️
nu = nu_cm*c; % cm^-1 to Hz
Tsurf1 = PCM_H4O_normTIPS.PCM_LBL_data.Tsurf_final;
Tsurf2 = PCM_H4O_adjTIPS.PCM_LBL_data.Tsurf_final;

B_temp1 = (2*h*nu.^3/c^2)./(exp(h*nu/(kB*Tsurf1)) - 1); % W/m2/sr/Hz
B1 = B_temp1*c*pi*(1/1000); %W/m2/cm^-1

B_temp2 = (2*h*nu.^3/c^2)./(exp(h*nu/(kB*Tsurf2)) - 1); % W/m2/sr/Hz
B2 = B_temp2*c*pi*(1/1000); %W/m2/cm^-1

figure()
plot(PCM_H4O_normTIPS.PCM_LBL_data.nu_lw,PCM_H4O_normTIPS.PCM_LBL_data.OLRnu,'linewidth',1.5,'Color',"#F6BE00")
hold on;
plot(PCM_H4O_adjTIPS.PCM_LBL_data.nu_lw,PCM_H4O_adjTIPS.PCM_LBL_data.OLRnu,'linewidth',1.5,'Color',"#ca0020")
plot(nu_cm,B1,'LineWidth',2.75,'Color',"#F6BE00")
plot(nu_cm,B2,'LineWidth',2.75,'Color',"#ca0020")

xlabel('Wavenumber [cm^{-1}]')
ylabel('OLR [W m^{-2} cm^{-1}] ')
legend('(T_{ref} / T)^{1.5}','(T_{ref} / T)^{2}',['BB(T=',num2str(Tsurf1,'%3.2f'),'K)'],['BB(T=',num2str(Tsurf2,'%3.2f'),'K)'],'NumColumns',2);
set(gca,'xLim',[0 2500])
set(gca,'fontsize',18,'linewidth',2.5)

%% Plot OLR 

figure()

subplot(2,1,1)
plot(PCM_H4O_normTIPS.PCM_LBL_data.nu_lw,PCM_H4O_normTIPS.PCM_LBL_data.OLRnu,'linewidth',1.5,'Color',"#F6BE00")
hold on;
plot(PCM_H4O_adjTIPS.PCM_LBL_data.nu_lw,PCM_H4O_adjTIPS.PCM_LBL_data.OLRnu,'linewidth',1.5,'Color',"#ca0020")
plot(nu_cm,B1,'LineWidth',2.75,'Color',"#F6BE00")
plot(nu_cm,B2,'LineWidth',2.75,'Color',"#ca0020")


xlabel('Wavenumber [cm^{-1}]')
ylabel('OLR [W m^{-2} cm^{-1}] ')
%title('Outgoing Longwave Radiation (OLR) vs. Wavenumber')
%subtitle('1 bar atmosphere of CO_2 and H_2O')
%legend('Q(T_{ref})/ Q(T) = (T_{ref}/T)^{1.5}','Q(T_{ref})/ Q(T) = (T_{ref}/T)^{2}'',['BB(T=',num2str(Tsurf1,%3.2f),'K)'],['BB(T=',num2str(Tsurf2,'%3.2f'),'K)'],'NumColumns',2)
set(gca,'xLim',[0 2500])
set(gca,'fontsize',18,'linewidth',2.5)
%nu_vec2 = linspace(1,2000,numel(gas_data_array.sigma_CO2_lw));


hf2 = subplot(2,1,2);
plot(table2array(gas_data_array.nu_lw),gas_data_array.sigma_CO2_lw,'Color','#911eb4','linewidth',1.5);
hold on;
plot(table2array(gas_data_array.nu_lw),gas_data_array.sigma_H2O_lw,'Color','#469990','linewidth',1.5);
plot(table2array(gas_data_array.nu_lw),gas_data_array.sigma_HCN_lw,'Color','b','linewidth',1.5)

set(hf2,'YScale','log')
set(hf2,'YLim',[10e-26 10e-16])
set(hf2,'xLim',[0 2500])
ylabel(hf2,'Absorption Cross Section [cm^2/molecule]');
xlabel(hf2,'Wavenumber [cm^{-1}]');
legend('CO_2','H_2O','HCN','location','best','LineWidth',2.5)
set(hf2,'fontsize',18 ,'linewidth',2.5)
title('Absorption Cross-Section vs. Wavenumber')
%title('H_2O Absorption Spectra using Air-Broadened and Estimated CO_2-Broadened Coefficnets')


%% Plot OLR 

figure();

c = 3e10; %  cm/s 
h = 6.626e-27;
kB = 1.38e-16; 

%plot(PCM_LBL_data.nu_lw,PCM_LBL_data.OLRnu)
hold on;
plot(PCM_LBL_data2.PCM_LBL_data.nu_lw,PCM_LBL_data2.PCM_LBL_data.OLRnu,'linewidth',1.25,'Color',color_vec(1))
%plot(PCM_LBL_data3.PCM_LBL_data.nu_lw,PCM_LBL_data3.PCM_LBL_data.OLRnu,'linewidth',1.25)
plot(PCM_LBL_data4.PCM_LBL_data.nu_lw,PCM_LBL_data4.PCM_LBL_data.OLRnu,'linewidth',1.25,'Color',color_vec(2))
plot(PCM_LBL_data5.PCM_LBL_data.nu_lw,PCM_LBL_data5.PCM_LBL_data.OLRnu,'linewidth',1.25,'Color',color_vec(3))
plot(PCM_LBL_data6.PCM_LBL_data.nu_lw,PCM_LBL_data6.PCM_LBL_data.OLRnu,'linewidth',1.25,'Color',color_vec(4))
plot(PCM_LBL_data7.PCM_LBL_data.nu_lw,PCM_LBL_data7.PCM_LBL_data.OLRnu,'linewidth',1.25,'Color',color_vec(5))

mat_dir = dir('/Users/jasonjorge/Thesis/Mars_Code/PCM_LBL/example_run/*.mat');

Tsurf_vec = [];

for u = 1:numel(mat_dir)

    load(mat_dir(u).name)
    mat_dir(u).name;
    Tsurf = PCM_LBL_data.Tsurf(end);
    Tsurf_vec = [Tsurf,Tsurf_vec];

end

nu_cm = (0:0.1:2500); %cm^-1 
nu = nu_cm*c; % cm^-1 to Hz

Tsurf_vec = sort(Tsurf_vec);

color_vec = ["#0571b0","#92c5de","#FFEF00","#f4a582","#ca0020"];

for l = 1:numel(Tsurf_vec)
 
        B_temp = (2*h*nu.^3/c^2)./(exp(h*nu/(kB*Tsurf_vec(l))) - 1); % W/m2/sr/Hz
        B = B_temp*c*pi*(1/1000); %W/m2/cm^-1
        plot(nu_cm,B,'color',color_vec(l),'LineWidth',2.75)

end

xlabel('Wavenumber [cm^{-1}]')
ylabel('OLR [W m^{-2} cm^{-1}] ')
title('Outgoing Longwave Radiation (OLR) vs. Wavenumber')
subtitle('1 bar atmosphere of CO_2, H_2O, and varying amounts of H_2O_2')
%legend('0 ppm','0.1 ppm','1 ppm','10 ppm',['BB(T=',num2str(Tsurf_vec(1),'%3.2f'),')'],...
  %  ['BB(T=',num2str(Tsurf_vec(2),'%3.2f'),')'],['BB(T=',num2str(Tsurf_vec(3),'%3.2f'),')'],...
   % ['BB(T=',num2str(Tsurf_vec(4),'%3.2f'),')'],'NumColumns',2)
legend('Air Broadening, CO_2-H_2O only','Air Broadening, 1 ppm H_2O_2','CO_2 Broadening, 1 ppm H_2O_2','Air Broadening, 10 ppm H_2O_2','CO_2 Broadening, 10 ppm H_2O_2',...
    ['BB(T=',num2str(Tsurf_vec(1),'%3.2f'),')'],['BB(T=',num2str(Tsurf_vec(2),'%3.2f'),')'],['BB(T=',num2str(Tsurf_vec(3),'%3.2f'),')'],['BB(T=',num2str(Tsurf_vec(4),'%3.2f'),')'],['BB(T=',num2str(Tsurf_vec(5),'%3.2f'),')'],'NumColumns',2)
ax = gca;
set(gca,'fontsize', 32,'linewidth',3)
%set(ax.Children(1),'LineStyle','--')
%set(ax.Children(2),'LineStyle',':')
%set(ax.Children(3),'LineStyle','-.')
%set(ax.Children(4),'LineStyle','-')

%% Plot Temperature vs Pressure 

figure();
hold on;
plot(PCM_CO2H2O.PCM_LBL_data.Tlev,PCM_CO2H2O.PCM_LBL_data.plev,'linewidth',2.25)
plot(PCM_NH3.PCM_LBL_data.Tlev,PCM_NH3.PCM_LBL_data.plev,'linewidth',2.25)
plot(PCM_CH6.PCM_LBL_data.Tlev,PCM_NH3.PCM_LBL_data.plev,'linewidth',2.25)
plot(PCM_SO2.PCM_LBL_data.Tlev,PCM_SO2.PCM_LBL_data.plev,'linewidth',2.25)
plot(PCM_HNO.PCM_LBL_data.Tlev,PCM_HNO.PCM_LBL_data.plev,'linewidth',2.25)
plot(PCM_H4O.PCM_LBL_data.Tlev,PCM_H4O.PCM_LBL_data.plev,'linewidth',2.25)

set(gca,'YScale','log')
set(gca,'YDir','reverse')
set(gca,'YLim',[5 1e5])
title('Temperature-Pressure Profiles at 1 bar')
%subtitle('1 bar atmosphere of CO_2 and H_2O, and varying amounts of H_2O_2')
xlabel('Temperature [K]')
ylabel('Pressure [Pa]')
legend('CO_2-H_2O only','NH_3','C_2H_4','SO_2','HNO_3','H_2O_2')
set(gca,'fontsize', 32)

%% Load HITRAN data with CO2 and Air HW 

filename = '/Users/jasonjorge/Downloads/SO2_CO2vsAir.txt';
filename2 = '/Users/jasonjorge/Downloads/N2O_CO2vsAir.txt';
delimiterIn = ' ';
headerlines = 1;
A = importdata(filename,delimiterIn);
B = importdata(filename2,delimiterIn,headerlines);

SO2_air_hw = A(:,6);
SO2_CO2_hw = A(:,13);
nu = A(:,3);

N2O_air_hw = cellfun(@str2num,B.textdata(2:end,6));
N2O_CO2_hw = cellfun(@str2num,B.textdata(2:end,13));
nu2 = cellfun(@str2num,B.textdata(2:end,3));

%% Plot CO2 and Air HW Comparisons

figure()
plot(nu,SO2_air_hw)
hold on;
plot(nu,SO2_CO2_hw)
legend('Air HW','CO2 HW')
title('SO2')
xlabel('Wavenumber (cm^-1)')
ylabel('Lorentzian HWHW (cm^-1 atm^-1)')


figure()
plot(nu2,N2O_air_hw)
hold on;
plot(nu2,N2O_CO2_hw)
legend('Air HW','CO2 HW')
title('N2O')
xlabel('Wavenumber (cm^-1)')
ylabel('Lorentzian HWHW (cm^-1 atm^-1)')

%% Plot Temperature and Pressure Grids

figure()

scatter(PCM_LBL_data.T_k_grid,PCM_LBL_data.p_k_grid,'LineWidth',1.5)
hold on;

plot(PCM_LBL_data.Tlev,PCM_LBL_data.plev,'LineWidth',2.5)

xlabel('Temperature [K]')
ylabel('Pressure [Pa]')
title('Absorption Cross Section Temperature-Pressure Grid')
legend('Absorption Cross Section Grid Point','Atmospheric Temperature Profile')

set(gca,'YScale','log')
set(gca,'YDir','reverse')
set(gca,'YLim',[1 2e5])
set(gca,'fontsize', 32)
