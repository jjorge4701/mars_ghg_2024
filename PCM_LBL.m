% PCM_LBL.m

% This script runs PCM_LBL, a 1D radiative-convective code designed to simulate the climates of diverse planetary atmospheres written by R. Wordsworth (2017-2021) with additional input from F. Ding.
% It assumes an atmosphere composed of CO2, H2O, and one additional greenhouse gas. Other parameters are set to reflect an early Martian atmosphere. 
% Inputs: 
%   - p_array: maximum surface pressure [Pa]
%   - gas_array: gases of model atmosphere (H40 = H2O2, CH6 = C2H4, HNO = HNO3, HCO = H2CO, CH8 = C2H6)
%   - conc_array: desired concentrations of each gas in atmosphere [mol/mol] in scientific notation
%   - conc_ppm_array: strings of concentrations
%
% Outputs:
%   - PCM_LBL_data: MAT file saved with gas name, gas concentration, and surface pressyre
%

%% Inputs
clear all

p_array = {'1e5'}; % [Pa], one pressure at a time!

gas_array = {'H4O','SO2','NH3','CH6','HNO','O3_','N2O','CO_','CH4','NO2', 'HBr', 'OCS', 'HCO', 'HCN','H2S','CH8'}; 

conc_array = {'1e-9','5e-9','1e-8','2.5-8','5e-8','7.5e-8','1e-7','5e-7','1e-6','5e-6','1e-5','2.5e-5','5e-5','7.5e-5','1e-4','1.5e-4','2e-4','2.5e-4','3.5e-4','5e-4'}; % [mol/mol], gas concentrations (ppm in scientific not.)
conc_ppm_array =  {'0_001ppm','0_005ppm','0_01ppm','0_025ppm','0_05ppm','0_075ppm','0_1ppm','0_5ppm','1ppm','5ppm','10ppm','25ppm','50ppm','75ppm','100ppm','150ppm','200ppm','250ppm','350ppm','500ppm'}; %ppm name array for gas concentrations

filepath = '/Users/jasonjorge/Thesis/Mars_Code/PCM_LBL/example_run/input.nml'; %filepath to input.nml file
filepath_2 = '/Users/jasonjorge/Thesis/Mars_Code/PCM_LBL/example_run/results'; %filepath to results folder
filepath_3 = '/Users/jasonjorge/Thesis/Mars_Code/PCM_LBL/example_run/PCM_LBL.e'; %filepath to PCM_LBL executable 

%% Edit Input File for PCM_LBL, Run PCM_LBL

pstr_int = 'ps                ='; % maximum surface pressure string in input.nml file
gstr_int = 'gas_name_MAX      = '; % gas name string
cstr_int = 'gas_molarconc_MAX = '; % gas concentration string 
sig_int =  'calc_sigma       = ';  % calc sigma string 

for k = 1:numel(p_array) %loop over maximum surface pressures, gases, and concentrations
    for l = 1:numel(gas_array)
        for j = 1:numel(conc_array)
          
            fid = fopen(filepath,'r+'); %open input.nml

            lines = textscan(fid,'%s','delimiter','\n'); %read input.nml into lines variable
            lines = lines{1};

            p_idx = find(contains(lines,pstr_int)); %find location of strings of interest 
            g_idx = find(contains(lines,gstr_int));
            c_idx = find(contains(lines,cstr_int));
            sig_idx = find(contains(lines,sig_int));

            lines{p_idx} = [pstr_int,p_array{k},',']; %change pressure in lines 
            lines{g_idx} = [gstr_int,'''CO2'',''',gas_array{l},''',''H2O'',']; %change gas
            lines{c_idx} = [cstr_int,'1,',conc_array{j},',-1/']; %change concentration
            
            % change calc_sigma variable in input file to match what cross-section data is stored for each gas 
            % at new maximum surface pressures, all cross-sections need to be calculated

            % !!uncomment below if stored CO2/H2O sigma data!! 

            % if j == 1
            %     lines{sig_idx} = [sig_int,'F,T,F,']; %if on first concentration of gas, calculate sigma
            % else
            %     lines{sig_idx} = [sig_int,'F,F,F']; %no need to recalculate sigma for other concentrations of same gas
            % end
            
            % !!uncomment below if no stored CO2/H2O sigma data!!

             % if l == 1 && j == 1 
             %     lines{sig_idx} = [sig_int,'T,T,T,']; %if on first concentration of first gas, calculate sigma of CO2, H2O, and new gas
             % elseif l == 1 && j > 1 
             %     lines{sig_idx} = [sig_int,'F,F,F,']; %no need to recalculate sigma for other concentrations of first gas
             % elseif l > 1 && j == 1 
             %     lines{sig_idx} = [sig_int,'F,T,F,']; %if on first concentration of second gas onward, calculate sigma for new gas
             % elseif l > 1 && j > 1
             %     lines{sig_idx} = [sig_int,'F,F,F,']; %no need to recalculate sigma for other concentrations of same gas
             % end

            fid2 = fopen('input.nml', 'w'); %open input.nml file for writing

            for i = 1:numel(lines) %rewrite new information in lines into input.nml

                fprintf(fid2,'%s\n', lines{i});

            end

            fclose(fid2); 

            system(filepath_3) %run PCM_LBL with new gas and concentration

            variable_array = [dir('/Users/jasonjorge/Thesis/Mars_Code/PCM_LBL/example_run/results/*.out')]; %make variable array from PCM_LBL results folder

            for h = 1:numel(variable_array)

                Fid = fopen(['/Users/jasonjorge/Thesis/Mars_Code/PCM_LBL/example_run/results/',variable_array(h).name],'r'); %open each variable in results folder
                data_array = textscan(Fid,'%s'); %scan results into data_array
                data_array = str2double(data_array{1}); %convert strings to numbers
                data_array = data_array(~isnan(data_array)); %filter out nans
                PCM_LBL_data.(variable_array(h).name(1:end-4)) = data_array; %save data_array of each variable into structure 'PCM_LBL_data'
                fclose(Fid);
                
            end
            
            save([conc_ppm_array{j},'-',gas_array{l},'-',p_array{k}],"PCM_LBL_data") %save PCM_LBL_data as mat file with concentration, gas, and surface pressure
        
        end % end concentration loop 
    end % end gas loop
end % end pressure loop
 
%% Load PCM_LBL MAT files into Data Structure

clear data

mat_dir = dir('/Users/jasonjorge/Thesis/Mars_Code/PCM_LBL/example_run/*.mat'); %get directory of all MAT files when simulations done

for g = 1:numel(mat_dir)
    
    name_array = regexp(mat_dir(g).name,'\w*','match'); %split name of mat file

    fieldname_str = [name_array{2},'_',name_array{1},'_',name_array{3}]; %create fieldname string of gas, conconcentration, and pressure

    data.(fieldname_str) = importdata(mat_dir(g).name); %import data from mat files into data structure 
   
end

%% Extract Surface Temps and Concentrations from Data Structure

fieldnames_array = fieldnames(data); %generate fieldnames of data array with all mat files
temp_conc_array = cell(numel(conc_ppm_array),numel(gas_array)); %create empty cell structure for temperature vs. concentration for each gas

for f = 1:numel(fieldnames_array)

    clear col_num
    clear row_num

        ppm_cell = extractBetween(fieldnames_array{f},"_","ppm",'Boundaries','Inclusive'); %extract ppm string from fieldname
        ppm_str = char(ppm_cell);
        ppm_str_edit = ppm_str(2:end);
        
        if numel(strfind(ppm_str_edit,'_')) > 1 %remove extra _ if needed

            ppm_str_edit2 = ppm_str_edit(2:end);
            row_num = find(strcmp(conc_ppm_array,ppm_str_edit2)); %find row associated with gas concentration
        else
            if strcmp(ppm_str_edit(1),'_')
                new_ppm_str = ppm_str_edit(2:end);
                row_num = find(strcmp(conc_ppm_array,new_ppm_str)); %find row associated with gas concentration
            else
                row_num = find(strcmp(conc_ppm_array,ppm_str_edit));
            end
        end

        col_num = find(strcmp(gas_array,fieldnames_array{f}(1:3))); %find column associated with gas

        temp_conc_array{row_num,col_num} = data.(fieldnames_array{f}).Tsurf_final; %add final Tsurf to appropriate spot on table


end

gas_name_array = {'H_2O_2','SO_2','NH_3','C_2H_4','HNO_3','O_3','N_2O','CO','CH_4','NO_2', 'HBr', 'OCS', 'H_2CO', 'HCN','H_2S','C_2H_6'};
Tsurf_data_table = cell2table(temp_conc_array,'VariableNames',gas_name_array,'RowNames',conc_ppm_array); %turn cell array into table






