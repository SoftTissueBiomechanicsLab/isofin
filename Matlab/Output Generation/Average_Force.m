%% Plot force vs displacment and stress vs stretch 
clearvars; %close all; clc;
warning off

% This averaging script was used to average force data for figure 9
%% Add paths

addpath('..\Formulation\German Formulation\C_files');
addpath('..\Formulation\German Formulation\nurbs_toolbox');
addpath('..\Formulation\\German Formulation\With Torsion\Functions_RN Parallel\3D_New');
addpath('..\Input File Generation\Functions');

%% Pick fiber number, and use matching beam order and power
type=2;% 1 for 3D straight and 2 for 3D undulated and 3 for 2D straight
fiber_num=302;
deformation_type = 'BIAX';
titles = {'Straight Network','Undulated Network'};
network_cases = [3 18 19 22];
data_cell = cell(length(network_cases),1);
loop_ind = 0;
New_Case=1;
max_disp = 0;
for Case = network_cases
    loop_ind = loop_ind + 1;
    
    if type == 1
        savename = ['forcedisp/network_',num2str(fiber_num),'_',num2str(Case),'_',deformation_type];
    else
        savename = ['forcedisp/network_',num2str(fiber_num),'_',num2str(Case),'_',num2str(New_Case),'_',deformation_type];
    end

    load(savename,'displacement','force','force_y')
    force = [force' force_y'];
    
    if max_disp < max(displacement)
        max_disp = max(displacement);
    end
    stress = force/1e6;
    data_cell{loop_ind} = [displacement' stress];

end
%%
[displacement,stress,stress_std] = average_data(data_cell,max_disp);


pos_stress = stress + stress_std;
neg_stress = stress-stress_std;

figure
hold on
a1 = area(displacement*100,pos_stress(:,1));
a2 = area(displacement*100,neg_stress(:,1));
p1 = plot(displacement*100,stress(:,1));
p1.Color = [0.1 0.1 0.5];
p1.LineWidth = 2;
a1.FaceAlpha = 0.5;
a2.FaceColor = [1 1 1];
xlabel('Strain (%)')
ylabel('\sigma_X/E')
legend('Standard Deviation','','Mean Curve','location','northwest')
title([titles{type},': 300 Fibers Biaxial Loading'])
ylim([0 3.5e-5])
figure
hold on
a1 = area(displacement*100,pos_stress(:,2));
a2 = area(displacement*100,neg_stress(:,2));
p1 = plot(displacement*100,stress(:,2));
p1.Color = [0.1 0.1 0.5];
p1.LineWidth = 2;
a1.FaceAlpha = 0.5;
a2.FaceColor = [1 1 1];
xlabel('Strain (%)')
ylabel('\sigma_Y/E')
legend('Standard Deviation','','Mean Curve','location','northwest')
title([titles{type},': 300 Fibers Biaxial Loading'])
ylim([0 3.5e-5])


%% 
function [displacement,force,force_std] = average_data(data_cell,max_disp)
    displacement = linspace(0,max_disp,20)';
    force_array_x = zeros(length(displacement),length(data_cell));
    force_array_y = zeros(size(force_array_x));
    for ii = 1:length(data_cell)
        [~,unique_ind] = unique(data_cell{ii}(:,1));
        force_array_x(:,ii) = interp1(data_cell{ii}(unique_ind,1),data_cell{ii}(unique_ind,2),displacement);
        force_array_y(:,ii) = interp1(data_cell{ii}(unique_ind,1),data_cell{ii}(unique_ind,3),displacement);

    end
    force = [mean(force_array_x,2) mean(force_array_y,2)];
    force_std = [std(force_array_x,[],2) std(force_array_y,[],2)];
end

