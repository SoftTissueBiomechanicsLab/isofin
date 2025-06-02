%% Run batches of simulations
clc; clearvars; close all;

%% Set simulations to run - set parameters to match input files
% 
% % S for straight, U for undulated
% network_type = {'S','U'};
% net_type = 1;
% % Number of fibers in network
% Network_cell = {[2 19 4 6 16 13] [1 24 13] [7 25] [12 14]};
% network_ind = 0;
% for num_fibers = [303 353 403 453]
%     network_ind = network_ind + 1;
%     % Set deformation mode
%     deformation_type = 'PS';
%     % Set network deformation direction
%     bc_rotation = 1;
%     % Networks to run
%     network_numbers = sort(Network_cell{network_ind});%1:25;
%     % Power and order from Input_File_Generation
%     power = 5; order = 4;
%     % Filename to copy input files as - set in ReadInputFile()
%     copyname = 'Input.txt';
%     % Set study name to save logfiles
%     study = ['ps_',num2str(num_fibers)];
%     logpath = ['..\Results\Outputs\voronoi\logfiles\',study,'\'];
%     if ~exist(logpath,'dir'); mkdir(logpath); end
% 
%     batch_check = [logpath,'batch_check'];
% 
% 
%     run_check = cell(length(network_numbers),1);
%     for ii = 1:length(network_numbers)        
% 
%         input_name = ['input_S_',deformation_type,'_',num2str(num_fibers),...
%             '_Case',num2str(network_numbers(ii)),'_P',num2str(power),'_MS',num2str(order),...
%             '_',num2str(bc_rotation),'.txt'];
% 
% 
%         input_path = '..\Input_files\voronoi\';
%         copy_input = [input_path,copyname];
%         if isfile(copy_input); delete(copy_input); end
%         cd_cmd = ['cd ',input_path];
%         copy_cmd = ['copy ', input_name, ' ', copyname];
%         system([cd_cmd, ' && ', copy_cmd]);
% 
%         % Run the simulations
%         logname = ['logfile_',network_type{net_type},'_',deformation_type,'_',num2str(num_fibers),...
%             '_Case',num2str(network_numbers(ii)),'_P',num2str(power),'_MS',num2str(order),...
%             '_',num2str(bc_rotation)];
%         logwogpog = [logpath,logname];
%         [~,logfile] = system('Analysis_3D_Batch');
%         save(logwogpog,"logfile")
%         % 1 means the run failed, 0 means success
%         run_check{ii} = contains(logfile,'STOPPING');
% 
%     end
% 
%     save(batch_check,'run_check')
% end

%% Set simulations to run - set parameters to match input files

% S for straight, U for undulated
network_type = {'S','U'};
net_type = 1;
    % Number of fibers in network
for num_fibers =  [313]
    % for new_Case = [9 10]
       
        % Set deformation mode
        deformation_type = 'BIAX';
        % Set network deformation direction
        bc_rotation = 1;
        % Networks to run
        network_numbers = [2 3 4 11 14];
        % Power and order from Input_File_Generation
        power = 5; order = 4;
        % Filename to copy input files as - set in ReadInputFile()
        copyname = 'Input.txt';
        % Set study name to save logfiles
        study = ['biax_',network_type{net_type},'_',num2str(num_fibers)];
        logpath = ['..\Results\Outputs\voronoi\logfiles\',study,'\'];
        if ~exist(logpath,'dir'); mkdir(logpath); end
        
        batch_check = [logpath,'batch_check'];
        
        
        run_check = cell(length(network_numbers),1);
        for ii = 1:length(network_numbers)
            % if net_type == 1
                input_name = ['input_',network_type{net_type},'_',deformation_type,'_',num2str(num_fibers),...
                    '_Case',num2str(network_numbers(ii)),'_P',num2str(power),'_MS',num2str(order),...
                    '_',num2str(bc_rotation),'.txt'];
            % else
                % input_name = ['input_',network_type{net_type},'_',deformation_type,'_',num2str(num_fibers),...
                %     '_Case',num2str(network_numbers(ii)),'_',num2str(new_Case),'_P',num2str(power),'_MS',num2str(order),...
                %     '_',num2str(bc_rotation),'.txt'];
            % end
        
            input_path = '..\Input_files\voronoi\';
            copy_input = [input_path,copyname];
            if isfile(copy_input); delete(copy_input); end
            cd_cmd = ['cd ',input_path];
            copy_cmd = ['copy ', input_name, ' ', copyname];
            [~,copy_result] = system([cd_cmd, ' && ', copy_cmd]);
             % system([cd_cmd, ' && ', copy_cmd]);

            % Run the simulations
            logname = ['logfile_',network_type{net_type},'_',deformation_type,'_',num2str(num_fibers),...
                '_Case',num2str(network_numbers(ii)),'_P',num2str(power),'_MS',num2str(order),...
                '_',num2str(bc_rotation)];
            logwogpog = [logpath,logname];
            fprintf('Running %s_%d Network %d: ',network_type{net_type},num_fibers,network_numbers(ii))
            [~,logfile] = system('Analysis_3D_Batch');
            save(logwogpog,"logfile")
            % 1 means the run failed, 0 means success
            run_check{ii} = contains(logfile,'STOPPING');
            if run_check{ii}
                fprintf('Failed\n')
            else
                fprintf('Converged... hopefully\n')
            end

        end

        save(batch_check,'run_check')
    % end
end


%% Set simulations to run - set parameters to match input files

% S for straight, U for undulated
network_type = {'S','U'};
net_type = 2;
    % Number of fibers in network
for num_fibers =  [303]
    for new_Case = [9 10]
       
        % Set deformation mode
        deformation_type = 'BIAX';
        % Set network deformation direction
        bc_rotation = 1;
        % Networks to run
        network_numbers = [2 3 4 11 14];
        % Power and order from Input_File_Generation
        power = 5; order = 4;
        % Filename to copy input files as - set in ReadInputFile()
        copyname = 'Input.txt';
        % Set study name to save logfiles
        study = ['biax_',network_type{net_type},'_',num2str(num_fibers),'_',num2str(new_Case)];
        logpath = ['..\Results\Outputs\voronoi\logfiles\',study,'\'];
        if ~exist(logpath,'dir'); mkdir(logpath); end
        
        batch_check = [logpath,'batch_check'];
        
        
        run_check = cell(length(network_numbers),1);
        for ii = 1:length(network_numbers)
            % if net_type == 1
            %     input_name = ['input_',network_type{net_type},'_',deformation_type,'_',num2str(num_fibers),...
            %         '_Case',num2str(network_numbers(ii)),'_P',num2str(power),'_MS',num2str(order),...
            %         '_',num2str(bc_rotation),'.txt'];
            % else
                input_name = ['input_',network_type{net_type},'_',deformation_type,'_',num2str(num_fibers),...
                    '_Case',num2str(network_numbers(ii)),'_',num2str(new_Case),'_P',num2str(power),'_MS',num2str(order),...
                    '_',num2str(bc_rotation),'.txt'];
            % end
        
            input_path = '..\Input_files\voronoi\';
            copy_input = [input_path,copyname];
            if isfile(copy_input); delete(copy_input); end
            cd_cmd = ['cd ',input_path];
            copy_cmd = ['copy ', input_name, ' ', copyname];
            [~,copy_result] = system([cd_cmd, ' && ', copy_cmd]);
        
            % Run the simulations
            logname = ['logfile_',network_type{net_type},'_',deformation_type,'_',num2str(num_fibers),...
                '_Case',num2str(network_numbers(ii)),'_P',num2str(power),'_MS',num2str(order),...
                '_',num2str(bc_rotation)];
            logwogpog = [logpath,logname];
            fprintf('Running %s_%d_%d Network %d: ',network_type{net_type},num_fibers,new_Case,network_numbers(ii))
            [~,logfile] = system('Analysis_3D_Batch');
            save(logwogpog,"logfile")
            % 1 means the run failed, 0 means success
            run_check{ii} = contains(logfile,'STOPPING');
            if run_check{ii}
                fprintf('Failed\n')
            else
                fprintf('Converged... hopefully\n')
            end

        end

        save(batch_check,'run_check')
    end
end


%% Set simulations to run - set parameters to match input files

% S for straight, U for undulated
network_type = {'S','U'};
net_type = 1;
    % Number of fibers in network
for num_fibers =  [323 333 343]
    % for new_Case = [9 10]
       
        % Set deformation mode
        deformation_type = 'BIAX';
        % Set network deformation direction
        bc_rotation = 1;
        % Networks to run
        network_numbers = [2 3 4 11 14];
        % Power and order from Input_File_Generation
        power = 5; order = 4;
        % Filename to copy input files as - set in ReadInputFile()
        copyname = 'Input.txt';
        % Set study name to save logfiles
        study = ['biax_',network_type{net_type},'_',num2str(num_fibers)];
        logpath = ['..\Results\Outputs\voronoi\logfiles\',study,'\'];
        if ~exist(logpath,'dir'); mkdir(logpath); end
        
        batch_check = [logpath,'batch_check'];
        
        
        run_check = cell(length(network_numbers),1);
        for ii = 1:length(network_numbers)
            % if net_type == 1
                input_name = ['input_',network_type{net_type},'_',deformation_type,'_',num2str(num_fibers),...
                    '_Case',num2str(network_numbers(ii)),'_P',num2str(power),'_MS',num2str(order),...
                    '_',num2str(bc_rotation),'.txt'];
            % else
                % input_name = ['input_',network_type{net_type},'_',deformation_type,'_',num2str(num_fibers),...
                %     '_Case',num2str(network_numbers(ii)),'_',num2str(new_Case),'_P',num2str(power),'_MS',num2str(order),...
                %     '_',num2str(bc_rotation),'.txt'];
            % end
        
            input_path = '..\Input_files\voronoi\';
            copy_input = [input_path,copyname];
            if isfile(copy_input); delete(copy_input); end
            cd_cmd = ['cd ',input_path];
            copy_cmd = ['copy ', input_name, ' ', copyname];
            [~,copy_result] = system([cd_cmd, ' && ', copy_cmd]);
             % system([cd_cmd, ' && ', copy_cmd]);

            % Run the simulations
            logname = ['logfile_',network_type{net_type},'_',deformation_type,'_',num2str(num_fibers),...
                '_Case',num2str(network_numbers(ii)),'_P',num2str(power),'_MS',num2str(order),...
                '_',num2str(bc_rotation)];
            logwogpog = [logpath,logname];
            fprintf('Running %s_%d Network %d: ',network_type{net_type},num_fibers,network_numbers(ii))
            [~,logfile] = system('Analysis_3D_Batch');
            save(logwogpog,"logfile")
            % 1 means the run failed, 0 means success
            run_check{ii} = contains(logfile,'STOPPING');
            if run_check{ii}
                fprintf('Failed\n')
            else
                fprintf('Converged... hopefully\n')
            end

        end

        save(batch_check,'run_check')
    % end
end

%% Set simulations to run - set parameters to match input files
% 
% % S for straight, U for undulated
% network_type = {'S','U'};
% 
%     % Number of fibers in network
% for num_fibers =  [353 403 453]
%     for net_type = 1:2
%         if net_type == 1 && num_fibers <= 353
%             fprintf('Skipping %s_%d\n', network_type{net_type},num_fibers)
%             continue
%         end
%         % Set deformation mode
%         deformation_type = 'BIAX';
%         % Set network deformation direction
%         bc_rotation = 1;
%         % Networks to run
%         network_numbers = 1:25;
%         % Power and order from Input_File_Generation
%         power = 5; order = 4;
%         % Filename to copy input files as - set in ReadInputFile()
%         copyname = 'Input.txt';
%         % Set study name to save logfiles
%         study = ['biax_',network_type{net_type},'_',num2str(num_fibers)];
%         logpath = ['..\Results\Outputs\voronoi\logfiles\',study,'\'];
%         if ~exist(logpath,'dir'); mkdir(logpath); end
% 
%         batch_check = [logpath,'batch_check'];
% 
% 
%         run_check = cell(length(network_numbers),1);
%         for ii = 1:length(network_numbers)
%             if net_type == 1
%                 input_name = ['input_',network_type{net_type},'_',deformation_type,'_',num2str(num_fibers),...
%                     '_Case',num2str(network_numbers(ii)),'_P',num2str(power),'_MS',num2str(order),...
%                     '_',num2str(bc_rotation),'.txt'];
%             else
%                 input_name = ['input_',network_type{net_type},'_',deformation_type,'_',num2str(num_fibers),...
%                     '_Case',num2str(network_numbers(ii)),'_1_P',num2str(power),'_MS',num2str(order),...
%                     '_',num2str(bc_rotation),'.txt'];
%             end
% 
%             input_path = '..\Input_files\voronoi\';
%             copy_input = [input_path,copyname];
%             if isfile(copy_input); delete(copy_input); end
%             cd_cmd = ['cd ',input_path];
%             copy_cmd = ['copy ', input_name, ' ', copyname];
%             [~,copy_result] = system([cd_cmd, ' && ', copy_cmd]);
% 
%             % Run the simulations
%             logname = ['logfile_',network_type{net_type},'_',deformation_type,'_',num2str(num_fibers),...
%                 '_Case',num2str(network_numbers(ii)),'_P',num2str(power),'_MS',num2str(order),...
%                 '_',num2str(bc_rotation)];
%             logwogpog = [logpath,logname];
%             fprintf('Running %s_%d Network %d: ',network_type{net_type},num_fibers,ii)
%             [~,logfile] = system('Analysis_3D_Batch');
%             save(logwogpog,"logfile")
%             % 1 means the run failed, 0 means success
%             run_check{ii} = contains(logfile,'STOPPING');
%             if run_check{ii}
%                 fprintf('Failed\n')
%             else
%                 fprintf('Converged... hopefully\n')
%             end
% 
%         end
% 
%         save(batch_check,'run_check')
%     end
% end

%% Set simulations to run - set parameters to match input files
% 
% % S for straight, U for undulated
% network_type = {'S','U'};
% net_type = 2;
% % Number of fibers in network
% num_fibers = 253;%[145 185 205 225 245 265 305]
% % Set deformation mode
% deformation_type = 'PS';
% % Set network deformation direction
% bc_rotation = 1;
% % Networks to run
% network_numbers = 1:10;
% % Power and order from Input_File_Generation
% power = 5; order = 4;
% % Filename to copy input files as - set in ReadInputFile()
% copyname = 'Input.txt';
% % Set study name to save logfiles
% study = 'und_rad_ps';
% logpath = ['..\Results\Outputs\voronoi\logfiles\',study,'\'];
% if ~exist(logpath,'dir'); mkdir(logpath); end
% 
% batch_check = [logpath,'batch_check'];
% 
% 
% run_check = cell(length(network_numbers),1);
% for ii = 1:length(network_numbers)
%     input_name = ['input_',network_type{net_type},'_',deformation_type,'_',num2str(num_fibers),...
%         '_Case',num2str(network_numbers(ii)),'_1_P',num2str(power),'_MS',num2str(order),...
%         '_',num2str(bc_rotation),'.txt'];
% 
%     input_path = '..\Input_files\voronoi\';
%     copy_input = [input_path,copyname];
%     if isfile(copy_input); delete(copy_input); end
%     cd_cmd = ['cd ',input_path];
%     copy_cmd = ['copy ', input_name, ' ', copyname];
%     system([cd_cmd, ' && ', copy_cmd]);
% 
%     % Run the simulations
%     logname = ['logfile_',network_type{net_type},'_',deformation_type,'_',num2str(num_fibers),...
%         '_Case',num2str(network_numbers(ii)),'_P',num2str(power),'_MS',num2str(order),...
%         '_',num2str(bc_rotation)];
%     logwogpog = [logpath,logname];
%     [~,logfile] = system('Analysis_3D_Batch');
%     save(logwogpog,"logfile")
%     % 1 means the run failed, 0 means success
%     run_check{ii} = contains(logfile,'STOPPING');
% 
% end
% 
% save(batch_check,'run_check')



%% Set simulations to run - set parameters to match input files
% 
% % S for straight, U for undulated
% network_type = {'S','U'};
% net_type = 2;
% % Number of fibers in network
% num_fibers = 252;%[145 185 205 225 245 265 305]
% % Set deformation mode
% deformation_type = 'PS';
% % Set network deformation direction
% for bc_rotation = 1:2
% % Networks to run
% network_numbers = 1:25;
% % Power and order from Input_File_Generation
% power = 5; order = 4;
% % Filename to copy input files as - set in ReadInputFile()
% copyname = 'Input.txt';
% % Set study name to save logfiles
% study = ['ps_',num2str(bc_rotation)];
% logpath = ['..\Results\Outputs\voronoi\logfiles\',study,'\'];
% if ~exist(logpath,'dir'); mkdir(logpath); end
% 
% batch_check = [logpath,'batch_check'];
% 
% 
% run_check = cell(length(network_numbers),1);
% for ii = 1:length(network_numbers)
%     input_name = ['input_',network_type{net_type},'_',deformation_type,'_',num2str(num_fibers),...
%         '_Case',num2str(network_numbers(ii)),'_1_P',num2str(power),'_MS',num2str(order),...
%         '_',num2str(bc_rotation),'.txt'];
% 
%     input_path = '..\Input_files\voronoi\';
%     copy_input = [input_path,copyname];
%     if isfile(copy_input); delete(copy_input); end
%     cd_cmd = ['cd ',input_path];
%     copy_cmd = ['copy ', input_name, ' ', copyname];
%     system([cd_cmd, ' && ', copy_cmd]);
% 
%     % Run the simulations
%     logname = ['logfile_',network_type{net_type},'_',deformation_type,'_',num2str(num_fibers),...
%         '_Case',num2str(network_numbers(ii)),'_P',num2str(power),'_MS',num2str(order),...
%         '_',num2str(bc_rotation)];
%     logwogpog = [logpath,logname];
%     [~,logfile] = system('Analysis_3D_Batch');
%     save(logwogpog,"logfile")
%     % 1 means the run failed, 0 means success
%     run_check{ii} = contains(logfile,'STOPPING');
% 
% end
% 
% save(batch_check,'run_check')
% end
% 
% 
% 
% 
%% Set simulations to run - set parameters to match input files
% 
% % S for straight, U for undulated
% network_type = {'S','U'};
% for net_type = 1:2
% % Number of fibers in network
% for num_fibers = [145 185 205 225 245 265 305]
% % Set deformation mode
% deformation_type = 'BIAX';
% % Set network deformation direction
% bc_rotation = 1;
% % Networks to run
% network_numbers = 1:25;
% % Power and order from Input_File_Generation
% power = 5; order = 4;
% % Filename to copy input files as - set in ReadInputFile()
% copyname = 'Input.txt';
% % Set study name to save logfiles
% study = ['biax_',network_type{net_type},'_',num2str(num_fibers)];
% logpath = ['..\Results\Outputs\voronoi\logfiles\',study,'\'];
% if ~exist(logpath,'dir'); mkdir(logpath); end
% 
% batch_check = [logpath,'batch_check'];
% 
% 
% run_check = cell(length(network_numbers),1);
% for ii = 1:length(network_numbers)
%     if net_type == 1
%         input_name = ['input_',network_type{net_type},'_',deformation_type,'_',num2str(num_fibers),...
%             '_Case',num2str(network_numbers(ii)),'_P',num2str(power),'_MS',num2str(order),...
%             '_',num2str(bc_rotation),'.txt'];
%     else
%         input_name = ['input_',network_type{net_type},'_',deformation_type,'_',num2str(num_fibers),...
%             '_Case',num2str(network_numbers(ii)),'_1_P',num2str(power),'_MS',num2str(order),...
%             '_',num2str(bc_rotation),'.txt'];
%     end
% 
%     input_path = '..\Input_files\voronoi\';
%     copy_input = [input_path,copyname];
%     if isfile(copy_input); delete(copy_input); end
%     cd_cmd = ['cd ',input_path];
%     copy_cmd = ['copy ', input_name, ' ', copyname];
%     system([cd_cmd, ' && ', copy_cmd]);
% 
%     % Run the simulations
%     logname = ['logfile_',network_type{net_type},'_',deformation_type,'_',num2str(num_fibers),...
%         '_Case',num2str(network_numbers(ii)),'_P',num2str(power),'_MS',num2str(order),...
%         '_',num2str(bc_rotation)];
%     logwogpog = [logpath,logname];
%     [~,logfile] = system('Analysis_3D_Batch');
%     save(logwogpog,"logfile")
%     % 1 means the run failed, 0 means success
%     run_check{ii} = contains(logfile,'STOPPING');
% 
% end
% 
% save(batch_check,'run_check')
% end
% end

%%
% %% Set simulations to run - set parameters to match input files
% 
% % S for straight, U for undulated
% network_type = 'S';
% % Number of fibers in network
% num_fibers = 220;
% % Set deformation mode
% deformation_type = 'PS';
% % Set network deformation direction
% bc_rotation = 1;
% % Networks to run
% network_numbers = 1:20;
% % Power and order from Input_File_Generation
% power = 5; order = 4;
% % Filename to copy input files as - set in ReadInputFile()
% copyname = 'Input.txt';
% % Set study name to save logfiles
% study = 'ps_odd';
% logpath = ['..\Results\Outputs\voronoi\logfiles\',study,'\'];
% if ~exist(logpath,'dir'); mkdir(logpath); end
% 
% batch_check = '..\Results\Outputs\voronoi\logfiles\ps_odd';
% 
% 
% run_check = cell(length(network_numbers),1);
% for ii = 1:length(network_numbers)
%     input_name = ['input_',network_type,'_',deformation_type,'_',num2str(num_fibers),...
%         '_Case',num2str(network_numbers(ii)),'_P',num2str(power),'_MS',num2str(order),...
%         '_',num2str(bc_rotation),'.txt'];
% 
%     input_path = '..\Input_files\voronoi\';
%     copy_input = [input_path,copyname];
%     if isfile(copy_input); delete(copy_input); end
%     cd_cmd = ['cd ',input_path];
%     copy_cmd = ['copy ', input_name, ' ', copyname];
%     system([cd_cmd, ' && ', copy_cmd]);
% 
%     % Run the simulations
%     logname = ['logfile_',network_type,'_',deformation_type,'_',num2str(num_fibers),...
%         '_Case',num2str(network_numbers(ii)),'_P',num2str(power),'_MS',num2str(order),...
%         '_',num2str(bc_rotation)];
%     logwogpog = [logpath,logname];
%     [~,logfile] = system('Analysis_3D_Batch_PS');
%     save(logwogpog,"logfile")
%     % 1 means the run failed, 0 means success
%     run_check{ii} = contains(logfile,'STOPPING');
% 
% end
% 
% save(batch_check,'run_check')


%% Set simulations to run - set parameters to match input files

% % S for straight, U for undulated
% network_type = 'S';
% % Number of fibers in network
% num_fibers = 600;
% % Set deformation mode
% deformation_type = 'UAX';
% % Set network deformation direction
% bc_rotation = 5;
% % Networks to run
% network_numbers = 1:20;
% % Power and order from Input_File_Generation
% power = 5; order = 4;
% % Filename to copy input files as - set in ReadInputFile()
% copyname = 'Input.txt';
% % Set name to save batch runs
% batch_check = '..\Results\Outputs\voronoi\logfiles\pure_shear';
% 
% 
% run_check = cell(length(network_numbers),1);
% for ii = 1:length(network_numbers)
%     input_name = ['input_',network_type,'_',deformation_type,'_',num2str(num_fibers),...
%         '_Case',num2str(network_numbers(ii)),'_P',num2str(power),'_MS',num2str(order),...
%         '_',num2str(bc_rotation),'.txt'];
% 
%     input_path = '..\Input_files\voronoi\';
%     copy_input = [input_path,copyname];
%     if isfile(copy_input); delete(copy_input); end
%     cd_cmd = ['cd ',input_path];
%     copy_cmd = ['copy ', input_name, ' ', copyname];
%     system([cd_cmd, ' && ', copy_cmd]);
% 
%     % Run the simulations
%     logname = ['input_',network_type,'_',deformation_type,'_',num2str(num_fibers),...
%         '_Case',num2str(network_numbers(ii)),'_P',num2str(power),'_MS',num2str(order),...
%         '_',num2str(bc_rotation)];
%     logpath = ['..\Results\Outputs\voronoi\logfiles\',logname];
% 
%     [~,logfile] = system('Analysis_3D_Batch_Pure_Shear');
%     save(logpath,"logfile")
% 
% end

% save(batch_check,'run_check')