% clear
% clc
% % This is the later version of EnergyCalculationMLAElongation that includes
% % the maximum length of the PF from the full range trials April 2017
% 
% direc = '/home/lauren/Desktop/MotionDataAus_April2017/';
% 
% list_files = dir([direc '*processed*']);
% 
% num_files = length(list_files);
% 
% % --------- To reprocess all of the files -----------------------
% ind1 = 1;
% ind2 = 1;
% % root_files = cell(num_files,1);
% for i = 1:num_files
%     subj_num = str2double(list_files(i).name(2:3));
%     subj_str = sprintf('S%.2i',subj_num);
%     rawdatadir = ['/media/lauren/Elements/AustraliaCollection/' subj_str '/SelectedTrials/'];
%     
%     
%     raw_data_files{i} = [rawdatadir list_files(i).name(1:end-20) '.mat'];
%     
%     if isempty(regexp(list_files(i).name(1:end-20),'fullrange')) % normal file
%         root_files{ind1} = [list_files(i).name(1:end-20) '.mat'];
%         ind1 = ind1 + 1;
%     else % full range file
%         root_filesFR{ind2} = [list_files(i).name(1:end-20) '.mat'];
%         ind2 = ind2 + 1;
%     end
%     
%     
%     try % basically I'm too lazy to move all the files over so if it isn't in this direcory, it's in the other one
%         load(raw_data_files{i})
%     catch
%         raw_data_files{i} = ['/media/lauren/Elements/AustraliaCollection/' subj_str '/Mocap/' list_files(i).name(1:end-20) '.mat'];
%     end
% end
% root_files = root_files';
% root_filesFR = root_filesFR';
% % MocapDataTransform(root_files,direc)
% 
% % MocapDataTransform(raw_data_files,direc)
% 
% clearvars -except root_files num_files direc root_filesFR
% % -----------------------------------------------------------------
% %%
% % load the offsets
% [subject_info,subj_head] = xlsread([direc 'Subject_list.xlsx']);
% close all;
% save_impulse_flag = 0;
% data_structFR = struct('subjects','','weight',[]);
% 
% for i = 1:9
%     data_structFR(i).weight = subject_info(i,1);
% end
% 
% 
% nFRfiles = length(root_filesFR);
% 
% 
% count = 1; % index for the structure so there is only one row per subject
% % the subjects will be replaced if a different trial for the same subject
% % has a greater maximum pf length
% 
% 
% % determine all the parameters for the full range trial for scaling
% % purposes
% for i = 1:nFRfiles
%     load([direc root_filesFR{i}(1:end-4) '_processedmotion.mat']);
%     if ~isfield(marker_data,'ME_') || ~isfield(marker_data,'LE_')  || ~isfield(marker_data,'MM_') || ~isfield(marker_data,'LM_')
%         marker_data = addMEandLEtoMarkerData(direc,[root_filesFR{i}(1:end-4) '_processedmotion.mat']);
%     end
%     model = createFootModel(marker_data);
%     
%     subj_num = str2double(root_filesFR{i}(2:3));
%     subj_num_list(i) = subj_num;
%     
%     k_ = strfind( root_filesFR{i},'_');
%     k_ = [0 k_];
%     % list the subjects
%     subj_name = root_filesFR{i}(k_(1)+1:k_(2)-1);
%     
%     % check to see if it's already in the structure; if it is, then
%     % calculate the max pf length and see if it's bigger. If it is, replace
%     % the trial with the info from that one
%     
%     subj_match = findInStruct(data_structFR,'subjects',subj_name);
%     disp(num2str(max(model.elongation)))
%     if ~isempty(subj_match) % there is a match, so subject is already in the structure
%         temp_max = max(model.elongation);
%         if temp_max > data_structFR(subj_match).max_pf
%             subj_ind = subj_match;
%             disp('Replaced by:')
%             data_structFR(subj_ind).subjects = subj_name;
%             data_structFR(subj_ind).MLA      = model.MLA;%-model.MLA(impulse_ends(1));
%             data_structFR(subj_ind).length_pf= model.elongation;
%             data_structFR(subj_ind).max_pf   = max(data_structFR(subj_ind).length_pf);
%             data_structFR(subj_ind).arch_height{1}   = model.marker_data.CST;
%             
%         end
%     else % place it in a new row, given by the index count
%         subj_ind = count;
%         data_structFR(subj_ind).subjects = subj_name;
%         data_structFR(subj_ind).MLA      = model.MLA;%-model.MLA(impulse_ends(1));
%         data_structFR(subj_ind).length_pf= model.elongation;
%         data_structFR(subj_ind).max_pf   = max(data_structFR(subj_ind).length_pf);
%         data_structFR(subj_ind).arch_height{1}   = model.marker_data.CST;
%         count = count + 1;
%         
%     end
%     
%     
%     
%     
%     disp(['The trial ' root_filesFR{i}(1:end-4) ' is added to the data structure'])
%     
%     
% end
% %% Make sure that the longest pf length is taken in all the trials
% 
% 
% for i =  1:length(root_files)
%     
%     load([direc root_files{i}(1:end-4) '_processedmotion.mat']);
%     if ~isfield(marker_data,'ME_') || ~isfield(marker_data,'LE_') || ~isfield(marker_data,'MM_') || ~isfield(marker_data,'LM_')
%         marker_data = addMEandLEtoMarkerData(direc,[root_files{i}(1:end-4) '_processedmotion.mat']);
%     end
%     model = createFootModel(marker_data);
%     
%     subj_num = str2double(root_files{i}(2:3));
%     subj_num_list(i) = subj_num;
%     
%     k_ = strfind( root_files{i},'_');
%     k_ = [0 k_];
%     % list the subjects
%     subj_name = root_files{i}(k_(1)+1:k_(2)-1);
%     
%     % check to see if it's already in the structure; if it is, then
%     % calculate the max pf length and see if it's bigger. If it is, replace
%     % the trial with the info from that one
%     
%     subj_match = findInStruct(data_structFR,'subjects',subj_name);
%     
%     if ~isempty(subj_match) % there is a match, so subject is already in the structure
%         temp_max = max(model.elongation);
%         if temp_max > data_structFR(subj_match).max_pf
%             subj_ind = subj_match;
%             
%             data_structFR(subj_ind).subjects = subj_name;
%             data_structFR(subj_ind).MLA      = model.MLA;%-model.MLA(impulse_ends(1));
%             data_structFR(subj_ind).length_pf= model.elongation;
%             data_structFR(subj_ind).max_pf   = max(model.elongation);
%             data_structFR(subj_ind).arch_height{1}   = model.marker_data.CST;
%             disp(num2str(max(model.elongation)))
%             
%             disp(['The trial ' root_files{i}(1:end-4) ' is replacing part of the data structure'])
%         end
%         
%         
%     end
%     
% end
% 





%% for all of the other trials, use the force to determine the frames where the impulse occurs
% then compute all parameters needed for analysis
clearvars -except root_files data_structFR direc subject_info save_impulse_flag
i = 1;
for j =  1:length(root_files)
    if isempty(regexp(root_files{j},'statcal'))
        load([direc root_files{j}(1:end-4) '_processedmotion.mat']);
        if ~isfield(marker_data,'ME_') || ~isfield(marker_data,'LE_')  || ~isfield(marker_data,'MM_') || ~isfield(marker_data,'LM_')
            marker_data = addMEandLEtoMarkerData(direc,[root_filesFR{j}(1:end-4) '_processedmotion.mat']);
        end
        subj_num = str2double(root_files{j}(2:3));
        
        raw_force = force_data.Force(1:3,:);
        raw_force = raw_force - repmat(subject_info(subj_num - 2,4:6)',1,length(raw_force));
        
        force_fr = force_data.Frequency;
        force_nfr = force_data.NrOfSamples;
        
        mot_fr = marker_data.FrameRate;
        mot_nfr = marker_data.nFrames;
        raw_force_norm = raw_force(3,:);
        %     for j = 1:force_nfr
        %         raw_force_norm(j) = norm(raw_force(:,j));
        %     end
        % the force is typically a multiple of the motion capture
        sample_factor = force_fr/mot_fr;
        
        if round(sample_factor) ~= sample_factor
            error('Error : The force and motion frame rates are not a constant sampling factor apart. ')
        end
        
        % resampled force
        force_res_norm = raw_force_norm(1:sample_factor:end);
        force_res = raw_force(1:3,1:sample_factor:end);
        moment_res = force_data.Moment(1:3,1:sample_factor:end);
%         Determine the centre of pressures
        COP_temp = force_data.COP(:,1:sample_factor:end); 
%         COPy = force_data.COP(2,1:sample_factor:end); 
%         COPz = force_data.COP(3,1:sample_factor:end);
        % determine when the impulse occurs by the shape of the force curve
        if j == 14 % issue with the fourteenth trial
            force_res_norm(1:150) = 0;
        end
        % only do this if it's selected by the flag. Otherwise load the file
        if save_impulse_flag == 1
            aveforce(1) = mean(force_res_norm(1:5));
            for t = 3:length(force_res_norm)-2
                aveforce(t-1) = mean(force_res_norm(t-2:t+2));% compute the average
            end
            
            ave_slope = 0.4;%abs(mean(diff(aveforce(1:5))))*5;
            
            force_filt = LowPassButterworth(aveforce,1,40,force_fr);
            
            diff_force_raw = abs(diff(force_filt)); % compute the difference between data points (slope approximation)
            [val,pk_locs] = findpeaks(diff_force_raw); % find the peaks of the curve (indicate roughly where the impulses are)
            while_ind = 1;
            pk_locs_new = [];
            while (length(pk_locs_new) ~= 4 && while_ind < 100) % find four peaks
                pk_locs_new = pk_locs(find(val > 0.25 + .25 * while_ind)); % only take the peak locations with a magnitude bigger than this value
                while_ind = while_ind + 1;
            end
            if while_ind == 100
                warning(['Iteration limit reached. There were not four peaks correctly identified. Skipping trial ' root_files{j} '.'])
                continue
            end
            
            for k = 1:4 %(for each peak)
                diff_search = 1; inc = 0; %increment
                if rem(k,2) == 1 % odd peak, look backwards
                    while (diff_search > ave_slope) || (pk_locs_new(k)-inc == 0)
                        diff_search = abs(diff(aveforce([pk_locs_new(k)-inc,pk_locs_new(k)-1-inc])));
                        inc = inc + 1;
                    end
                    
                    if pk_locs_new(k)-inc-1 == 0
                        impulse_ends(k) = 1;
                    else
                        impulse_ends(k) = pk_locs_new(k)-inc-1;
                    end
                else % even peak (look forwards)
                    while (diff_search > ave_slope) || (pk_locs_new(k)+inc == length(diff_force_raw))
                        diff_search = abs(diff(aveforce([pk_locs_new(k)+inc,pk_locs_new(k)+1+inc])));
                        inc = inc + 1;
                    end
                    
                    if  (pk_locs_new(k)+inc+1 == length(diff_force_raw))
                        impulse_ends(k) = length(diff_force_raw);
                    else
                        impulse_ends(k) = pk_locs_new(k)+inc+1;
                    end
                end
            end
            
            % % --------- UNCOMMENT FOR FIGURES OF IMPULSE QUALITY---------------------
            figure(); plot(aveforce); hold on; plot(impulse_ends,aveforce(impulse_ends),'or');
            title(root_files{j}(1:end-4))
            drawnow
            pause(1)
            % % ----------------------------------------------------------------------
            
            % save the impulse files
            
            save([direc root_files{j}(1:end-4) '_impulse.mat'],'impulse_ends')
            disp(['The trial ' root_files{j}(1:end-4) '_impulse.mat is saved.'])
        else
            %     % load the impulse_ends variable
            load([direc root_files{j}(1:end-4) '_impulse.mat']);
        end
        
        % -------Now load the force and the elongation data and compute ---------
        % ----------the energy change over the cycle---------------------------
        
        %
        %     % load the motion and force data
        %     load([direc root_files{i}(1:end-4) '_processedmotion.mat']);
        
        k_ = strfind( root_files{j},'_');
        k_ = [0 k_];
        data_struct(i).subjects = root_files{j}(k_(1)+1:k_(2)-1);
        data_struct(i).trials = root_files{j}(k_(2)+1:k_(3)-1);
        data_struct(i).conditions = root_files{j}(k_(3)+1:k_(4)-1);
        
        subj_FR_ind = findInStruct(data_structFR,'subjects',root_files{j}(k_(1)+1:k_(2)-1));
        
        max_pf = data_structFR(subj_FR_ind).max_pf;
        
        n1 = length(impulse_ends(1):impulse_ends(2));
        n2 = length(impulse_ends(3):impulse_ends(4));
        
        model = createFootModel(marker_data);
        
        force{1}         = force_res(1:3,impulse_ends(1):impulse_ends(2));
        force{2}         = force_res(1:3,impulse_ends(3):impulse_ends(4));
        data_struct(i).force{1} = force{1}/(data_structFR(subj_FR_ind).weight*9.81);
        data_struct(i).force{2} = force{2}/(data_structFR(subj_FR_ind).weight*9.81);
        
        data_struct(i).force_val{1}     = -force_res(3,impulse_ends(1):impulse_ends(2))/data_structFR(subj_FR_ind).weight;
        data_struct(i).force_val{2}     = -force_res(3,impulse_ends(3):impulse_ends(4))/data_structFR(subj_FR_ind).weight;
        
        
        COP{1}     = COP_temp(:,impulse_ends(1):impulse_ends(2));
        COP{2}     = COP_temp(:,impulse_ends(3):impulse_ends(4));
        
%         COPy{1}     = COPy(impulse_ends(1):impulse_ends(2));
%         COPy{2}     = COPy(impulse_ends(3):impulse_ends(4));
%         
%         COPz{1}     = COPz(impulse_ends(1):impulse_ends(2));
%         COPz{2}     = COPz(impulse_ends(3):impulse_ends(4));
        
        
        data_struct(i).MLA{1}           = model.MLA(impulse_ends(1):impulse_ends(2))-model.MLA(impulse_ends(1));
        data_struct(i).MLA{2}           = model.MLA(impulse_ends(3):impulse_ends(4))-model.MLA(impulse_ends(3));
        
        data_struct(i).length_pf{1}    = model.elongation(impulse_ends(1):impulse_ends(2));
        data_struct(i).length_pf{2}    = model.elongation(impulse_ends(3):impulse_ends(4));
        
        data_struct(i).elongation{1}    = model.elongation(impulse_ends(1):impulse_ends(2))/max_pf;% -  repmat(model.elongation(impulse_ends(1)),1,n1);
        data_struct(i).elongation{2}    = model.elongation(impulse_ends(3):impulse_ends(4))/max_pf;% -  repmat(model.elongation(impulse_ends(3)),1,n2);
        
        lengthPF{1}    = model.elongation(impulse_ends(1):impulse_ends(2));
        lengthPF{2}    = model.elongation(impulse_ends(3):impulse_ends(4));
        
        arch_marker{1}   = model.marker_data.CST(1:3,impulse_ends(1):impulse_ends(2));% -  model.marker_data.CST(1:3,impulse_ends(1));
        arch_marker{2}   = model.marker_data.CST(1:3,impulse_ends(3):impulse_ends(4));% -  model.marker_data.CST(1:3,impulse_ends(3));
        
        ca_cst_proj{1} = model.sagittal_arch.ca_cst_proj(:,impulse_ends(1):impulse_ends(2));
        ca_cst_proj{2} = model.sagittal_arch.ca_cst_proj(:,impulse_ends(3):impulse_ends(4));
        
        mh1_cst_proj{1} = model.sagittal_arch.mh1_cst_proj(:,impulse_ends(1):impulse_ends(2));
        mh1_cst_proj{2} = model.sagittal_arch.mh1_cst_proj(:,impulse_ends(3):impulse_ends(4));
          
        mh1_ca_proj{1} = model.sagittal_arch.mh1_ca_proj(:,impulse_ends(1):impulse_ends(2));
        mh1_ca_proj{2} = model.sagittal_arch.mh1_ca_proj(:,impulse_ends(3):impulse_ends(4));
        
        sagittal_plane{1} = model.sagittal_arch.sag_plane_foot(:,impulse_ends(1):impulse_ends(2));
        sagittal_plane{2} = model.sagittal_arch.sag_plane_foot(:,impulse_ends(3):impulse_ends(4));
        
        foot_pose{1} = model.pose.foot(:,:,impulse_ends(1):impulse_ends(2));
        foot_pose{2} = model.pose.foot(:,:,impulse_ends(3):impulse_ends(4));
        
        cal_pose{1} = model.pose.cal(:,:,impulse_ends(1):impulse_ends(2));
        cal_pose{2} = model.pose.cal(:,:,impulse_ends(3):impulse_ends(4));
        
        cal_met_pose{1} = model.pose.cal_met(:,:,impulse_ends(1):impulse_ends(2));
        cal_met_pose{2} = model.pose.cal_met(:,:,impulse_ends(3):impulse_ends(4));
        
        met_pose{1} = model.pose.met(:,:,impulse_ends(1):impulse_ends(2));
        met_pose{2} = model.pose.met(:,:,impulse_ends(3):impulse_ends(4));
        
        data_struct(i).marker_data.ca_ma{1} = model.marker_data.CA_(:,impulse_ends(1):impulse_ends(2));
        data_struct(i).marker_data.ca_ma{2} = model.marker_data.CA_(:,impulse_ends(3):impulse_ends(4));
        data_struct(i).marker_data.cst_ma{1} = model.marker_data.CST(:,impulse_ends(1):impulse_ends(2));
        data_struct(i).marker_data.cst_ma{2} = model.marker_data.CST(:,impulse_ends(3):impulse_ends(4));
        data_struct(i).marker_data.mh1_ma{1} = model.marker_data.MH1(:,impulse_ends(1):impulse_ends(2));
        data_struct(i).marker_data.mh1_ma{2} = model.marker_data.MH1(:,impulse_ends(3):impulse_ends(4));
          
        data_struct(i).flex_ankle{1}   = model.sagittal_arch.flex_ankle(1,impulse_ends(1):impulse_ends(2));% -  model.marker_data.CST(1:3,impulse_ends(1));
        data_struct(i).flex_ankle{2}   = model.sagittal_arch.flex_ankle(1,impulse_ends(3):impulse_ends(4));% -  model.marker_data.CST(1:3,impulse_ends(3));

        dot_angle = @(v1,v2) acosd( dot(v1,v2)/ (norm(v1) * norm(v2)));
        
        
        
        force_local = cell(0,0);
        F_pf= cell(0,0);
        
        
        
        
        
        for k = 1:2
            
            n_pts = length(data_struct(i).elongation{k});
            
            
            
            for kk = 1:n_pts
                inv_foot = invTranspose(foot_pose{k}(:,:,kk));
                % determine the arch height in the foot co-ordinate system
                arch_marker_foot = inv_foot * [arch_marker{k}(1:3,kk);1];
                arch_height(1:3,kk) = arch_marker_foot(1:3); % in the y co-ordinate in the foot co-ordinate system
                
                
                % determine helical axis parameters
                [rot_temp,norm_temp,~,ax_pt_temp] = helical(cal_met_pose{k}(:,:,kk));
                %convert the helical parameters to be in global 
                norm_glob_temp = cal_pose{k}(:,:,kk)*[norm_temp;0];
                ax_pt_glob_temp = cal_pose{k}(:,:,kk)*[ax_pt_temp;1];
                
                helical_norm_ax = norm_glob_temp(1:3)/norm(norm_glob_temp(1:3));
                helical_rot = rot_temp;
                helical_ax_pt = ax_pt_glob_temp(1:3);
                % determine the force vector resolved in the sagittal plane of
                % the foot - using a static model to determine the joint forces
                proj_perp = dot(force{k}(:,kk),sagittal_plane{k}(:,kk)) * sagittal_plane{k}(:,kk); % project the vector in plane normal direcion
                force_proj = force{k}(:,kk) - proj_perp; % proj_perp + plane_vec = v_proj as they are the two components in that plane; solve for plane_vec
                
                % convert the projected force into the sagittal co-ordinate
                % system to get two components to input into the two
                % dimensional model
                force_local_temp = inv_foot * [force_proj;0];
                force_local{k}(:,kk) = force_local_temp(1:3)/(data_structFR(subj_FR_ind).weight*9.81); % negative to be in the right direction
                theta = dot_angle(-ca_cst_proj{k}(:,kk),mh1_ca_proj{k}(:,kk));
                phi = dot_angle(-mh1_cst_proj{k}(:,kk),-mh1_ca_proj{k}(:,kk));
           %--------------- determine the moment around the helical axis
                % start by determining the point of intersection of the
                % force vector and the plane made by the foot
                pN = foot_pose{k}(1:3,2,kk); % plane normal
                D = -dot(pN,foot_pose{k}(1:3,4,kk));
%               
                inter_pt = COP{k}(:,kk) + force{k}(:,kk) * (dot(pN,COP{k}(:,kk)) + D) / dot(pN,-force{k}(:,kk));
                % this is the raised COP
                
                markerPlottingVerification
                pause(0.1)
                
                
                
                %    --------- calculate the plantar fascia force from the model
                %     determine lengths of truss sides and put them in the foot co-ordinate
                %     system
                l_mh_temp = inv_foot * [mh1_cst_proj{k}(:,kk);0];
                l_mh = l_mh_temp(1:3);
                
                l_ca_temp = inv_foot * [ca_cst_proj{k}(:,kk);0];
                l_ca = l_ca_temp(1:3);
                
%                 force_loc_temp =  inv_foot * [force_local{k}(:,kk);0];
                F = force_local{k}(:,kk);% force_loc_temp(1:3);
                
                
                A_ext = [1 0 0; 0 1 1; 0 0 (abs(l_mh(1))+abs(l_ca(1)))];
                B_ext = [-F(1);-F(2) ; (-F(1)* abs(l_mh(2))-F(2)*abs(l_mh(1)))];
                f_ext(1:3) = A_ext\B_ext;
                
                F_M = f_ext(1:2);
                F_C= [0 ; f_ext(3)];
                
                
                F_CA_mag = f_ext(3)/-cosd(theta);
                F_CA = F_CA_mag * [-cosd(theta); -sind(theta)];
                
                F_pf = F_CA(2);
                
                F_MH_mag = 1/cosd(phi) *(F(1) - F_CA(1));
                F_MH = F_MH_mag * [cosd(phi); sind(phi)];
                
                
                
                data_struct(i).force_localx{k}(:,kk) =  force_local{k}(1,kk);
                data_struct(i).force_localy{k}(:,kk) =  force_local{k}(2,kk);
                data_struct(i).theta{k}(kk) = theta;
                data_struct(i).phi{k}(kk) = phi;
                data_struct(i).F_pf{k}(kk) =  F_pf;
                data_struct(i).F_Mx{k}(kk) =  F_M(1);
                data_struct(i).F_My{k}(kk) =  F_M(2);
                data_struct(i).F_C{k}(kk) =  F_C(2);
                data_struct(i).F_CAx{k}(kk) =  F_CA(1);
                data_struct(i).F_CAy{k}(kk) =  F_CA(2);
                data_struct(i).F_MHx{k}(kk) =  F_MH(1);
                data_struct(i).F_MHy{k}(kk) =  F_MH(2);
                
                data_struct(i).helical_norm_ax{k}(:,kk) = helical_norm_ax;
                data_struct(i).helical_rot{k}(kk) = helical_rot;
                data_struct(i).helical_ax_pt{k}(:,kk) =  helical_ax_pt ;
                
                if kk > 2 % start at the second so I don't have to make a second loop to differentiate
                    data_struct(i).arch_vel{k}(1:3,kk-1) = diff(arch_height(1:3,[kk-1,kk]),1,2) * mot_fr;
                    data_struct(i).power{k}(kk-1)  = dot(F,data_struct(i).arch_vel{k}(:,kk-1)); % local force and local arch velocity (not that it matters for arch vel)
                    data_struct(i).elong_speed{k}(kk-1) = diff(lengthPF{k}([kk-1,kk])) * mot_fr;
                    data_struct(i).elong_power{k}(kk-1) = data_struct(i).elong_speed{k}(:,kk-1) * data_struct(i).F_pf{k}(kk-1);
                    
                    data_struct(i).omega{k}(kk-1) = diff(data_struct(i).helical_rot{k}([kk-1,kk])) * mot_fr;
                end
                
                
            end
            data_struct(i).arch_height{k} = arch_height(2,:);% in the local co-ordinate system it's in the y direction for height
            clearvars('arch_height')
%             data_struct(i).helical_vals{k} = helical_vals;
            
            
            data_struct(i).arch_energy{k}  = cumtrapz(data_struct(i).elongation{k},data_struct(i).F_pf{k}(1,:));
            
            data_struct(i).work{k} = cumtrapz(0:1/mot_fr:(n_pts-2)/mot_fr,data_struct(i).power{k});
            data_struct(i).work{k}(n_pts-1) = data_struct(i).work{k}(n_pts - 2); % add an extra frame for consistency
            % determine the max energy,
            data_struct(i).energy_max{k} = max(data_struct(i).arch_energy{k});
            data_struct(i).energy_final{k} = (data_struct(i).arch_energy{k}(end));
            data_struct(i).energy_return{k} = data_struct(i).energy_final{k} - data_struct(i).energy_max{k};
            data_struct(i).energy_ratio{k} = (data_struct(i).arch_energy{k}(end))/data_struct(i).energy_max{k};
            
            data_struct(i).work_max{k} = min(data_struct(i).work{k});
            data_struct(i).work_final{k} = data_struct(i).work{k}(end);
            data_struct(i).work_return{k} = data_struct(i).work_final{k} - data_struct(i).work_max{k};
            data_struct(i).work_ratio{k} = data_struct(i).work_final{k}/data_struct(i).work_max{k};
            
            data_struct(i).mla_max{k} = max(data_struct(i).MLA{k});
            data_struct(i).length_max{k} = max(data_struct(i).length_pf{k});
            
            n_pts_new = kk;
            n_pts_mid = floor(n_pts_new/2);
            
            % take the mean of 5 data pts at the 50 % of the trial
            data_struct(i).energy_50{k} = mean(data_struct(i).arch_energy{k}(n_pts_mid-2:n_pts_mid+2));
            data_struct(i).work_50{k} = mean(data_struct(i).work{k}(n_pts_mid-2:n_pts_mid+2));
            
        end
      
        
        disp(['The trial ' root_files{j}(1:end-4) ' is added to the data structure'])
        
        i = i+1;
    end
end

save([direc 'energy_datastruct.mat'],'data_struct')






