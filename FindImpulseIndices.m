
% function to determine the points to track using the force impulses
% will save the four points marking the start and end of each impulse.

clear
clc
%%
for k = 4
    direct = sprintf('/media/lauren/Elements/AustraliaCollection/S%0.2d/SelectedTrials/',k);
    list_files = dir(direct);
    cd(direct)
    file_names = {'dyn05_neut'};
    rf_file = cell(0);
    mot_file = cell(0);
    
    % organize the file names so that there is a motion set and an rf set
    for i = 1:length(file_names)
        for j = 1:length(list_files)
            
            match = strfind(list_files(j).name,file_names{i});
            if numel(match) == 1
                if numel(strfind(list_files(j).name,'image')) == 1
                    rf_file{end+1} = list_files(j).name;
                elseif numel(strfind(list_files(j).name,'motion')) == 1
                    mot_file{end+1} = list_files(j).name;
                end
            end
        end
    end
    
    
    
    
    for i = 1:length(rf_file)
        
        ind_mot = regexp(mot_file{i},'_[12345]_')+1; % find the trial number
        ind_rf = regexp(rf_file{i},'_[12345]_')+1; % find the trial number
        if ~strcmp(mot_file{i}(ind_mot),rf_file{i}(ind_rf))
            disp(mot_file{i});
            disp(rf_file{i});
            error('Motion and RF files do not match. ')
        end
        
        
        raw_mot_data = load(mot_file{i});
        raw_RF_data = load(rf_file{i});
        
        raw_force = raw_mot_data.force_data.Force(1:3,:);
        raw_force = raw_force - repmat(mean(raw_force(1:3,1:10),2),1,length(raw_force));
        force_fr = raw_mot_data.force_data.Frequency;
        force_nfr = raw_mot_data.force_data.NrOfSamples;
        for j = 1:length(raw_force)
            raw_force_norm(j) = norm(raw_force(:,j));
        end
        
        t_rf = 0:1/raw_RF_data.RFheader.framerate:(raw_RF_data.RFheader.nframes-1)/raw_RF_data.RFheader.framerate;
        t_force = 0:1/force_fr:(force_nfr-1)/force_fr;
        % resampled force
        force_res = interp1(t_force,raw_force_norm,t_rf,'nearest');
        
        % determine when the impulse occurs by the shape of the force curve
        
        aveforce(1) = mean(force_res(1:5));
        for t = 3:length(force_res)-2
            aveforce(t-1) = mean(force_res(t-2:t+2));% compute the average
            diff_force_raw(t-2) = aveforce(t-1)-aveforce(t-2); % compute the difference
            %
        end
        nan_ind = isnan(diff_force_raw);
        diff_force_raw(nan_ind) = [];
        diff_force = LowPassButterworth(diff_force_raw,4,2,raw_RF_data.RFheader.framerate);
        %         the largest peaks indicate when the impulse is starting
        impulse_index = find(abs(diff_force) < 8);
        diff_force(impulse_index) = 0;
        figure(); plot(diff_force);
        [~,impulse_ends] = findpeaks(abs(diff_force));
        % it finds the middle of the increasing/decreasing portion of the
        % wave so expand a threshold amount to include the rest of the wave
        
        %     thresh = 30; %for static trials
        thresh = 4;
        
        % catch if we exceed the minimum/maximum indices
        if impulse_ends(1) < 0; impulse_ends(1) = 1; end;
        
        t = 1;
        while t <= length(impulse_ends)
            if rem(t,2) == 1 % if it's odd, we subtract from impulse ends to widen the waveform
                impulse_ends(t) = impulse_ends(t)-thresh/2;
            elseif rem(t,2) == 0
                impulse_ends(t) = impulse_ends(t)+thresh+2;
            end
            if impulse_ends(t) > length(force_res) % if the indices are higher than the force
                impulse_ends(t:end) = [];
            end
            t = t+1;
        end
        
        
        
        disp(impulse_ends)
        figure()
        plot(force_res)
        hold on
        plot(impulse_ends,force_res(impulse_ends),'og')
        title(mot_file{i}(1:ind_mot))
        
        drawnow
        save([direct, mot_file{i}(1:ind_mot) '_impulse.mat'],'impulse_ends')
        
    end
    
end