function marker_data = addMEandLEtoMarkerData(direc,file)
% add the medial and lateral epicondyles using the static trials
% input file  and directory are needed to load and save files in
% appropriate places
%   outputs the new marker structure if output variable is
%   selected
subj = file(1:3);

load([direc subj '_statcal_processedmotion.mat'])
% load the static file, make the cluster coordinate system, determine
% marker position in TC coordinate system
lat_epi = marker_data.LE_(:,10);
med_epi = marker_data.ME_(:,10);
lat_mal = marker_data.LM_(:,10);
med_mal = marker_data.MM_(:,10);


knee_centre = (med_epi + lat_epi)/2;
ankle_centre = (med_mal + lat_mal)/2;

% T_GTA = makeCoordSysYZ([ankle_centre,knee_centre,lat_epi],ankle_centre);
T_GTC = makeCoordSysYZ([marker_data.TB1(:,10),marker_data.TB2(:,10),marker_data.TB3(:,10)],marker_data.TB1(:,10));

% T_TCTA = invTranspose(T_GTC) * T_GTA;

lat_epi_TC = invTranspose(T_GTC) * [lat_epi;1]; % in the tibia cluster coordinate system
lat_epi_TC = lat_epi_TC(1:3,:);

med_epi_TC = invTranspose(T_GTC) * [med_epi;1]; % in the tibia cluster coordinate system
med_epi_TC = med_epi_TC(1:3,:);

lat_mal_TC = invTranspose(T_GTC) * [lat_mal;1]; % in the tibia cluster coordinate system
lat_mal_TC = lat_mal_TC(1:3,:);

med_mal_TC = invTranspose(T_GTC) * [med_mal;1]; % in the tibia cluster coordinate system
med_mal_TC = med_mal_TC(1:3,:);

% now load the trial missing the ME and LE markers. Remake the global-
% cluster co-ordinate system and transform the markers back into global and
% save them
load([direc file])
flags = zeros(1,4); % make a set of flags to know which markers need to be redone
for i = 1:marker_data.nFrames
    T_GTC = makeCoordSysYZ([marker_data.TB1(:,i),marker_data.TB2(:,i),marker_data.TB3(:,i)],marker_data.TB1(:,i));
    
    % replace the markers that are missing
    if ~isfield(marker_data,'LE_')|| flags(1) == 1
        flags(1) = 1;
        lat_epi = T_GTC * [lat_epi_TC;1]; % in the tibia cluster coordinate system
        marker_data.LE_(:,i) = lat_epi(1:3,:);
    end
    if ~isfield(marker_data,'ME_')|| flags(2) == 1
        flags(2) = 1;
        med_epi = T_GTC * [med_epi_TC;1]; % in the tibia cluster coordinate system
        marker_data.ME_(:,i) = med_epi(1:3,:);
    end
    if ~isfield(marker_data,'LM_')|| flags(3) == 1
        flags(3) = 1;
        lat_mal = T_GTC * [lat_mal_TC;1]; % in the tibia cluster coordinate system
        marker_data.LM_(:,i) = lat_mal(1:3,:);
    end
    if ~isfield(marker_data,'MM_') || flags(4) == 1
        flags(4) = 1;
        med_mal = T_GTC * [med_mal_TC;1]; % in the tibia cluster coordinate system
        marker_data.MM_(:,i) = med_mal(1:3,:);
    end
    
end

save([direc file],'marker_data','force_data')
disp(['Resaved file : ' file 'with medial and lateral epicondyle positions']);
if nargout == 1
    vargout = marker_data;
end
% plot checking if needed
% plot3quick([marker_data.LE_(:,10),marker_data.ME_(:,10),marker_data.TB1(:,10),marker_data.TB2(:,10),marker_data.TB3(:,10),marker_data.MM_(:,10),marker_data.LM_(:,10),marker_data.MH1(:,10),marker_data.MH5(:,10)],'k','o')
% axis equal