% c = get(gca,'ColorOrder');
% % close all
% figure;
% ind = 1;
% for i = 1%:108
%     if rem(i,12) == 1
% %         subplot(3,3,ind)
%         ind = ind+1;
%         hold on;
%     end
%     for k = 1
% %     % by trial type
% %     if rem(i,12) >= 1 && rem(i,12) <= 3 
% %         cc = c(1,:);
% %     elseif rem(i,12) >= 4 && rem(i,12) <= 6
% %         cc = c(2,:);
% %     elseif rem(i,12) >= 7 && rem(i,12) <= 9
% %         cc = c(3,:);
% %     elseif rem(i,12) == 10 || rem(i,12) == 11 || rem(i,12) == 0
% %         cc = c(4,:);
% %     end
%         
%         % by toe angle
%         
%     if rem(i,12) == 1 || rem(i,12) == 4 || rem(i,12) == 7|| rem(i,12) == 10
%         cc = c(1,:);
%     elseif rem(i,12) == 2 || rem(i,12) == 5 || rem(i,12) == 8 || rem(i,12) == 11
%         cc = c(2,:);
%         continue
%     elseif rem(i,12) == 3 || rem(i,12) == 6 || rem(i,12) == 9 || rem(i,12) == 0
%         cc = c(3,:);
%         continue
% %     elseif rem(i,12) == 10 || rem(i,12) == 11 || rem(i,12) == 0
% %         cc = c(4,:);
%     end
%     n = length(data_struct(i).moment_arm{k});
%     for j = 1:3:75
%         
%     plotvector3(data_struct(i).helical_ax_pt{k}(:,j),data_struct(i).hel_ax_stable{k}*500,cc);
%     plotvector3(data_struct(i).inter_pt{k}(:,j),data_struct(i).moment_arm{k}(:,j),c(4,:));
%     end
%     if rem(i,12) == 0
%         axis equal
%         view([90 90])
%     end
%     end
% end

%% PLOT THE 9 GRAPHS - each subject independently
figure;
c = get(gca,'ColorOrder');

for iii = 1
    figure
    ind = 1;
    for i = 1:108
        if rem(i,12) == 1
            subplot(3,3,ind)
            ind = ind+1;
            hold on;
        end
%         if isempty(regexp(data_struct(i).trials, 'dyn'))
%             continue
%         end
%         
        for k = 1
            %     % by trial type
            %     if rem(i,12) >= 1 && rem(i,12) <= 3
            %         cc = c(1,:);
            %     elseif rem(i,12) >= 4 && rem(i,12) <= 6
            %         cc = c(2,:);
            %     elseif rem(i,12) >= 7 && rem(i,12) <= 9
            %         cc = c(3,:);
            %     elseif rem(i,12) == 10 || rem(i,12) == 11 || rem(i,12) == 0
            %         cc = c(4,:);
            %     end
            
            % by toe angle
            
            if rem(i,12) == 1 || rem(i,12) == 4 || rem(i,12) == 7|| rem(i,12) == 10
                cc = c(1,:);
            elseif rem(i,12) == 2 || rem(i,12) == 5 || rem(i,12) == 8 || rem(i,12) == 11
                cc = c(2,:);
                continue
            elseif rem(i,12) == 3 || rem(i,12) == 6 || rem(i,12) == 9 || rem(i,12) == 0
                cc = c(3,:);
                %     elseif rem(i,12) == 10 || rem(i,12) == 11 || rem(i,12) == 0
                %         cc = c(4,:);
            end
            var_invest = 'helical_ax_local';
            plot(normalise(data_struct(i).(var_invest){k}(iii,:)'),'color',cc);
            title(sprintf('Subject %i',round(i/12)))
            ylabel(var_invest)
            xlabel('Frame')
%             ylim([-50000 50000])
%             xlim([0 100])
            
        end
    end
end


%% run pca on helical axes

k = 1; %impulse
for i = 1:108
    n_pts = length(data_struct(i).helical_ax_local{k});
    mid_pt_d45 = round(n_pts/2);
    
    pca_array(i,:) = data_struct(i).helical_ax_local{k}(:,mid_pt_d45);
    
    cond_d45(i) = ~isempty(regexp(data_struct(i).conditions, 'd45'));
    cond_neut(i) = ~isempty(regexp(data_struct(i).conditions, 'neut'));
    cond_p45(i) = ~isempty(regexp(data_struct(i).conditions, 'p45'));
    
     
    cond_slow(i) = ~isempty(regexp(data_struct(i).trials, 'stat'));
    cond_fast(i) = ~isempty(regexp(data_struct(i).trials, 'dyn'));
    
    
     
    cond_light(i) = ~isempty(regexp(data_struct(i).trials, '05'));
    cond_heavy(i) = ~isempty(regexp(data_struct(i).trials, '1'));
    
    
    
    
end

d45 = find(cond_d45);
neut = find(cond_neut);
p45 = find(cond_p45);

slow = find(cond_slow);
fast = find(cond_fast);

light = find(cond_light);
heavy = find(cond_heavy);

v1 = d45;
v2 = p45;


[coeff,score,latent,tsquared,explained,~] = pca(pca_array);
comp = 1; % which component
figure;
for t = 1:5
% t = 1; % which level of std

sc = [-2,-1,0,1,2]';
mu = mean(pca_array([v1 v2],:));
score_lvl = score([v1 v2],comp) + sc(t) * std(score([v1 v2],comp));
overall_wave = score_lvl*coeff(:,comp)'+repmat(mu,size([v1 v2]',1),1);

hold all

xlabel('X (anterior - posterior)')

ylabel('Z (prox- distal)');
zlabel('Y (medial - lateral)')
h = plotvector3([0 0 0],mean(overall_wave(:,:)),[1-t*0.15 1-t*0.15 1]);
end
axis equal
%%
c = [c;c];
high_var = find(tsquared>10);
subj_grid = reshape(1:108,12,9)';
ind = 1;
figure;
for i = [1,2,3,7,8,9]
    hold on;
    plot(score(subj_grid(i,:),1),score(subj_grid(i,:),2),'x','color',c(ind,:))
%     plot(score(neut,1),score(neut,2),'kx')
    
    
    
    ind = ind + 1;
    
end
%%
%---- by toe angle
figure;
plot(score(v1,1),score(v1,2),'x','color',c(1,:));
hold on
% plot(score(neut,1),score(neut,2),'x','color',c(2,:));

plot(score(v2,1),score(v2,2),'x','color',c(3,:));

plot(score(high_var,1),score(high_var,2),'x','color','r')

xlabel('PC1')
ylabel('PC2')

[h,p] = ttest2(score(v1,:),score(v2,:));
fprintf('P45 - d45 %1.6f \n',p)
%%
c = get(gca,'ColorOrder');
k = 2;
figure(10);hold on
a1 = 10;
a2 = 12;
for jj = 1:2
    

    j = [a1 a2];
    filename1 = dir(['/home/lauren/Desktop/MotionDataAus_April2017/' data_struct(j(jj)).subjects '_' data_struct(j(jj)).trials '_' data_struct(j(jj)).conditions '_*_processedmotion.mat']);
    filename2 = dir(['/home/lauren/Desktop/MotionDataAus_April2017/' data_struct(j(jj)).subjects '_' data_struct(j(jj)).trials '_' data_struct(j(jj)).conditions '_*_impulse.mat']);
    load(['/home/lauren/Desktop/MotionDataAus_April2017/' filename1.name])
    load(['/home/lauren/Desktop/MotionDataAus_April2017/' filename2.name])
    kk = round((impulse_ends(2)-impulse_ends(1))/2);
    
field_list = fields(marker_data);

for iii = 1:length(field_list)
    if length(field_list{iii}) == 3 % i.e. it's a marker
        
        m = marker_data.(field_list{iii});
      
        scatter3(m(1,impulse_ends(1)+kk),m(2,impulse_ends(1)+kk),m(3,impulse_ends(1)+kk),'MarkerEdgeColor',c(jj,:))
    
    end
end
end


figure(10)
hold on
plotvector3(data_struct(a1).hel_pt_stable{k},data_struct(a1).hel_ax_stable{k}*-80,c(1,:));
plotvector3(data_struct(a1).hel_pt_stable{k},data_struct(a1).hel_ax_stable{k}*80,c(1,:));
% axis equal
% figure(2)
plotvector3(data_struct(a2).hel_pt_stable{k},data_struct(a2).hel_ax_stable{k}*-80,c(2,:));
plotvector3(data_struct(a2).hel_pt_stable{k},data_struct(a2).hel_ax_stable{k}*80,c(2,:));
axis equal
legend('p45','d45')

xlabel('X (anterior - posterior)')
ylabel('Y (medial - lateral)')
zlabel('Z (prox- distal)');
axis equal
%% make a power graph
figure;
hold on;
% h = plot(-normalise(data_struct(42).force_power{1, 1}'/1000)');

h = plot(-normalise(data_struct(42).work{1, 1}')/1000');
set(h,'Linewidth',5)
% set(h2, 'Alpha',0.4,'color',c(1,:))
xlabel('% Trial')
xlim([0 100])
ylabel('Work [J]')
set(gca,'Fontsize',20)
set(h,'color',c(1,:))

%% make an iv file with the helical axes

mot_cal_pose_temp = data_struct(a1).cal_pose{k};
n_pts = size(mot_cal_pose_temp,3);
mid_pt_d45 = round(n_pts/2);

mot_cal_d45 = mot_cal_pose_temp(:,:,mid_pt_d45); % DORSIFLEX POSE

mot_cal_pose_temp = data_struct(a2).cal_pose{k};
n_pts = size(mot_cal_pose_temp,3);
mid_pt_p45 = round(n_pts/2);

mot_cal_p45 = mot_cal_pose_temp(:,:,mid_pt_p45); % PLANTAR FLEX POSE
  
O_cal = load('/home/lauren/Documents/VirtualDesktopTemp/ISB/bones/calc_marker.stack')';
CST = load('/home/lauren/Documents/VirtualDesktopTemp/ISB/bones/CST.stack')';
CPT = load('/home/lauren/Documents/VirtualDesktopTemp/ISB/bones/CPT.stack')';

CIC = (CST + CPT)/2;
x_cal = CIC - O_cal;
temp_cal = O_cal - CST;



y_cal = cross(temp_cal,x_cal);
z_cal = cross(x_cal,y_cal);

% normalise the vectors
x_cal = x_cal/norm(x_cal);
y_cal = y_cal/norm(y_cal);
z_cal = z_cal/norm(z_cal);

% make the pose matrix
xray_cal(1:3,1:4) = [x_cal,y_cal,z_cal,O_cal];
xray_cal(4,1:4) = [0 0 0 1];

cal_d45_xray = invTranspose(mot_cal_d45) * xray_cal;
cal_p45_xray = invTranspose(mot_cal_p45) * xray_cal;
% 
% hel_pt_d45_xray = invTranspose(mot_cal_d45)*[data_struct(a1).hel_pt_stable{k};1];
% hel_pt_p45_xray = invTranspose(mot_cal_p45)*[data_struct(a2).hel_pt_stable{k};1];
% 
% hel_ax_d45_xray = invTranspose(mot_cal_d45)*[data_struct(a1).hel_ax_stable{k};0];
% hel_ax_p45_xray = invTranspose(mot_cal_p45)*[data_struct(a2).hel_ax_stable{k};0];


T_myles = (xray_cal);

hel_pt_d45_xray = T_myles * [data_struct(a1).helical_pt_local{k}(:,mid_pt_d45);1];
hel_pt_p45_xray = T_myles * [data_struct(a2).helical_pt_local{k}(:,mid_pt_p45);1];

hel_ax_d45_xray = T_myles * [data_struct(a1).helical_ax_local{k}(:,mid_pt_d45);0];
hel_ax_p45_xray = T_myles * [data_struct(a2).helical_ax_local{k}(:,mid_pt_p45);0];

hel_pt_d45_xray = hel_pt_d45_xray(1:3);
hel_pt_p45_xray = hel_pt_p45_xray(1:3);
hel_ax_d45_xray = hel_ax_d45_xray(1:3);
hel_ax_p45_xray = hel_ax_p45_xray(1:3);

% figure
% hold on
% plotvector3(hel_pt_d45_xray , hel_ax_d45_xray*80,c(1,:));
% plotvector3(hel_pt_p45_xray , hel_ax_p45_xray*80,c(2,:));
% xlabel('X')
% ylabel('Y')
% zlabel('Z')
close all
figure
plotPointsAndCoordSys([O_cal,CST,CPT],xray_cal)
plotPointsAndCoordSys([O_cal,CST,CPT],mot_cal_d45)
plotPointsAndCoordSys([O_cal,CST,CPT],mot_cal_p45)
plotPointsAndCoordSys([],eye(4))
plotvector3(data_struct(a1).hel_pt_stable{k},data_struct(a1).hel_ax_stable{k}*120,c(1,:));
plotvector3(data_struct(a2).hel_pt_stable{k},data_struct(a2).hel_ax_stable{k}*120,c(2,:));
plotvector3(hel_pt_d45_xray , hel_ax_d45_xray*60,c(1,:));
plotvector3(hel_pt_p45_xray , hel_ax_p45_xray*60,c(2,:));
xlabel('X')
ylabel('Y')
zlabel('Z')
%%
            ivstring = createInventorHeader();
            centroid = hel_pt_d45_xray;
            
            ivstring = [ivstring createInventorArrow(centroid-90*hel_ax_d45_xray,hel_ax_d45_xray,120, 2,c(1,:),0)];
            ivstring = [ivstring createInventorArrow(centroid-90*hel_ax_p45_xray,hel_ax_p45_xray,120, 2,c(2,:),0)];
%             ivstring = [ivstring createInventorArrow(centroid,eig{3},60, 2,[0 0 1],0)];
            
            % Write the extra components to a file
            fid = fopen(['/home/lauren/Documents/VirtualDesktopTemp/ISB/bones/helical_axes.iv'],'w'); 
            %remove the '_coord' to write to the same .iv files and add a
            %CreateInventorLink to the bone file before the ivstring, also
            %comment out the .pos sections below
            fprintf(fid,ivstring);
            fclose(fid);


%%

for i = 25:36
figure;
for k = 1:2
subplot(2,1,k)
plot(data_struct(i).helical_norm_ax{k}','DisplayName','data_struct(43).helical_norm_ax{1, 1}');
end
end




%%

figure; hold on; 

scatter3(pt_q(1),pt_q(2),pt_q(3))
scatter3(pt_p(1),pt_p(2),pt_p(3))
plotvector3(pt_q,pq,'k');
plotvector3(pt_q,hel_ax_stable*85,'r');
plotvector3(pt_q+helical_proj,moment_arm,'b')
axis equal

%% plot the calcaneal co-ordinate system and the global/local poitns

figure % in global

i = 1;
k = 1;
kk = 41;
plotPointsAndCoordSys1(data_struct(i).helical_ax_pt{k}(:,kk),data_struct(i).cal_pose{k}(:,:,kk))



figure; % local

plotPointsAndCoordSys1(data_struct(i).helical_pt_local{k}(:,kk),eye(4))

%% Plot the force arch height graphs





















