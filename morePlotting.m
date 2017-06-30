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

%% figure;
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
    if isempty(regexp(data_struct(i).trials, 'stat'))
        continue
    end
    
    for k = 1:2
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
    var_invest = 'moment_power';
    plot(data_struct(i).(var_invest){k}(iii,:)','color',cc);
    title(sprintf('Subject %i',round(i/12)))
    ylabel(var_invest)
    xlabel('Frame')
%         ylim([0 20])
    
    end
end
end

%% run pca on helical axes

k = 1; %impulse
for i = 1:108
    pca_array(i,:) = data_struct(i).hel_ax_stable{k};
    
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

v1 = light;
v2 = heavy;


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
ylabel('Y (medial - lateral)')
zlabel('Z (prox- distal)');
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


figure;
a1 = 9;
a2 = 7;
plotvector3(data_struct(a1).hel_pt_stable{k},data_struct(a1).hel_ax_stable{k}*30,c(5,:));
hold on

plotvector3(data_struct(a2).hel_pt_stable{k},data_struct(a2).hel_ax_stable{k}*30,c(7,:));
legend('p45','d45')

xlabel('X (anterior - posterior)')
ylabel('Y (medial - lateral)')
zlabel('Z (prox- distal)');
axis equal
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

