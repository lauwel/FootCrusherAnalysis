% close all
% clear
% 
% 
% 
% 
% %%
% load('/home/lauren/Desktop/Motion_Dec15/energy_export.mat')
% save_var = boxplot_data;
% %%
% % centred_data = [];%save_var(:,1:3);
% % p = size(save_var,2); % number of variables
% % % subj = unique(save_var(:,3))';
% % for i = 1:p
% % centred_data(:,end+1) = (save_var(    :,i)-mean(save_var(:,i)))/std(save_var(:,i));
% % end
% i = 1;
% g1 = save_var(:,[1]);
% g2 = save_var(:,[2]);
% % [h,p] = ttest2(g1(:),g2(:));
% % disp(sprintf('P45 - d45 %1.6f',p))
% [h,p] = ttest(g1(:),g2(:));
% disp(sprintf('P45 - d45 %1.6f',p))
% % disp(sprintf('05 bw - 1 bw %1.3f',p))
% %2 way repeated measures 2 way repeated measures anova- post hoc after
% 

%% determine the two way anova repeated measures indices and structure
clear
clc
% clearvars -except data_struct
c = colormap(parula(16));
set(gca,'ColorOrder',c)
c = get(gca,'ColorOrder');
load('/home/lauren/Desktop/MotionDataAus_Nov2017/energy_datastruct.mat')
close all
clearvars('anovaS','f1','f2','sub')
ind = 0;
imp = [1];
subj_exc = 12;%[3,4,5,6,7,8,9,10]; EXCLUDE SUBJECT
for i = 1:108
    for k = imp; %impulse
        ind = ind + 1;
        n_pts = length(data_struct(i).helical_ax_local{k});
        mid_pt = round(n_pts/2);
        
        anovaS.midtoe(ind) = data_struct(i).toe_mid{k};
        anovaS.ROMtoe(ind) = data_struct(i).ROMtoe{k};
        anovaS.minarch(ind) = data_struct(i).min_arch_height{k};
        anovaS.ROMarch(ind) = data_struct(i).ROM_arch{k};
        anovaS.maxMLA(ind) = data_struct(i).mla_max{k};
        anovaS.ROMphi(ind) =  data_struct(i).ROM_phi{k};
        
        anovaS.maxelong(ind) = data_struct(i).max_elong{k};
        anovaS.ROMelong(ind) = data_struct(i).ROM_elong{k};
        anovaS.workratio(ind) = data_struct(i).work_ratio{k};
        anovaS.workmax(ind) = data_struct(i).work_max{k};
        anovaS.workreturn(ind) = data_struct(i).work_return{k};
        anovaS.workfinal(ind) = data_struct(i).work_final{k};
        
        anovaS.midtrans(ind) = data_struct(i).mid_trans{k};
        anovaS.medlatarch(ind) = data_struct(i).ROM_medlatarch{k} ;
        
        pcaS.hel_ax_local(ind,:) = mean(data_struct(i).helical_ax_local{k}(:,mid_pt-2:mid_pt+2),2);
        pcaS.hel_ax_orient(ind,:) = data_struct(i).hel_ax_stable{k}';
        % make arrays with 1's where the condition is met in that row
        % all 1x108
        cond_d45(ind,1) = ~isempty(regexp(data_struct(i).conditions, 'd45', 'once'));
        cond_neut(ind,1) = ~isempty(regexp(data_struct(i).conditions, 'neut', 'once'));
        cond_p45(ind,1) = ~isempty(regexp(data_struct(i).conditions, 'p45', 'once'));
        
        cond_slow(ind,1) = ~isempty(regexp(data_struct(i).trials, 'stat', 'once'));
        cond_fast(ind,1) = ~isempty(regexp(data_struct(i).trials, 'dyn', 'once'));
        
        cond_light(ind,1) = ~isempty(regexp(data_struct(i).trials, '05', 'once'));
        cond_heavy(ind,1) = ~isempty(regexp(data_struct(i).trials, '1', 'once'));
    end
end

condtoe =  cond_d45 | cond_p45 ;%| cond_neut; %logical array with either of these conditions
condspeed = cond_slow | cond_fast ;%
condweight = cond_heavy ;

trial_ind = condtoe & condspeed & condweight;

subj_ind = find(trial_ind);
ind = 1;
for i = 1:length(subj_ind)
    subj_num = str2double(data_struct(subj_ind(i)).subjects(2:3));
    if subj_num  ~= subj_exc
        sub(ind,1) = subj_num;
        ind = ind + 1;
    else
       trial_ind(subj_ind(i)) = 0; 
    end
end
% mark each trial
toe(cond_d45) = 1;
toe(cond_neut) = 2;
toe(cond_p45) = 3;

speed(cond_slow) = 1;
speed(cond_fast) = 2;

weight(cond_light) = 1;
weight(cond_heavy) = 2;

% get all of the indices we need to add
% trial_ind = find(cond_all == 1);

f1 = toe(trial_ind)';
f2 = speed(trial_ind)';

toeDF = find(f1 == 1 & f2 == 2);
toePF = find(f1 == 3 & f2 == 2);
toeDS = find(f1 == 1 & f2 == 1);
toePS = find(f1 == 3 & f2 == 1);

%% do the two way rm anova

nimp = length(imp);

fie = fields(anovaS); % all the fields
% fie{1} = 'maxelong';
for i= 1:length(fie) % for each field conduct a two way rm anova
anova_table = anovaS.(fie{i})(trial_ind); % make a new table every time

disp(fie(i))

stats = rm_anova2(anova_table,sub,f1,f2,{'toes','speed'})

fprintf('Dorsiflexed average is %0.4f +- %0.4f \n',mean([anova_table(toeDS),anova_table(toeDF)]),std([anova_table(toeDS),anova_table(toeDF)]));
fprintf('Plantarflexed average is %0.4f +- %0.4f \n',mean([anova_table(toePS), anova_table(toePF)]),std([anova_table(toePS), anova_table(toePF)]));

fprintf('Fast average is %0.4f +- %0.4f \n',mean([anova_table(toePF),anova_table(toeDF)]),std([anova_table(toePF),anova_table(toeDF)]));
fprintf('Slow average is %0.4f +- %0.4f \n',mean([anova_table(toePS), anova_table(toeDS)]),std([anova_table(toePS), anova_table(toeDS)]));

fprintf('Slow plant average is %0.4f +- %0.4f \n',mean([anova_table(toePS)]),std([anova_table(toePS)]));
fprintf('Fast plant average is %0.4f +- %0.4f \n',mean([anova_table(toePF)]),std([anova_table(toePF)]));

fprintf('Slow dors average is %0.4f +- %0.4f \n',mean([anova_table(toeDS)]),std([anova_table(toeDS)]));
fprintf('Fast dors average is %0.4f +- %0.4f \n',mean([anova_table(toeDF)]),std([anova_table(toeDF)]));



% fprintf('Average  between dorsi + plantar is %0.4f +- %0.4f \n',mean([anova_table(toeDS)-anova_table(toePS),anova_table(toeDF)-anova_table(toePF)]),std([anova_table(toeDS)-anova_table(toePS),anova_table(toeDF)-anova_table(toePF)]));


% figure;
% boxplot([anova_table(toeDS)' anova_table(toePS)' anova_table(toeDF)' anova_table(toePF)'], {'DorsiSlow' 'PlantarSlow' 'DorsiFast'  'PlantarFast'})
% title(fie(i))
% 
% drawnow
% 
% figure;
% plot(unique(sub),[anova_table(toeDS)' anova_table(toePS)' anova_table(toeDF)' anova_table(toePF)'],'x')

% title(fie(i))
% legend({'DorsiSlow' 'PlantarSlow' 'DorsiFast'  'PlantarFast'})
% drawnow

end
%% make the paired subject plot
leg_str = {'3','4','5','6','7','8','9','10','11'};
var_name = 'workratio';

anova_table = anovaS.(var_name)(trial_ind); % make a new table every time

% make the paired plots
figure;
subplot(1,4,1)
plot([anova_table(toeDS);anova_table(toeDF)])
title('Slow-fast, dorsiflexed')
legend(leg_str)
subplot(1,4,2)
plot([anova_table(toePS);anova_table(toePF)])
title('Slow-fast, plantarflexed')
legend(leg_str)
subplot(1,4,3)
plot([anova_table(toeDS);anova_table(toePS)])
title('Dors-Plant, slow')
legend(leg_str)
subplot(1,4,4)
plot([anova_table(toeDF);anova_table(toePF)])
title('Dors-Plant, fast')
legend(leg_str)
%make an overall linear model
set(gca,'ColorOrder',c)
%%
mdl = fitlm(anovaS.midtoe(trial_ind),anovaS.maxelong(trial_ind));
figure;
plot(anovaS.midtoe(trial_ind),anovaS.maxelong(trial_ind),'x')
hold on;
plot(anovaS.midtoe(trial_ind),mdl.Fitted)


% show all the subject specific relationships
c = [c;c];
figure;
ns = 9; % number of subjects
nt = 12; % number of trials
ind = 1;
for i = 1:ns
    % make it so that the trials relevant to that subject are kept and then
    % put into an empty array so none of the other subject's values show up
    % for that round
    if length(find((subj_exc-2) == i)) == 1
        continue
    end
    trial_ind_temp = trial_ind((i-1)*nt+1:(i*nt));
    trial_ind_sub = logical(zeros(ns*nt,1));
    trial_ind_sub((i-1)*nt+1:(i*nt)) = trial_ind_temp;
    mdl = fitlm(anovaS.midtoe(trial_ind_sub),anovaS.maxelong(trial_ind_sub));
    slopes(ind) = mdl.Coefficients{2,1};
    %save subject data in a specific cell array
    toe_subj{ind} = anovaS.midtoe(trial_ind_sub);
    elong_subj{ind} = anovaS.maxelong(trial_ind_sub);
    
    hold on;
    plot(anovaS.midtoe(trial_ind_sub),anovaS.maxelong(trial_ind_sub),'.','color',c(i,:),'markersize',15);
    plot(anovaS.midtoe(trial_ind_sub),mdl.Fitted,'color',c(i,:));
    ind = ind+1;
end


% figure;
%% compare all of the slopes for each subject
ns = length(unique(sub));
ntr = length(unique(f1))*length(unique(f2)); % number of different trial types
for i = 1:ns
    
    for j = i+1:ns
        tbl = table([elong_subj{i}'; elong_subj{j}'],[toe_subj{i}'; toe_subj{j}'],[repmat(i,ntr,1);repmat(j,ntr,1)],'VariableNames',{'Elongation','ToeDorsiflexion','Subject'});
        %         tbl.Subject = categorical(tbl.Subject);
        
        %         mdl = fitlm(tbl,'Elongation~ToeDorsiflexion+Subject'
        sprintf('Subject %i and %i',i+2,j+2);
        mdl = fitlm(tbl,'interactions','ResponseVar','Elongation',...
            'PredictorVars',{'ToeDorsiflexion','Subject'},...
            'CategoricalVar',{'Subject'});
        
        sig_array(i,j) = mdl.Coefficients{4,4};
        
    end
    
end

%% now compare them all to the average 
clearvars('sig_array')
for i = 1:ns
    j = 1:ns;
    si = find(j == i);
    j(si) = [];% remove the subject being compared from the average
    
    tbl = table([elong_subj{i}, [elong_subj{j}]]',[toe_subj{i}, [toe_subj{j}]]',[repmat(i,ntr,1);repmat(0,ntr*length(j),1)],'VariableNames',{'Elongation','ToeDorsiflexion','Subject'});
    
    sprintf('Subject %i',i)
    mdl = fitlm(tbl,'interactions','ResponseVar','Elongation',...
        'PredictorVars',{'ToeDorsiflexion','Subject'},...
        'CategoricalVar',{'Subject'});
    
    sig_array(i) = mdl.Coefficients{4,4};
    
end


%% save the neutral values of the trial
% taking the mean of the first 5 values of the neutral trial, first impulse,
% % 1 bw, slow trial


% for i = 1:length(subj_ind)
% neut_norm_vals(i,1) = mean(data_struct(subj_ind(i)).F2Ps{imp}(1:5));
% % std(data_struct(subj_ind(i)).F2Ps{imp}(1:5))
% end
% 
% save('/home/lauren/Desktop/MotionDataAus_April2017/neut_toe_angles.mat','neut_norm_vals')


%% run pca on helical axes

% run pca on the helical axis for the trials indicated

pca_array = pcaS.hel_ax_orient(trial_ind,:); % global
% pca_array = pcaS.hel_ax_local(trial_ind,:); % local
[coeff,score,latent,tsquared,explained,~] = pca(pca_array);

disp('Explained variance')

explained

% plot levels of the wave

    % toe cond
v1 = [toeDF toeDS];
v2 = [toePF toePS];
%     % speed
% v1 = [toeDF; toePF];
% v2 = [toeDS; toePS];

comp = 2; % which component


sc = [-2,-1,0,1,2]';
mu = mean(pca_array([v1, v2],:));
muv1 = (pca_array([v1],:));
muv2 = (pca_array([v2],:));

figure;
% get projections in all of the planes of the foot
for i = 1:18
    %local
% muv1_sag = [muv1(i,1:2),0];
% muv2_sag = [muv2(i,1:2),0];
% muv1_fro = [0,muv1(i,2:3)];
% muv2_fro = [0,muv2(i,2:3)];
% muv1_tra = [muv1(i,1),0,muv1(i,3)];
% muv2_tra = [muv2(i,1),0,muv2(i,3)];

%global
muv1_tra = [muv1(i,1:2),0];
muv2_tra = [muv2(i,1:2),0];
muv1_fro = [0,muv1(i,2:3)];
muv2_fro = [0,muv2(i,2:3)];
muv1_sag = [muv1(i,1),0,muv1(i,3)];
muv2_sag = [muv2(i,1),0,muv2(i,3)];

sag_ang(i) = acosd(dot(muv1_sag,muv2_sag)/(norm(muv1_sag)*norm(muv2_sag)));
fro_ang(i) = acosd(dot(muv1_fro,muv2_fro)/(norm(muv1_fro)*norm(muv2_fro)));
trans_ang(i) = acosd(dot(muv1_tra,muv2_tra)/(norm(muv1_tra)*norm(muv2_tra)));
% fprintf('sag = %0.2f, fro = %0.2f, trans = %0.2f \n',sag_ang(i),fro_ang(i),trans_ang(i))

% plot each angle view
% subplot(1,3,1)
% h = plotvector3([0 0 0],muv1(i,:),'k');
% hold on
% h = plotvector3([0 0 0],muv2(i,:),'r');
% view(0,90) % x y - sagittal
% hold off;
% axis equal
% 
% subplot(1,3,2)
% h = plotvector3([0 0 0],muv1(i,:),'k');
% hold on
% h = plotvector3([0 0 0],muv2(i,:),'r');
% hold off;
% view(90,0) % y z - frontal
% axis equal
% 
% subplot(1,3,3)
% h = plotvector3([0 0 0],muv1(i,:),'k');
% hold on
% h = plotvector3([0 0 0],muv2(i,:),'r');
% hold off;
% view(0,0) % x z - transverse
% axis equal
end

fprintf('sagittal %0.2f +- %0.2f \n',mean(sag_ang),std(sag_ang))
fprintf('frontal %0.2f +- %0.2f \n',mean(fro_ang),std(fro_ang))
fprintf('transverse %0.2f +- %0.2f \n',mean(trans_ang),std(trans_ang))
% figure;
% h = plotvector3([0 0 0],muv1,'k');
% hold on
% h = plotvector3([0 0 0],muv2,'r');


for comp = 1:2
    figure;
    for t = 1:5
        % t = 1; % which level of std
        
        
        score_lvl = score([v1; v2],comp) + sc(t) * std(score([v1; v2],comp));
        overall_wave = score_lvl*coeff(:,comp)'+repmat(mu,size([v1; v2],1),1);
        hold all
        
        % note : using the global axes these are the orientations
        xlabel('X (anterior - posterior)')
        ylabel('Y (medial - lateral)');
        zlabel('Z (prox- distal)')
        
        h = plotvector3([0 0 0],mean(overall_wave(:,:)),[1-t*0.15 1-t*0.15 1]);
    end
    
    axis equal
end

[h,p] = ttest2(score(v1,:),score(v2,:));
fprintf('P45 - d45 %1.6f \n',p)
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
c = get(gca,'ColorOrder');
%---- by toe angle
figure;
plot(score(v1,1),score(v1,2),'x','color',c(1,:));
hold on
% plot(score(neut,1),score(neut,2),'x','color',c(2,:));

plot(score(v2,1),score(v2,2),'x','color',c(3,:));

% plot(score(high_var,1),score(high_var,2),'x','color','r')

xlabel('PC1')
ylabel('PC2')

[h,p] = ttest2(score(v1,:),score(v2,:));
fprintf('P45 - d45 %1.6f \n',p)
%----------------------------------------------------------------------
%%
% save_var(:,4) = save_var(:,4)*75;
[loading, score,eig,tsquare,expl,~] = pca(centred_data);
% plot of each component
figure()
hold on
n = length(loading);
% n = 3; % look at the contributions from 3 components
% each column in loading is the contribution to one principal component
% from each variable
bar(loading(:,1:n)')%,'color',[1-i*0.15 i*0.10 1-i*0.1],'linewidth',4)

ylabel('Contribution')
xlabel('PC')
% legend('PC1','PC2','PC3','PC4','PC5','PC6')
% legend('strain','MLA','pron','rot')
legend(save_var_headers(5:end))
% determine the indices of the group
p45 = find(save_var(:,2) == 1);
neut = find(save_var(:,2) == 2);
d45 = find(save_var(:,2) == 3);
% this indexing won't work if there are different numbers of subjects
leg_str = cell(0,0);
for i = subj
subj_ind(i,:) = find(save_var(:,3) == i);
end
%%
% plot the components

pca_a = 1;
for i = 2:3
pca_b = i;

cmap = colormap('parula');


figure()
hold on;
for i = subj
h = plot(score(subj_ind(i,:),pca_a),score(subj_ind(i,:),pca_b),'x');
set(h,'markerfacecolor',cmap(i*5,:),'markersize',15)
leg_str{end+1} = num2str(i);
end
legend(leg_str)
xlabel(['PCA ' num2str(pca_a)])
ylabel(['PCA ' num2str(pca_b)])
axis equal

figure()
hold on;
plot(score(p45,pca_a),score(p45,pca_b),'rx','markersize',15)
plot(score(neut,pca_a),score(neut,pca_b),'bx','markersize',15)
plot(score(d45,pca_a),score(d45,pca_b),'gx','markersize',15)
legend('p45','neut','d45')

xlabel(['PCA ' num2str(pca_a)])
ylabel(['PCA ' num2str(pca_b)])
axis equal



end

%% PLS


[XL,YL,XS,YS,beta,pctvar,mse,stats] = plsregress(centred_data(:,2:end),centred_data(:,1),2);

figure()
bar(XL')
legend(save_var_headers(6:end))
expl
pctvar

figure()
% plot(XS(subj_ind,1),XS(subj_ind,2),'o')
leg_str = cell(0,0);
hold on
for i = subj
h = plot(XS(subj_ind(i,:),1),'.');
set(h,'markerfacecolor',cmap(i*5,:),'markersize',15)
leg_str{end+1} = num2str(i);
end

legend(leg_str)
%%
var1 = 4;
for i = [1,2,3,5]
var2 = i;
mod = fitlm(centred_data(:,var1),centred_data(:,var2),'linear')
hf = figure();
plot(save_var(:,4+var2),save_var(:,4+var1),'x');
xlabel(save_var_headers(4+var2))
ylabel(save_var_headers(4+var1))
% h_a = get(hf,'Children');
% axis_pos = get(h_a,'Position');
text(0,0,['R^2 = ' num2str(mod.Rsquared.Adjusted)])
disp(mod.Rsquared)
end

%% Multiple Regression
xvals = [save_var(:,6) save_var(:,7) save_var(:,8) save_var(:,6).*save_var(:,7) save_var(:,7).*save_var(:,8) save_var(:,6).*save_var(:,8)];
yvals = save_var(:,5);
mod = fitlm(xvals,yvals,'linear')

%% Run a t-test to see if impulses are significantly different


% imp1 = load('tracking_final_mid50_impulse1.mat');
% imp2 = load('tracking_final_mid50_impulse2.mat');
imp1 = load('/home/lauren/tracking_final_range.mat');
% for i = 5:9
% % paired ttest
%     [h,p] = ttest(imp1.save_var(:,i),imp2.save_var(:,i));
%     if h == 1
%         disp([imp1.save_var_headers{i} ' have unequal means between impulses (p = ' num2str(p) '). '])
%     else
%         disp([imp1.save_var_headers{i} ' have equal means between impulses (p = ' num2str(p) '). '])
%     end
% end

p45 = find(save_var(:,2) == 1);
neut = find(save_var(:,2) == 2);
d45 = find(save_var(:,2) == 3);
stat05 = find(save_var(:,1) == 1);
stat1 = find(save_var(:,1) == 2);
dyn05 = find(save_var(:,1) == 3);
dyn1 = find(save_var(:,1) == 4);
close all
for i = 5:9
figure()
boxplot(imp1.save_var(:,i),imp1.save_var(:,1))
title([imp1.save_var_headers{i} '- classified by ' imp1.save_var_headers{1}])

pause(0.1)

drawnow
figure()
boxplot(imp1.save_var(:,i),imp1.save_var(:,2))
title([imp1.save_var_headers{i} '- classified by ' imp1.save_var_headers{2}])
drawnow
pause(0.1)
figure()
boxplot(imp1.save_var(:,i),imp1.save_var(:,3))
title([imp1.save_var_headers{i} '- classified by ' imp1.save_var_headers{3}])
drawnow

pause(0.1)
imp1.save_var_headers{i} 
[h,p] = ttest2(save_var(p45,i),imp1.save_var(neut,i));
disp(sprintf('P45 - neut %1.3f',p))

[h,p] = ttest2(imp1.save_var(d45,i),imp1.save_var(neut,i));
disp(sprintf('d45 - neut %1.3f',p))

[h,p] = ttest2(imp1.save_var(p45,i),imp1.save_var(d45,i));
disp(sprintf('P45 - d45 %1.3f',p))

[h,p] = ttest2(imp1.save_var(stat05,i),imp1.save_var(stat1,i));
disp(sprintf('stat05- stat1 %1.5f',p))

[h,p] = ttest2(imp1.save_var([stat05;stat1],i),imp1.save_var([dyn05;dyn1],i));
disp(sprintf('slow- fast %1.5f',p))
[h,p] = ttest2(imp1.save_var([stat05;dyn05],i),imp1.save_var([stat1;dyn1],i));
disp(sprintf('.5BW vs 1BW %1.5f',p))

end
