% function markerPlottingVerification
figure(1);
c = get(gca,'ColorOrder');
% plot the plane of the foot
x = 0:10:400;
y = 0:10:400;
[X,Y] = meshgrid(x,y);
C = ones(length(x),length(y));
Z = (-pN(1)*X -pN(2) * Y -D)/pN(3);
h = mesh(X,Y,Z,C);
set(h, 'FaceAlpha', 0.2)
axis equal
hold on

scatter3(COP{k}(1,kk),COP{k}(2,kk),COP{k}(3,kk),'r')
plot3(inter_pt(1),inter_pt(2),inter_pt(3),'.k','MarkerSize',10)
scatter3(foot_pose{k}(1,4,kk),foot_pose{k}(2,4,kk),foot_pose{k}(3,4,kk),'b')

force_ch = -force{k}(:,kk);
plotvector3(COP{k}(:,kk),force_ch,c(1,:));

% plot the markers for a specific frame
field_list = fields(marker_data);

for iii = 1:length(field_list)
    if length(field_list{iii}) == 3 % i.e. it's a marker
        
        m = marker_data.(field_list{iii});
        if k == 1
        scatter3(m(1,impulse_ends(1)+kk),m(2,impulse_ends(1)+kk),m(3,impulse_ends(1)+kk),'MarkerEdgeColor',c(2,:))
        elseif k == 2 
        scatter3(m(1,impulse_ends(3)+kk),m(2,impulse_ends(3)+kk),m(3,impulse_ends(3)+kk),'MarkerEdgeColor',c(2,:))
        end
    end
end

% plot the helical axis parameters

plotvector3(data_struct(i).helical_ax_pt{k}(:,kk),data_struct(i).helical_norm_ax{k}(:,kk)*100,c(4,:));
plotvector3(data_struct(i).helical_ax_pt{k}(:,kk),-data_struct(i).helical_norm_ax{k}(:,kk)*100,c(4,:));

% 
% 
% scatter3(pt_q(1),pt_q(2),pt_q(3));
% scatter3(pt_p(1),pt_p(2),pt_p(3));
% plotvector3(pt_q,pq,'k');
% plotvector3(pt_q,hel_ax_stable*85,'r');
% plotvector3(pt_q+helical_proj,moment_arm,'b');
axis equal
hold off;

view([-52 42])
drawnow 