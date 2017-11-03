%%
%  load('/home/lauren/Desktop/MotionDataAus_April2017/energy_datastruct.mat')
 %%
 
close all
ind_s = 1; ind_f = 1; 
for i = 1:108
    
     if strcmp(data_struct(i).trials(1:4),'stat')
        [val, peak_locs] = findpeaks(-data_struct(i).force{1}(3,:));
        ind_big = find(peak_locs > 50);
        max_ind_slow(ind_s) = peak_locs(ind_big(1));% find the index at the end of the increased loading
        max_mag_slow(ind_s) = val(ind_big(1)); %get the force 
        
        % find the end of the impulse peak by flipping the dataset and
        % looking at the beginning
        
         [val, peak_locs] = findpeaks(-fliplr(data_struct(i).force{1}(3,:)));
         mmax = max(-data_struct(i).force{1}(3,:));
        ind_big = find(val > (mmax-0.05));
        max_ind_slow_end(ind_s) = length(data_struct(i).force{1}(3,:)) - peak_locs(ind_big(1))+1;
        max_mag_slow_end(ind_s) = val(ind_big(1));
        
        
        figure;
        plot(-data_struct(i).force{1}(3,:)); hold on;
        plot(max_ind_slow_end(ind_s), max_mag_slow_end(ind_s) ,'ro')
%         
        ind_s = ind_s + 1;
     else strcmp(data_struct(i).trials(1:3),'dyn')
         [val, peak_locs] = max(-data_struct(i).force{1}(3,:));
         max_ind_fast(ind_f) = peak_locs(1);
         max_mag_fast(ind_f) = val(1);
         
         %find the other end
         n = length(data_struct(i).force{1}(3,:));
         [val, peak_locs] = findpeaks(-fliplr(data_struct(i).force{1}(3,:)));
         mmax = (-data_struct(i).force{1}(3,round(n/2)));
        ind_big = find(val > (mmax-0.15));
        max_ind_fast_end(ind_f) =  n - peak_locs(ind_big(1))+1;
        max_mag_fast_end(ind_f) = val(ind_big(1));
        
%         figure;
%         plot(-data_struct(i).force{1}(3,:)); hold on;
%          plot(max_ind_fast_end(ind_f), max_mag_fast_end(ind_f) ,'ro')
         ind_f = ind_f + 1;
%         
     end
end
%%
% get the index transformed to the time domain
slow_t2pk = (max_ind_slow)*1/185;
fast_t2pk = (max_ind_fast)*1/185 ;  
slow_t2pk_end = (max_ind_slow_end)*1/185;
fast_t2pk_end = (max_ind_fast_end)*1/185;  

% length of the held compresion
slow_length = slow_t2pk_end - slow_t2pk;
fast_length = fast_t2pk_end - fast_t2pk;

mean(slow_length)
std(slow_length)

mean(fast_length)
std(fast_length)

% average time to peak
slow_t2pk_ave = mean(slow_t2pk);
fast_t2pk_ave = mean(fast_t2pk);

% [a b] = ttest2(max_val_slow,max_val_fast)

loading_rate_slow = max_mag_slow./slow_t2pk;
loading_rate_fast = max_mag_fast./fast_t2pk;


loading_rate_slow_end = max_mag_slow./slow_t2pk_end;
loading_rate_fast_end = max_mag_fast./fast_t2pk_end;

ind_1BW = logical(repmat([0 0 0 1 1 1],1,9));

sim_fz_slow = 0.09*loading_rate_slow(ind_1BW)+0.520;
sim_fz_fast = 0.09*loading_rate_fast(ind_1BW)+0.520;

diff_fz_slow = max_mag_slow(ind_1BW)-(sim_fz_slow);
diff_fz_fast = max_mag_fast(ind_1BW)-(sim_fz_fast);
%%
close all
figure; hold on;
plot(max_mag_slow(ind_1BW));plot(sim_fz_slow);
ylabel('Force, slow trials, in z [BW]')
legend('Measured','Simulated')
ylim([0.6,1.5])
figure; hold on;
plot(max_mag_fast(ind_1BW));plot(sim_fz_fast);
ylabel('Force, fast trials, in z [BW]')
legend('Measured','Simulated')
ylim([0.6,1.5])
figure; hold on;
plot(diff_fz_slow)
text(17, 0.26,sprintf('Mean difference = %0.2f BW',mean(diff_fz_slow)))
ylabel('Force DIFFERENCE, slow trials, in z [BW]')
figure; hold on;
plot(diff_fz_fast)
text(17, -0.14,sprintf('Mean difference = %0.2f BW',mean(diff_fz_fast)))
ylabel('Force DIFFERENCE, fast trials, in z [BW]')

fem_time = [1.15,1.36,1.73,2.11,2.36,2.33,2.54,2.28]./[7.77,11.5,14.6,16.9,19.1,19.6,23.7,22.3]
male_time = [1.23,1.42,1.62,2.10,2.45,2.35,2.46]./[8.2,11,14.6,16,18.32,18.9,22.8]