function create_FSF_script (subject_name, session_date, stimuli_folders, non_WA_ind, WA_ind)
% session_dir
% subject_name
% session_date
% stimuli_folders
% non_WA_ind
% WA_ind

session_dir = fullfile('/data/jag/MELA/MOUNT_SINAI', subject_name, session_date);
output_dir = fullfile('/data/jag/MELA/Matlab/gkaguirrelab/SC-7T/Analysis', subject_name, session_date, 'FSF_files');
bold_dirs = find_bold(session_dir);
fid = fopen(fullfile(output_dir,[subject_name '_' session_date '_' 'FSF_shell.txt']), 'wt' );
fprintf(fid, ['#### Analysis for ' subject_name ' ' session_dir '\n\n\n\n']);

%%%% NON WRAP AROUND

for ii = 1:length(non_WA_ind)
    current_series = bold_dirs{non_WA_ind(ii)};
    current_stimulus = stimuli_folders{non_WA_ind(ii)};
    
    fprintf (fid, 'cp Templ_6freq_AT.fsf Run_%d_raw.fsf\n\n\n', non_WA_ind(ii));
    fprintf (fid, 'sed -i -- ''s/Series_005_bold_1.6_P2_mb5_L_minus_M_A_run1/%s/g'' ./Run_%d_raw.fsf\n', current_series, non_WA_ind(ii));
    fprintf (fid, ' sed -i -- ''s/HERO_gka1-HCLV_Photo_7T-01/%s/g'' ./Run_%d_raw.fsf\n\n\n', current_stimulus, non_WA_ind(ii));
    fprintf (fid, 'cp Run_%d_raw.fsf Run_%d_5mm.fsf\n\n', non_WA_ind(ii), non_WA_ind(ii));
    fprintf (fid, 'sed -i -- ''s/wdrf.tf/s5.wdrf.tf/g'' ./Run_%d_5mm.fsf\n\n\n', non_WA_ind(ii));
    fprintf (fid, 'feat Run_%d_5mm.fsf\n',non_WA_ind(ii));
    fprintf (fid, 'feat Run_%d_raw.fsf\n\n\n\n\n', non_WA_ind(ii));
end

if WA_ind == 0
    return
else
    for ii = 1:length(WA_ind)
    current_series = bold_dirs{WA_ind(ii)};
    current_stimulus = stimuli_folders{WA_ind(ii)};
    
    fprintf( fid, 'cp Templ_6freq_AT_WA.fsf Run_%d_raw.fsf\n\n\n', WA_ind(ii));
    fprintf( fid, 'sed -i -- ''s/Series_005_bold_1.6_P2_mb5_L_minus_M_A_run1/%s/g'' ./Run_%d_raw.fsf\n', current_series, WA_ind(ii));
    fprintf (fid, ' sed -i -- ''s/HERO_gka1-HCLV_Photo_7T-01/%s/g'' ./Run_%d_raw.fsf\n\n\n', current_stimulus, WA_ind(ii));
    fprintf (fid, 'cp Run_%d_raw.fsf Run_%d_5mm.fsf\n\n', WA_ind(ii), WA_ind(ii));
    fprintf (fid, 'sed -i -- ''s/wdrf.tf/s5.wdrf.tf/g'' ./Run_%d_5mm.fsf\n\n\n', WA_ind(ii));
    fprintf (fid, 'feat Run_%d_5mm.fsf\n',WA_ind(ii));
    fprintf (fid, 'feat Run_%d_raw.fsf\n\n\n\n\n', WA_ind(ii));
    end
end
fclose(fid);





% cp Templ6freq_AT.fsf TTF_run09_raw.fsf
% 
% 
% 
% sed -i -- 's/Series_005_bold_1.6_P2_mb5_L_minus_M_A_run1/CURRENT_SERIES/g' ./TTF_run09_raw.fsf
% sed -i -- 's/HERO_gka1-HCLV_Photo_7T-01/CURRENT_STIMULI/g' ./TTF_run09_raw.fsf
% 
% 
% 
% cp TTF_run01_raw.fsf TTF_run09_5mm.fsf
% 
% sed -i -- 's/wdrf.tf/s5.wdrf.tf/g' ./TTF_run09_5mm.fsf
% 
% 
% feat TTF_run09_5mm.fsf
% feat TTF_run09_raw.fsf
