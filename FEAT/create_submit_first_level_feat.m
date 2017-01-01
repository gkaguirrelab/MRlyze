function create_submit_first_level_feat (session_dir,subj_name,numRuns)

% creates a simple script that calls feat for each run first level stat fsf file.

out_dir = fullfile (session_dir, 'shell_scripts');
FSF_dir = fullfile(session_dir,'first_level_stat');
if ~isdir(out_dir)
    mkdir(out_dir);
end
fname = fullfile(out_dir, [subj_name 'submit_first_level_feat.sh']);
fid = fopen(fname,'w');
fprintf(fid,'#!/bin/bash\n');
for rr = 1:numRuns
 fprintf(fid, ['feat ' FSF_dir '/Run_%02d_5mm.fsf\n'],rr);
 fprintf(fid, ['feat ' FSF_dir '/Run_%02d_raw.fsf\n'],rr);
end