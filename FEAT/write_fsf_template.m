function write_fsf_template (templateDir,templateName, fsf, EVs, Contrasts, Ftests)

% requires fsf struct (create_fsf_struct), EVs struct (create_EVs_struct),
% Contrasts struct and Ftests struct.

% June 2016 - written (GF)

%% Open fsf file 
% Open the .fsf file based on information in the structure
fid = fopen(fullfile(templateDir,templateName),'wt');

%% Part 1 - 'static' parameters
fsf_str = sprintf([ ...
       '\n# FEAT version number' ...
       '\nset fmri(version) %1.2f\n'], fsf.version);
 fprintf(fid,'%s', fsf_str);
 
 fsf_str = sprintf([ ...
       '\n# Are we in MELODIC?' ...
       '\nset fmri(inmelodic) %d\n',fsf.inmelodic]);
 fprintf(fid,'%s',fsf_str);
 
 fsf_str = sprintf([...
       '\n# Analysis level' ...
       '\n# 1 : First-level analysis' ...
       '\n# 2 : Higher-level analysis' ...
       '\nset fmri(level) %d\n'], fsf.level);
 fprintf(fid,'%s', fsf_str);
 
 fsf_str = sprintf([...
       '\n# Which stages to run' ...
       '\n# 0 : No first-level analysis (registration and/or group stats only)' ...
       '\n# 7 : Full first-level analysis' ...
       '\n# 1 : Pre-processing' ...
       '\n# 2 : Statistics' ...
       '\nset fmri(analysis) %d\n'], fsf.analysis);
 fprintf(fid,'%s', fsf_str);
 
 fsf_str = sprintf([...
       '\n# Use relative filenames' ...
       '\nset fmri(relative_yn) %d\n'], fsf.relative_yn);
 fprintf(fid,'%s', fsf_str);
 
 fsf_str = sprintf([...
       '\n# Balloon help' ...
       '\nset fmri(help_yn) %d\n'], fsf.help_yn);
 fprintf(fid,'%s', fsf_str);
 
 fsf_str = sprintf([...
       '\n# Run Featwatcher' ...
       '\nset fmri(featwatcher_yn) %d\n'], fsf.featwatcher_yn);
 fprintf(fid,'%s', fsf_str);
 
 fsf_str = sprintf([...
       '\n# Cleanup first-level standard-space images' ...
       '\nset fmri(sscleanup_yn) %d\n'], fsf.sscleanup_yn);
 fprintf(fid,'%s', fsf_str);
 
 fsf_str = sprintf([...
       '\n# Output directory' ...
       '\nset fmri(outputdir) "%s"\n'], fsf.outputdir);
 fprintf(fid,'%s', fsf_str);
 
 fsf_str = sprintf([...
       '\n# TR(s)' ...
       '\nset fmri(tr) %s\n'], fsf.tr);
 fprintf(fid,'%s', fsf_str);
 
 fsf_str = sprintf([...
       '\n# Total volumes' ...
       '\nset fmri(npts) %s\n'], fsf.npts);
 fprintf(fid,'%s', fsf_str);
 
 fsf_str = sprintf([...
       '\n# Delete volumes' ...
       '\nset fmri(ndelete) %d\n'], fsf.ndelete);
 fprintf(fid,'%s', fsf_str);
 
   fsf_str = sprintf([...
       '\n# Perfusion tag/control order' ...
       '\nset fmri(tagfirst) %d\n'], fsf.tagfirst);
   fprintf(fid,'%s', fsf_str);

   
 fsf_str = sprintf([...
       '\n# Number of first-level analyses' ...
       '\nset fmri(multiple) %d\n'], fsf.multiple);
 fprintf(fid,'%s', fsf_str);
 
 fsf_str = sprintf([...
       '\n# Higher-level input type' ...
       '\n# 1 : Inputs are lower-level FEAT directories' ...
       '\n# 2 : Inputs are cope images from FEAT directories' ...
       '\nset fmri(inputtype) %d\n'], fsf.inputtype);
 fprintf(fid,'%s', fsf_str);
 
 fsf_str = sprintf([...
       '\n# Carry out pre-stats processing?' ...
       '\nset fmri(filtering_yn) %d\n'], fsf.filtering_yn);
 fprintf(fid,'%s', fsf_str);
 
 fsf_str = sprintf([...
       '\n# Brain/background threshold, %%' ...
       '\nset fmri(brain_thresh) %d\n'], fsf.brain_thresh);
 fprintf(fid,'%s', fsf_str);
 
 fsf_str = sprintf([...
     '\n# Critical z for design efficiency calculation' ...
     '\nset fmri(critical_z) %1.2f\n'], fsf.critical_z);
 fprintf(fid,'%s', fsf_str);
 
 fsf_str = sprintf([...
     '\n# Noise level' ...
     '\nset fmri(noise) %1.6f\n'], fsf.noise);
 fprintf(fid,'%s', fsf_str);
 
 fsf_str = sprintf([...
     '\n# Noise AR(1)' ...
     '\nset fmri(noisear) %1.6f\n'], fsf.noisear);
 fprintf(fid,'%s', fsf_str);
 
 fsf_str = sprintf([...
       '\n# Motion correction' ...
       '\n# 0 : None' ...
       '\n# 1 : MCFLIRT' ...
       '\nset fmri(mc) %d\n'], fsf.mc);
 fprintf(fid,'%s', fsf_str);
 
 fsf_str = sprintf([...
       '\n# Spin-history (currently obsolete)' ...
       '\nset fmri(sh_yn) %d\n'], fsf.sh_yn);
 fprintf(fid,'%s', fsf_str);
 
 fsf_str = sprintf([... 
       '\n# B0 fieldmap unwarping?' ...
       '\nset fmri(regunwarp_yn) %d\n'], fsf.regunwarp_yn);
 fprintf(fid,'%s', fsf_str);
 
 fsf_str = sprintf([...
     '\n# EPI dwell time (ms)' ...
     '\nset fmri(dwell) %d\n'], fsf.dwell);
 fprintf(fid,'%s', fsf_str);
 
 fsf_str = sprintf([...
     '\n# EPI TE (ms)' ...
     '\nset fmri(te) %d\n'], fsf.te);
 fprintf(fid,'%s', fsf_str);
 
 fsf_str = sprintf([...
     '\n# % Signal loss threshold' ...
     '\nset fmri(signallossthresh) %d\n'], fsf.signallossthresh);
 fprintf(fid,'%s', fsf_str);
 
 fsf_str = sprintf([...
     '\n# Unwarp direction' ...
     '\nset fmri(unwarp_dir) %s\n'], fsf.unwarp_dir);
 fprintf(fid,'%s', fsf_str);
 
 fsf_str = sprintf([...
     '\n# Slice timing correction' ...
     '\n# 0 : None' ...
     '\n# 1 : Regular up (0, 1, 2, 3, ...)' ...
     '\n# 2 : Regular down' ...
     '\n# 3 : Use slice order file' ...
     '\n# 4 : Use slice timings file' ...
     '\n# 5 : Interleaved (0, 2, 4 ... 1, 3, 5 ... )' ...
     '\nset fmri(st) %d\n'], fsf.st);
 fprintf(fid,'%s', fsf_str);
 
 fsf_str = sprintf([...
       '\n# Slice timings file' ...
       '\nset fmri(st_file) "%s"\n'], fsf.st_file);
 fprintf(fid,'%s', fsf_str);
 
 fsf_str = sprintf([...
       '\n# BET brain extraction' ...
       '\nset fmri(bet_yn) %d\n'], fsf.bet_yn);
 fprintf(fid,'%s', fsf_str);
 
 fsf_str = sprintf([...
       '\n# Spatial smoothing FWHM (mm)' ...
       '\nset fmri(smooth) %d\n'], fsf.smooth);
 fprintf(fid,'%s', fsf_str);
 
 fsf_str = sprintf([...
       '\n# Intensity normalization' ...
       '\nset fmri(norm_yn) %d\n'], fsf.norm_yn);
 fprintf(fid,'%s', fsf_str);
 
 fsf_str = sprintf([...
       '\n# Perfusion subtraction' ...
       '\nset fmri(perfsub_yn) %d\n'], fsf.perfsub_yn);
   fprintf(fid,'%s', fsf_str);
   
   fsf_str = sprintf([...
       '\n# Highpass temporal filtering' ...
       '\nset fmri(temphp_yn) %d\n'], fsf.temphp_yn);
 fprintf(fid,'%s', fsf_str);
 
 fsf_str = sprintf([...
       '\n# Lowpass temporal filtering' ...
       '\nset fmri(templp_yn) %d\n'], fsf.templp_yn);
 fprintf(fid,'%s', fsf_str);
 
 fsf_str = sprintf([...
       '\n# MELODIC ICA data exploration' ...
       '\nset fmri(melodic_yn) %d\n'], fsf.melodic_yn);
 fprintf(fid,'%s', fsf_str);
 
 fsf_str = sprintf([...
       '\n# Carry out main stats?' ...
       '\nset fmri(stats_yn) %d\n'], fsf.stats_yn);
 fprintf(fid,'%s', fsf_str);
 
 fsf_str = sprintf([...
       '\n# Carry out prewhitening?' ...
       '\nset fmri(prewhiten_yn) %d\n'], fsf.prewhiten_yn);
 fprintf(fid,'%s', fsf_str);
 
 
fsf_str = sprintf([...
       '\n# Add motion parameters to model' ...
       '\n# 0 : No' ...
       '\n# 1 : Yes' ...
       '\nset fmri(motionevs) %d\n' ...
       '\nset fmri(motionevsbeta) "%s"\n' ...
       '\nset fmri(scriptevsbeta) "%s"\n'...
       ], fsf.motionevs, fsf.motionevbeta, fsf.scriptevsbeta);
   fprintf(fid,'%s', fsf_str);
 
   fsf_str = sprintf([...
       '\n# Robust outlier detection in FLAME?' ...
       '\nset fmri(robust_yn) %d\n'], fsf.robust_yn);
   fprintf(fid,'%s', fsf_str);
   
   fsf_str = sprintf([...
       '\n# Higher-level modelling' ...
       '\n# 3 : Fixed effects' ...
       '\n# 0 : Mixed Effects: Simple OLS' ...
       '\n# 2 : Mixed Effects: FLAME 1' ...
       '\n# 1 : Mixed Effects: FLAME 1+2' ...
       '\nset fmri(mixed_yn) %d\n'], fsf.mixed_yn);
 fprintf(fid,'%s', fsf_str);
 
 fsf_str = sprintf([...
       '\n# Number of EVs' ...
       '\nset fmri(evs_orig) %d' ...
       '\nset fmri(evs_real) %d\n' ...
       '\nset fmri(evs_vox) %d\n'], fsf.evs_orig,fsf.evs_real, fsf.evs_vox);
 fprintf(fid,'%s', fsf_str);
 
 fsf_str = sprintf([...
       '\n# Number of contrasts' ...
       '\nset fmri(ncon_orig) %d' ...
       '\nset fmri(ncon_real) %d\n'], fsf.ncon_orig, fsf.ncon_real);
 fprintf(fid,'%s', fsf_str);
 
 fsf_str = sprintf([...
       '\n# Number of F-tests' ...
       '\nset fmri(nftests_orig) %d' ...
       '\nset fmri(nftests_real) %d\n'], fsf.nftests_orig, fsf.nftests_real);
 fprintf(fid,'%s', fsf_str);
 
 fsf_str = sprintf([...
       '\n# Add constant column to design matrix? (obsolete)' ...
       '\nset fmri(constcol) %d\n'], fsf.constcol);
 fprintf(fid,'%s', fsf_str);
 
 fsf_str = sprintf([...
       '\n# Carry out post-stats steps?' ...
       '\nset fmri(poststats_yn) %d\n'], fsf.poststats_yn);
 fprintf(fid,'%s', fsf_str);
 
 fsf_str = sprintf([...
       '\n# Pre-threshold masking?' ...
       '\nset fmri(threshmask) "%s"\n'], fsf.threshmask);
 fprintf(fid,'%s', fsf_str);
 
 fsf_str = sprintf([...
       '\n# Thresholding' ...
       '\n# 0 : None' ...
       '\n# 1 : Uncorrected' ...
       '\n# 2 : Voxel' ...
       '\n# 3 : Cluster' ...
       '\nset fmri(thresh) %d\n'], fsf.thresh);
 fprintf(fid,'%s', fsf_str);
 
 fsf_str = sprintf([...
       '\n# P threshold' ...
       '\nset fmri(prob_thresh) %1.4f\n'], fsf.prob_thresh);
 fprintf(fid,'%s', fsf_str);
 
 fsf_str = sprintf([...
       '\n# Z threshold' ...
       '\nset fmri(z_thresh) %1.2f\n'], fsf.z_thresh);
 fprintf(fid,'%s', fsf_str);
 
 fsf_str = sprintf([...
       '\n# Z min/max for colour rendering' ...
       '\n# 0 : Use actual Z min/max' ...
       '\n# 1 : Use preset Z min/max' ...
       '\nset fmri(zdisplay) %d\n'], fsf.zdisplay);
 fprintf(fid,'%s', fsf_str);
 
 fsf_str = sprintf([...
       '\n# Z min in colour rendering' ...
       '\nset fmri(zmin) %d\n'], fsf.zmin);
 fprintf(fid,'%s', fsf_str);
 
 fsf_str = sprintf([...
       '\n# Z max in colour rendering' ...
       '\nset fmri(zmax) %d\n'], fsf.zmax);
 fprintf(fid,'%s', fsf_str);
 
 fsf_str = sprintf([...
       '\n# Colour rendering type' ...
       '\n# 0 : Solid blobs' ...
       '\n# 1 : Transparent blobs' ...
       '\nset fmri(rendertype) %d\n'], fsf.rendertype);
 fprintf(fid,'%s', fsf_str);
 
 fsf_str = sprintf([...
       '\n# Background image for higher-level stats overlays' ...
       '\n# 1 : Mean highres' ...
       '\n# 2 : First highres' ...
       '\n# 3 : Mean functional' ...
       '\n# 4 : First functional' ...
       '\n# 5 : Standard space template' ...
       '\nset fmri(bgimage) %d\n'], fsf.bgimage);
 fprintf(fid,'%s', fsf_str);
 
 fsf_str = sprintf([...
       '\n# Create time series plots' ...
       '\nset fmri(tsplot_yn) %d\n'], fsf.tsplot_yn);
   fprintf(fid,'%s', fsf_str);
   
   fsf_str = sprintf([...
       '\n# Registration to initial structural' ...
       '\nset fmri(reginitial_highres_yn) %d\n'], fsf.reginitial_highres_yn);
 fprintf(fid,'%s', fsf_str);
 
 fsf_str = sprintf([...
       '\n# Search space for registration to initial structural' ...
       '\n# 0   : No search' ...
       '\n# 90  : Normal search' ...
       '\n# 180 : Full search' ...
       '\nset fmri(reginitial_highres_search) %d\n'], fsf.reginitial_highres_search);
 fprintf(fid,'%s', fsf_str);
 
 fsf_str = sprintf([...
       '\n# Degrees of Freedom for registration to initial structural' ...
       '\nset fmri(reginitial_highres_dof) %d\n'], fsf.reginitial_highres_dof);
 fprintf(fid,'%s', fsf_str);
 
 fsf_str = sprintf([...
       '\n# Registration to main structural' ...
       '\nset fmri(reghighres_yn) %d\n'], fsf.reghighres_yn);
 fprintf(fid,'%s', fsf_str);
 
 fsf_str = sprintf([...
       '\n# Search space for registration to main structural' ...
       '\n# 0   : No search' ...
       '\n# 90  : Normal search' ...
       '\n# 180 : Full search' ...
       '\nset fmri(reghighres_search) %d\n'], fsf.reghighres_search);
   fprintf(fid,'%s', fsf_str);
   
   fsf_str = sprintf([...
       '\n# Degrees of Freedom for registration to main structural' ...
       '\nset fmri(reghighres_dof) %d\n'], fsf.reghighres_dof);
   fprintf(fid,'%s', fsf_str);
   
   fsf_str = sprintf([...
       '\n# Registration to standard image?' ...
       '\nset fmri(regstandard_yn) %d\n'], fsf.regstandard_yn);
   fprintf(fid,'%s', fsf_str);
   
   fsf_str = sprintf([...
       '\n# Use alternate reference images?' ...
       '\nset fmri(alternateReference_yn) %d\n'], fsf.alternateReference_yn);
   fprintf(fid,'%s', fsf_str);
   
   fsf_str = sprintf([...
       '\n# Standard image' ...
       '\nset fmri(regstandard) "%s"\n'], fsf.regstandard);
   fprintf(fid,'%s', fsf_str);
   
   fsf_str = sprintf([...
       '\n# Search space for registration to standard space' ...
       '\n# 0   : No search' ...
       '\n# 90  : Normal search' ...
       '\n# 180 : Full search' ...
       '\nset fmri(regstandard_search) %d\n'], fsf.regstandard_search);
   fprintf(fid,'%s', fsf_str);
   
   fsf_str = sprintf([...
       '\n# Degrees of Freedom for registration to standard space' ...
       '\nset fmri(regstandard_dof) %d\n'], fsf.regstandard_dof);
   fprintf(fid,'%s', fsf_str);
   
   fsf_str = sprintf([...
       '\n# Do nonlinear registration to standard space?' ...
       '\nset fmri(regstandard_nonlinear_yn) %d\n'], fsf.regstandard_nonlinear_yn);
   fprintf(fid,'%s', fsf_str);
   
   fsf_str = sprintf([...
       '\n# Control nonlinear warp field resolution' ...
       '\nset fmri(regstandard_nonlinear_warpres) %d\n'],...
       fsf.regstandard_nonlinear_warpres);
   fprintf(fid,'%s', fsf_str);
   
   fsf_str = sprintf([...
       '\n# High pass filter cutoff' ...
       '\nset fmri(paradigm_hp) %d\n'], fsf.paradigm_hp);
   fprintf(fid,'%s', fsf_str);
   
   fsf_str = sprintf([...
       '\n# Total voxels' ...
       '\nset fmri(totalVoxels) %s\n'], fsf.totalVoxels);
   fprintf(fid,'%s', fsf_str);
   
   fsf_str = sprintf([...
       '\n# Number of lower-level copes feeding into higher-level analysis' ...
       '\nset fmri(ncopeinputs) %d\n'], fsf.ncopeinputs);
   fprintf(fid,'%s', fsf_str);
   
   nfiles = length({fsf.feat_files});
   for ifile = 1:nfiles
       fsf_str = sprintf([...
           '\n# 4D AVW data or FEAT directory (%d)' ...
           '\nset feat_files(%d) "%s"\n'], ifile, ifile, fsf.feat_files);
       fprintf(fid,'%s', fsf_str);
   end
   
   fsf_str = sprintf([...
       '\n# Add confound EVs text file' ...
       '\nset fmri(confoundevs) %d\n'], fsf.confoundevs);
   fprintf(fid,'%s', fsf_str);
   
   
   %% Part 2 -  EVs
   
   fsf_str = sprintf([...
       '\n'...
       '\n'...
       '\n']);
   fprintf(fid,'%s', fsf_str);
   
   for ct = 1:fsf.evs_orig
       
       fsf_str = sprintf([...
           '# EV %d title\n'...
           'set fmri(evtitle%d) "%s"\n'...
           '\n'...
           ], ct,ct,EVs.title{ct});
       fprintf(fid,'%s', fsf_str);
       
       fsf_str = sprintf([...
           '# Basic waveform shape (EV %d)\n' ...
           '# 0 : Square\n' ...
           '# 1 : Sinusoid\n' ...
           '# 2 : Custom (1 entry per volume)\n' ...
           '# 3 : Custom (3 column format)\n' ...
           '# 4 : Interaction\n' ...
           '# 10 : Empty (all zeros)\n' ...
           'set fmri(shape%d) %d\n'...
           '\n'...
           ], ct,ct,EVs.shape(ct));
       fprintf(fid,'%s', fsf_str);
       
       fsf_str = sprintf([...
           '# Convolution (EV %d)\n' ...
           '# 0 : None\n' ...
           '# 1 : Gaussian\n' ...
           '# 2 : Gamma\n' ...
           '# 3 : Double-Gamma HRF\n' ...
           '# 4 : Gamma basis functions\n' ...
           '# 5 : Sine basis functions\n' ...
           '# 6 : FIR basis functions\n' ...
           'set fmri(convolve%d) %d\n' ...
           '\n' ...
           ], ct,ct,EVs.convolve(ct));
       fprintf(fid,'%s', fsf_str);
       
       fsf_str = sprintf([...
           '# Convolve phase (EV %d)\n' ...
           'set fmri(convolve_phase%d) %d\n' ...
           '\n' ...
           ], ct,ct,EVs.convolve_phase(ct));
       fprintf(fid,'%s', fsf_str);
       
       fsf_str = sprintf([...
           '# Apply temporal filtering (EV %d)\n' ...
           'set fmri(tempfilt_yn%d) %d\n' ...
           '\n' ...
           ], ct,ct,EVs.tempfilt_yn(ct));
       fprintf(fid,'%s', fsf_str);
       
       fsf_str = sprintf([...
           '# Add temporal derivative (EV %d)\n' ...
           'set fmri(deriv_yn%d) %d\n' ...
           '\n' ...
           ], ct,ct,EVs.deriv_yn(ct));
       fprintf(fid,'%s', fsf_str);
       
       fsf_str = sprintf([...
           '# Custom EV file (EV %d)\n' ...
           'set fmri(custom%d) "%s"\n' ...
           '\n' ...
           ], ct,ct,EVs.custom{ct});
       fprintf(fid,'%s', fsf_str);
       
       fsf_str = sprintf([...
           '# Orthogonalise EV %d wrt EV 0\n' ...
           'set fmri(ortho%d.0) %d\n' ...
           '\n' ...
           ], ct,ct,EVs.ortho(ct,1));
       fprintf(fid,'%s', fsf_str);
       
       for nn = 1:(fsf.evs_real)
           fsf_str = sprintf([...
               '# Orthogonalise EV %d wrt EV %d\n' ...
               'set fmri(ortho%d.%d) %d\n' ...
               '\n' ...
               ], ct,nn,ct,nn, EVs.ortho(ct,nn+1));
           fprintf(fid,'%s', fsf_str);
           
       end
   end

   %% Part 3 - Contrasts and F tests
 %%%%% Here starts the part with contrasts and F-tests on the text file.
 
 fsf_str = sprintf([...
       '# Contrast & F-tests mode\n'...
    '# real : control real EVs\n'...
    '# orig : control original EVs\n'...
    'set fmri(con_mode_old) %s\n'...
    'set fmri(con_mode) %s\n' ...
    '\n'], ...
    Contrasts.con_mode_old, Contrasts.con_mode);
   fprintf(fid,'%s', fsf_str);
   
    % portion for 'real EVs'
   for ct = 1:fsf.ncon_real
       fsf_str = sprintf([...
           '# Display images for contrast_real %d\n' ...
           'set fmri(conpic_real.%d) %d\n' ...
           '\n'] , ...
           ct, ct, Contrasts.conpic_real(ct));
       fprintf(fid,'%s', fsf_str);
       
       fsf_str = sprintf([...
           '# Title for contrast_real %d\n' ...
           'set fmri(conname_real.%d) %s\n' ...
           '\n'] , ...
           ct, ct, Contrasts.conname_real{ct});
       fprintf(fid,'%s', fsf_str);
       
       for nn = 1:fsf.ncon_real
           fsf_str = sprintf([...
               '# Real contrast_real vector %d element %d\n' ...
               'set fmri(con_real%d.%d) %d\n' ...
               '\n'] , ...
               ct, nn, ct, nn, Contrasts.con_real(ct,nn));
           fprintf(fid,'%s', fsf_str);
       end
       
       for ff = 1:fsf.nftests_real
           fsf_str = sprintf([...
               '# F-test %d element %d\n' ...
               'set fmri(ftest_real%d.%d) %d\n' ...
               '\n'] , ...
               ff, ct, ff, ct, Ftests.ftest_real(ff,ct));
           fprintf(fid,'%s', fsf_str);
       end
   end
   
   
   % portion for 'original EVs'
   for ct = 1:fsf.ncon_real
       fsf_str = sprintf([...
           '# Display images for contrast_orig %d\n' ...
           'set fmri(conpic_orig.%d) %d\n' ...
           '\n'] , ...
           ct, ct, Contrasts.conpic_orig(ct));
       fprintf(fid,'%s', fsf_str);
       
       fsf_str = sprintf([...
           '# Title for contrast_orig %d\n' ...
           'set fmri(conname_orig.%d) %s\n' ...
           '\n'] , ...
           ct, ct, Contrasts.conname_orig{ct});
       fprintf(fid,'%s', fsf_str);
       
       for nn = 1:fsf.ncon_orig
           fsf_str = sprintf([...
               '# Real contrast_orig vector %d element %d\n' ...
               'set fmri(con_orig%d.%d) %d\n' ...
               '\n'] , ...
               ct, nn, ct, nn, Contrasts.con_orig(ct,nn));
           fprintf(fid,'%s', fsf_str);
       end
       
       for ff = 1:fsf.nftests_orig
           fsf_str = sprintf([...
               '# F-test %d element %d\n' ...
               'set fmri(ftest_orig%d.%d) %d\n' ...
               '\n'] , ...
               ff, ct, ff, ct, Ftests.ftest_orig(ff,ct));
           fprintf(fid,'%s', fsf_str);
       end
   end
   
   % portion for contrast masking
   fsf_str = sprintf([...
       '# Contrast masking - use >0 instead of thresholding?\n' ...
       'set fmri(conmask_zerothresh_yn) %d\n'...
       '\n'] , ...
       Contrasts.zerothresh_yn);
   fprintf(fid,'%s', fsf_str);
   
   for mm = 1:(fsf.ncon_real +fsf.nftests_real)
       for kk = 2:(fsf.ncon_real +fsf.nftests_real)
           fsf_str = sprintf([...
               '# Mask real contrast/F-test %d with real contrast/F-test %d?\n'...
               'set fmri(conmask%d_%d) %d\n'...
               '\n'] , ...
               mm, kk, mm, kk, Contrasts.conmask(mm,kk));
           fprintf(fid,'%s', fsf_str);
       end
   end
   fsf_str = sprintf([...
       '# Do contrast masking at all?\n'...
       'set fmri(conmask1_1) %d\n'...
       '\n'] , ...
       Contrasts.conmask(1,1));
   fprintf(fid,'%s', fsf_str);
   
 %% Part 4 - Static parameters (not on GUI)
 
    fsf_str = sprintf([...
       '\n#############################################################' ...
       '\n# Now options that don''t appear in the GUI\n' ...
       '\n# Alternative (to BETting) mask image' ...
       '\nset fmri(alternative_mask) "%s"\n' ...
       '\n# Initial structural space registration initialisation transform' ...
       '\nset fmri(init_initial_highres) "%s"\n' ...
       '\n# Structural space registration initialisation transform' ...
       '\nset fmri(init_highres) "%s"\n' ...
       '\n# Standard space reigstration initialisation transform' ...
       '\nset fmri(init_standard) "%s"\n' ...
       '\n# For full FEAT analysis: overwrite .feat output dir?' ...
       '\nset fmri(overwrite_yn) %d\n'],...
       fsf.alternative_mask,...
       fsf.init_initial_highres,fsf.init_highres,...
       fsf.init_standard,fsf.overwrite_yn);
   fprintf(fid,'%s', fsf_str);
 
 
 %% close fsf file
 fclose(fid);
 return
 
 
 
 