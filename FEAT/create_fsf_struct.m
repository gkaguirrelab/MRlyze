function create_fsf_struct

%% all fields needed in write fsf
% fsf.templateDir
% fsf.templateName
% 
% fsf.version
% fsf.inmelodic
% fsf.level
% fsf.analysis
% fsf.relative_yn
% fsf.help_yn
% fsf.featwatcher_yn
% fsf.sscleanup_yn
% fsf.outputdir
% fsf.tr
% fsf.npts
% fsf.ndelete
% fsf.multiple
% fsf.inputtype

%% create an empty fsf_struct (to be filled in manually)
fsf.version = 0;
fsf.inmelodic = 0;
fsf.level  = 0;
fsf.analysis = 0;
fsf.relative_yn = 0;
fsf.help_yn = 0;
fsf.featwatcher_yn = 0;
fsf.sscleanup_yn = 0;
fsf.outputdir = 0;
fsf.tr = 0;
fsf.npts = 0;
fsf.ndelete = 0;
fsf.tagfirst = 0;
fsf.multiple = 0;
fsf.inputtype = 0;
fsf.filtering_yn = 0;
fsf.brain_thresh = 0;
fsf.critical_z = 0;
fsf.noise = 0;
fsf.noisear = 0;
fsf.mc = 0;
fsf.sh_yn = 0;
fsf.regunwarp_yn = 0;
fsf.dwell = 0;
fsf.te = 0;
fsf.signallossthresh = 10;
fsf.unwarp_dir = 0;
fsf.st = 0;

%% create a precompiled fsf_struct (default values + fields to design template)

 % FEAT version number
 fsf.version = 6.00;
 
 %Are we in MELODIC?' ...
 fsf.inmelodic = 0;
 
 % Analysis level
 % 1 : First-level analysis
 % 2 : Higher-level analysis
 fsf.level = 1;
 
 % Which stages to run
 % 0 : No first-level analysis (registration and/or group stats only)
 % 7 : Full first-level analysis
 % 1 : Preprocessing
 % 2 : Statistics
 fsf.analysis = 7;
 
 % Use relative filenames
 fsf.relative_yn = 0;
 
 % Balloon help
 fsf.help_yn = 1;
 
 % Run Featwatcher
 fsf.featwatcher_yn = 1;
 
 % Cleanup first-level standard-space images
 fsf.sscleanup_yn = 0;
 
 % Output directory
 fsf.outputdir = ''; %%%%%%%%%%%%%%%%%
 
 % TR(s)
 fsf.tr = []; %%%%%%%%%%%%%%%%%%%%%
 
 % Total volumes
 fsf.npts = []; %%%%%%%%%%%%%%%%%%%%%%%%%
 
 % Delete volumes
 fsf.ndelete = 0; 
 
 %  Perfusion tag/control order
 fsf.tagfirst = 1;
 
 % Number of first-level analyses
 fsf.multiple = 1;  % Not a flag - actual num of analyses
 
 % Higher-level input type
 % 1 : Inputs are lower-level FEAT directories
 % 2 : Inputs are cope images from FEAT directories
 fsf.inputtype = 2;
 
 % Carry out pre-stats processing?
 fsf.filtering_yn = 1;
 
 % Brain/background threshold % 
 fsf.brain_thresh = 10;
 
 % Critical z for design efficiency calculation
 fsf.critical_z = 5.3;
 
 % Noise level
 fsf.noise = 0.66;
 
 % Noise AR(1)
 fsf.noisear = 0.34;
 
 % Motion correction
 % 0 : None
 % 1 : MCFLIRT
 fsf.mc = 0;
 
 % Spin-history (currently obsolete)
 fsf.sh_yn = 0;
 
 % B0 fieldmap unwarping?
 fsf.regunwarp_yn = 0;
 
 % EPI dwell time (ms)
 fsf.dwell = 0.7;
 
 % EPI TE (ms)
 fsf.te = 35;
 
 % Signal loss threshold
 fsf.signallossthresh = 10;
 
 % Unwarp direction
 fsf.unwarp_dir = 'y-';
 
 % Slice timing correction
 % 0 : None
 % 1 : Regular up (0, 1, 2, 3, ...)
 % 2 : Regular down
 % 3 : Use slice order file
 % 4 : Use slice timings file
 % 5 : Interleaved (0, 2, 4 ... 1, 3, 5 ... )
 fsf.st = 0;
  
 % Slice timings file
 fsf.st_file='';
 
 
 
 
 
 
 
 
 
 
 
 
 % Post-stats-only directory copying
 % 0 : Overwrite original post-stats results
 % 1 : Copy original FEAT directory for new Contrasts, Thresholding, Rendering
 fsf.newdir_yn=0;
 

 
 % Slice timings file
 fsf.st_file='';
 

 

 
 % BET brain extraction
 fsf.bet_yn=1;
 
 % Spatial smoothing FWHM (mm)
 fsf.smooth=5;
 
 % Intensity normalization
 fsf.norm_yn=0;
 
 % Highpass temporal filtering
 fsf.temphp_yn=0;
 
 % Lowpass temporal filtering
 fsf.templp_yn=0;
 
 % MELODIC ICA data exploration
 fsf.melodic_yn=0;
 
 % Carry out main stats?
 fsf.stats_yn=1;
 
 % Carry out prewhitening?
 fsf.prewhiten_yn=1;
 
 % Higher-level modelling
 % 3 : Fixed effects
 % 0 : Mixed Effects: Simple OLS
 % 2 : Mixed Effects: FLAME (stage 1 only)
 % 1 : Mixed Effects: FLAME (full)
 fsf.mixed_yn=1;
 
 % Number of EVs
 fsf.evs_orig=0;
 fsf.evs_real=0;
 
 % Number of contrasts
 fsf.ncon_orig=0;
 fsf.ncon_real=0;
 
 % Number of F-tests
 fsf.nftests_orig=0;
 fsf.nftests_real=0;
 
 % Add constant column to design matrix? (obsolete)
 fsf.constcol=0;
 
 % Carry out post-stats steps?
 fsf.poststats_yn=1;
 
 % Pre-threshold masking?
 fsf.threshmask='';
 
 % Thresholding
 % 0 : None
 % 1 : Uncorrected
 % 2 : Voxel
 % 3 : Cluster
 fsf.thresh=3;

 % P threshold
 fsf.prob_thresh=0.05;
 
% Z threshold
 fsf.z_thresh=2.3;

 % Z min/max for colour rendering
 % 0 : Use actual Z min/max
 % 1 : Use preset Z min/max
 fsf.zdisplay=0;
 
 % Z min in colour rendering
 fsf.zmin=1;
 
 % Z max in colour rendering
 fsf.zmax=15;
 
 % Colour rendering type
 % 0 : Solid blobs
 % 1 : Transparent blobs
 fsf.rendertype=1;
 
 % Background image for higher-level stats overlays
 % 1 : Mean highres
 % 2 : First highres
 % 3 : Mean functional
 % 4 : First functional
 % 5 : Standard space template
 fsf.bgimage=1;
 
 % Registration?
 fsf.reg_yn=1;
 

 
 % Dwell/Asymmetry ratio
 fsf.dwellasym=-1;
 
 % Registration to initial structural
 fsf.reginitial_highres_yn=0;
 
 % Search space for registration to initial structural
 % 0   : No search
 % 90  : Normal search
 % 180 : Full search
 fsf.reginitial_highres_search=90;
 
 % Degrees of Freedom for registration to initial structural
 fsf.reginitial_highres_dof=12;
 
 % Registration to main structural
 fsf.reghighres_yn=0;

 % Search space for registration to main structural
 % 0   : No search
 % 90  : Normal search
 % 180 : Full search
 fsf.reghighres_search=90;
 
 % Degrees of Freedom for registration to main structural
 fsf.reghighres_dof=12;
 
 % Registration to standard image?
 fsf.regstandard_yn=1;
 
 % Standard image
 fsf.regstandard='/usr/local/fsl/etc/standard/avg152T1_brain.hdr';
 
 % Search space for registration to standard space
 % 0   : No search
 % 90  : Normal search
 % 180 : Full search
 fsf.regstandard_search=90;
 
 % Degrees of Freedom for registration to standard space
 fsf.regstandard_dof=12;
 
 % Do nonlinear registration to standard space?
 fsf.regstandard_nonlinear_yn=0;
 
 % High pass filter cutoff
 fsf.paradigm_hp=-1;
 
 % Number of lower-level copes feeding into higher-level analysis
 fsf.ncopeinputs=0;
 
 % 4D AVW data or FEAT directory (1)
 fsf.feat_files={};