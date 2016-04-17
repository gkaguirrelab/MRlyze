function project_2_decimate(in_name,out_name,sind)

% Projects a surface file to a decimated surface, within subject
%
%   Usage:
%   project_2_decimate(in_name,out_name)
%
%   Written by Andrew S Bock Aug 2015

%% Project in_val to decimated surface
[~,~,ext] = fileparts(in_name);
if strcmp(ext,'.gz')
    in = load_nifti(in_name);
    out = in;
    out.vol = nan(length(sind),1);
    out.dim(2) = length(sind);
    if length(size(in.vol)) == 4
        out.vol = in.vol(sind,1,1,:);
    else
        out.vol = in.vol(sind);
    end
    save_nifti(out,out_name);
elseif strcmp(ext,'.mgh') || strcmp(ext,'.mgz')
    [in.vol,M,p] = load_mgh(in_name);
    out.vol = nan(length(sind),1);
    out.vol = in.vol(sind);
    save_mgh(out.vol,out_name,M,p);
end