function ProgressBar_parfor(varargin)
%% Progress bar for parfor loops
%
%   Usage:
%
%   tstart = clock;
%   ProgressBar_parfor(tstart);
%   parfor i = 1:<number of loops>
%       some stuff
%       ProgressBar_parfor(tstart,i,<number of loops>);
%   end
%   ProgressBar_parfor(tstart); % Clean up files and display total loop time
%
%   Parfor Example:
%   poolobj = gcp; % start parpool
%   tstart = clock;
%   nloops = 100; % number of loops
%   ProgressBar_parfor(tstart);
%   parfor v = 1:nloops
%       <some stuff>
%       ProgressBar_parfor(tstart,v,nloops);
%   end
%   ProgressBar_parfor(tstart,'clean'); % Display total loop time
%   delete(poolobj); % close parpool
%
%   Written by Andrew S Bock Mar 2014
fileid = num2str(varargin{1});
fileid(isspace(fileid)) = [];
if nargin == 1;
    if exist(fileid,'file'); delete(fileid);end
    f = fopen(fileid,'a');
    fclose(f);
elseif nargin == 3;
    f = fopen(fileid,'a');
    fprintf(f, [num2str(1) '\n']);
    fclose(f);
    f = fopen(fileid, 'r');
    prog = fscanf(f, '%d');
    fclose(f);
    perc = 100*sum(prog(1:end))/varargin{3};
    tstart = varargin{1};
    if mod(perc,1) < 0.1 % every 10% (or so) give an update
        %disp(['Working on input: ' num2str(varargin{2}) ' out of ' num2str(varargin{3})]);
        disp(['Iterations completed: ' num2str(perc) '%']);
        tmp = clock;
        elapsed = etime(tmp,tstart);
        TL = ((100/perc*elapsed) - elapsed)/60; % Time left in minutes
        disp(['Time left: ' num2str(TL) ' minutes'])
        %msg = ['Time left: ' num2str(TL) ' minutes'];
        %reverseStr = repmat(sprintf('\b'), 1, length(msg));
        % fprintf([reverseStr, msg]);
    end
else
    tstart = varargin{1};
    tmp = clock;
    TT = datevec(datenum(tmp)-datenum(tstart));
    fprintf('\n\n\nLoop took: %dh %dm %ds\n\n\n',TT(4),TT(5),round(TT(6)));
    delete(fileid);
end