function printDuplicates

% Prints the name and location of duplicate functions on MATLAB's current 
%   search path
%
%   Written by Andrew S Bock Jun 2016

%% Find and print duplicate functions
P=path;
P=strsplit(P, pathsep());
P=cellfun(@(x) what(x),P,'UniformOutput',false);
P=vertcat(P{:});
Q=arrayfun(@(x) x.m,P,'UniformOutput',false); % Q is a cell of cells of strings
Q=vertcat(Q{:});
R=arrayfun(@(x) repmat({x.path},size(x.m)),P,'UniformOutput',false); % R is a cell of cell of strings
R=vertcat(R{:});
[C]=unique(Q);
for c=1:numel(C)
    ind=strcmpi(C{c},Q);
   if sum(ind)>1
       fprintf('duplicate %s at paths\n\t',C{c});
       fprintf('%s\n\t',R{ind});
       fprintf('\n');
   end
end
