function checks = TypeChecker()
% Return a structure with handles to functions that handle type checks
% and possibly range checks.

%   Copyright 2007-2011 The MathWorks, Inc.

    checks = struct('isStructWithFields', @isStructWithFields, ...
                    'isIntegerScalar', @isIntegerScalar, ...
                    'isRealScalar', @isRealScalar, ...
                    'isClusterObject', @isClusterObject, ...
                    'isIndependentJobObject', @isIndependentJobObject, ...
                    'isCommunicatingJobObject', @isCommunicatingJobObject);
end % End of pTypeChecker.

function valid = isStructWithFields(value, varargin)
% Verify that value is a struct and that it has the fieldnames specified  
% in varargin and only those fieldnames.
% Varargin must list at least one fieldname.
    wantedFields = varargin;
    valid = isstruct(value);
    if ~valid
        return;
    end
    actualFields = fieldnames(value);
    % Verify that the actualFields and wantedFields contain the same 
    % elements
    valid = isempty(setxor(actualFields, wantedFields));
end % End of isStructWithFields.

function valid = isIntegerScalar(value, lowerBound, upperBound)
%valid = isIntegerScalar(value) Return true if and only if value is a  
% finite, scalar integer in the specified range.
    valid = isreal(value) && isscalar(value) ...
            && (value >= lowerBound ) && (value <= upperBound) ...
            && (value  == floor(value)) && isfinite(value);
end % End of isIntegerScalar.


function valid = isRealScalar(value, lowerBound, upperBound)
%valid = isRealScalar(value) Return true if and only if value is a 
% finite, real scalar in the specified range.
    valid = isreal(value) && isscalar(value) && (value >= lowerBound) ...
            && (value <= upperBound) && isfinite(value);
end % End of isRealScalar.

function valid = isClusterObject(value)
%valid = isClusterObject(value) Return true if and only if value is a 
% single cluster object handle.
    valid = isscalar(value) && isa(value, 'parallel.Cluster');
end

function valid = isIndependentJobObject(value)
%valid = isIndependentJobObject(value) Return true if and only if value is a 
% single independent job object handle.
    valid = isscalar(value) && isa(value,'parallel.IndependentJob');
end

function valid = isCommunicatingJobObject(value)
%valid = isCommunicatingJobObject(value) Return true if and only if value is a 
% single communicating job object handle.
    valid = isscalar(value) && isa(value,'parallel.CommunicatingJob');
end
