% Kronecker Delta
% Richard Rieber
% January 26, 2010
% rrieber@gmail.com
%
% Purpose:  This function calculates the kronecker delta from two given
% subscripts, j and k. If j == k, it returns 1. If j != k, it returns 0
%
% Inputs:  o j - the first subscript
%          o k - the second subscript
%
% Outputs: o d - the kronecker delta, 1 if j == k and - if j != k
%

function [d] = KronD(j,k)

if nargin < 2
    error('Too few inputs.  See help KronD')
elseif nargin > 2
    error('Too many inputs.  See help KronD')
end

if j == k
    d = 1;
else
    d = 0;
end
