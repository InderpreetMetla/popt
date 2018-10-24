%******************************************************************%
% This code adds the required paths to the MATLAB path directory
% in order to run POPT
%
% Distributed under the MIT License
%
% Copyright (c) 2018 Inderpreet Metla
%******************************************************************%
directory_current = pwd;
if directory_current(end-4:end) ~= '\POPT'
    error('You must be in the ''<path>\POPT'' folder')
end
addpath(genpath(directory_current))
prompt = 'Would you like to save this to your path permanently? Y/N [Y]: ';
RESULT = input(prompt,'s');
if isempty(RESULT)
    RESULT = 'Y';
end

if RESULT == 'Y'
    savepath
    disp('The following folders have been added to the path permanantely:')
    disp('          <path>\POPT\examples')
    disp('          <path>\POPT\nlp')
    disp('          <path>\POPT\src')
else
    disp('The following folders have been added to the path for this session:')
    disp('          <path>\POPT\examples')
    disp('          <path>\POPT\nlp')
    disp('          <path>\POPT\src')
end