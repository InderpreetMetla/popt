function inputValidator(problem)
%******************************************************************%
% This function valiidates that the user has set up the problem 
% in the correct format and all the required fields have been 
% provided
% 
% Inputs:
%
% - problem struct
%
% Distributed under the MIT License
%
% Copyright (c) 2018 Inderpreet Metla
%******************************************************************%

% Check problem struct fields
if ~isfield(problem,'bounds')
    error('Error: bounds Field Has Not Been Provided in problem Struct')
elseif ~isfield(problem,'name')
    error('Error: problem.name Must Be Specified. Please name the problem.')
elseif ~isfield(problem,'funcs')
    error('Error: funcs Field Has Not Been Provided in problem Struct')
elseif ~isfield(problem,'guess')
    error('Error: guess Field Has Not Been Provided in problem Struct')  
end

nphases = length(problem.bounds.phase);
bounds  = problem.bounds.phase;
guess   = problem.guess.phase;

% Check if all phases have a guess
if length(guess)<nphases
    error(strcat('Error: The guess Struct Has Not Been',...
        ' Specified For Every Phase'));
end

for iphase = 1:nphases
    
    %********************************************************%
    % Check Bounds
    %********************************************************%
    
    %****************************%
    % States Bounds
    %****************************%
    if isfield(bounds(iphase),'state')
        % Call subfunction #1 to check bounds
        boundsinputCheck1(bounds(iphase).state,iphase,' States ');
    end
     
    if isfield(bounds(iphase),'initialstate')
        % Call subfunction #2 to check bounds
        boundsinputCheck2(bounds(iphase).initialstate,...
            bounds(iphase).state,iphase,' Initial State ');
    end
      
    if isfield(bounds(iphase),'finalstate')
        % Call subfunction #2 to check bounds
        boundsinputCheck2(bounds(iphase).finalstate,...
            bounds(iphase).state,iphase,' Final State ');
    end
    
    %****************************%
    % Control Bounds
    %****************************%
    if isfield(bounds(iphase),'control')
        % Call subfunction #1 to check bounds
        boundsinputCheck1(bounds(iphase).control,iphase,' Controls ');
    end
    
    %****************************%
    % Path Constaints Bounds
    %****************************%
    if isfield(bounds(iphase),'path')
        % Call subfunction #1 to check bounds
        boundsinputCheck1(bounds(iphase).path,iphase,' Path Constraints ');
    end
    
    %****************************%
    % Boundary Constaints Bounds
    %****************************%
    if isfield(bounds(iphase),'boundary')
        % Call subfunction #1 to check bounds
        boundsinputCheck1(bounds(iphase).boundary,...
            iphase,' Boundary Constraints ');
    end
    
    %*****************************************%
    % Check that upper bounds > lower bounds
    %*****************************************%
    
    % Check input on time bounds
    if ~isequal(size(bounds(iphase).initialtime.lb),[1, 1]) || ...
            ~isequal(size(bounds(iphase).initialtime.ub),[1, 1]) || ...
            ~isequal(size(bounds(iphase).finaltime.lb),[1, 1]) || ...
            ~isequal(size(bounds(iphase).finaltime.ub),[1, 1])
        error(strcat('Error: Upper and Lower Bounds on Time',...
            ' must be Scalars in phase %i.'),iphase);
    end
    
    % Check input on state bounds
    if (any(bounds(iphase).initialstate.ub < bounds(iphase).state.lb)...
       || any(bounds(iphase).finalstate.ub < bounds(iphase).state.lb)...
       || any(bounds(iphase).state.ub < bounds(iphase).initialstate.lb)...
       || any(bounds(iphase).state.ub < bounds(iphase).finalstate.lb)...
       || any(bounds(iphase).state.ub < bounds(iphase).state.lb))
        error(strcat('Error: Inconsistent Bounds on States ',...
            ' in phase %i (i.e. max < min)'),iphase);
    end
    % Check input on control bounds
    if isfield(bounds(iphase),'control')
        if any(bounds(iphase).control.ub < bounds(iphase).control.lb)
            error(strcat('Error: Inconsistent Bounds on Controls in',...
                ' phase %i (i.e. max < min)'),iphase);
        end
    end
    % Check input on path constraint bounds
    if isfield(bounds(iphase),'path')
        if any(bounds(iphase).path.ub < bounds(iphase).path.lb)
            error(strcat('Error: Inconsistent Bounds on Path Constraints',...
                ' in phase %i (i.e. max < min)'),iphase);
        end
    end
    % Check input on boundary constraint bounds
    if isfield(bounds(iphase),'boundary')
        if any(bounds(iphase).boundary.ub < bounds(iphase).boundary.lb)
            error(strcat('Error: Inconsistent Bounds on Boundary',...
                ' Constraints in phase %i (i.e. max < min)'),iphase);
        end
    end
    
    %****************************%
    % Linkage Constraint Bounds
    %****************************%
    if isfield(problem.bounds,'link')
        linkB   = problem.bounds.link;
        nlinks  = length(linkB);
        for ilink=1:nlinks
            % Check inputs on link constraint bounds
            if ~isfield(linkB(ilink),'lb') || ~isfield(linkB(ilink),'ub')
                error(strcat('Error: Both Upper and Lower Bounds Need',...
                    ' to be Specified for State Linkage Constraints ',...
                    ' across Link %i'),ilink);
            else
                if ~isequal(size(linkB(ilink).lb),size(linkB(ilink).ub))
                    error(strcat('Error: Upper and Lower Linkage Bounds ',...
                        ' Need to be the Same Size across',...
                        ' Link %i'),ilink);
                end
                if ~isequal(size(linkB(ilink).lb,1), 1)
                    error(strcat('Error: Upper and Lower Linkage Bounds ',...
                        ' across Link %i Should be Single',...
                        ' Row Vectors'),ilink);
                end
                if ~isequal(size(linkB(ilink).lb,2)-1,...
                        size(bounds(1).state.lb,2))
                    error(strcat('Error: Upper and Lower Linkage Bounds ',...
                        ' across Linkage %i Should be Have',...
                        ' %i Columns (as there are %i States',...
                        ' across Phases %i and %i Plus Time)'),ilink,...
                        size(bounds(1).state.lb,2)+1,...
                        size(bounds(1).state.lb,2),...
                        ilink,ilink+1);
                end
                if any(linkB(ilink).lb - linkB(ilink).ub > 0)
                    error(strcat('Error: Inconsistent Linkage Bounds',...
                        ' i.e. (ub < lb) across Link %i'),ilink);
                end
                if isempty(linkB(ilink).lb) || isempty(linkB(ilink).ub)
                    error(strcat('Error: Linkage Bounds are Empty',...
                        ' Across Link %i'),ilink);
                end
            end
        end
    end
    
end

%******************************%
% More Linkage Constaint Checks
%******************************%
if isfield(problem.bounds,'link')
    linkB   = problem.bounds.link;
    nlinks  = length(linkB);
    if ~isequal(nlinks,nphases-1)
        error(strcat('Error: %i Phases Have Been Specified ',...
            ' and Bounds Have Been Provided for %i Linkages.',...
            ' There should be %i Linkages',...
            ' for this Problem.'),nphases, nlinks, nphases-1);
    end
end

%********************************************************%
% Check Guess
%********************************************************%
if ~isfield(guess(iphase),'time')
    error(strcat('Error: Field ''time'' has been ommitted',...
        ' in iphase = %i from ''problem.guess.phase(iphase)'''),iphase);
end

if ~isfield(guess(iphase),'state')
    error(strcat('Error: Field ''state'' has been ommitted',...
        ' in iphase = %i from ''problem.guess.phase(iphase)'''),iphase);
end

end

function boundsinputCheck1(Struct,iphase,name)
%******************************************************************%
% Validates input for the state, control and constraints
% Does NOT validate for initial and final state
%
% Inputs:
%
% - Struct:         The actual struct that is being error validated
% - iphase:         Phase Number
% - name  :         string value of name for the Struct
%******************************************************************%

if ~isfield(Struct,'lb')
    error(strcat('Error: Lower Bounds have not been',...
        ' provided on ', name, ' in phase %i.\nIf',...
        ' there are no ', name, ' in this phase, then',...
        ' lower and upper bounds should be set to [].'),iphase);
elseif ~isfield(Struct,'ub')
    error(strcat('Error: Upper Bounds have not been',...
        ' provided on ', name, ' in phase %i.\nIf',...
        ' there are no ', name, ' in this phase, then',...
        ' lower and upper bounds should be set to [].'),iphase);
else
    if ~isempty(Struct.lb) && ~isempty(Struct.ub)
        if ~isequal(size(Struct.lb),size(Struct.ub))
            if name == ' States '
                error(strcat('Error: Upper and Lower Bound ',...
                    ' Matrices for ', name, ' must be same size',...
                    ' in phase %i.'),iphase);
            else
                error(strcat('Error: Upper and Lower Bound ',...
                    ' Vectors for ', name, ' must be same size',...
                    ' in phase %i.'),iphase);
            end
        end
    elseif (~isempty(Struct.lb) && isempty(Struct.ub)) ||...
            (isempty(Struct.lb) && ~isempty(Struct.ub))
        error(strcat('Error: Either the Upper or Lower Bounds for',...
            name,' is empty whilst the other is non-empty in phase %i.'),...
            iphase);
    end
end         

end

function boundsinputCheck2(Struct1,Struct2,iphase,name)
%******************************************************************%
% Validates input for the initial and final state
%
% Inputs:
%
% - problem struct
%******************************************************************%

if ~isfield(Struct1,'lb')
    error(strcat('Error: Lower Bounds have not been',...
        ' provided on ', name, ' Vector in phase %i.\nIf',...
        ' there are no States in this phase, then', name,...
        ' lower and upper bounds should be set to [].'),iphase);
    
elseif ~isfield(Struct1,'ub')
    error(strcat('Error: Upper Bounds have not been',...
        ' provided on ', name, ' Vector in phase %i.\nIf',...
        ' there are no States in this phase, then', name,...
        ' lower and upper bounds should be set to [].'),iphase);
else
    if ~isequal(size(Struct1.lb),size(Struct1.ub))
        error(strcat('Error: Sizes of ', name, ' Lower ',...
            ' and Upper Bound Vectors do not match each other',...
            ' in phase %i.'),iphase);
    else
        if ~isequal(size(Struct1.lb),size(Struct2.lb))
            error(strcat('Error: Size of ', name, ' Vector',...
                ' does not match the size of the State Vector',...
                ' in phase %i.'),iphase);
        end
    end
end

end