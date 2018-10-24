function problem = NumofDecVarsandCsts(problem)
%******************************************************************%
% Get sizes in each phase of the problem for the: states, controls,
% path constraints and boundary constraints          
%
% Also gets the Bounds on the Decision Variables and Constraints
% 
% Inputs:
%  - problem struct
%
% Outputs:
%   - problem struct with:
%       - .nphases field :  Number of phases for the problem
%       - .sizes field   :  Contains a matrix of the sizes of the 
%                           problem problem. each row is the number 
%                           of states, controls, path constraints 
%                           and boundary constraints in the iphase
%                           (in that order).
%
% Distributed under the MIT License
%
% Copyright (c) 2018 Inderpreet Metla
%******************************************************************%

% Number of phases
nphases = length(problem.bounds.phase);
problem.nphases = nphases;
% Rename structs for better readability
bounds  = problem.bounds.phase;
% Number of structures to get sizes of
nstructs = 4;
% Pre-allocate zero matrix for sizes of each struct
SizeofDecVarsandCsts = zeros(nstructs,nphases);

for iphase=1:nphases
    %SizeofDecVarsandCsts(:,iphase) = subFunction(bounds,iphase);
    if isfield(bounds(iphase),'state')
        nStates = subFunction(bounds(iphase).state,iphase,'State');
    else
        nStates = 0;
    end
    
    if isfield(bounds(iphase),'control')
        nControls = subFunction(bounds(iphase).control,iphase,'Control');
    else
        nControls = 0;
    end
    
    if isfield(bounds(iphase),'path')
        nPathCsts = subFunction(bounds(iphase).path,iphase,'Path');
    else
        nPathCsts = 0;
    end
    
    if isfield(bounds(iphase),'boundary')
        nBndCsts = subFunction(bounds(iphase).boundary,iphase,'Boundary');
    else
        nBndCsts = 0;
    end
    
    SizeofDecVarsandCsts(:,iphase) = ...
        [nStates; nControls; nPathCsts; nBndCsts];
end
% Assign a field to the matrix of sizes
problem.sizes = SizeofDecVarsandCsts;
end

function Number = subFunction(Struct,iphase,name)
%******************************************************************%
% This function calculates the number of states, controls, 
% path constraints or boundary constraints in each phase 
% of the problem depending on what the input name is
% 
% Inputs:
%  - Struct:    bounds(iphase).X struct where X = state, controls,
%               path or boundary
%  - iphase:    Phase Number
%  - name  :    String associated with X. Used in an if condition. 
%
% Outputs:
%   - Number:   Size of Struct .lb and .ub (number of columns)  
%******************************************************************%

if (isfield(Struct,'lb') && ...
        isfield(Struct,'ub'))
    if isempty(Struct.lb) && isempty(Struct.ub)
        Number = 0;
    else
        Number = size(Struct.lb,2);
        % Check for correct number of rows
        if ~isequal(size(Struct.lb),[1,Number])
            if name == 'Path' || name == 'Boundary'
            error(strcat('Error: ', name, ' Bound Vectors must',...
                ' have size [1,#',name,'Constraints] in phase %i.'),...
                iphase);
            else
            error(strcat('Error: ', name, ' Bound Vectors must',...
                ' have size [1,#',name,'s] in phase %i.'),iphase);    
            end
        end
    end
end

end

% function Number = subFunction(bounds,iphase)
% %******************************************************************%
% % This function calculates the number of states, controls, 
% % path constraints and boundary constraints in each phase 
% % of the problem
% % 
% % Inputs:
% %  - bounds:    User supplied problem bounds struct
% %  - iphase:    Phase Number
% %
% % Outputs:
% %   - Number:   Outputs a column vector where each row is the number 
% %               of states, controls, path constraints and boundary
% %               constraints in the iphase (in that order). 
% %******************************************************************%
% 
% %  States
% if isfield(bounds(iphase),'state')
%     if (isfield(bounds(iphase).state,'lb') && ...
%             isfield(bounds(iphase).state,'ub'))
%         stateslb = bounds(iphase).state.lb;
%         statesub = bounds(iphase).state.ub;
%         if isempty(stateslb) && isempty(statesub)
%             nStates = 0;
%         else
%             nStates = size(stateslb,2);
%             if ~isequal(size(stateslb),[1,nStates])
%                 error(strcat('Error: State Bound Vectors must',...
%                     ' have size [1,#states] in phase %i.'),iphase);
%             end
%         end
%     end
% else
%     nStates = 0;
% end
% 
% % Controls
% if isfield(bounds(iphase),'control')
%     if (isfield(bounds(iphase).control,'lb') && ...
%             isfield(bounds(iphase).control,'ub'))
%         controlslb = bounds(iphase).control.lb;
%         controlsub = bounds(iphase).control.ub;
%         if isempty(controlslb) && isempty(controlsub)
%             nControls = 0;
%         else
%             nControls = size(controlslb,2);
%             if ~isequal(size(controlslb),[1,nControls])
%                 error(strcat('Error: Control Bound Vectors must',...
%                     ' have size [1,#controls] in phase %i.'),iphase);
%             end
%         end
%     end
% else
%     nControls = 0;
% end
% 
% % Path Constraints
% if isfield(bounds(iphase),'path')
%     if (isfield(bounds(iphase).path,'lb') && ...
%             isfield(bounds(iphase).path,'ub'))
%         pathCstlb = bounds(iphase).path.lb;
%         pathCstub = bounds(iphase).path.ub;
%         if isempty(pathCstlb) && isempty(pathCstub)
%             nPathCsts = 0;
%         else
%             nPathCsts = size(pathCstlb,2);
%             if ~isequal(size(pathCstlb),[1,nPathCsts])
%                 error(strcat('Error: Path Constraint Bound',...
%                     ' Vectors must have size [1,#pathCsts] in',...
%                     ' phase %i.'),iphase);
%             end
%         end
%     end
% else
%     nPathCsts = 0;
% end
% 
% % Boundary Constraints
% if isfield(bounds(iphase),'boundary')
%     if (isfield(bounds(iphase).boundary,'lb') && ...
%             isfield(bounds(iphase).boundary,'ub'))
%         bndCstlb = bounds(iphase).boundary.lb;
%         bndCstub = bounds(iphase).boundary.ub;
%         if isempty(bndCstlb) && isempty(bndCstub)
%             nBndCsts = 0;
%         else
%             nBndCsts = size(bndCstlb,2);
%             if ~isequal(size(bndCstlb),[1, nBndCsts])
%                 error(strcat('Error: Boundary Constraint Bound',...
%                     ' Vectors must have size [1,#bndCsts] in',...
%                     ' phase %i.'),iphase);
%             end
%         end
%     end
% else
%     nBndCsts = 0;
% end
% % Output vector
% Number = [nStates nControls nPathCsts nBndCsts]';
% 
% end

