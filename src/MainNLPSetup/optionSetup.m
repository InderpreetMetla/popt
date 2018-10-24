function problem = optionSetup(problem)
%******************************************************************%
% This function assigns some default settings for IPOPT.
% These can be overwritten by the user in the main setup script.
%
% Inputs:
% - problem: Main problem information file.
%
% Outputs:
% - problem: Updated problem struct with options info for ipopt
%
% Distributed under the MIT License
%
% Copyright (c) 2018 Inderpreet Metla
%******************************************************************%

% Define variable name for readability
if isfield(problem,'options')
    options = problem.options;
end

% Lower and Upper bounds on Decision Variables
options.lb = problem.DecVarLB.*problem.scaling.decvar_scales + ...
    problem.scaling.decvar_shifts;
options.ub = problem.DecVarUB.*problem.scaling.decvar_scales + ...
    problem.scaling.decvar_shifts;

% Lower and Upper bounds on Constraints
options.cl = problem.scaling.W*problem.CstLB;
options.cu = problem.scaling.W*problem.CstUB;

% Auxdata struct
if isfield(problem,'auxdata')
    options.auxdata = problem.auxdata;
end

% -------------------------------------------------------- %
% Set Default Settings
% -------------------------------------------------------- %

if isfield(options,'ipopt')
    
    if ~isfield(options.ipopt,'max_iter')
        options.ipopt.max_iter = 3000;
    end
    
    if isfield(options.ipopt,'linear_solver')
        %         if ~isequal(options.ipopt.linear_solver,'mumps')
        %             error('Error: IPOPT Solver must be mumps')
        %         end
        if isequal(options.ipopt.linear_solver,'mumps')
            options.ipopt.linear_solver = 'mumps';
        elseif isequal(options.ipopt.linear_solver,'ma57')
            options.ipopt.linear_solver = 'ma57';
        else
            error('Error: IPOPT Solver must be mumps or ma57')
        end
    else
        options.ipopt.linear_solver = 'mumps';
    end
    
    if ~isfield(options.ipopt,'tol')
        options.ipopt.tol = 1e-7;
    end
    
    if ~isfield(options.ipopt,'print_level')
        options.ipopt.print_level = 5;
    end
    
    if ~isfield(options.ipopt,'constr_viol_tol')
        options.ipopt.constr_viol_tol = 0.0001;
    end
    
    if ~isfield(options.ipopt,'nlp_scaling_method')
        options.ipopt.nlp_scaling_method = 'gradient-based';
    end
    
    if ~isfield(options.ipopt,'nlp_scaling_max_gradient')
        options.ipopt.nlp_scaling_max_gradient  = 100;
    end
    
    if ~isfield(options.ipopt,'nlp_scaling_min_value')
        options.ipopt.nlp_scaling_min_value = 1e-8;
    end
    
    if ~isfield(options.ipopt,'mu_strategy')
        options.ipopt.mu_strategy = 'adaptive';
    end
    
    if problem.grid.iter == 0
        if ~isfield(options.ipopt,'print_user_options')
            options.ipopt.print_user_options = 'yes';
        end
    else
        options.ipopt.print_user_options = 'no';
    end
    
    if ~isfield(options.ipopt,'output_file')
        options.ipopt.output_file = strcat(problem.name,...
            '-ipopt-problem-info');
    end
else
    options.ipopt.max_iter                  = 3000;
    options.ipopt.linear_solver             = 'mumps';
    options.ipopt.tol                       = 1e-7;
    options.ipopt.print_level               = 5;
    options.ipopt.constr_viol_tol           = 0.0001;
    options.ipopt.nlp_scaling_method        = 'gradient-based';
    options.ipopt.nlp_scaling_max_gradient  = 100;
    options.ipopt.nlp_scaling_min_value     = 1e-8;
    options.ipopt.mu_strategy               = 'adaptive';
    if problem.grid.iter == 0
        options.ipopt.print_user_options    = 'yes';
    else
        options.ipopt.print_user_options    = 'no';
    end
    options.ipopt.output_file = strcat(problem.name,...
        '-ipopt-problem-info');
end

% Derivative Options
if isfield(problem,'derivatives')
    if isfield(problem.derivatives,'order')
        if isequal(problem.derivatives.order,2)
            options.ipopt.hessian_approximation = 'exact';
        else
            options.ipopt.hessian_approximation = 'limited-memory';
        end
    else
        problem.derivatives.order = 1;
        options.ipopt.hessian_approximation = 'limited-memory';
    end
    if ~isfield(problem.derivatives,'first')
        if ~isfield(problem.derivatives,'method') || ...
                isequal(problem.derivatives.method,'fd') || ...
                isequal(problem.derivatives.method,'FD') || ...
                isequal(problem.derivatives.method,'forward-difference')
            problem.derivatives.first.stepsize.gradient = 3e-06;
            problem.derivatives.first.stepsize.jacobian = 3e-08;
        elseif isequal(problem.derivatives.method,'cs') || ...
                isequal(problem.derivatives.method,'CS')
            problem.derivatives.first.stepsize.gradient = 1e-20;
            problem.derivatives.first.stepsize.jacobian = 1e-20;
        end
    else
        if isfield(problem.derivatives.first,'stepsize')
            if ~isfield(problem.derivatives,'method') || ...
                    isequal(problem.derivatives.method,'fd') || ...
                    isequal(problem.derivatives.method,'FD') || ...
                    isequal(problem.derivatives.method,'forward-difference')
                if ~isfield(problem.derivatives.first.stepsize,'gradient')
                    problem.derivatives.first.stepsize.gradient = 3e-06;
                end
                if ~isfield(problem.derivatives.first.stepsize,'jacobian')
                    problem.derivatives.first.stepsize.jacobian = 3e-08;
                end
            elseif isequal(problem.derivatives.method,'cs') || ...
                    isequal(problem.derivatives.method,'CS')
                if ~isfield(problem.derivatives.first.stepsize,'gradient')
                    problem.derivatives.first.stepsize.gradient = 1e-20;
                end
                if ~isfield(problem.derivatives.first.stepsize,'jacobian')
                    problem.derivatives.first.stepsize.jacobian = 1e-20;
                end
            end
        else
            if ~isfield(problem.derivatives,'method') || ...
                    isequal(problem.derivatives.method,'fd') || ...
                    isequal(problem.derivatives.method,'FD') || ...
                    isequal(problem.derivatives.method,'forward-difference')
                problem.derivatives.first.stepsize.gradient = 3-06;
                problem.derivatives.first.stepsize.jacobian = 3e-08;
            elseif isequal(problem.derivatives.method,'cs') || ...
                    isequal(problem.derivatives.method,'CS')
                problem.derivatives.first.stepsize.gradient = 1e-20;
                problem.derivatives.first.stepsize.jacobian = 1e-20;
            end
        end
    end
    if problem.derivatives.order == 2 && ...
            (~isfield(problem.derivatives,'second') || ...
            ~isfield(problem.derivatives.second,'stepsize'))
        problem.derivatives.second.stepsize = 1.75e-05;
    end
else
    problem.derivatives.method = 'fd';
    problem.derivatives.order = 1;
    problem.derivatives.first.stepsize.gradient = 3e-06;
    problem.derivatives.first.stepsize.jacobian = 3e-08;
    options.ipopt.hessian_approximation = 'limited-memory';
end


if isequal(options.ipopt.nlp_scaling_method,'none')
    fields = {'nlp_scaling_max_gradient'
        'nlp_scaling_min_value'};
    options.ipopt=rmfield(options.ipopt,fields);
end

% Output struct
problem.options = options;
end
