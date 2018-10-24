function outputdisplay(PrintFlag)
%**********************************************************************%
% This function prints some information at the end of the solution
%
% Distributed under the MIT License
%
% Copyright (c) 2018 Inderpreet Metla
%**********************************************************************%
if PrintFlag == 1
    disp(' ')
    disp('********************************************************************************')
    disp('*                               SOLUTION DATA                                  *')
    disp('********************************************************************************')
    disp('*            solution.phase.time       =    Time Solution                      *')
    disp('*            solution.phase.state      =    State Solution                     *')
    disp('*            solution.phase.control    =    Control Solution                   *')
    disp('*            solution.cost             =    Total Cost                         *')
    disp('*            solution.bndcost          =    Mayer Cost                         *')
    disp('*            solution.pathcost         =    Lagrange Cost                      *')
    disp('*            solution.status           =    Ipopt Termination Status           *')
    disp('*            solution.nlp_iters        =    Number of Ipopt Iterations         *')
    disp('*            solution.f_evals          =    Number of Function Evals           *')
    disp('*            solution.cpu_time.grids   =    Ipopt CPU Time Per Grid Refinement *')
    disp('*            solution.cpu_time.total   =    Total Ipopt CPU Time               *')
    disp('*            solution.grid_iters       =    Number of Grid Iterations          *')
    disp('*            solution.grid_analysis    =    Data from Each Solution Grid       *')
    disp('********************************************************************************')
    disp('*                               PLOTTING DATA                                  *')
    disp('********************************************************************************')
    disp('*            plotdata.phase.time       =    Time Plotting Data                 *')
    disp('*            plotdata.phase.state      =    State Plotting Data                *')
    disp('*            plotdata.phase.control    =    Control Plotting Data              *')
    disp('********************************************************************************')
elseif PrintFlag == 2
    disp(' ')
    disp('********************************************************************************')
    %disp('*                                                                              *')
    disp('*                               SOLUTION DATA                                  *')
    disp('********************************************************************************')
    disp('*            solution.phase.time       =    Time Solution                      *')
    disp('*            solution.phase.state      =    State Solution                     *')
    disp('*            solution.phase.control    =    Control Solution                   *')
    disp('*            solution.cost             =    Total Cost                         *')
    disp('*            solution.bndcost          =    Mayer Cost                         *')
    disp('*            solution.pathcost         =    Lagrange Cost                      *')
    disp('*            solution.status           =    Ipopt Termination Status           *')
    disp('*            solution.nlp_iters        =    Number of Ipopt Iterations         *')
    disp('*            solution.f_evals          =    Number of Function Evals           *')
    disp('*            solution.cpu_time.grids   =    Ipopt CPU Time Per Grid Refinement *')
    disp('*            solution.cpu_time.total   =    Total Ipopt CPU Time               *')
    disp('*            solution.grid_iters       =    Number of Grid Iterations          *')
    disp('*            solution.grid_analysis    =    Data from Each Solution Grid       *')
    disp('********************************************************************************')
end