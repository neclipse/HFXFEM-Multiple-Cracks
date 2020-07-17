classdef Postprocessor < matlab.mixin.SetGet
properties
    IInc                % The index of increment when the postprocessor is created
    Inc                 % The value of increment when the postprocessor is created
    RLOADX              % Converged internal load vector is actually the reaction force vector.
    RLOADY              % No need for enriched nodes
    XN                  % X-coordinates of all nodes in global order
    YN                  % Y-coordinates of all nodes in global order
    UX                  % Total X-displacement of all nodes in global order
    UY                  % Total Y-displacement of all nodes in global order
    P                   % Total Pore pressure of all nodes in global order
    ESnN                % Averaged strain at all nodes in global order
    SN                  % Averaged total stress at all nodes in global order
    SPN                 % Averaged effective stress at all nodes in global order
    EnrichItems         % Practical results information from Domain.EnrichItem (Deep copy)
end

methods
    painter(obj,varargin); 
    h=plot3d(obj,variable,varargin);                      % varargin is to specify which component of the variable
    h=plotcontour(obj,variable,meshstructure,varargin);   % varargin is to specify which component of the variable
    h=plotsurf(obj,variable,meshstructure,varargin);      % varargin is to specify which component of the variable
    [h,x,y]=plotxy(obj,variable,indices,varargin);              % varargin is to specify which component of the variable
    plotcrack( obj,mesh,varargin );
end

end
