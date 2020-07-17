function matprop = materialProp( mt,et,E,nu,c0,phi,psi,H )
%MATERIALPROP Calculate the material properties for one elements
%   mt is the material property type set in Domain
%   et is the element type set in Domain
    G=0.5*E/(1+nu);                             % Shear modulus,in the unit of Mpa
    lamda=2*nu*G/(1-2*nu);                      % The first Lame constant, in the unit of MPa
    xi=3/(sqrt(9+12*tan(phi)^2));               % parameter for Drucker-Prager yield surface
    eta=tan(phi)*xi;                            % parameter for Drucker-Prager yield surface
    etabar=3*tan(psi)/(sqrt(9+12*tan(psi)^2));  % parameter for non-associative Drucker-Prager flow rule
%% plane strain or axisymmetric
if et==2|| et==3                                
    % -Linear elastic
    matprop.lamda=lamda;                        
    matprop.G=G; 
    % -Drucker-Prager
    if mt==2                
    matprop.c0=c0;
    matprop.xi=xi;          
    matprop.eta=eta;        
    matprop.etabar=etabar;  
    matprop.H=H;                             
    end
%% plane stress
elseif et==1 
    % -Linear elastic
    matprop.lamda=3367.0;    
    matprop.gb=1893.9;  
    % -Linear elastic
    if mt==2                
    matprop.theta=33;
    matprop.c=1000;
    end
end
end


