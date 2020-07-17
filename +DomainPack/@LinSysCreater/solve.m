function solve(obj,newmark,varargin)
%% Choose solver depending on the size of the linear equation system obj.Drow
% Direct solver based on factorization or elimination technique
if isempty(varargin)
    mode=1;
else
    mode=varargin{1};    
end
if mode==2
    LHS=obj.LHSnew;     % linear because Kc is left out in the LHSnew
else
    LHS=obj.LHS;        % nonlinear because of Kc in the LHS
end
if obj.Drow<=100000 
    if isempty(obj.EnrichItems)
        % The standard matrix is already well Banded
        obj.Unknowns=LHS\obj.RHS;
    else
        % The enriched matrix may disturb the banded structure
        % Use the Sparse reverse Cuthill-McKee ordering to reduce the
        % bandwidth of the LHS in order to improve the solution speed. (03032019)
        r=symrcm(LHS);
        rcmlhs=LHS(r,r);
        rcmrhs=obj.RHS(r);
        rcmunknowns=rcmlhs\rcmrhs;
        [~,rI]=sort(r);
        obj.Unknowns=rcmunknowns(rI);   %recover the orginal order
    end
else
% Iterative solver GMRES with ilu preconditioner, ILU(0) or milu
    setup.type='nofill';
    [L,U]=ilu(LHS,setup);
    xiter=gmres(LHS,obj.RHS,10,1e-10,20,L,U);
    obj.Unknowns=xiter;
end
% update all field parameters using the newly solved [u(n+1) and p(n+1)]
obj.upconf(newmark);
end