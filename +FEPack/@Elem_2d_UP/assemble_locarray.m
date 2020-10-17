function assemble_locarray( obj )
% Assemble the locarray from individual JacobianMat within JacobianMatDict
% into the comprehensive elemental JacobianMat.
obj.JacobianMat.LocarrayEnr=[obj.JacobianMatDict.LocarrayEnr];
obj.JacobianMat.LocarrayUEnr=[obj.JacobianMatDict.LocarrayUEnr];
obj.JacobianMat.LocarrayPEnr=[obj.JacobianMatDict.LocarrayPEnr];
end

