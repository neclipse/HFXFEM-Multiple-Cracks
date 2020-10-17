function  obj=updatedofarray_enriched( obj )
%UPDATEDOFS Method of Domain class to update both the standard and the
%enriched dofs of all nodes within the domain.
%  Enriched dofs
totnenrdof=obj.NoDofs+obj.NoEnrDofs;
totuenrdof=obj.NoUDofs+obj.NoUenrDofs;
totpenrdof=obj.NoPDofs+obj.NoPenrDofs;
% loop over all EnrichItems
BCTableqEn=[];
% BCTableqEn=cell(length(obj.EnrichItems));
% First loop to update doflocarray by the storage order in obj.EnrichItem
for iEnrich=1:length(obj.EnrichItems)
   enrichitem=obj.EnrichItems{iEnrich};
   % update the dofs indices for all the interacted nodes
   [totnenrdof, totuenrdof, totpenrdof]=updatedofarray(enrichitem,totnenrdof, totuenrdof, totpenrdof);
   % update the crack specific locarrayenr within the enriched elements
   id=enrichitem.Id;
   for ielem=1:length(enrichitem.NewElems)
       obj.ElemDict(enrichitem.NewElems(ielem)).givelocarray_enriched(id);
   end
   if ~isempty(enrichitem.Qtable)
      BCTableqEn{iEnrich}=enrichitem.Qtable;
   end
end
% Second loop to concatenate the crack-specific locarrayenr stored in
% JacobianMatDict together and store in JacobianMat. 10/12/20
for iEnrich=1:length(obj.EnrichItems)
    enrichitem=obj.EnrichItems{iEnrich};
    for ielem=1:length(enrichitem.NewElems)
       obj.ElemDict(enrichitem.NewElems(ielem)).assemble_locarry;
   end
end

% BCTableqEn=BCTableqEn(~isempty(BCTableqEn));
obj.NoEnrDofs=totnenrdof-obj.NoDofs;                % number of total enriched dofs
obj.NoUenrDofs=totuenrdof-obj.NoUDofs;              % number of total enriched displacement dofs
obj.NoPenrDofs=totpenrdof-obj.NoPDofs;              % number of total enriched pressure dofs

%% add the enriched dofs of the stdnodes in the psddofs
obj=obj.updatelinsys;                               % update the linsystem here
obj.LinSysCrt.BCTableqEn=BCTableqEn;
obj.LinSysCrt.upbctableen;                          % update the dirichlet dofs from enrichitems
obj.PsdEnrDofs=obj.LinSysCrt.BCTableEn(:,1);        % all the enriched dofs that need to be dropped
end

