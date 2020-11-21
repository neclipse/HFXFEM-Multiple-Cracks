function realgrow(obj)
obj.NewElems=[];
obj.NewNodes=[];
%% 2. if obj.GrowCheck.Growflag, then calculate grow direction
if obj.Isactive
    if obj.Growcheck.Growflag
        % theta is the crack grow turning direction from the current tip
        theta=obj.Growcheck.growdirection;
        % 3. calculate obj.Growincrement (embeded in obj.findnextelem)
%         elemahead1=obj.NextElem.Ind;
        [elemahead2,pnts,~,Phi]=obj.findnextelem(theta);
%         if elemahead1~=elemahead2
%             warning('The grow direction may turn too much.')
%         end
%         map=[2,1];
%         oldtip_temp=pnts(map(obj.Itip),:);
%         dist=sqrt((oldtip_temp(1)-obj.Xct)^2+(oldtip_temp(2)-obj.Yct)^2);
%         if dist/min(obj.INTELEM.Length)>0.1
%             warning('The current INTELEM may need get updated');
%         end
%         obj.NextElem.Seeds{obj.Id}=seeds; % seeds generated in
%         elem.subdomain
%         obj.NextElem.LocalInt{obj.Id}=localpnts;
%         obj.NextElem.GlobalInt{obj.Id}=pnts;
        %% 4. Really advance the crack front
        % a. The increment is inserted to Mygeo, the tip
        % element, tip nodes are accordingly updated inside
        % Mygeo, the level sets are also to be updated
        newpoint=[size(obj.Mygeo.Segments,1)+1,pnts(obj.Itip,:)];
        % newpoint is obtained this way because pnts has been ordered in
        % mygeo.intersection.
        if obj.Itip==1
            obj.Mygeo.Segments=[obj.Mygeo.Segments;newpoint];
        elseif obj.Itip==2
            obj.Mygeo.Segments=[newpoint;obj.Mygeo.Segments];
        end
        %% 5. update related info, especially about geometry and the lsv
        % here we can simply use obj.Mygeo.initiate_2nd to
        % finish the update although it may not be the optimal
        % solution when crack is simply exented.
        % elements
        obj.Mygeo.Intelements=[obj.Mygeo.Intelements;elemahead2];
        % The Rtipelements were corrected on 07092020 as previous code did
        % not prevent repeated values in mygeo.initiate_2nd.
        obj.Mygeo.Rtipelements(obj.Itip)=elemahead2;
        % Bodyelements are correctly determined this way
        % althought it is not used later on. 07082019. Given correct
        % Rtipelements and Intelements. 07092020
        obj.Mygeo.Bodyelements=setdiff(obj.Mygeo.Intelements,...
            obj.Mygeo.Rtipelements);
        obj.NewElems=[obj.NewElems,elemahead2];
        % nodes
        % body nodes, rtipnodes are unnecessarily to be updated
        % so they are not updated by 07/08/2019
        % Updated on 07/09/2020 as crack.plotme will use these info.
        obj.Mygeo.Rtipnodes(4*obj.Itip-3:4*obj.Itip)=obj.Elemdict(elemahead2).NodList';
        [AllNodes,ia,~]=unique([obj.Mygeo.Nodes;obj.NextElem.NodList']);
        obj.NewNodes=setdiff(AllNodes,obj.Mygeo.Nodes);
        % (Optional)
        % The function to update old tip element was not
        % finished by 06/21/2019 because the
        % elem.crtstif_enriched will also be changed and it may
        % induce more instabilities, I guess. So it was
        % postponed until further investigation validates it.
        %                    oldnodestoupdate=setdiff(obj.NextElem.NodList,obj.NewNodes);
        %                    obj.UpdateElems=obj.Interactedelem;
        obj.Mygeo.Nodes=AllNodes;
        obj.Mygeo.Bodynodes=setdiff(AllNodes,obj.Mygeo.Rtipnodes);
        % Important level set information will be used in enfs.
        temp=[obj.Mygeo.Phi(:,2);Phi];
        AllPhi=temp(ia);
        %                    for inold=1:length(oldnodestoupdate)
        %                        AllPhi(AllNodes==oldnodestoupdate(inold))=...
        %                        Phi(obj.NextElem.NodList==oldnodestoupdate(inold));
        %                    end
        obj.Mygeo.Phi=[AllNodes,AllPhi];
        % To update the obj.Mygeo.Stdnodes;
        obj.Mygeo.Toedge2(obj.Itip);
        % update other info of obj.Mygeo
%         obj.Mygeo.findblending(2);        %  update the blending standard elements
%         elements if the isolation is only needed at the crack mouth.
        obj.Mygeo.callength;
        obj.Mygeo.calelen;
        %% 6. Update the crack tip
        obj.update;
    end
end
end