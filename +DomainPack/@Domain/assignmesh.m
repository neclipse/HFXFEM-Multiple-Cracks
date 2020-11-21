function obj=assignmesh(obj) 
   Mesher=obj.Preprocess.Mesher;
% NodeDict
if obj.MatType==4
   nodedict(1:Mesher.Totnodes)=FEPack.Node_UP();
   for inod=1:Mesher.Totnodes
       x=Mesher.VX(inod);
       y=Mesher.VY(inod);
       type=Mesher.VT(inod);
       nodedict(inod)=FEPack.Node_UP(inod,type,x,y); 
   end
else
   nodedict(1:Mesher.Totnodes)=FEPack.Node();
   for inod=1:Mesher.Totnodes
       x=Mesher.VX(inod);
       y=Mesher.VY(inod);
       type=Mesher.VT(inod);
       nodedict(inod)=FEPack.Node(inod,type,x,y); 
   end
end
   obj.NodeDict=nodedict;
% ElemDict
gaussdictm=obj.crtgpd;                                    % create an object array of gaussian points (value type) 
if obj.ElemType == 2
    % Bug: CHANGE COLON TO COMMA AS COLON WOULD GIVE A DICTIONARY OF ALL
    % REPEATED HANDLES
    elemdictionary(1,Mesher.Totelem)=FEPack.Elem_2d_EP(); % element class should have default constructor
    for ielem=1:Mesher.Totelem
        elemdictionary(ielem)=FEPack.Elem_2d_EP(ielem,obj.ElemType,obj.MatType, Mesher.EToV(ielem,:), nodedict(Mesher.EToV(ielem,:)),gaussdictm);
    end
elseif obj.ElemType == 4
    elemdictionary(1,Mesher.Totelem)=FEPack.Elem_2d_UP(); % element class should have default constructor
    for ielem=1:Mesher.Totelem
        elemdictionary(ielem)=FEPack.Elem_2d_UP(ielem,obj.ElemType,obj.MatType, Mesher.EToV(ielem,:), nodedict(Mesher.EToV(ielem,:)),gaussdictm);
    end
end
   obj.ElemDict=elemdictionary; %IMPORTANT: IT IS NOT ALLOWED TO CREATE AN OBJECT ARRAY AND ASSIGN IT TO A PROPERTY OF ANOTHER OBJECT
end
