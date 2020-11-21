function elems=findelems( obj,plist,varargin )
%FINDRTIPELEMS Find the elements according to the points in plist
%   Detailed explanation goes here
% Method of opengeo objects Varargin sets the search options, 
%1. By default, no varargin, option ='initial', search enveloping element 
% and element with the point on its edge, also check if the element is cut 
% by the segment (after extension to edge)
%2. option = 'edge', search elements with point on its edge
%3. option= 'in_edge', search elements with point inside or on their edges
if isempty(varargin) || strcmp(varargin{1},'initial')
    option=1;
elseif strcmp(varargin{1},'edge')
    option=2;
elseif strcmp(varargin{1},'in_edge')
    option=3;
end
    elems=zeros(4*size(plist,1),1);
    iloc=1;
    for ip=1:size(plist,1)
        star_n=plist(ip,:);
        [elemnums,~,~]=obj.Mesh.locate(star_n(1),star_n(2));
        for ielem=1:length(elemnums)
            em=elemnums(ielem);                                             % elem index
            elem=obj.Elemdict(em);                                          % elem obj
            [flagie,~,flage,~,~] = isinside_vec(elem,star_n);
            if option ==1
                if flagie
                    %check if the one of the two tips are exactly on the
                    %element edge, if so, one should not count the one not
                    %cut by the crack
                    if ip==1 && size(plist,1)>1
                        star_nnext=plist(ip+1,:);
                    elseif ip==size(plist,1) && size(plist,1)>1
                        star_nnext=plist(ip-1,:);
                    else
                        elems(iloc)=em;
                        iloc=iloc+1;
                        continue;
                    end
                    xy1=[star_n,star_nnext];
                    xy2=[elem.X(1),elem.Y(1),elem.X(3),elem.Y(3);];
                    out=lineSegmentIntersect(xy1,xy2);
                    if out.intAdjacencyMatrix==1                            % the diagnal of the this element intersects the crack
                        elems(iloc)=em;
                        iloc=iloc+1;
                    end
                end
            elseif option==2
                if flage
                    elems(iloc)=em;
                    iloc=iloc+1;
                end
            elseif option==3              
                if flagie
                    elems(iloc)=em;
                    iloc=iloc+1;
                end
            end
        end
    end
    elems=elems(1:iloc-1);                                                  % only counts the filled elements (dump the zeros)
    elems=unique(elems);
end


