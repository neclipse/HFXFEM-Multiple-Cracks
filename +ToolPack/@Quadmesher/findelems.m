function elems=findelems( obj,plist,elemdict,varargin )
%FINDRTIPELEMS Find the elements according to the points in plist
%   Detailed explanation goes here
% Method of Mesh objects
% Varargin sets the search options,
% 1. By default, no varargin, option ='inside', search enveloping element
% 2. option = 'edge', search elements with point on one of their edge
% 3. option= 'in_edge', search elements with point inside or on their edges
if isempty(varargin) || strcmp(varargin{1},'inside')
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
        [elemnums,~,~]=obj.locate(star_n(1),star_n(2));
        for ielem=1:length(elemnums)
            em=elemnums(ielem);                                             % elem index
            elem=elemdict(em);                                          % elem obj
            [flagie,flagi,flage,~] = isinside_vec(elem,star_n);             
            if option ==1
                if flagi
                    elems(iloc)=em;
                    iloc=iloc+1;
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
    elems=elems(1:iloc-1);
    elems=unique(elems);
end

