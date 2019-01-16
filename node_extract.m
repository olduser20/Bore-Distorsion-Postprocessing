function coord_s = node_extract( coord,disp )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    n=1;
    for i=1:length(coord)
        for ii=1:length(disp)
          if coord(i,1)== disp (ii,1)
             coord_s(n,1:4)= coord(i,1:4);
             n=n+1;
          end
        end
    end
%     coord_s=coord(any(coord(:,1)==disp(:,1),2)==1,1:4);
    
    
end

