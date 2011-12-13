function [patch] = getImgPatch(Im, bb)

patch = zeros(bb(4), bb(3), 3);

for row = 1:bb(4)
    for col = 1:bb(3)
        posx = bb(1) + col - 1;
        posy = bb(2) + row - 1;
        
        posx = max(min(posx, size(Im, 2)), 1);
        posy = max(min(posy, size(Im, 1)), 1);
        
        %patch(row, col, :) = Im(posy, posx, :); 
        patch(row, col, 1) = Im(posy, posx, 1); 
        patch(row, col, 2) = Im(posy, posx, 2); 
        patch(row, col, 3) = Im(posy, posx, 3); 
    end
end

end