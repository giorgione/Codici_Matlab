function ol = getOverlap(A, B)
% inputs = [x, y, width, height]

A1 = A + [0, 0, A(1:2) - 1];
B1 = B + [0, 0, B(1:2) - 1];

obox = [max(A1(1), B1(1)), max(A1(2), B1(2)), min(A1(3), B1(3)), min(A1(4), B1(4))];

oa = max(0, obox(3)-obox(1)+1) * max(0, obox(4)-obox(2)+1);
ol = oa / (A(3) * A(4) + B(3) * B(4) - oa);
end