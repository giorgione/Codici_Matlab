function KLT2 = preprocessKLT(KLT)
featcnt = sum(sum(KLT.val ~= 0));

KLT2.x = sparse(zeros(featcnt, size(KLT.x, 2)));
KLT2.y = sparse(zeros(featcnt, size(KLT.x, 2)));
KLT2.val = zeros(featcnt, 1);

idx = 1:size(KLT.x, 1);

KLT2.x(idx, 1) = KLT.x(:, 1);
KLT2.y(idx, 1) = KLT.y(:, 1);
KLT2.val(idx, 1) = KLT.val(:, 1);

lastidx = idx(end);
for i =2:size(KLT.x, 2)
    tempidx = (KLT.val(:,i) ~= 0);
    newidx = (lastidx+1):(lastidx+ sum(tempidx));
    lastidx = lastidx+ sum(tempidx);
    idx(tempidx) = newidx;
    
    KLT2.x(idx, i) = KLT.x(:, i);
    KLT2.y(idx, i) = KLT.y(:, i);
    KLT2.val(newidx) = KLT.val(tempidx, i);
end

KLT2.vidx = 1:length(KLT2.val);

end