function [result] = Evaluation(Beaconlocation,Distance,h,X)
w = zeros(length(h),1); 
for i = 1:length(h)
    w(i,1) = h(i) / sum(h);
end 
[m, n] = size(Beaconlocation);
result = sum(w .* abs(sqrt(sum(((repmat(X,m,1) - Beaconlocation).^2),2)) - Distance));
