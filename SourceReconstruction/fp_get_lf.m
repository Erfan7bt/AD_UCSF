function L = fp_get_lf(L1)
% keyboard
%L1 = inverse.MEG.L;
ns = size(L1,2);
L = zeros(size(L1,1),size(L1,2),2);
for is=1:ns
%    clear L2
%    L2 = L1{is}./10^-12;   
    %remove radial orientation
    clear u s
    [u, s, ~] = svd(squeeze(L1(:,is,:)));
    %     L(:,is,:) = u(:,1:2)*s(1:2,1:2);
    %     L(:,is) = u(:,:)*s(:,1);
    L(:,is,:) = u(:,:)*s(:,1:2);
end
