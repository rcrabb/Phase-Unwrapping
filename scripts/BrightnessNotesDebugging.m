%     %test our projection %%%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Zproj = Rcalc.*((sum(A.^2,3)./A(:,:,3).^2)).^-.5;
%     Xproj = Zproj.*A(:,:,1)./A(:,:,3);
%     Yproj = Zproj.*A(:,:,2)./A(:,:,3);
%     Rtest = dR(nanmap);
%     figure;hist(Rtest(sample(inliers)),100);
% %     figure;plot3(XYZ(1,:),XYZ(2,:),XYZ(3,:),'r.');
% %     hold;plot3(XYZ(1,sample),XYZ(2,sample),XYZ(3,sample),'b.');
% %     plot3(XYZ(1,sample(inliers)),XYZ(2,sample(inliers)),XYZ(3,sample(inliers)),'g.');
% %     plot3(Xproj(nanmap),Yproj(nanmap),Zproj(nanmap),'m.');
%     above = dR(nanmap);
%     below = above(sample)<=0;
%     above = above(sample)>0;
%     figure;plot3(XYZ(1,below),XYZ(2,below),XYZ(3,below),'r.',...
%         XYZ(1,above),XYZ(2,above),XYZ(3,above),'b.');
%     hold;
%     %plot3(XYZ(1,sample),XYZ(2,sample),XYZ(3,sample),'b.');
%     plot3(XYZ(1,sample(inliers)),XYZ(2,sample(inliers)),XYZ(3,sample(inliers)),'g.');
%     plot3(Xproj(nanmap),Yproj(nanmap),Zproj(nanmap),'m.');
%     
%     [coefsT, PT, inliersT] = ransacfitplane(XYZ(:,sample), 30, 1);
%     normT = coefsT/norm(coefsT(1:3))*sign(coefsT(3));
%     CosBT = abs(sum(A.*repmat(reshape(normT(1:3),[1 1 3]),[M N 1]),3));
%     RcalcT = -normT(4)./CosBT;
%     dR = RcalcT-Rcalc;
%     dR = dR(nanmap)
%     figure;hist(dR(sample(inliers)),100);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

figure;hist(Rtest(sample(inliers)),100);
   XYZsample = XYZ(:,sample);
    figure;plot3(XYZ(1,:),XYZ(2,:),XYZ(3,:),'r.');
    hold;
    plot3(XYZ(1,sample),XYZ(2,sample),XYZ(3,sample),'b.');
    plot3(XYZ(1,sample(inliers)),XYZ(2,sample(inliers)),XYZ(3,sample(inliers)),'m.');
    
    
    plot3(XYZsample(1,:),XYZsample(2,:),XYZsample(3,:),'b.');
    plot3(XYZsample(1,inliers),XYZsample(2,inliers),XYZsample(3,inliers),'g.');
    il = zeros(size(sample));
  
    Zproj = Rcalc.*((sum(A.^2,3)./A(:,:,3).^2)).^-.5;
    Xproj = Zproj.*A(:,:,1)./A(:,:,3);
    Yproj = Zproj.*A(:,:,2)./A(:,:,3);
    Rtest = dR(nanmap);
    proj = Xproj;
    proj(:,:,2) = Yproj;
    proj(:,:,3) = Zproj;
    XYZproj = [Xproj(nanmap)';Yproj(nanmap)';Zproj(nanmap)'];
    
    
    plot3(Xproj(nanmap),Yproj(nanmap),Zproj(nanmap),'k.');
    plot3(XYZproj(1,sample),XYZproj(2,sample),XYZproj(3,sample),'m.');
    above = dR(nanmap);
    below = above(sample)<=0;
    above = above(sample)>0;
    figure;plot3(XYZ(1,below),XYZ(2,below),XYZ(3,below),'r.',...
        XYZ(1,above),XYZ(2,above),XYZ(3,above),'b.');
    hold;
    %plot3(XYZ(1,sample),XYZ(2,sample),XYZ(3,sample),'b.');
    plot3(XYZ(1,sample(inliers)),XYZ(2,sample(inliers)),XYZ(3,sample(inliers)),'g.');
    plot3(Xproj(nanmap),Yproj(nanmap),Zproj(nanmap),'m.');
    
    
    hist(sum(XYZ(:,sample).^2,1)-sum(XYZproj(:,sample).^2),100)
    
    hist(sum(XYZ(:,inliers).^2,1)-sum(XYZproj(:,inliers).^2),100)