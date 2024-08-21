function [out] = importImages2(images)
    N=size(images,3);
    L=size(images,1);
    for I=1:N
        curImg=images(:,:,I);
        if L>256
            curImg=fadeBorderCos(curImg,10);
        else
            curImg=fadeBorderCos(curImg,5);
        end
        images(:,:,I)=curImg;
    end
    out=images;
end