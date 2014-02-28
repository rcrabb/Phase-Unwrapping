function zseqmask(grayscale, masklayer, saturationlayer)

grayscale = double(grayscale);
masklayer = double(masklayer);

if nargin < 3
    saturationlayer = .8*ones(size(grayscale));
end

if max(grayscale(:)) > 1
    grayscale = (grayscale)/(max(grayscale(:)));%(grayscale+max(grayscale(:)))/(max(grayscale(:))+2);
% elseif max(grayscale(:)) > 1
%     grayscale = grayscale/255;
end

masklayer = masklayer / (max(masklayer(:))+1);

hsvimg = masklayer;
hsvimg(:,:,2) = saturationlayer;
hsvimg(:,:,3) = grayscale;

rgbimg = hsv2rgb(hsvimg);

zseq(uint8(rgbimg*255));