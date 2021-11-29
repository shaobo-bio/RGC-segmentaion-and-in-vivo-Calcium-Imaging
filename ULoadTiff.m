function [Img,Planes]= ULoadTiff(filein)
% Function for loading in tif file from computer
% file in gives location of file for reading
% Reads in for 8 or 16 bit - for other formats saves as double

ImgProps    = imfinfo(filein);
Planes      = length(ImgProps);
ImgSize     = [ImgProps(1).Height, ImgProps(1).Width];
if(ImgProps(1).BitDepth == 8)
        Img     = uint8(zeros(ImgSize(1),ImgSize(2),Planes));
elseif(ImgProps(1).BitDepth == 16)
        Img     = uint16(zeros(ImgSize(1),ImgSize(2),Planes));
else
        disp('Warning, input image is neither 8 or 16 bit')
        Img     = zeros(ImgSize(1),ImgSize(2),Slice,TimeZ);
end

for i = 1:Planes
    Img(:,:,i)  = imread(filein,i);

end
