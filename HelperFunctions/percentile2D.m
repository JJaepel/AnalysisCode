function A = percentile2D(Mat,Perc,Dim)
if length(size(Mat))~=2
    error('MAt needs to be 2D')
end
Mat=sort(Mat,Dim);
if Dim == 1
    Mat=Mat';
end
A=Mat(:,round((size(Mat,2)*(Perc/100))));