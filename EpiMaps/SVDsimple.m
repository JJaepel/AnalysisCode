function [mov, mixedfilters, percent] = SVDsimple(Movie, nPCs)
%quick and dirty code for viewing and removing bad PCs from data
%AG 2020

% this line is optional if you want to hav this run faster; simply
% downsizes movie
%Movie=imresize(Movie,0.5);


if nargin < 2 || isempty(nPCs)
    nPCs=300;
end

sz=size(Movie);
npix = sz(1)*sz(2);
Movie = double(reshape(Movie, npix, sz(3)));


% Create covariance matrix
% Perform SVD on temporal covariance
c1 = (Movie'*Movie)/npix; movtm = mean(Movie,1); % Average over space
covmat = c1 - movtm'*movtm; clear c1
covtrace = trace(covmat) / npix;

%-----------------------
% Perform SVD
[mixedsig, CovEvals] = eigs(covmat, nPCs, 'LM');  % pca_mixedsig are the temporal signals, mixedsig
CovEvals = diag(CovEvals);
mixedsig = mixedsig' * sz(3);
CovEvals = CovEvals / npix;

%Report percent variance per PC
percentvar = 100*sum(CovEvals)/covtrace;
fprintf([' First ',num2str(nPCs),' PCs contain ',num2str(percentvar,3),'%% of the variance.\n'])
percent=CovEvals/covtrace;
figure; plot(percent(1:30))

% Re-load movie data to view PCs
Sinv = inv(diag(CovEvals.^(1/2)));
movtm = mean(Movie,1); % Average over space
movuse = Movie - ones(npix,1) * movtm;
mixedfilters = reshape(movuse * mixedsig' * Sinv, npix, nPCs);
mixedfilters = reshape(mixedfilters, sz(1),sz(2),nPCs);

implay(mixedfilters);

badPCs = input('bad PCs= (example "[1:2]")'); 
PCuse=setdiff([1:nPCs],badPCs); 
mixedfilters2 = reshape(mixedfilters(:,:,PCuse),npix,length(PCuse));  
mov = mixedfilters2 * diag(CovEvals(PCuse).^(1/2)) * mixedsig(PCuse,:);  
mov = zscore(reshape(mov,npix*sz(3),1));
mov = reshape(mov, sz); 
end