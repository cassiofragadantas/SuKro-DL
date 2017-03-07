function DEMO_image_denoise
% This DEMO is based on the image denoising demo available at the KSVDBOX
% package. 
%
%  Reference:
%  [1] M. Elad and M. Aharon, "Image Denoising via Sparse and Redundant
%      representations over Learned Dictionaries", the IEEE Trans. on Image
%      Processing, Vol. 15, no. 12, pp. 3736-3745, December 2006.
%
% The dictionary learning algorithm is exchanged by the SuKro (Sum of Kroneckers)
% dictionary learning technique (sum_separable_dict_learn.m)
%
%  [2] C.F. Dantas, M.N. da Costa and R.R. Lopes, "Learning Dictionaries as
%       a sum of Kronecker products"

%  DEMO_image_denoise reads an image, adds random white noise and denoises it
%  The input and output PSNR are compared, and the
%  trained dictionary is displayed.
%
%  To run the demo, type DEMO_image_denoise from the Matlab prompt.
%
%  See also image_denoise for the denoising algorithm.


%  ORIGINAL DEMO BY:
%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  August 2009
%
%  ADAPTED BY:
%  Cassio Fraga Dantas
%  DSPCom Laboratory - Unicamp
%  Campinas, Brasil
%  cassio@decom.fee.unicamp.br
%
%  September 2016

%% Downloading OMP and KSVD Toolboxes if necessary
FS=filesep;
pathstr = fileparts(which('DEMO_image_denoise'));
addpath([pathstr,filesep,'..']); addpath([pathstr,filesep,'..',filesep,'misc']);

KSVD_path = [pwd,FS,'toolbox'];

if ~exist('toolbox','dir')
    fprintf('\n ATTENTION: This is the first time you run this demonstration.');
    fprintf('\n Before starting the simulations, we will install some toolboxes.');
    fprintf('\n\n ****************** SETUP: Toolboxes installation ******************');
    fprintf('\n\n The following toolboxes will be downloaded:');
    fprintf('\n OMPbox version 10');
    fprintf('\n KSVDbox version 13');
    fprintf('\n\n IMPORTANT: To successfully install the toolboxes');
    fprintf('\n you will need to have MEX setup to compile C files.');
    fprintf('\n\n If this is not already setup, please type "n" to exit and then ');
    fprintf('\n run "mex -setup" or type "help mex" in the MATLAB command prompt.');
    fprintf('\n\n IMPORTANT: You must have an internet connection.');
    fprintf('\n\n ******************************************************************');

    install_ack = input('\n\n Do you wish to continue: (y/n)? ','s');

    if strcmp(install_ack,'"n"'),
      install_ack = 'n';
    end

    if install_ack == 'n',
      return;
    else
      fprintf('\n\n Installation now beginning...');

    end

    fprintf('\n ******************************************************************');
    fprintf('\n\n Initialising OMPbox and KSVDBox Setup');
    
    try
        if exist([KSVD_path, FS, 'ompbox10.zip'],'file'),
            omp_zip=[KSVD_path, FS, 'ompbox10.zip'];
        else
            omp_zip='http://www.cs.technion.ac.il/%7Eronrubin/Software/ompbox10.zip';
            fprintf('\n\n Downloading OMP toolbox, please be patient\n\n');
        end
        unzip(omp_zip,[KSVD_path, FS, 'ompbox']);
        
        cd([KSVD_path, FS, 'ompbox', FS, 'private']);
        make;
        cd(pathstr);
        
        if exist([KSVD_path, FS, 'ksvdbox13.zip'],'file'),
            KSVD_zip=[KSVD_path, FS, 'ksvdbox13.zip'];
        else
            KSVD_zip='http://www.cs.technion.ac.il/%7Eronrubin/Software/ksvdbox13.zip';
            fprintf('\n\n Downloading KSVD toolbox, please be patient\n\n');
        end
        unzip(KSVD_zip,[KSVD_path, FS, 'ksvdbox']);
        cd([KSVD_path, FS, 'ksvdbox', FS, 'private']);
        make;
        cd(pathstr);
        movefile('image_denoise.m',[KSVD_path, FS, 'ksvdbox']);
        fprintf('\n KSVDBox and OMPBox Installation Successful\n');
        fprintf('\n ******************************************************************');
        fprintf('\n\n >>> Now, please RERUN THIS SCRIPT to perform the demonstration <<< \n\n');
        return
    catch
        fprintf('\n KSVDBox and OMPBox Installation Failed');
        cd(pathstr);
    end
end
KSVD_p=genpath(KSVD_path);
addpath(KSVD_p);

%% Prompt user %%
disp(' ');
disp('***********************  SuKro Dictionary Denoising Demo  ***********************');
disp('*                                                                               *');
disp('* This demo reads an image, adds random Gaussian noise, and denoises the image  *');
disp('* using a dictionary which is the sum of a few separable terms. The function    *');
disp('* displays the original, noisy, and denoised images, and shows the resulting    *');
disp('* trained dictionary.                                                           *');
disp('*                                                                               *');
disp('*********************************************************************************');

dirname = fullfile(pathstr, 'images', '*.png');
imglist = dir(dirname);

disp('  Available test images:');
disp(' ');
for k = 1:length(imglist)
  fprintf('  %d. %s\n', k, imglist(k).name);
end
fprintf('  %d. All images', length(imglist)+1);
disp(' ');

imnum = 0;
while (~isnumeric(imnum) || (rem(imnum,1)~=0) || imnum<1 || imnum>length(imglist)+1)
  imnum = input(sprintf('  Image to denoise (%d-%d): ', 1, length(imglist)), 's');
  imnum = sscanf(imnum, '%d');
end

total_images = 1;
if imnum == length(imglist)+1, total_images = length(imglist), end

fprintf('\n\n  Choose the experiment type:\n');
fprintf('\n  1. Single-run (a few minutes)\n  2. Complete   (a few hours)\n');
exp_type = 0;
while (~isnumeric(exp_type) || (rem(exp_type,1)~=0) || exp_type<1 || exp_type>2)
  exp_type = input(sprintf('  Experiment type (1 or 2): '), 's');
  exp_type = sscanf(exp_type, '%d');
end

for image = 1:total_images 

if image > 1, imnum = image; end
if imnum == 6, imnum = 1; end
imgname = fullfile(pathstr, 'images', imglist(imnum).name);

%% Simulation parameters %%
sigma = 20; fprintf('\nSIGMA = %d\n',sigma)

samples_training = 40000; % Number of training samples
params.iternum = 100;
params.blocksize = 8;
params.dictsize = 256;

params.sigma = sigma;
params.maxval = 255;
params.memusage = 'high';

if (exp_type == 1)  % Single-run experiment
    alpha = 250;
else                % Complete experiment
    alpha = logspace(log10(100),log10(600),50); % Nuclear norm regularization coefficient
end

im = imread(imgname);
im = double(im);

% results to store
SNR = zeros(size(alpha));
nuclear_norm_Dmod = zeros(size(alpha));
rank_Dmod = zeros(size(alpha));


%% Run Simulations

for ktrain = 1:length(samples_training)
for k = 1:length(alpha)
%% Generate noisy image %%
fprintf('\n>>>>>>>>>> Running chosen parameters - %d of %d <<<<<<<<<<\n',k,length(alpha));
disp('Generating noisy image...');

n = randn(size(im)) * sigma;
imnoise = im + n;

%% Denoise %%
params.x = imnoise;
params.trainnum = samples_training(ktrain);
params.alpha = alpha(k);

[imout, dict, dict_unnorm] = image_denoise(params);

%% Show results (single-run experiment) %%
if (exp_type == 1)  % Single-run experiment
    dictimg = showdict(dict,[1 1]*params.blocksize,round(sqrt(params.dictsize)),round(sqrt(params.dictsize)),'lines','highcontrast');
    figure; ax = subplot(1,1,1); imshow(imresize(dictimg,2,'nearest'));
    title(ax,'Trained dictionary');

    figure; ax = subplot(1,1,1); imshow(im/params.maxval);
    title(ax,'Original image'); drawnow

    figure; ax = subplot(1,1,1); imshow(imnoise/params.maxval); 
    title(ax,sprintf('Noisy image, PSNR = %.2fdB', 20*log10(params.maxval * sqrt(numel(im)) / norm(im(:)-imnoise(:))) ));

    figure; ax = subplot(1,1,1); imshow(imout/params.maxval);
    title(ax,sprintf('Denoised image, PSNR: %.2fdB', 20*log10(params.maxval * sqrt(numel(im)) / norm(im(:)-imout(:))) ));
end

SNR(k) = 20*log10(params.maxval * sqrt(numel(im)) / norm(im(:)-imout(:)));

%% Saving results %%
reord_dict = reord(dict_unnorm);
nuclear_norm_Dmod(k) = sum(svd(reord_dict));
rank_Dmod(k) = rank(reord_dict,norm(dict)*2e-7);
params_light = rmfield(params,'x');
save(strcat('new_SNR_SuKro_sigma',num2str(sigma),'_',imglist(imnum).name(1:end-4),'_trainnum',num2str(params.trainnum)),'params_light','alpha','SNR','nuclear_norm_Dmod','rank_Dmod');

end
end

%% Show results (complete experiment) %%
if (exp_type == 2)  % Complete experiment
    figure, hold on
    ylabel('Recovered Image PSNR [dB]');
    xlabel('Number of Separable Terms');
    title([imglist(imnum).name, '(input PSNR = ' num2str(10*log10(255^2/(sigma^2))) ')']);

    % Averaging results with same number of separable terms 
    SNR_mean = zeros(1,params.blocksize^2);
    for kk = 1:params.blocksize^2
        if ~isempty(SNR(rank_Dmod==kk))
            SNR_mean(kk) = mean(SNR(rank_Dmod==kk));
        end
    end
    SNR_mean = SNR_mean(SNR_mean~=0);
    rank_Dmod_mean =  unique(rank_Dmod,'first');
    
    plot(rank_Dmod_mean,SNR_mean,'o');
    
    % ksvd and odct results
    pos = find([10 20  50 75 100] == sigma);
    if (pos) && (params.blocksize == 8) && (params.dictsize == 256) && (params.trainnum == 40000)
        % Pre-calculated results for a specific parameter setting
        SNR_KSVD = [34.42 33.64 35.98 35.47 34.80 ;... % sigma 10
                    30.93 30.44 33.37 32.42 32.31 ;... % sigma 20
                    25.43 25.95 28.05 27.79 28.07 ;... % sigma 50
                    22.94 23.98 25.21 25.75 25.73 ;... % sigma 75
                    21.86 22.81 23.65 24.43 24.21];   % sigma 100

        SNR_ODCT = [33.97 33.44 35.41 35.28 34.61 ;... % sigma 10
                    29.95 29.91 32.17 32.00 31.87 ;... % sigma 20
                    24.75 25.57 27.41 27.44 27.65 ;... % sigma 50
                    22.83 23.85 25.10 25.63 25.57 ;... % sigma 75
                    21.89 22.79 23.78 24.42 24.22];   % sigma 100

        plot([0 64], SNR_KSVD(pos,imnum)*ones(1,2),'k'); %ksvd
        plot([0 64], SNR_ODCT(pos,imnum)*ones(1,2),'k--'); %odct
    end
    grid on;
end

end
