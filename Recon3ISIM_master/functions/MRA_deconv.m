function [f] = MRA_deconv(g,lambda,pixelsize,NA,MRAreg,Mode)
%% 2D-MRA deconvolution
        w2=size(g,1);
        w1=size(g,2);
        [X,Y]=meshgrid(linspace(0,w1-1,w1),linspace(0,w2-1,w2));
        scale=2*pi*NA/lambda*pixelsize;
        R=sqrt(min(X,abs(X-w1)).^2+min(Y,abs(Y-w2)).^2);
        PSF=abs(2*besselj(1,scale*R+eps,1)./(scale*R+eps)).^2;
        PSF = PSF/(sum(sum(PSF)));
        PSF=fftshift(PSF);
%% Pre-process
    g=single(g);
    [d1,d2]=size(PSF);
    PSF=PSF/sum(sum(PSF));
    d=min(d1,d2);
    center=[fix(d/2)+1,fix(d/2)+1];
%% Definition of the waveform
    wav1=@(U)FrameletDec(U,2);
    iwav1=@(C)FrameletRec(C,2);
    wav2 = @(X)CurveletDec(X,1,2); 
    iwav2 =@(C)real(CurveletRec(C,1));
    if Mode==2
        Maxiteration=100;
    else
    Maxiteration=50;
    end
    m=max(max(g));
    mi=min(min(g));
    g=(g-mi)/(m-mi);

    [f,fun_all]=MRAFISTA_Core(g,PSF,center,wav1,iwav1,wav2,iwav2,MRAreg,MRAreg,Maxiteration,0,0);
    if Mode==3
    f = deconvlucy(f,PSF,2);
    else
      f = deconvlucy(f,PSF,5);   
    end
    m=max(max(f));
    mi=min(min(f));
    f=2^16*(f-mi)/(m-mi);

end
function [f,fun_all]=MRAFISTA_Core(g,PSF,center,W1,WT1,W2,WT2,lambda1,lambda2,Maxiteration,flag1,flag11)
% Assigning parameters according to pars and/or default values
    [m,n]=size(g);
    PSF=padPSF(PSF,m,n);
    trans=@(X) 1/sqrt(m*n)*fft2(X);
    itrans=@(X) sqrt(m*n)*ifft2(X);     
    Sbig=fft2(circshift(PSF,1-center));

% computing the two dimensional transform of Bobs
    Btrans=trans(g);
%The Lipschitz constant of the gradient of ||A(X)-Bobs||^2
    L=2*max(max(abs(Sbig).^2));
% initialization
    f_iter=g;
    Y=f_iter;
    t_new=1;

%Rationalize the soft-thresholding weights
        interval=10;

fun_all=norm(Sbig.*trans(f_iter)-Btrans,'fro')^2;
for i=1:Maxiteration
    % Store the old value of the iterate and the t-constant
    f_old=f_iter;
    t_old=t_new;
    % Gradient step
    D=Sbig.*trans(Y)-Btrans;
    Y=Y-2/L*itrans(conj(Sbig).*D);
    Y=real(Y);  
    f_iter=Y;

 %% Transform1 soft-thresholding
 if i==1||mod(i,interval)==0
    % Framelet transform
    WX=W1(Y);
    % Soft thresholding 
    [~,m]=size(WX);
    for ii=1:m
            tmp1=WX{1,ii};
            [dx,dy]=size(tmp1);
            for jj=1:dx
                for kk=1:dy
                tmp=tmp1{jj,kk};
    D=abs(tmp)-lambda1/(L);
    tmp=sign(tmp).*((D>0).*D);
    WX{1,ii}{jj,kk}=tmp;
                end
            end       
    end
    % The new iterate inverse wavelet transform of WY
    f_iter=WT1(WX);
%% Transform2 soft-thresholding
    % Curvelet transform
    WY=W2(f_iter);
    % Soft thresholding 
    [m,n]=size(WY);
    for ii=1:m
        for jj=1:n
            tmp1=WY{ii,jj};
            z=length(tmp1);
            for kk=1:z
                tmp=tmp1{1,kk};
    D=abs(tmp)-lambda2/(L);
    tmp=sign(tmp).*((D>0).*D);
    WY{ii,jj}{1,kk}=tmp;
            end
        end
    end
    % The new iterate inverse wavelet transform of WY
    f_iter=WT2(WY);
    %updating t and Y
    t_new=(1+sqrt(1+4*t_old^2))/2;
    Y=f_iter+(t_old-1)/t_new*(f_iter-f_old);

 end
end
fun_all=[fun_all,norm(Sbig.*trans(f_iter)-Btrans,'fro')^2];
f_iter(f_iter<0)=0;
f=f_iter;
end

%% 2D Curvelet
function [wl,wr] = fdct_wrapping_window(x)
wr = zeros(size(x));
wl = zeros(size(x));
x(abs(x) < 2^-52) = 0;
wr((x > 0) & (x < 1)) = exp(1-1./(1-exp(1-1./x((x > 0) & (x < 1)))));
wr(x <= 0) = 1;
wl((x > 0) & (x < 1)) = exp(1-1./(1-exp(1-1./(1-x((x > 0) & (x < 1))))));
wl(x >= 1) = 1;
normalization = sqrt(wl.^2 + wr.^2);
wr = wr ./ normalization;
wl = wl ./ normalization;
end

function C = CurveletDec(x, is_real, finest, nbscales, nbangles_coarse)
X = fftshift(fft2(ifftshift(x)))/sqrt(prod(size(x)));
[N1,N2] = size(X);
if nargin < 2, is_real = 0; end;
if nargin < 3, finest = 2; end;
if nargin < 4, nbscales = ceil(log2(min(N1,N2)) - 3); end;
if nargin < 5, nbangles_coarse = 16; end;

% Initialization: data structure
nbangles = [1, nbangles_coarse .* 2.^(ceil((nbscales-(nbscales:-1:2))/2))];
if finest == 2, nbangles(nbscales) = 1; end;
C = cell(1,nbscales);
for j = 1:nbscales
    C{j} = cell(1,nbangles(j));
end;

% Loop: pyramidal scale decomposition
M1 = N1/3;
M2 = N2/3;
if finest == 1,

    % Initialization: smooth periodic extension of high frequencies
    bigN1 = 2*floor(2*M1)+1;
    bigN2 = 2*floor(2*M2)+1;
    equiv_index_1 = 1+mod(floor(N1/2)-floor(2*M1)+(1:bigN1)-1,N1);
    equiv_index_2 = 1+mod(floor(N2/2)-floor(2*M2)+(1:bigN2)-1,N2);
    X = X(equiv_index_1,equiv_index_2);
        % Invariant: equiv_index_1(floor(2*M1)+1) == (N1 + 2 - mod(N1,2))/2
        % is the center in frequency. Same for M2, N2.
    window_length_1 = floor(2*M1) - floor(M1) - 1 - (mod(N1,3)==0);
    window_length_2 = floor(2*M2) - floor(M2) - 1 - (mod(N2,3)==0);
        % Invariant: floor(M1) + floor(2*M1) == N1 - (mod(M1,3)~=0)
        % Same for M2, N2.
    coord_1 = 0:(1/window_length_1):1;
    coord_2 = 0:(1/window_length_2):1;
    [wl_1,wr_1] = fdct_wrapping_window(coord_1);
    [wl_2,wr_2] = fdct_wrapping_window(coord_2);
    lowpass_1 = [wl_1, ones(1,2*floor(M1)+1), wr_1];
    if mod(N1,3)==0, lowpass_1 = [0, lowpass_1, 0]; end;
    lowpass_2 = [wl_2, ones(1,2*floor(M2)+1), wr_2];
    if mod(N2,3)==0, lowpass_2 = [0, lowpass_2, 0]; end;
    lowpass = lowpass_1'*lowpass_2;
    Xlow = X .* lowpass;

    scales = nbscales:-1:2;

else
    
    M1 = M1/2;
    M2 = M2/2;
    window_length_1 = floor(2*M1) - floor(M1) - 1;
    window_length_2 = floor(2*M2) - floor(M2) - 1;
    coord_1 = 0:(1/window_length_1):1;
    coord_2 = 0:(1/window_length_2):1;
    [wl_1,wr_1] = fdct_wrapping_window(coord_1);
    [wl_2,wr_2] = fdct_wrapping_window(coord_2);
    lowpass_1 = [wl_1, ones(1,2*floor(M1)+1), wr_1];
    lowpass_2 = [wl_2, ones(1,2*floor(M2)+1), wr_2];
    lowpass = lowpass_1'*lowpass_2;
    hipass = sqrt(1 - lowpass.^2);
    Xlow_index_1 = ((-floor(2*M1)):floor(2*M1)) + ceil((N1+1)/2);
    Xlow_index_2 = ((-floor(2*M2)):floor(2*M2)) + ceil((N2+1)/2);
    Xlow = X(Xlow_index_1, Xlow_index_2) .* lowpass;
    Xhi = X;
    Xhi(Xlow_index_1, Xlow_index_2) = Xhi(Xlow_index_1, Xlow_index_2) .* hipass;
    C{nbscales}{1} = fftshift(ifft2(ifftshift(Xhi)))*sqrt(prod(size(Xhi)));
    if is_real, C{nbscales}{1} = real(C{nbscales}{1}); end;
    
    scales = (nbscales-1):-1:2;

end;
for j = scales,

    M1 = M1/2;
    M2 = M2/2;
    window_length_1 = floor(2*M1) - floor(M1) - 1;
    window_length_2 = floor(2*M2) - floor(M2) - 1;
    coord_1 = 0:(1/window_length_1):1;
    coord_2 = 0:(1/window_length_2):1;
    [wl_1,wr_1] = fdct_wrapping_window(coord_1);
    [wl_2,wr_2] = fdct_wrapping_window(coord_2);
    lowpass_1 = [wl_1, ones(1,2*floor(M1)+1), wr_1];
    lowpass_2 = [wl_2, ones(1,2*floor(M2)+1), wr_2];
    lowpass = lowpass_1'*lowpass_2;
    hipass = sqrt(1 - lowpass.^2);
    Xhi = Xlow;                 % size is 2*floor(4*M1)+1 - by - 2*floor(4*M2)+1
    Xlow_index_1 = ((-floor(2*M1)):floor(2*M1)) + floor(4*M1) + 1;
    Xlow_index_2 = ((-floor(2*M2)):floor(2*M2)) + floor(4*M2) + 1;
    Xlow = Xlow(Xlow_index_1, Xlow_index_2);
    Xhi(Xlow_index_1, Xlow_index_2) = Xlow .* hipass;
    Xlow = Xlow .* lowpass;     % size is 2*floor(2*M1)+1 - by - 2*floor(2*M2)+1
    
    % Loop: angular decomposition
    l = 0;
    nbquadrants = 2 + 2*(~is_real);
    nbangles_perquad = nbangles(j)/4;
    for quadrant = 1:nbquadrants
        M_horiz = M2 * (mod(quadrant,2)==1) + M1 * (mod(quadrant,2)==0);
        M_vert = M1 * (mod(quadrant,2)==1) + M2 * (mod(quadrant,2)==0);
        if mod(nbangles_perquad,2),
            wedge_ticks_left = round((0:(1/(2*nbangles_perquad)):.5)*2*floor(4*M_horiz) + 1);
            wedge_ticks_right = 2*floor(4*M_horiz) + 2 - wedge_ticks_left;
            wedge_ticks = [wedge_ticks_left, wedge_ticks_right(end:-1:1)];
        else
            wedge_ticks_left = round((0:(1/(2*nbangles_perquad)):.5)*2*floor(4*M_horiz) + 1);
            wedge_ticks_right = 2*floor(4*M_horiz) + 2 - wedge_ticks_left;
            wedge_ticks = [wedge_ticks_left, wedge_ticks_right((end-1):-1:1)];
        end;
        wedge_endpoints = wedge_ticks(2:2:(end-1));         % integers
        wedge_midpoints = (wedge_endpoints(1:(end-1)) + wedge_endpoints(2:end))/2;
                % integers or half-integers
        
        % Left corner wedge
        l = l+1;
        first_wedge_endpoint_vert = round(2*floor(4*M_vert)/(2*nbangles_perquad) + 1);
        length_corner_wedge = floor(4*M_vert) - floor(M_vert) + ceil(first_wedge_endpoint_vert/4);
        Y_corner = 1:length_corner_wedge;
        [XX,YY] = meshgrid(1:(2*floor(4*M_horiz)+1),Y_corner);
        width_wedge = wedge_endpoints(2) + wedge_endpoints(1) - 1;
        slope_wedge = (floor(4*M_horiz) + 1 - wedge_endpoints(1))/floor(4*M_vert);
        left_line = round(2 - wedge_endpoints(1) + slope_wedge*(Y_corner - 1));
                                                            % integers
        [wrapped_data, wrapped_XX, wrapped_YY] = deal(zeros(length_corner_wedge,width_wedge));
        first_row = floor(4*M_vert)+2-ceil((length_corner_wedge+1)/2)+...
            mod(length_corner_wedge+1,2)*(quadrant-2 == mod(quadrant-2,2));
        first_col = floor(4*M_horiz)+2-ceil((width_wedge+1)/2)+...
            mod(width_wedge+1,2)*(quadrant-3 == mod(quadrant-3,2));
                % Coordinates of the top-left corner of the wedge wrapped
                % around the origin. Some subtleties when the wedge is
                % even-sized because of the forthcoming 90 degrees rotation
        for row = Y_corner
            cols = left_line(row) + mod((0:(width_wedge-1))-(left_line(row)-first_col),width_wedge);
            admissible_cols = round(1/2*(cols+1+abs(cols-1)));
            new_row = 1 + mod(row - first_row, length_corner_wedge);
            wrapped_data(new_row,:) = Xhi(row,admissible_cols) .* (cols > 0);
            wrapped_XX(new_row,:) = XX(row,admissible_cols);
            wrapped_YY(new_row,:) = YY(row,admissible_cols);
        end;
        slope_wedge_right = (floor(4*M_horiz)+1 - wedge_midpoints(1))/floor(4*M_vert);
        mid_line_right = wedge_midpoints(1) + slope_wedge_right*(wrapped_YY - 1);
                % not integers in general
        coord_right = 1/2 + floor(4*M_vert)/(wedge_endpoints(2) - wedge_endpoints(1)) * ...
            (wrapped_XX - mid_line_right)./(floor(4*M_vert)+1 - wrapped_YY);
        C2 = 1/(1/(2*(floor(4*M_horiz))/(wedge_endpoints(1) - 1) - 1) + 1/(2*(floor(4*M_vert))/(first_wedge_endpoint_vert - 1) - 1));
        C1 = C2 / (2*(floor(4*M_vert))/(first_wedge_endpoint_vert - 1) - 1);
        wrapped_XX((wrapped_XX - 1)/floor(4*M_horiz) + (wrapped_YY-1)/floor(4*M_vert) == 2) = ...
            wrapped_XX((wrapped_XX - 1)/floor(4*M_horiz) + (wrapped_YY-1)/floor(4*M_vert) == 2) + 1;
        coord_corner = C1 + C2 * ((wrapped_XX - 1)/(floor(4*M_horiz)) - (wrapped_YY - 1)/(floor(4*M_vert))) ./ ...
            (2-((wrapped_XX - 1)/(floor(4*M_horiz)) + (wrapped_YY - 1)/(floor(4*M_vert))));
        wl_left = fdct_wrapping_window(coord_corner);
        [wl_right,wr_right] = fdct_wrapping_window(coord_right);
        wrapped_data = wrapped_data .* (wl_left .* wr_right);

        switch is_real
            case 0
                wrapped_data = rot90(wrapped_data,-(quadrant-1));
                C{j}{l} = fftshift(ifft2(ifftshift(wrapped_data)))*sqrt(prod(size(wrapped_data)));
            case 1
                wrapped_data = rot90(wrapped_data,-(quadrant-1));
                x = fftshift(ifft2(ifftshift(wrapped_data)))*sqrt(prod(size(wrapped_data)));
                C{j}{l} = sqrt(2)*real(x);
                C{j}{l+nbangles(j)/2} = sqrt(2)*imag(x);
        end;
                
        % Regular wedges
        length_wedge = floor(4*M_vert) - floor(M_vert);
        Y = 1:length_wedge;
        first_row = floor(4*M_vert)+2-ceil((length_wedge+1)/2)+...
            mod(length_wedge+1,2)*(quadrant-2 == mod(quadrant-2,2));
        for subl = 2:(nbangles_perquad-1);
            l = l+1;
            width_wedge = wedge_endpoints(subl+1) - wedge_endpoints(subl-1) + 1;
            slope_wedge = ((floor(4*M_horiz)+1) - wedge_endpoints(subl))/floor(4*M_vert);
            left_line = round(wedge_endpoints(subl-1) + slope_wedge*(Y - 1));
            [wrapped_data, wrapped_XX, wrapped_YY] = deal(zeros(length_wedge,width_wedge));
            first_col = floor(4*M_horiz)+2-ceil((width_wedge+1)/2)+...
                mod(width_wedge+1,2)*(quadrant-3 == mod(quadrant-3,2));
            for row = Y
                cols = left_line(row) + mod((0:(width_wedge-1))-(left_line(row)-first_col),width_wedge);
                new_row = 1 + mod(row - first_row, length_wedge);
                wrapped_data(new_row,:) = Xhi(row,cols);
                wrapped_XX(new_row,:) = XX(row,cols);
                wrapped_YY(new_row,:) = YY(row,cols);             
            end;
            slope_wedge_left = ((floor(4*M_horiz)+1) - wedge_midpoints(subl-1))/floor(4*M_vert);
            mid_line_left = wedge_midpoints(subl-1) + slope_wedge_left*(wrapped_YY - 1);
            coord_left = 1/2 + floor(4*M_vert)/(wedge_endpoints(subl) - wedge_endpoints(subl-1)) * ...
                (wrapped_XX - mid_line_left)./(floor(4*M_vert)+1 - wrapped_YY);
            slope_wedge_right = ((floor(4*M_horiz)+1) - wedge_midpoints(subl))/floor(4*M_vert);
            mid_line_right = wedge_midpoints(subl) + slope_wedge_right*(wrapped_YY - 1);
            coord_right = 1/2 + floor(4*M_vert)/(wedge_endpoints(subl+1) - wedge_endpoints(subl)) * ...
                (wrapped_XX - mid_line_right)./(floor(4*M_vert)+1 - wrapped_YY);
            wl_left = fdct_wrapping_window(coord_left);
            [wl_right,wr_right] = fdct_wrapping_window(coord_right);
            wrapped_data = wrapped_data .* (wl_left .* wr_right);
            switch is_real
                case 0
                    wrapped_data = rot90(wrapped_data,-(quadrant-1));
                    C{j}{l} = fftshift(ifft2(ifftshift(wrapped_data)))*sqrt(prod(size(wrapped_data)));
                case 1
                    wrapped_data = rot90(wrapped_data,-(quadrant-1));
                    x = fftshift(ifft2(ifftshift(wrapped_data)))*sqrt(prod(size(wrapped_data)));
                    C{j}{l} = sqrt(2)*real(x);
                    C{j}{l+nbangles(j)/2} = sqrt(2)*imag(x);
            end;
        end;

        % Right corner wedge
        l = l+1;
        width_wedge = 4*floor(4*M_horiz) + 3 - wedge_endpoints(end) - wedge_endpoints(end-1);
        slope_wedge = ((floor(4*M_horiz)+1) - wedge_endpoints(end))/floor(4*M_vert);
        left_line = round(wedge_endpoints(end-1) + slope_wedge*(Y_corner - 1));
        [wrapped_data, wrapped_XX, wrapped_YY] = deal(zeros(length_corner_wedge,width_wedge));
        first_row = floor(4*M_vert)+2-ceil((length_corner_wedge+1)/2)+...
            mod(length_corner_wedge+1,2)*(quadrant-2 == mod(quadrant-2,2));
        first_col = floor(4*M_horiz)+2-ceil((width_wedge+1)/2)+...
            mod(width_wedge+1,2)*(quadrant-3 == mod(quadrant-3,2));
        for row = Y_corner
            cols = left_line(row) + mod((0:(width_wedge-1))-(left_line(row)-first_col),width_wedge);
            admissible_cols = round(1/2*(cols+2*floor(4*M_horiz)+1-abs(cols-(2*floor(4*M_horiz)+1))));
            new_row = 1 + mod(row - first_row, length_corner_wedge);
            wrapped_data(new_row,:) = Xhi(row,admissible_cols) .* (cols <= (2*floor(4*M_horiz)+1));
            wrapped_XX(new_row,:) = XX(row,admissible_cols);
            wrapped_YY(new_row,:) = YY(row,admissible_cols);
        end;
        slope_wedge_left = ((floor(4*M_horiz)+1) - wedge_midpoints(end))/floor(4*M_vert);
        mid_line_left = wedge_midpoints(end) + slope_wedge_left*(wrapped_YY - 1);
        coord_left = 1/2 + floor(4*M_vert)/(wedge_endpoints(end) - wedge_endpoints(end-1)) * ...
            (wrapped_XX - mid_line_left)./(floor(4*M_vert) + 1 - wrapped_YY);
        C2 = -1/(2*(floor(4*M_horiz))/(wedge_endpoints(end) - 1) - 1 + 1/(2*(floor(4*M_vert))/(first_wedge_endpoint_vert - 1) - 1));
        C1 = -C2 * (2*(floor(4*M_horiz))/(wedge_endpoints(end) - 1) - 1);
        wrapped_XX((wrapped_XX - 1)/floor(4*M_horiz) == (wrapped_YY - 1)/floor(4*M_vert)) = ...
            wrapped_XX((wrapped_XX - 1)/floor(4*M_horiz) == (wrapped_YY - 1)/floor(4*M_vert)) - 1;
        coord_corner = C1 + C2 * (2-((wrapped_XX - 1)/(floor(4*M_horiz)) + (wrapped_YY - 1)/(floor(4*M_vert)))) ./ ...
            ((wrapped_XX - 1)/(floor(4*M_horiz)) - (wrapped_YY - 1)/(floor(4*M_vert)));
        wl_left = fdct_wrapping_window(coord_left);
        [wl_right,wr_right] = fdct_wrapping_window(coord_corner);

        wrapped_data = wrapped_data .* (wl_left .* wr_right);
        switch is_real
            case 0
                wrapped_data = rot90(wrapped_data,-(quadrant-1));
                C{j}{l} = fftshift(ifft2(ifftshift(wrapped_data)))*sqrt(prod(size(wrapped_data)));
            case 1
                wrapped_data = rot90(wrapped_data,-(quadrant-1));
                x = fftshift(ifft2(ifftshift(wrapped_data)))*sqrt(prod(size(wrapped_data)));
                C{j}{l} = sqrt(2)*real(x);
                C{j}{l+nbangles(j)/2} = sqrt(2)*imag(x);
        end;

        if quadrant < nbquadrants, Xhi = rot90(Xhi); end;
    end;
end;

% Coarsest wavelet level
C{1}{1} = fftshift(ifft2(ifftshift(Xlow)))*sqrt(prod(size(Xlow)));
if is_real == 1,
    C{1}{1} = real(C{1}{1});
end;
end

function x = CurveletRec(C, is_real, M, N)
% Initialization
nbscales = length(C);
nbangles_coarse = length(C{2});
nbangles = [1, nbangles_coarse .* 2.^(ceil((nbscales-(nbscales:-1:2))/2))];
if length(C{end}) == 1, finest = 2; else finest = 1; end;
if finest == 2, nbangles(nbscales) = 1; end;
if nargin < 2, is_real = 0; end;
if nargin < 4,
    if finest == 1, error('Syntax: IFCT_wrapping(C,M,N) where the matrix to be recovered is M-by-N'); end;
    [N1,N2] = size(C{end}{1});
else
    N1 = M;
    N2 = N;
end;

M1 = N1/3;
M2 = N2/3;

if finest == 1;
    
    bigN1 = 2*floor(2*M1)+1;
    bigN2 = 2*floor(2*M2)+1;
    X = zeros(bigN1,bigN2);

    % Initialization: preparing the lowpass filter at finest scale
    window_length_1 = floor(2*M1) - floor(M1) - 1 - (mod(N1,3)==0);
    window_length_2 = floor(2*M2) - floor(M2) - 1 - (mod(N2,3)==0);
    coord_1 = 0:(1/window_length_1):1;
    coord_2 = 0:(1/window_length_2):1;
    [wl_1,wr_1] = fdct_wrapping_window(coord_1);
    [wl_2,wr_2] = fdct_wrapping_window(coord_2);
    lowpass_1 = [wl_1, ones(1,2*floor(M1)+1), wr_1];
    if mod(N1,3)==0, lowpass_1 = [0, lowpass_1, 0]; end;
    lowpass_2 = [wl_2, ones(1,2*floor(M2)+1), wr_2];
    if mod(N2,3)==0, lowpass_2 = [0, lowpass_2, 0]; end;
    lowpass = lowpass_1'*lowpass_2;

    scales = nbscales:-1:2;
   
else

    M1 = M1/2;
    M2 = M2/2;
    
    bigN1 = 2*floor(2*M1)+1;
    bigN2 = 2*floor(2*M2)+1;
    X = zeros(bigN1,bigN2);
    
    window_length_1 = floor(2*M1) - floor(M1) - 1;
    window_length_2 = floor(2*M2) - floor(M2) - 1;
    coord_1 = 0:(1/window_length_1):1;
    coord_2 = 0:(1/window_length_2):1;
    [wl_1,wr_1] = fdct_wrapping_window(coord_1);
    [wl_2,wr_2] = fdct_wrapping_window(coord_2);
    lowpass_1 = [wl_1, ones(1,2*floor(M1)+1), wr_1];
    lowpass_2 = [wl_2, ones(1,2*floor(M2)+1), wr_2];
    lowpass = lowpass_1'*lowpass_2;
    hipass_finest = sqrt(1 - lowpass.^2);
    
    scales = (nbscales-1):-1:2;
    
end;

% Loop: pyramidal reconstruction

Xj_topleft_1 = 1;
Xj_topleft_2 = 1;
for j = scales,

    M1 = M1/2;
    M2 = M2/2;
    window_length_1 = floor(2*M1) - floor(M1) - 1;
    window_length_2 = floor(2*M2) - floor(M2) - 1;
    coord_1 = 0:(1/window_length_1):1;
    coord_2 = 0:(1/window_length_2):1;
    [wl_1,wr_1] = fdct_wrapping_window(coord_1);
    [wl_2,wr_2] = fdct_wrapping_window(coord_2);
    lowpass_1 = [wl_1, ones(1,2*floor(M1)+1), wr_1];
    lowpass_2 = [wl_2, ones(1,2*floor(M2)+1), wr_2];
    lowpass_next = lowpass_1'*lowpass_2;
    hipass = sqrt(1 - lowpass_next.^2);
    Xj = zeros(2*floor(4*M1)+1,2*floor(4*M2)+1);
    
    % Loop: angles
    l = 0;
    nbquadrants = 2 + 2*(~is_real);
    nbangles_perquad = nbangles(j)/4;
    for quadrant = 1:nbquadrants
        
        M_horiz = M2 * (mod(quadrant,2)==1) + M1 * (mod(quadrant,2)==0);
        M_vert = M1 * (mod(quadrant,2)==1) + M2 * (mod(quadrant,2)==0);
        if mod(nbangles_perquad,2),
            wedge_ticks_left = round((0:(1/(2*nbangles_perquad)):.5)*2*floor(4*M_horiz) + 1);
            wedge_ticks_right = 2*floor(4*M_horiz) + 2 - wedge_ticks_left;
            wedge_ticks = [wedge_ticks_left, wedge_ticks_right(end:-1:1)];
        else
            wedge_ticks_left = round((0:(1/(2*nbangles_perquad)):.5)*2*floor(4*M_horiz) + 1);
            wedge_ticks_right = 2*floor(4*M_horiz) + 2 - wedge_ticks_left;
            wedge_ticks = [wedge_ticks_left, wedge_ticks_right((end-1):-1:1)];
        end;
        wedge_endpoints = wedge_ticks(2:2:(end-1));         % integers
        wedge_midpoints = (wedge_endpoints(1:(end-1)) + wedge_endpoints(2:end))/2;
        
        % Left corner wedge
        
        l = l+1;
        first_wedge_endpoint_vert = round(2*floor(4*M_vert)/(2*nbangles_perquad) + 1);
        length_corner_wedge = floor(4*M_vert) - floor(M_vert) + ceil(first_wedge_endpoint_vert/4);
        Y_corner = 1:length_corner_wedge;
        [XX,YY] = meshgrid(1:(2*floor(4*M_horiz)+1),Y_corner);
        width_wedge = wedge_endpoints(2) + wedge_endpoints(1) - 1;
        slope_wedge = (floor(4*M_horiz) + 1 - wedge_endpoints(1))/floor(4*M_vert);
        left_line = round(2 - wedge_endpoints(1) + slope_wedge*(Y_corner - 1));
        [wrapped_XX, wrapped_YY] = deal(zeros(length_corner_wedge,width_wedge));
        first_row = floor(4*M_vert)+2-ceil((length_corner_wedge+1)/2)+...
            mod(length_corner_wedge+1,2)*(quadrant-2 == mod(quadrant-2,2));
        first_col = floor(4*M_horiz)+2-ceil((width_wedge+1)/2)+...
            mod(width_wedge+1,2)*(quadrant-3 == mod(quadrant-3,2));
        
        for row = Y_corner
            cols = left_line(row) + mod((0:(width_wedge-1))-(left_line(row)-first_col),width_wedge);
            new_row = 1 + mod(row - first_row, length_corner_wedge);
            admissible_cols = round(1/2*(cols+1+abs(cols-1)));
            wrapped_XX(new_row,:) = XX(row,admissible_cols);
            wrapped_YY(new_row,:) = YY(row,admissible_cols);
        end;

        slope_wedge_right = (floor(4*M_horiz)+1 - wedge_midpoints(1))/floor(4*M_vert);
        mid_line_right = wedge_midpoints(1) + slope_wedge_right*(wrapped_YY - 1);
                                                            % not integers
                                                            % in general
        coord_right = 1/2 + floor(4*M_vert)/(wedge_endpoints(2) - wedge_endpoints(1)) * ...
            (wrapped_XX - mid_line_right)./(floor(4*M_vert)+1 - wrapped_YY);
        C2 = 1/(1/(2*(floor(4*M_horiz))/(wedge_endpoints(1) - 1) - 1) + 1/(2*(floor(4*M_vert))/(first_wedge_endpoint_vert - 1) - 1));
        C1 = C2 / (2*(floor(4*M_vert))/(first_wedge_endpoint_vert - 1) - 1);
        wrapped_XX((wrapped_XX - 1)/floor(4*M_horiz) + (wrapped_YY-1)/floor(4*M_vert) == 2) = ...
            wrapped_XX((wrapped_XX - 1)/floor(4*M_horiz) + (wrapped_YY-1)/floor(4*M_vert) == 2) + 1;
        coord_corner = C1 + C2 * ((wrapped_XX - 1)/(floor(4*M_horiz)) - (wrapped_YY - 1)/(floor(4*M_vert))) ./ ...
            (2-((wrapped_XX - 1)/(floor(4*M_horiz)) + (wrapped_YY - 1)/(floor(4*M_vert))));
        wl_left = fdct_wrapping_window(coord_corner);
        [wl_right,wr_right] = fdct_wrapping_window(coord_right);
        switch is_real
         case 0
          wrapped_data = fftshift(fft2(ifftshift(C{j}{l})))/sqrt(prod(size(C{j}{l})));
          wrapped_data = rot90(wrapped_data,(quadrant-1));
         case 1
          x = C{j}{l} + sqrt(-1)*C{j}{l+nbangles(j)/2};
          wrapped_data = fftshift(fft2(ifftshift(x)))/sqrt(prod(size(x)))/sqrt(2);
          wrapped_data = rot90(wrapped_data,(quadrant-1));
        end;
        wrapped_data = wrapped_data .* (wl_left .* wr_right);
 
        % Unwrapping data
        for row = Y_corner
            cols = left_line(row) + mod((0:(width_wedge-1))-(left_line(row)-first_col),width_wedge);
            admissible_cols = round(1/2*(cols+1+abs(cols-1)));
            new_row = 1 + mod(row - first_row, length_corner_wedge);
            Xj(row,admissible_cols) = Xj(row,admissible_cols) + wrapped_data(new_row,:);
                                % We use the following property: in an assignment
                                % A(B) = C where B and C are vectors, if
                                % some value x repeats in B, then the
                                % last occurrence of x is the one
                                % corresponding to the eventual assignment.
        end;

        % Regular wedges
        length_wedge = floor(4*M_vert) - floor(M_vert);
        Y = 1:length_wedge;
        first_row = floor(4*M_vert)+2-ceil((length_wedge+1)/2)+...
            mod(length_wedge+1,2)*(quadrant-2 == mod(quadrant-2,2));
        for subl = 2:(nbangles_perquad-1);
            l = l+1;
            width_wedge = wedge_endpoints(subl+1) - wedge_endpoints(subl-1) + 1;
            slope_wedge = ((floor(4*M_horiz)+1) - wedge_endpoints(subl))/floor(4*M_vert);
            left_line = round(wedge_endpoints(subl-1) + slope_wedge*(Y - 1));
            [wrapped_XX, wrapped_YY] = deal(zeros(length_wedge,width_wedge));
            first_col = floor(4*M_horiz)+2-ceil((width_wedge+1)/2)+...
                mod(width_wedge+1,2)*(quadrant-3 == mod(quadrant-3,2));
            for row = Y
                cols = left_line(row) + mod((0:(width_wedge-1))-(left_line(row)-first_col),width_wedge);
                new_row = 1 + mod(row - first_row, length_wedge);
                wrapped_XX(new_row,:) = XX(row,cols);
                wrapped_YY(new_row,:) = YY(row,cols);
            end;
            slope_wedge_left = ((floor(4*M_horiz)+1) - wedge_midpoints(subl-1))/floor(4*M_vert);
            mid_line_left = wedge_midpoints(subl-1) + slope_wedge_left*(wrapped_YY - 1);
            coord_left = 1/2 + floor(4*M_vert)/(wedge_endpoints(subl) - wedge_endpoints(subl-1)) * ...
                (wrapped_XX - mid_line_left)./(floor(4*M_vert)+1 - wrapped_YY);
            slope_wedge_right = ((floor(4*M_horiz)+1) - wedge_midpoints(subl))/floor(4*M_vert);
            mid_line_right = wedge_midpoints(subl) + slope_wedge_right*(wrapped_YY - 1);
            coord_right = 1/2 + floor(4*M_vert)/(wedge_endpoints(subl+1) - wedge_endpoints(subl)) * ...
                (wrapped_XX - mid_line_right)./(floor(4*M_vert)+1 - wrapped_YY);
            wl_left = fdct_wrapping_window(coord_left);
            [wl_right,wr_right] = fdct_wrapping_window(coord_right);
            switch is_real
             case 0
              wrapped_data = fftshift(fft2(ifftshift(C{j}{l})))/sqrt(prod(size(C{j}{l})));
              wrapped_data = rot90(wrapped_data,(quadrant-1));
             case 1
              x = C{j}{l} + sqrt(-1)*C{j}{l+nbangles(j)/2};
              wrapped_data = fftshift(fft2(ifftshift(x)))/sqrt(prod(size(x)))/sqrt(2);
              wrapped_data = rot90(wrapped_data,(quadrant-1));
            end;
            wrapped_data = wrapped_data .* (wl_left .* wr_right);
            
            % Unwrapping data
            for row = Y
                cols = left_line(row) + mod((0:(width_wedge-1))-(left_line(row)-first_col),width_wedge);
                new_row = 1 + mod(row - first_row, length_wedge);
                Xj(row,cols) = Xj(row,cols) + wrapped_data(new_row,:);
            end;

        end;    % for subl
        
        % Right corner wedge
        l = l+1;
        width_wedge = 4*floor(4*M_horiz) + 3 - wedge_endpoints(end) - wedge_endpoints(end-1);
        slope_wedge = ((floor(4*M_horiz)+1) - wedge_endpoints(end))/floor(4*M_vert);
        left_line = round(wedge_endpoints(end-1) + slope_wedge*(Y_corner - 1));
        [wrapped_XX, wrapped_YY] = deal(zeros(length_corner_wedge,width_wedge));
        first_row = floor(4*M_vert)+2-ceil((length_corner_wedge+1)/2)+...
            mod(length_corner_wedge+1,2)*(quadrant-2 == mod(quadrant-2,2));
        first_col = floor(4*M_horiz)+2-ceil((width_wedge+1)/2)+...
            mod(width_wedge+1,2)*(quadrant-3 == mod(quadrant-3,2));
        
        for row = Y_corner
            cols = left_line(row) + mod((0:(width_wedge-1))-(left_line(row)-first_col),width_wedge);
            admissible_cols = round(1/2*(cols+2*floor(4*M_horiz)+1-abs(cols-(2*floor(4*M_horiz)+1))));
            new_row = 1 + mod(row - first_row, length_corner_wedge);
            wrapped_XX(new_row,:) = XX(row,admissible_cols);
            wrapped_YY(new_row,:) = YY(row,admissible_cols);        
        end;
        YY = Y_corner'*ones(1,width_wedge);
        slope_wedge_left = ((floor(4*M_horiz)+1) - wedge_midpoints(end))/floor(4*M_vert);
        mid_line_left = wedge_midpoints(end) + slope_wedge_left*(wrapped_YY - 1);
        coord_left = 1/2 + floor(4*M_vert)/(wedge_endpoints(end) - wedge_endpoints(end-1)) * ...
            (wrapped_XX - mid_line_left)./(floor(4*M_vert)+1 - wrapped_YY);
        C2 = -1/(2*(floor(4*M_horiz))/(wedge_endpoints(end) - 1) - 1 + 1/(2*(floor(4*M_vert))/(first_wedge_endpoint_vert - 1) - 1));
        C1 = -C2 * (2*(floor(4*M_horiz))/(wedge_endpoints(end) - 1) - 1);
        wrapped_XX((wrapped_XX - 1)/floor(4*M_horiz) == (wrapped_YY-1)/floor(4*M_vert)) = ...
            wrapped_XX((wrapped_XX - 1)/floor(4*M_horiz) == (wrapped_YY-1)/floor(4*M_vert)) - 1;
        coord_corner = C1 + C2 * (2-((wrapped_XX - 1)/(floor(4*M_horiz)) + (wrapped_YY - 1)/(floor(4*M_vert)))) ./ ...
            ((wrapped_XX - 1)/(floor(4*M_horiz)) - (wrapped_YY - 1)/(floor(4*M_vert)));
        wl_left = fdct_wrapping_window(coord_left);
        [wl_right,wr_right] = fdct_wrapping_window(coord_corner);
        switch is_real
         case 0
          wrapped_data = fftshift(fft2(ifftshift(C{j}{l})))/sqrt(prod(size(C{j}{l})));
          wrapped_data = rot90(wrapped_data,(quadrant-1));
         case 1
          x = C{j}{l} + sqrt(-1)*C{j}{l+nbangles(j)/2};
          wrapped_data = fftshift(fft2(ifftshift(x)))/sqrt(prod(size(x)))/sqrt(2);
          wrapped_data = rot90(wrapped_data,(quadrant-1));
        end;
        wrapped_data = wrapped_data .* (wl_left .* wr_right);
        
         % Unwrapping data
        for row = Y_corner
            cols = left_line(row) + mod((0:(width_wedge-1))-(left_line(row)-first_col),width_wedge);
            admissible_cols = round(1/2*(cols+2*floor(4*M_horiz)+1-abs(cols-(2*floor(4*M_horiz)+1))));
            new_row = 1 + mod(row - first_row, length_corner_wedge);
            Xj(row,fliplr(admissible_cols)) = Xj(row,fliplr(admissible_cols)) + wrapped_data(new_row,end:-1:1);
                                % We use the following property: in an assignment
                                % A(B) = C where B and C are vectors, if
                                % some value x repeats in B, then the
                                % last occurrence of x is the one
                                % corresponding to the eventual assignment.
        end;

        Xj = rot90(Xj);
        
    end;    % for quadrant
    
    Xj = Xj .* lowpass;
    Xj_index1 = ((-floor(2*M1)):floor(2*M1)) + floor(4*M1) + 1;
    Xj_index2 = ((-floor(2*M2)):floor(2*M2)) + floor(4*M2) + 1;
    Xj(Xj_index1, Xj_index2) = Xj(Xj_index1, Xj_index2) .* hipass;
    
    loc_1 = Xj_topleft_1 + (0:(2*floor(4*M1)));
    loc_2 = Xj_topleft_2 + (0:(2*floor(4*M2)));
    X(loc_1,loc_2) = X(loc_1,loc_2) + Xj;

    % Preparing for loop reentry or exit
    
    Xj_topleft_1 = Xj_topleft_1 + floor(4*M1) - floor(2*M1);
    Xj_topleft_2 = Xj_topleft_2 + floor(4*M2) - floor(2*M2);
    
    lowpass = lowpass_next;
    
end;    % for j

if is_real
    Y = X;
    X = rot90(X,2);
    X = X + conj(Y);
end
    
% Coarsest wavelet level
M1 = M1/2;
M2 = M2/2;
Xj = fftshift(fft2(ifftshift(C{1}{1})))/sqrt(prod(size(C{1}{1})));
loc_1 = Xj_topleft_1 + (0:(2*floor(4*M1)));
loc_2 = Xj_topleft_2 + (0:(2*floor(4*M2)));
X(loc_1,loc_2) = X(loc_1,loc_2) + Xj .* lowpass;

% Finest level
M1 = N1/3;
M2 = N2/3;
if finest == 1,

    % Folding back onto N1-by-N2 matrix
    shift_1 = floor(2*M1)-floor(N1/2);
    shift_2 = floor(2*M2)-floor(N2/2);
    Y = X(:,(1:N2)+shift_2);
    Y(:,N2-shift_2+(1:shift_2)) = Y(:,N2-shift_2+(1:shift_2)) + X(:,1:shift_2);
    Y(:,1:shift_2) = Y(:,1:shift_2) + X(:,N2+shift_2+(1:shift_2));
    X = Y((1:N1)+shift_1,:);
    X(N1-shift_1+(1:shift_1),:) = X(N1-shift_1+(1:shift_1),:) + Y(1:shift_1,:);
    X(1:shift_1,:) = X(1:shift_1,:) + Y(N1+shift_1+(1:shift_1),:);
    
else
    
    % Extension to a N1-by-N2 matrix
    Y = fftshift(fft2(ifftshift(C{nbscales}{1})))/sqrt(prod(size(C{nbscales}{1})));
    X_topleft_1 = ceil((N1+1)/2) - floor(M1);
    X_topleft_2 = ceil((N2+1)/2) - floor(M2);
    loc_1 = X_topleft_1 + (0:(2*floor(M1)));
    loc_2 = X_topleft_2 + (0:(2*floor(M2)));
    Y(loc_1,loc_2) = Y(loc_1,loc_2) .* hipass_finest + X;
    X = Y;
    
end;

x = fftshift(ifft2(ifftshift(X)))*sqrt(prod(size(X)));
if is_real, x = real(x); end;
end

%% 2D Framelet
function C = ConvSymAsym(A,M,L)
[m,n] = size(A);
nM = length(M);
step = 2^(L-1);
ker = zeros(step*(nM-1)+1,1);
ker(1:step:step*(nM-1)+1,1) = M;
lker = floor(length(ker)/2);
C = imfilter(A,ker,'circular');
end
function Dec = FrameletDec(A,L)
    D{1} = [1 2 1]/4;
    D{2} = [1 0 -1]/4*sqrt(2);
    D{3} = [-1 2 -1]/4;
    kDec = A;
for k = 1:L
    Dec{k} = FraDec(kDec,D,k);
    kDec = Dec{k}{1,1};
end
end
function Dec = FraDec(A,D,L)
nD = length(D);
for i = 1:nD
    M1 = D{i};
    tempi = ConvSymAsym(A,M1,L);
    for j = 1:nD
        M2 = D{j};        
        tempj = ConvSymAsym(tempi',M2,L);
        Dec{i,j} = single(tempj');
    end
end
end
function Rec = FrameletRec(C,L)
    R{1} = [1 2 1]/4;
    R{2} = [-1 0 1]/4*sqrt(2);
    R{3} = [-1 2 -1]/4;
for k = L:-1:2
    C{k-1}{1,1} = FraRec(C{k},R,k);
end
Rec = FraRec(C{1},R,1);
end
function Rec = FraRec(C,R,L)
nR = length(R);
ImSize = size(C{1,1});
Rec = zeros(ImSize);

for i = 1:nR
    temp = zeros(ImSize);
    for j = 1:nR
        M2 = R{j};
        temp = temp + (ConvSymAsym((C{i,j})',M2,L))';
    end
    M1 = R{i};
    Rec = Rec + ConvSymAsym(temp,M1,L);
end
end

%% Basic operation
function [filename] = GenerateTimeFolder(filename)
aa=clock;
for i=1:5
    filename=[filename,num2str(aa(i)),'_'];
end
end
function P = padPSF(PSF, m, n)
if nargin == 2
  if length(m) == 1
    n = m;
  else
    n = m(2); 
    m = m(1);
  end
end

%
% Pad the PSF with zeros.
%
P = zeros(m, n);
P(1:size(PSF,1), 1:size(PSF,2)) = PSF;
end
function PSF = PSF_Generator(lambada,pixelsize,NA,w,factor)
[X,Y]=meshgrid(linspace(0,w-1,w),linspace(0,w-1,w));
scale=2*pi*NA/lambada*pixelsize;
scale=scale*factor;
R=sqrt(min(X,abs(X-w)).^2+min(Y,abs(Y-w)).^2);
PSF=abs(2*besselj(1,scale*R+eps,1)./(scale*R+eps)).^2;
PSF = PSF/(sum(sum(PSF)));
PSF=fftshift(PSF);
end
function raw_data = readMTiffn(raw_file,n)
  info = imfinfo(raw_file);
  frames = numel(info);
  if n==8
  raw_data = zeros(info(1).Height, info(1).Width, frames, 'uint8');
  for k = 1:frames
      raw_data(:,:,k) = im2uint8(imread(raw_file, k));
  end
  elseif n==16
  raw_data = zeros(info(1).Height, info(1).Width, frames, 'uint16');
  for k = 1:frames
      raw_data(:,:,k) = im2uint16(imread(raw_file, k));
  end      
  else
        raw_data = zeros(info(1).Height, info(1).Width, frames, 'uint32');
  for k = 1:frames
      raw_data(:,:,k) = im2uint32(imread(raw_file, k));
  end  
end
end
function [] = writeMTiffn(f_stack,filename,n)
stackfilename = filename;
[~,~,z]=size(f_stack);
% m=max(max(max(f_stack)));
% mi=min(min(min(f_stack)));
for k = 1:z
    X=f_stack(:,:,k);
    X=double(X);
    m=max(max(X));
    mi=min(min(X));
    if n==8
    X=(X-mi)*(2^8-1)/(m-mi);
    imwrite(uint8(X), stackfilename, 'WriteMode','append') % 写入stack图像
    elseif n==16
    X=(X-mi)*(2^16-1)/(m-mi);
    imwrite(uint16(X), stackfilename, 'WriteMode','append') % 写入stack图像
    else
            X=(X-mi)*(2^32-1)/(m-mi);
    imwrite(uint32(X), stackfilename, 'WriteMode','append') % 写入stack图像
    end
end

end
