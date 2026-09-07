%
clear all;
close all;
clc;
tic;  % 
load(fullfile('data', 'ori2.mat'));
data=ori2;



drr_plot3d(data,[500,100,80],0.1); % 
missing_ratio =0.1;
mask =rand(size(data))>missing_ratio;
data_missing= data .* mask;
drr_plot3d(data_missing,[100,100,80],0.1);

[hk,mk, n]=size(data);
fresult=zeros(size(data));
fZ = dct(data_missing,[],1);
drr_plot3d(fZ,[100,100,80],0.1);  % 
tic;
% s_cplot(squeeze(fZ(150,:,:)));
for hi = 1 : 80%
%         hi=140
hi
            alpha=5;
            k=5;
            dt=0.005;
            flow=0;
            fhigh=130;
            DATA=squeeze(fZ(hi,:,:));
            DATA=permute(DATA,[2,1]);
            [nt,nx] = size(DATA);
            nf = 2*2^nextpow2(nt);
            DATA_FX_f = zeros(nf,nx);
            sing = [];        
            % First and last samples of the DFT.
            ilow  = floor(flow*dt*nf)+1;
            if ilow<1
                     ilow=1;
            end
            ihigh = floor(fhigh*dt*nf)+1;
            if ihigh > floor(nf/2)+1
                        ihigh=floor(nf/2)+1;
            end
            % Transform to FX
             DATA_FX_tmp = fft(DATA,nf,1);
            samp = 1:1:nf;
            f = (samp-1)/(nf*dt);
            f = f(1:ihigh);
            % Size of Hankel Matrix
            Lcol = floor(nx/2)+1;
            Lrow = nx-Lcol+1;
            % Form level-1 block Hankel matrix 
             for j= ilow:ihigh
                 M = zeros(Lrow,Lcol);
                 for lc = 1:Lcol
                         M(:,lc)  = DATA_FX_tmp(j,lc:lc+Lrow-1);
                 end
            datam=M;

            %%
            %
            [m,n]=size(datam);
            [initial_u, initial_d, initial_v] = svd(datam);
            U = initial_u(:, 1:k);
            V = initial_v(:, 1:k);
            Lambda = diag(initial_d);
            Lambda = Lambda(1:k);
            U = U * diag(Lambda);
            for c=1:3%
            RL=datam-U * (V');
            r=reshape(RL, [], 1);
            seigma= 1.4862*median(abs(r - median(r)));
            W=customFunction(RL, alpha,seigma);
            for cout=1:3%
            for i = 1:m
                ui=U(i,:)';
                mi=datam(i,:)';
                wi=diag(W(i,:));
                A = V' * wi' * V;
                b = V' * wi' * mi;
                [Q, R] = qr(A);
                y = Q' * b;
                u_i = backSubstitution(R, y);
                U(i,:)=u_i';
            end
            for q = 1:n
                vj=V(q,:)';
                mj=datam(:,q);
                wj=diag(W(:,q));
                A = U' * wj' * U;
                b = U' * wj' * mj;
                [Q, R] = qr(A);
                y = Q' * b;
                v_j = backSubstitution(R, y);
                V(q,:)=v_j';
            end
            end
            end
            result=U*V';
                Count = zeros(nx,1);
               tmp2 = zeros(nx,1);

                 for ic = 1:Lcol;
                  for ir = 1:Lrow;
                   Count(ir+ic-1,1) = Count(ir+ic-1,1)+1;
                   tmp2(ir+ic-1,1)  = tmp2(ir+ic-1,1) + result(ir,ic);%Mout涔熸槸476*476
                   end;
                  end

                   tmp2 = tmp2./Count;

                 DATA_FX_f(j,:) = tmp2;
            end
                    % Honor symmetries
               for k2=nf/2+2:nf
                DATA_FX_f(k2,:) = conj(DATA_FX_f(nf-k2+2,:));
               end
            % Back to TX (the output)
               DATA_f = real(ifft(DATA_FX_f,[],1));
               DATA_f = DATA_f(1:nt,:);
               DATA_f=permute(DATA_f,[2,1]);
               fresult(hi,:,:)=DATA_f;   
end
elapsedTime = toc;
disp(['时间: ', num2str(elapsedTime), ' s']);
fresult(1:140,:,:)=fZ(1:140,:,:);
construct_result=idct(fresult,[],1);
drr_plot3d(construct_result,[500,100,80],0.2); % 
drr_plot3d(data-construct_result,[500,100,80],0.2);% 



%%
function resultMatrix = customFunction(RL, alpha,seigma)
    resultMatrix = zeros(size(RL));
    [m,n]=size(RL);
    for i = 1:m
        for j=1:n
            x=RL(i,j)/seigma;
            if abs(x) <= alpha
                resultMatrix(i,j) = (1 - (abs(x) / alpha)^2)^2;
            else
                resultMatrix(i,j) = 0;
            end
        end
    end
end

function x = backSubstitution(R, y)
    n = length(y);
    x = zeros(n, 1);
    for i = n:-1:1
        x(i) = y(i) / R(i, i);
        y(1:i-1) = y(1:i-1) - R(1:i-1, i) * x(i);
    end
end