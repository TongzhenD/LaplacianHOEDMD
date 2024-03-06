% This code demonstrates the application of Laplacian HOEDMD method to 
% conduct the experiment of Robust network exploration in task-fMRI data 
% in the main text. Other experiments are similar.

clc
clear
% get the Laplacian observations snapshot matrix of subject 102109
path='..\tfMRI_language_data\102109\102109_language_data_filtered.mat';
math_path='..\tfMRI_language_data\102109\EVs\math.txt';
story_path='..\tfMRI_language_data\102109\EVs\story.txt';
[n,new_W1,ltrial1] = GetLapofSubject(path,math_path,story_path);

% get the Laplacian observations snapshot matrix of subject 101915
path='..\tfMRI_language_data\101915\101915_language_data_filtered.mat';
math_path='..\tfMRI_language_data\101915\EVs\math.txt';
story_path='..\tfMRI_language_data\101915\EVs\story.txt';
[~,new_W2,ltrial2] = GetLapofSubject(path,math_path,story_path);
 
% get the Laplacian observations snapshot matrix of subject 1020085
path='..\tfMRI_language_data\102008\102008_language_data_filtered.mat';
math_path='..\tfMRI_language_data\102008\EVs\math.txt';
story_path='..\tfMRI_language_data\102008\EVs\story.txt';
[~,new_W3,ltrial3] = GetLapofSubject(path,math_path,story_path);

% get the Laplacian observations snapshot matrix of subject 100206
path='..\tfMRI_language_data\100206\100206_language_data_filtered.mat';
math_path='..\tfMRI_language_data\100206\EVs\math.txt';
story_path='..\tfMRI_language_data\100206\EVs\story.txt';
[~,new_W4,ltrial4] = GetLapofSubject(path,math_path,story_path);

new_W=cat(2,new_W1,new_W2,new_W3,new_W4); % Psi0
ltrial=cat(1,ltrial1,ltrial2,ltrial3,ltrial4); % list of lengths for each trajetory

M=size(new_W,1);

p = 2;
S=min(M*p,sum(ltrial)-length(ltrial)*p);
multi_traje=zeros(M*(p+1),sum(ltrial)-p*length(ltrial)); % Psi
for j=1:length(ltrial)
    c=sum(ltrial(1:j));
    X=new_W(:,c-ltrial(j)+1:c);
    for i=1:p+1
        multi_traje((M*(i-1)+1):(M*i),c-ltrial(j)-(j-1)*p+1:c-j*p)=X(:,i:(i+ltrial(j)-p-1));
    end
end
% apply Laplacian HOEDMD based on Gram matrix
tic;
[U,Lambda,V] = hoedmd_gram(new_W,multi_traje,p,S,ltrial);
toc;

% get the Laplacian observations snapshot matrix of subject
function [n,new_W,ltrial] = GetLapofSubject(path,math_path,story_path)
    load(path);
    data=zscore(data,0,2); % zscore the filtered data
    n=size(data,1);
    % compute the Point Cloud Laplace opreator
    k=2;alpha=1;tau=k+2+alpha;
    window_size=5;
    sliding_step=1;
    T=floor((size(data,2)-window_size)/sliding_step)+1;
    W=zeros(n,n,T);
    for i=1:T
        W(:,:,i)=exp(squareform(pdist(data(:,sliding_step*(i-1)+1:sliding_step*(i-1)+window_size)))./(-0.2*tau^2));
    end
    
    % remove noisy parts at the beginning and end
    load(math_path);
    load(story_path);

    trail=[math(:,1:2);story(:,1:2)];
    trail=trail./0.72;
    trail(:,2)=sum(trail(:,1:2),2);
    
    trail(:,1)=ceil(trail(:,1))+1;
    trail(:,2)=floor(trail(:,2))-1;
    
    [~,index]=sort(trail(:,1),'ascend');
    trail=trail(index,:);
    n_trail=size(trail,1)-1;
    ltrial=trail(1:n_trail,2)-trail(1:n_trail,1)+1;
    
    extr_W=[];
    for i=1:length(ltrial)
        temp=W(:,:,trail(i,1):trail(i,2));
        extr_W=cat(3,extr_W,temp);
    end
    T=size(extr_W,3);
    
    % take out the elements of the lower triangle to generate the Laplacian observations
    index_Matrix=tril(ones(n),-1);
    new_W=zeros((n*(n-1))/2,T);
    for i =1:T
        temp=extr_W(:,:,i); % 
        new_W(:,i)=temp(index_Matrix==1);
    end
end