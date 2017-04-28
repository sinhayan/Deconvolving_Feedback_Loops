function  [Scoremat, fracrat]=Deconvolve_feedback_loops(R,minrem)
% R is the ratings matrix
%minrem is the minimum number of ratings of an item to be considered
%Scoremat is the 1-score as defined in paper for every user-item pair
%fracrat is the fraction of items undergoing feedback effects. 

   
    [Rcentered,V,Sobs,~,~,svdnum,~]=deconvoluteSVD(R,minrem);
   
    [Scoremat,~,~,~,~,~,fracrat]=deconvoluteprocess(Rcentered,V,Sobs,svdnum);
end


function [Rcentered,V,Sobs,idxj,idxi,svdnum,R]=deconvoluteSVD(R,minrem)

if size(R,1)*size(R,2)>25000000
    svdnum=3000;
else
    svdnum=size(R,2);
end

%remove colums with less than minrem ratings
idxR=sum(R);
R(:,idxR==0)=[];

[~,n]=size(R);
[~,j,~]=find(R);
NRj=histc(j,1:n);

% NR=R>0;
% NRj=sum(NR);
idxj=NRj<minrem;
R(:,idxj)=[];

%remove empty rows
% NR=R~=0;
[m,~]=size(R);
[i,~,~]=find(R);
NR=histc(i,1:m);
% NR=sum(NR,2);
idxi=NR==0;
R(idxi,:)=[];


%normalize
[i,j,valint]=find(R);
[m,n]=size(R);

NR=histc(i,1:m);

NRat=sum(R,2);
NRR=NRat./NR;

NRR=full(NRR);
valcen=NRR(i);
val=valint-valcen;

clear valcen valint
Rcentered=sparse(i,j,val,m,n);


Rcentered = normRat(Rcentered);
clear i j val NR NRR NRat
disp( 'matrix normalized' )


RcenteredT=Rcentered';
f = @(x) RcenteredT*(RcenteredT'*x);
[V Sobssq]=eigs(f,n,min(svdnum,n),'LA',struct('issym',1,'disp',0));
Sobs=sqrt(Sobssq);



end

function Rn = normRat(R)
itemnorms = (sum(R.^2));
Rn = R*sparse(1:length(itemnorms),1:length(itemnorms),1./sqrt(itemnorms),length(itemnorms),length(itemnorms));
end

function [Scoremat,RconvmatN,RcenteredN,Rconv,ratingsort,timessort,fracrat]=deconvolute_process(Rcentered,V,Sobs,svdnum)

close all
plotfig=1; %to plot figures



threshslope=5;
[~,n]=size(Rcentered);
SobsI=spdiags(1./diag(Sobs),0,min(svdnum,n),min(svdnum,n));
U=Rcentered*(V*SobsI);

alpha=1;
[~,~,singvalobs]=find(Sobs);
idx= singvalobs<1e-10;
singvaltrue=(-1./(2*alpha*singvalobs) + sqrt((1./(4*alpha^2*singvalobs.^2))+1/alpha));
singvaltrue(idx)=0;

%calulate deconvolution
Strue=sparse(1:min(size(Sobs)),1:min(size(Sobs)),singvaltrue,size(Sobs,1),size(Sobs,2));
Vconv=Strue*V';
Vconv=Vconv';
[i,j,valR]=find(Rcentered);
Rconv=zeros(length(i),1);


lump=100000;


for kk=1:ceil(length(i)/lump)
    jj=(kk-1)*lump+1:min(kk*lump,length(i));
    val=sum(U(i(jj),:).*Vconv(j(jj),:),2);
    Rconv(jj)=val;
    disp(kk);
end

clear U Sobs Vconv V RcenteredT
disp('Values calculated')

if plotfig==1
    % plot original
    p = randperm(length(i));
    lump=1000000;
    figure;
    plot(zeros(length(min(valR):0.01:max(valR)),1),min(valR):0.01:max(valR),'k-');hold on
    plot(min(Rconv):0.01:max(Rconv),zeros(length(min(Rconv):0.01:max(Rconv)),1),'k-');
    for kk=1:ceil(length(i)/lump)
        jj=(kk-1)*lump+1:min(kk*lump,length(i));
        idx=p(jj);
        plot(Rconv(idx),valR(idx),'c.');
        hold on
    end
    xlabel('Deconvoluted','FontSize',40)
    ylabel('Observed','FontSize',40)
end

scaling=max(abs(valR))/max(abs(Rconv));
Rconv=scaling*Rconv;



%% sort by movie or user
disp('calculating ransac')



RconvmatN=zeros(length(Rconv),1);
RcenteredN=zeros(length(Rconv),1);

[~,index,~]=unique(j);
%calculate randac


for kk=1:length(index)
    
    if kk==1
        iint=1;
    elseif kk>1
        iint=index(kk-1)+1;
    end
    ient=index(kk);
    
    iintient=iint:ient;
    RX=Rconv(iintient);
    RY=valR(iintient);
    
    if length(RY)>1;
        
        [lineeq(1),lineeq(2)] = ransacfast([RX,RY]',2,50,0.01,0.3);
        
        
        if lineeq(1)==0
            [lineeq(1),lineeq(2)] = ransacfast([RX,RY]',2,100,0.1,0.3);
        end
        
        if lineeq(1)==0
            [lineeq(1),lineeq(2)] = ransacfast([RX,RY]',2,100,1,0.3);
        end
        
        if lineeq(1)>0
            theta=pi/2-atan(lineeq(1));
        else
            theta=-pi/2-atan(lineeq(1));
        end
        
        if abs(lineeq(1))>threshslope
            theta=0;
            lineeq(2)=0;
        end
        
        
        Rot = [cos(theta) -sin(theta); sin(theta) cos(theta)];
        RYT=RY-lineeq(2);
        
        NewXY= Rot*[RX';RYT'];
        RXn=NewXY(1,:);
        RYn=NewXY(2,:)+lineeq(2);
        RconvmatN(iintient)=RXn;
        RcenteredN(iintient)=RYn;
        
        
        if mod(kk,100)==0
            disp(kk)
        end
    else
        RconvmatN(iintient)=RX;
        RcenteredN(iintient)=RY;
    end
end

%% plot transformed ratings

if plotfig==1
    p = randperm(length(i));
    lump=1000000;
    figure;
    plot(zeros(length(min(valR):0.01:max(valR)),1),min(valR):0.01:max(valR),'k-');hold on
    plot(min(Rconv):0.01:max(Rconv),zeros(length(min(Rconv):0.01:max(Rconv)),1),'k-');
    for kk=1:ceil(length(i)/lump)
        jj=(kk-1)*lump+1:min(kk*lump,length(i));
        idx=p(jj);
        Rp=RconvmatN(idx);
        Rc=RcenteredN(idx);
        plot(Rp,Rc,'c.');
        hold on
        
    end
    xlabel('Deconvoluted','FontSize',40)
    ylabel('Observed','FontSize',40)
    
end


%% assign score
scale=max(abs(RcenteredN))/max(abs(RconvmatN));
RconvmatN=scale*RconvmatN;



slope=1;
score=sqrt(RconvmatN.^2-(RcenteredN.^2/slope.^2));
score=real(score);
clear dist p


Scoremat=sparse(i,j,score,size(Rcentered,1),size(Rcentered,2));
allscore=sum(Scoremat);
Ntimes=histc(j,1:n);
allscore=allscore./Ntimes';

[~,ratingsort]=sort(allscore);
timessort=Ntimes(ratingsort);
fracrat=1-sum(score==0)/length(j);

end


function [bestParameter1,bestParameter2] = ransacfast(data,num,iter,threshDist,inlierRatio)
% data: a 2xn dataset with #n data points
% num: the minimum number of points. For line fitting problem, num=2
% iter: the number of iterations
% threshDist: the threshold of the distances between points and the fitting line
% inlierRatio: the threshold of the numer of inliers

%% Plot the data points
%figure;plot(data(1,:),data(2,:),'o');hold on;
number = size(data,2); % Total number of points
bestInNum = 0; % Best fitting line with largest number of inliers
bestParameter1=0;bestParameter2=0; % parameters for best fitting line
for i=1:iter
    %% Randomly select 2 points
    idx = randperm(number);
    idx=idx(1:num);
    sample = data(:,idx);
    %% Compute the distances between all points with the fitting line
    kLine = sample(:,2)-sample(:,1);
    kLineNorm = kLine/norm(kLine);
    normVector = [-kLineNorm(2),kLineNorm(1)];
    distance = normVector*(data - repmat(sample(:,1),1,number));
    %% Compute the inliers with distances smaller than the threshold
    inlierIdx = find(abs(distance)<=threshDist);
    inlierNum = length(inlierIdx);
    %% Update the number of inliers and fitting model if better model is found
    if inlierNum>=round(inlierRatio*number) && inlierNum>bestInNum
        bestInNum = inlierNum;
        parameter1 = (sample(2,2)-sample(2,1))/(sample(1,2)-sample(1,1));
        parameter2 = sample(2,1)-parameter1*sample(1,1);
        bestParameter1=parameter1; bestParameter2=parameter2;
    end
end

%  %% Plot the best fitting line
%  yAxis = min(data(2,:)):0.01:max(data(2,:));
%  xAxis = (yAxis-bestParameter2)/bestParameter1;
%  plot(xAxis,yAxis,'r-','LineWidth',2);
%  ylim([min(data(2,:)),max(data(2,:))]);
end


