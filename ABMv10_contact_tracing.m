clear all;
clc
close all;
rng(1)


%% load population
load("popGeneration.mat")


%% Default parameters 
% Population parameters   set by users or import
nOfIniE = 2;
nOfIniP = 1;
nOfIniAsym = 2;              
nOfIniSym = 2;

% Simulation parameters
nDays = 15;

% Basic disease transmission parameters
R0 = 2.5;           %basic reproduction number
mu = 1/2.9;         %2.9 a rate of moving from I to R
r = 0.55;          %ratio of transmission rate for asymptomatic over symptomatic cases
epsilon = 1/2.9;   %a rate of moving from E to P
gamma = 1/2.3;     %a rate of moving from P to I
p = 30.8/100;      %Asymptomatic rate
beta = R0*mu/(mu*r/gamma+p*r+1-p);   %a rate of moving from S to E  0.6644


meanNOfContactPerDay = 10;
varNOfContactPerDay = 0.6;


% Defined by users:
nOfIniS = popSize - nOfIniE - nOfIniP - nOfIniAsym - nOfIniSym;
nOfIniInfected = nOfIniE + nOfIniP + nOfIniAsym + nOfIniSym;    %infected
nOfIniR = popSize - nOfIniS - nOfIniInfected;

%% Generation Cases
idxOfIniS = uid(1:nOfIniS);
idxOfIniE = uid((nOfIniS + 1) : (nOfIniS + nOfIniE));
idxOfIniP = uid((nOfIniS + nOfIniE + 1) : (nOfIniS + nOfIniE + nOfIniP));
idxOfIniAsym = uid((nOfIniS + nOfIniE + nOfIniP + 1) : (nOfIniS + nOfIniE + nOfIniP + nOfIniAsym));
idxOfIniSym = uid((nOfIniS + nOfIniE + nOfIniP + nOfIniAsym + 1) : (nOfIniS + nOfIniInfected));
idxOfIniRe = uid((nOfIniS + nOfIniInfected + 1): popSize);
iniDurE2P = truncatedPoisson(epsilon,1,nOfIniE);
iniDurP2AsymOrSym = truncatedPoisson(gamma,1,nOfIniE+nOfIniP);
iniDurAsymOrSym2Re = truncatedPoisson(mu,1,nOfIniInfected);

People(idxOfIniS,indSCol) = 1;
People = fIniInfect(People,idxOfIniE,idxOfIniP,idxOfIniAsym,idxOfIniSym,iniDurE2P,iniDurP2AsymOrSym,iniDurAsymOrSym2Re,0,indECol,indPCol,indAsymCol,indSymCol,durE2PCol,durP2AsymOrSymCol,durAsymOrSym2ReCol,timeOfECol,timeOfPCol,timeOfAsymOrSymCol,timeOfReCol);
People(idxOfIniRe,indReCol) = 1;

[nS,nE,nP,nAsym,nSym,nRe,nQ,nC] = countNumOfAll(People,indSCol,indECol,indPCol,indAsymCol,indSymCol,indReCol,indQCol,indCCol);

nEPerDay = 0;

colorOfS = [55,126,184]/256;
colorOfE = [152,78,163]/256;
colorOfP = [228,26,28]/256;
colorOfAsym = [228,26,28]/256;
colorOfSym = [228,26,28]/256;
colorOfRe = [77,175,74]/256;


listOfBackwardTracing = [];
recordOfTransmission =[];

%% Paramaters about Contact Tracing
minDelayC = 3.89;  maxDelayC = 8.93;
minDelayT = 1;     maxDelayT = 5.53;
kappaC = 9.5e-2;    thldC = 20.19;
kappaT = 4.2e-3;    thldT = 19.85;
p2 = 24.4/100;   %probability of back-traced
q = 11.4/100;    %quarantined to unquarantined
startDayOfContactTracing = 5;


%% Run
% figure
% mapshow(S, 'FaceColor', 'white');
% mapshow(coordinate_points_set(idxOfIniS), 'FaceColor', colorOfS);
% mapshow(coordinate_points_set(idxOfIniE), 'FaceColor', colorOfE);
% mapshow(coordinate_points_set(idxOfIniP), 'FaceColor', colorOfP);
% mapshow(coordinate_points_set(idxOfIniAsym), 'FaceColor', colorOfAsym);
% mapshow(coordinate_points_set(idxOfIniSym), 'FaceColor', colorOfSym);
% mapshow(coordinate_points_set(idxOfIniRe), 'FaceColor', colorOfRe);
countOfRecord = 1;

for day = 1:nDays
    tic
    nEDay = 0;

    idxS = fGetIdxS(People,indSCol);
    [idxE,idxP,idxAsym,idxSym] = fGetIdxOfAllInfected(People,indECol,indPCol,indAsymCol,indSymCol);
    idxRe = fGetIdxRe(People,indReCol);
    
    [idxE2P,idxP2AsymOrSym,idxAsym2Re,idxSym2Re] = check4StatesChange(People,idxE,idxP,idxAsym,idxSym,durE2PCol,durP2AsymOrSymCol,durAsymOrSym2ReCol);
    

    idxOfKeepE = setdiff(idxE,idxE2P);
    idxOfKeepP = setdiff(idxP,idxP2AsymOrSym);
    idxOfKeepAsym = setdiff(idxAsym,idxAsym2Re);
    idxOfKeepSym = setdiff(idxSym,idxSym2Re);
    
    numOfContact = round(lognrnd(log(meanNOfContactPerDay),varNOfContactPerDay,nP(day)+nAsym(day)+nSym(day),1));
    countInfect = 1;
    for countSym = 1:length(idxSym)
        flagOfInfType = 3;
        sourceDist = People(idxSym(countSym),DistCol);
        nSymContact = numOfContact(countInfect);
        sourceUid = People(idxSym(countSym),uidCol);
        referRandNumOfDist = rand(nSymContact,1);

        [idxContact, uidContact, idxInfected, uidInfected]  = fGetIdxUidOfContactAndInfected(People,cumPostPOfContactRate,idxAnduidOfDists,referRandNumOfDist,nSymContact,sourceDist,flagOfInfType,r,beta,indSCol);

        symBackwardTracing = [uidInfected,repmat(sourceUid,length(uidInfected),1)];
        listOfBackwardTracing = [listOfBackwardTracing;symBackwardTracing];
        if isempty(recordOfTransmission)
            recordOfTransmission(countOfRecord).sourceUid = sourceUid;
            recordOfTransmission(countOfRecord).uidContact = uidContact;
            recordOfTransmission(countOfRecord).uidInfected = uidInfected;
            countOfRecord = countOfRecord + 1;
        else
            existSourceUid = cat(1,recordOfTransmission.sourceUid);
            if find(existSourceUid == sourceUid)
                idxSource = find(existSourceUid == sourceUid);
                allUidContact = [recordOfTransmission(idxSource).uidContact; uidContact];
                allUidContact = unique(allUidContact);
                recordOfTransmission(idxSource).uidContact = allUidContact;
                allUidInfected = [recordOfTransmission(idxSource).uidInfected;uidInfected];
                recordOfTransmission(idxSource).uidInfected = allUidInfected;
            else
                recordOfTransmission(countOfRecord).sourceUid = sourceUid;
                recordOfTransmission(countOfRecord).uidContact = uidContact;
                recordOfTransmission(countOfRecord).uidInfected = uidInfected;   
                countOfRecord = countOfRecord + 1;
            end
        end
        nEDay = nEDay + length(uidInfected);
        People = fInfect(People,uidInfected,length(uidInfected),day,epsilon,gamma,mu,uidCol,indSCol,indECol,durE2PCol,durP2AsymOrSymCol,durAsymOrSym2ReCol,timeOfECol,timeOfPCol,timeOfAsymOrSymCol,timeOfReCol);
        countInfect = countInfect + 1;
        
    end
  
    for countAsym = 1:length(idxAsym)
        flagOfInfType = 2;
        sourceDist = People(idxAsym(countAsym),DistCol);
        nAsymContact = numOfContact(countInfect);
        sourceUid = People(idxAsym(countAsym),uidCol);
        referRandNumOfDist = rand(nAsymContact,1);
        [idxContact, uidContact, idxInfected, uidInfected]  = fGetIdxUidOfContactAndInfected(People,cumPostPOfContactRate,idxAnduidOfDists,referRandNumOfDist,nAsymContact,sourceDist,flagOfInfType,r,beta,indSCol);
        asymBackwardTracing = [uidInfected,repmat(sourceUid,length(uidInfected),1)];
        listOfBackwardTracing = [listOfBackwardTracing; asymBackwardTracing];
        if isempty(recordOfTransmission)
            recordOfTransmission(countOfRecord).sourceUid = sourceUid;
            recordOfTransmission(countOfRecord).uidContact = uidContact;
            recordOfTransmission(countOfRecord).uidInfected = uidInfected;
            countOfRecord = countOfRecord + 1;
        else
            existSourceUid = cat(1,recordOfTransmission.sourceUid);
            if find(existSourceUid==sourceUid)
                idxSource = find(existSourceUid==sourceUid);
                allUidContact = [recordOfTransmission(idxSource).uidContact;uidContact];
                allUidContact = unique(allUidContact);
                recordOfTransmission(idxSource).uidContact = allUidContact;
                allUidInfected = [recordOfTransmission(idxSource).uidInfected;uidInfected];
                recordOfTransmission(idxSource).uidInfected = allUidInfected;
            else
                recordOfTransmission(countOfRecord).sourceUid = sourceUid;
                recordOfTransmission(countOfRecord).uidContact = uidContact;
                recordOfTransmission(countOfRecord).uidInfected = uidInfected;
                countOfRecord = countOfRecord + 1;
            end
        end
        nEDay = nEDay + length(uidInfected);
        People = fInfect(People,uidInfected,length(uidInfected),day,epsilon,gamma,mu,uidCol,indSCol,indECol,durE2PCol,durP2AsymOrSymCol,durAsymOrSym2ReCol,timeOfECol,timeOfPCol,timeOfAsymOrSymCol,timeOfReCol);
        countInfect = countInfect + 1;
    end

    for countP = 1:length(idxP)
        flagOfInfType = 1;
        sourceDist = People(idxP(countP),DistCol);
        nPContact = numOfContact(countInfect);
        sourceUid = People(idxP(countP),uidCol);
        referRandNumOfDist = rand(nPContact,1);
        [idxContact, uidContact, idxInfected, uidInfected]  = fGetIdxUidOfContactAndInfected(People,cumPostPOfContactRate,idxAnduidOfDists,referRandNumOfDist,nPContact,sourceDist,flagOfInfType,r,beta,indSCol);
        pBackwardTracing = [uidInfected,repmat(sourceUid,length(uidInfected),1)];
        listOfBackwardTracing = [listOfBackwardTracing; pBackwardTracing];
        if isempty(recordOfTransmission)
            recordOfTransmission(countOfRecord).sourceUid = sourceUid;
            recordOfTransmission(countOfRecord).uidContact = uidContact;
            recordOfTransmission(countOfRecord).uidInfected = uidInfected;
            countOfRecord = countOfRecord + 1;
        else
            existSourceUid = cat(1,recordOfTransmission.sourceUid);
            if find(existSourceUid==sourceUid)
                idxSource = find(existSourceUid==sourceUid);
                allUidContact = [recordOfTransmission(idxSource).uidContact;uidContact];
                allUidContact = unique(allUidContact);
                recordOfTransmission(idxSource).uidContact = allUidContact;
                allUidInfected = [recordOfTransmission(idxSource).uidInfected;uidInfected];
                recordOfTransmission(idxSource).uidInfected = allUidInfected;
            else
                recordOfTransmission(countOfRecord).sourceUid = sourceUid;
                recordOfTransmission(countOfRecord).uidContact = uidContact;
                recordOfTransmission(countOfRecord).uidInfected = uidInfected;
                countOfRecord = countOfRecord + 1;
            end
        end
        nEDay = nEDay + length(uidInfected);
        People = fInfect(People,uidInfected,length(uidInfected),day,epsilon,gamma,mu,uidCol,indSCol,indECol,durE2PCol,durP2AsymOrSymCol,durAsymOrSym2ReCol,timeOfECol,timeOfPCol,timeOfAsymOrSymCol,timeOfReCol);
        countInfect = countInfect + 1;
    end

    indOfAsymOrSym = rand(numel(idxP2AsymOrSym),1);
    indOfAsymOrSym(indOfAsymOrSym > p) = 1;
    indOfAsymOrSym(indOfAsymOrSym <= p) = 0;

    People(idxE2P,indECol) = 0;
    People(idxE2P,indPCol) = 1;
    People(idxP2AsymOrSym,indPCol) = 0;
    People(idxP2AsymOrSym,indAsymCol) = 1-indOfAsymOrSym;
    People(idxP2AsymOrSym,indSymCol) = indOfAsymOrSym;
    People(idxAsym2Re,indAsymCol) = 0;
    People(idxAsym2Re,indReCol) = 1;
    People(idxSym2Re,indSymCol) = 0;
    People(idxSym2Re,indReCol) = 1;


    People(idxOfKeepE,durE2PCol) = People(idxOfKeepE,durE2PCol)-1;
    People(idxOfKeepP,durP2AsymOrSymCol) = People(idxOfKeepP,durP2AsymOrSymCol)-1;
    People(idxOfKeepAsym,durAsymOrSym2ReCol)= People(idxOfKeepAsym,durAsymOrSym2ReCol)-1;
    People(idxOfKeepSym,durAsymOrSym2ReCol) = People(idxOfKeepSym,durAsymOrSym2ReCol)-1;

    [nSDay,nEDay,nPDay,nAsymDay,nSymDay,nReDay,nQDay,nCDay] = countNumOfAll(People,indSCol,indECol,indPCol,indAsymCol,indSymCol,indReCol,indQCol,indCCol);


    nS = [nS nSDay];  nE = [nE nEDay];  nP = [nP nPDay];
    nAsym =[nAsym nAsymDay];  nSym =[nSym nSymDay];  nRe = [nRe nReDay];
    
    nEPerDay = [nEPerDay; nEDay];
    

%     idxSNow = fGetIdxS(People,indSCol);
%     idxS2E = setdiff(idxS, idxSNow);
%     coordS2E = coordinate_points_set(idxS2E);
%     mapshow(coordS2E,'FaceColor', colorOfE); 
% 
%     coordE2P = coordinate_points_set(idxE2P);
%     mapshow(coordE2P,'FaceColor', colorOfP); 
% 
%     coordP2I = coordinate_points_set(idxP2AsymOrSym);
%     mapshow(coordP2I,'FaceColor', colorOfAsym); 
% 
%     coordAsym2Re = coordinate_points_set(idxAsym2Re);
%     mapshow(coordAsym2Re,'FaceColor', colorOfRe); 
% 
%     coordSym2Re = coordinate_points_set(idxSym2Re);
%     mapshow(coordSym2Re,'FaceColor', colorOfRe); 
%     pause(0.5)
%     toc     

end

nOfCumInfected = cumsum(nEPerDay) + nOfIniInfected;
nInfectious = nP + nAsym + nSym;


%% Definition of function

function [nS,nE,nP,nAsym,nSym,nRe,nQ,nC] = countNumOfAll(People,indSCol,indECol,indPCol,indAsymCol,indSymCol,indReCol,indQCol,indCCol)
    nS = sum(People(:,indSCol)); 
    nE = sum(People(:,indECol)); 
    nP = sum(People(:,indPCol)); 
    nAsym = sum(People(:,indAsymCol)); 
    nSym = sum(People(:,indSymCol)); 
    nRe = sum(People(:,indReCol)); 
    nQ = sum(People(:,indQCol)); 
    nC = sum(People(:,indCCol)); 
end

function idxS = fGetIdxS(People,indSCol)
    idxS =  find(People(:,indSCol)); 
end

function [idxE,idxP,idxAsym,idxSym] = fGetIdxOfAllInfected(People,indECol,indPCol,indAsymCol,indSymCol)
    idxE =  find(People(:,indECol)); 
    idxP =  find(People(:,indPCol)); 
    idxAsym =  find(People(:,indAsymCol)); 
    idxSym =  find(People(:,indSymCol)); 
end

function idxRe = fGetIdxRe(People,indReCol)
    idxRe =  find(People(:,indReCol)); 
end

function [idxE2P,idxP2AsymOrSym,idxAsym2Re,idxSym2Re] = check4StatesChange(People,idxE,idxP,idxAsym,idxSym,durE2PCol,durP2AsymOrSymCol,durAsymOrSym2ReCol)
    durE2P = People(idxE,durE2PCol);
    idxE2P = idxE(durE2P == 0);
    durP2AsymOrSym = People(idxP,durP2AsymOrSymCol);
    idxP2AsymOrSym = idxP(durP2AsymOrSym==0);
    durAsym2Re = People(idxAsym,durAsymOrSym2ReCol);
    idxAsym2Re = idxAsym(durAsym2Re==0);
    durSym2Re = People(idxSym,durAsymOrSym2ReCol);
    idxSym2Re = idxSym(durSym2Re==0);
end

function truncatedPeriod = truncatedPoisson(mu,minvalue,size)
    truncatedPeriod = poissrnd(1/mu,size,1);
    truncatedPeriod(truncatedPeriod < minvalue) = minvalue;
end

function People = fInfect(People,uidInfect,nInfect,nDay,epsilon,gamma,mu,uidCol,indSCol,indECol,durE2PCol,durP2AsymOrSymCol,durAsymOrSym2ReCol,timeOfECol,timeOfPCol,timeOfAsymOrSymCol,timeOfReCol)
    idxInfect = [];
    for countInfect = 1:nInfect
        idxInfect = [idxInfect find(People(:,uidCol) == uidInfect(countInfect))];
    end
    durE2P = truncatedPoisson(epsilon,1,nInfect);
    durP2AsymOrSym = truncatedPoisson(gamma,1,nInfect);
    durAsymOrSym2Re = truncatedPoisson(mu,1,nInfect);
    
    People(idxInfect,indSCol) = 0;
    People(idxInfect,indECol) = 1;
    People(idxInfect,durE2PCol) = durE2P;
    People(idxInfect,durP2AsymOrSymCol) = durP2AsymOrSym;
    People(idxInfect,durAsymOrSym2ReCol) = durAsymOrSym2Re;
    People(idxInfect,timeOfECol) = nDay;
    People(idxInfect,timeOfPCol) = nDay + durE2P;
    People(idxInfect,timeOfAsymOrSymCol) = nDay + durE2P + durP2AsymOrSym;
    People(idxInfect,timeOfReCol) = nDay + durE2P + durP2AsymOrSym + durAsymOrSym2Re;
end

function People = fIniInfect(People,idxOfIniE,idxOfIniP,idxOfIniAsym,idxOfIniSym,iniDurE2P,iniDurP2AsymOrSym,iniDurAsymOrSym2Re,nDay,indECol,indPCol,indAsymCol,indSymCol,durE2PCol,durP2AsymOrSymCol,durAsymOrSym2ReCol,timeOfECol,timeOfPCol,timeOfAsymOrSymCol,timeOfReCol)
    People(idxOfIniE,indECol) = 1;
    People(idxOfIniE,durE2PCol) = iniDurE2P;
    People(idxOfIniE,durP2AsymOrSymCol) = iniDurP2AsymOrSym(1:length(idxOfIniE));
    People(idxOfIniE,durAsymOrSym2ReCol) = iniDurAsymOrSym2Re(1:length(idxOfIniE));
    People(idxOfIniE,timeOfECol) = nDay;
    People(idxOfIniE,timeOfPCol) = nDay + iniDurE2P;
    People(idxOfIniE,timeOfAsymOrSymCol) = nDay + iniDurE2P + iniDurP2AsymOrSym(1:length(idxOfIniE));
    People(idxOfIniE,timeOfReCol) = nDay + iniDurE2P + iniDurP2AsymOrSym(1:length(idxOfIniE)) + iniDurAsymOrSym2Re(1:length(idxOfIniE));
    
    People(idxOfIniP,indPCol) = 1;
    People(idxOfIniP,durP2AsymOrSymCol) = iniDurP2AsymOrSym(length(idxOfIniE)+1:length(idxOfIniE)+length(idxOfIniP));
    People(idxOfIniP,durAsymOrSym2ReCol) = iniDurAsymOrSym2Re(length(idxOfIniE)+1:length(idxOfIniE)+length(idxOfIniP));
    People(idxOfIniP,timeOfPCol) = nDay;
    People(idxOfIniP,timeOfAsymOrSymCol) = nDay + iniDurP2AsymOrSym(length(idxOfIniE)+1:length(idxOfIniE)+length(idxOfIniP));
    People(idxOfIniP,timeOfReCol) = nDay + iniDurP2AsymOrSym(length(idxOfIniE)+1:length(idxOfIniE)+length(idxOfIniP)) + iniDurAsymOrSym2Re(length(idxOfIniE)+1:length(idxOfIniE)+length(idxOfIniP));
    
    People(idxOfIniAsym,indAsymCol) = 1;
    People(idxOfIniAsym,durAsymOrSym2ReCol) = iniDurAsymOrSym2Re(length(idxOfIniE)+length(idxOfIniP)+1:length(idxOfIniE)+length(idxOfIniP)+length(idxOfIniAsym));
    People(idxOfIniAsym,timeOfAsymOrSymCol) = nDay;
    People(idxOfIniAsym,timeOfReCol) = nDay + iniDurAsymOrSym2Re(length(idxOfIniE)+length(idxOfIniP)+1:length(idxOfIniE)+length(idxOfIniP)+length(idxOfIniAsym));
   
    People(idxOfIniSym,indSymCol) = 1;
    People(idxOfIniSym,durAsymOrSym2ReCol) = iniDurAsymOrSym2Re(length(idxOfIniE)+length(idxOfIniP)+length(idxOfIniAsym)+1:length(idxOfIniE)+length(idxOfIniP)+length(idxOfIniAsym)+length(idxOfIniSym));
    People(idxOfIniSym,timeOfAsymOrSymCol) = nDay;
    People(idxOfIniSym,timeOfReCol) = nDay + iniDurAsymOrSym2Re(length(idxOfIniE)+length(idxOfIniP)+length(idxOfIniAsym)+1:length(idxOfIniE)+length(idxOfIniP)+length(idxOfIniAsym)+length(idxOfIniSym));
end

function [idxContact, uidContact, idxInfected, uidInfected]  = fGetIdxUidOfContactAndInfected(People,cumPostPOfContactRate,idxAnduidOfDists,referRandNumOfDist,nContact,sourceDist,flagOfInfType,r,beta,indSCol)
    idxContact=zeros(nContact,1); 
    uidContact=zeros(nContact,1);
    idxInfected=[];
    uidInfected=[];
    for i=1:nContact
        greaterDist = find(cumPostPOfContactRate(:,sourceDist) >= referRandNumOfDist(i));
        contactDist = greaterDist(1);
        idxContactCandidate = idxAnduidOfDists(contactDist).index;
        uidContactCandidate = idxAnduidOfDists(contactDist).uid;
        idxRand = randperm(length(idxContactCandidate),1);
        idxContact(i) = idxContactCandidate(idxRand);
        uidContact(i) = uidContactCandidate(idxRand);
        if People(idxContact(i),indSCol)==1
            referRandNumOfInfect = rand();
            if flagOfInfType == 1    %P
                if referRandNumOfInfect < r * beta
                    idxInfected = [idxInfected; idxContact(i)];
                    uidInfected = [uidInfected; uidContact(i)];
                end
            elseif flagOfInfType == 2 %Asym
                if referRandNumOfInfect < r * beta
                    idxInfected = [idxInfected; idxContact(i)];
                    uidInfected = [uidInfected; uidContact(i)];
                end                
            else
                if referRandNumOfInfect < beta
                    idxInfected = [idxInfected; idxContact(i)];
                    uidInfected = [uidInfected; uidContact(i)];
                end
            end
        end
    end
    idxContact = unique(idxContact);
    uidContact = unique(uidContact);
    idxInfected = unique(idxInfected);
    uidInfected = unique(uidInfected);
end

function timeOfContracing = fGetContactTracingTime(kappaC,minDelayC,maxDelayC,nCaseC,thldC)
    timeOfContracing = minDelayC + (maxDelayC - minDelayC)*(1 - exp(-(nCaseC - thldC)/kappaC));
end

function timeOfTestDelay = fGetTestingDelay(kappaT,minDelayT,maxDelayT,nCaseT,thldT)
    timeOfTestDelay = minDelayT + (maxDelayT - minDelayT)*(1 - exp(-(nCaseT - thldT)/kappaT));
end