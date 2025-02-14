% Make Figure

clear all;
close all

% import data

comp = importdata('../../../MoleCompo.dat');
H = comp(:,1);
O = comp(:,2);
N = comp(:,3);
C = comp(:,4);
S = comp(:,5);

C = [0;C];
N = [0;N];
S = [0;S];
O = [0;O];
H = [0;H];

R = importdata('../../Data/Reaction_R.txt');
M = importdata('../../Data/Reaction_M.txt');
P = importdata('../../Data/Reaction_P.txt');
T = importdata('../../Data/Reaction_T.txt');

CBR = C(R(:,1)+1)+C(R(:,2)+1)-C(R(:,4)+1)-C(R(:,5)+1)-C(R(:,6)+1);

CBM = C(M(:,1)+1)+C(M(:,2)+1)-C(M(:,3)+1)-C(M(:,4)+1);

CBP = C(P(:,1)+1)-C(P(:,2)+1)-C(P(:,3)+1)-C(P(:,4)+1);

CBT = C(T(:,1)+1)-C(T(:,2)+1)-C(T(:,3)+1);

CB = sum(CBR) + sum(CBM) + sum(CBP) + sum(CBT);

SBR = S(R(:,1)+1)+S(R(:,2)+1)-S(R(:,4)+1)-S(R(:,5)+1)-S(R(:,6)+1);

SBM = S(M(:,1)+1)+S(M(:,2)+1)-S(M(:,3)+1)-S(M(:,4)+1);

SBP = S(P(:,1)+1)-S(P(:,2)+1)-S(P(:,3)+1)-S(P(:,4)+1);

SBT = S(T(:,1)+1)-S(T(:,2)+1)-S(T(:,3)+1);

SB = sum(SBR) + sum(SBM) + sum(SBP) + sum(SBT);


NBR = N(R(:,1)+1)+N(R(:,2)+1)-N(R(:,4)+1)-N(R(:,5)+1)-N(R(:,6)+1);

NBM = N(M(:,1)+1)+N(M(:,2)+1)-N(M(:,3)+1)-N(M(:,4)+1);

NBP = N(P(:,1)+1)-N(P(:,2)+1)-N(P(:,3)+1)-N(P(:,4)+1);

NBT = N(T(:,1)+1)-N(T(:,2)+1)-N(T(:,3)+1);

NB = sum(NBR) + sum(NBM) + sum(NBP) + sum(NBT);



OBR = O(R(:,1)+1)+O(R(:,2)+1)-O(R(:,4)+1)-O(R(:,5)+1)-O(R(:,6)+1);

OBM = O(M(:,1)+1)+O(M(:,2)+1)-O(M(:,3)+1)-O(M(:,4)+1);

OBP = O(P(:,1)+1)-O(P(:,2)+1)-O(P(:,3)+1)-O(P(:,4)+1);

OBT = O(T(:,1)+1)-O(T(:,2)+1)-O(T(:,3)+1);

OB = sum(OBR) + sum(OBM) + sum(OBP) + sum(OBT);


HBR = H(R(:,1)+1)+H(R(:,2)+1)-H(R(:,4)+1)-H(R(:,5)+1)-H(R(:,6)+1);

HBM = H(M(:,1)+1)+H(M(:,2)+1)-H(M(:,3)+1)-H(M(:,4)+1);

HBP = H(P(:,1)+1)-H(P(:,2)+1)-H(P(:,3)+1)-H(P(:,4)+1);

HBT = H(T(:,1)+1)-H(T(:,2)+1)-H(T(:,3)+1);

HB = sum(HBR) + sum(HBM) + sum(HBP) + sum(HBT);
