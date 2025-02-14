% Make Figure

clear all;
close all

% import data

IM=importdata('GlobalBalance.dat','\t',1);
B=IM.data;
std = B(:,1);
emi = B(:,2);
dep = B(:,5)+B(:,6);
esc = B(:,7);

redox = importdata('../../../MoleRedox.dat');
red = redox(std);

redin = red'*emi;
redout = red'*dep + red'*esc;