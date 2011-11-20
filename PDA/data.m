%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demo Probabilistic Data Association
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc
n=50;  
%Tempo Iniziale
T=1;    
%monte carlos run times                                               
MC_number=5;  
%STATO Iniziale: [x vx y vy]
target_position=[1500 500 1500  400];

PDAF(target_position,n,T,MC_number)