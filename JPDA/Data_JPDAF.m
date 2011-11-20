 

clc;clear all;close all;
%Numero di Evoluzioni  dello stato sistema che desidero stimare
n=50;

%Time
T=1;

%Monte Carlo 
MC_number=5;     

%Monte Carlo 
%Dimensione Stato
c=2;  

%Posizioni Iniziali
target_position=[1500 300 500 400; 
                 500 400 1500 300];                                              
JPDAF(target_position,n,T,MC_number,c);                                           