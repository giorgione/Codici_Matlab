%Test Modello Grafico
%Genero un grafo NON ORIENTATO fully connected
clc;clear; close all
   
A=[0 0 0 0 1 0 0 ;
   0 0 0 0 1 1 1 ;
   0 0 0 0 0 1 0 ;
   0 0 0 0 0 0 1 ;
   1 1 0 0 0 0 0 ;
   0 1 1 0 0 0 0 ;
   0 1 0 1 0 0 0 ];
   

G=GraficalModel(7,A);
G=G.setNodes(1:4,5:7);

G=G.DisegnaGrafo();
G=G.setAnimazione(1);
[Res G]=G.ForwardMsg(3)


