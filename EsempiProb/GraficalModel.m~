classdef GraficalModel
   properties
      AdjMat=0;
      Nodes=0;
      PotentialFun=0;
      FactorNodes=[];
      VariableNodes=[];
      ActiveNodes=[];
      MsgTable;
      %Handles della finestra contenete il Grafo
      MyHnodi=[];
      MyHarchi=[];
      animazione=0;
   end % properties
   methods
      function obj = GraficalModel(N,Adj)
         obj.Nodes=N;
         obj.AdjMat=Adj;
      end % set.Material
       
      function obj = setNodes(obj,Variable,Factor)
         
         obj.FactorNodes=Factor;
         obj.VariableNodes=Variable;
         obj.ActiveNodes=1:obj.Nodes;
         %Matrice ( FactorNodes x VariableNodes ) contenete i msg 
         %MsgTable(f,x)= u f->x
         %MsgTable(x,f)
         obj.MsgTable=cell(obj.Nodes);
         
      end % set.Material
      
      function obj = setAnimazione(obj,val)
         
          obj.animazione=val;
         
      end 
      
      
      
      function obj=DisegnaGrafo(obj)         
        clf
       
        if(isempty(obj.FactorNodes))
             names=cell(1,obj.Nodes);
            for i=1:obj.Nodes
                names(i)={[ 'X' num2str(i)]};
            end
            [x,y,h]=draw_layout(obj.AdjMat,names,zeros(1,obj.Nodes));
            obj.MyH=h;
        else
            
            nF=length(obj.FactorNodes)
            nV=length(obj.VariableNodes)
            namesF=cell(1,nF);
            namesV=cell(1,nV);
            for i=1:nV
                namesV(i)={[ 'X' num2str(i)]};
            end
            
            for i=1:nF 
                namesF(i)={[ 'F' num2str(obj.FactorNodes(i))]};
            end
            
            [x,y,h,Archi]=draw_layout(obj.AdjMat,[namesV namesF],[zeros(1,nV) ones(1,nF)]);
            obj.MyHnodi=h; 
            obj.MyHarchi=Archi; 
        end
        
      end
      
      
      %Ritorna l'indice del nodo padre Xi connesso alla foglia
      function [res index]=IsLeaf(obj,Xi)
        [j Xj]=find(obj.AdjMat(Xi,:) > 0);
        if(length(Xj) ==1)
            res=1;
            index=Xj;
        else
            res=0;
            index=0;
        end
      end
      
      %Ritorna i nodi vicini del Nodo Xi
      function [Ni]=Neighborhood(obj,Xi)
        [i Ni]=find( obj.AdjMat(Xi,:) > 0);
        
      end
      
      %Ritorna il Tipo del Nodo Xi
      function [type]=NodeType(obj,Xi)
        if(ismeber(Xi,obj.FactorNodes))
            type='F'; %nodo Variabile
        else
            type='V'; %nodo Fattore
        end
      end

      %MSG Node --> Factor  
      function [Res obj]=Ux_f(obj,X,F)
         
        if IsLeaf(obj,X)
            Val=1;
            %elimino il Nodo Var dalla lista dei Nodi Attivi
            obj.ActiveNodes=setdiff(obj.ActiveNodes,X);
            obj.MsgTable(X,F)={'1'};
            Res='1';
            if(obj.animazione)
                    %Coloro l'arco attivo di verde
                    set(obj.MyHnodi(X,2), 'FaceColor','g');
                    pause;
                    set(obj.MyHarchi(X,F), 'Color','g');
                    pause;
                end
        else
            %Recupero tutti i MSG ENTRANTI IN X dai Fattori: Ne(X)\F=
            Ne=obj.Neighborhood(X);
            Ne=setdiff(Ne,F);
            Res='';
            i=1;
            
            if(obj.animazione)
                %Coloro l'arco attivo di rosso
                set(obj.MyHarchi(X,F), 'Color','r');       
                set(obj.MyHnodi(F,2), 'FaceColor','r');
                pause;
                
            end
            
            
            for fl=Ne
                
                if(obj.animazione)
                    %Coloro l'arco attivo di verde
                    set(obj.MyHarchi(X,fl), 'Color','r');
                    pause;
                    set(obj.MyHnodi(fl,2), 'FaceColor','r');                    
                    pause;
                end
                
                [R obj]=Uf_x(obj,X,fl);
                if(i==1)
                    Res=[Res ' ' R];
                else
                     Res=[Res '*' R ];
                end
                i=i+1;
                
                
                if(obj.animazione)
                    %Coloro l'arco attivo di verde
                    set(obj.MyHnodi(fl,2), 'FaceColor','g');                    
                    pause;
                    set(obj.MyHarchi(X,fl), 'Color','g');
                    pause;
                    
                end
                
            end
            
            obj.MsgTable(X,F)={Res};
            if(obj.animazione)
                %Coloro l'arco attivo di rosso
                set(obj.MyHarchi(X,F), 'Color','g');       
                set(obj.MyHnodi(X,2), 'FaceColor','g');
                pause;                
            end
            
        end
      end
       
      %MSG Factor --> Node   : F --> X  
      function [Res obj]=Uf_x(obj,X,F)
         
        if IsLeaf(obj,F)
            obj.MsgTable(X,F)={[ 'f' num2str(F) '( x' num2str(X) ')']};
            %elimino il Nodo Fattore dalla lista dei Nodi Attivi
            obj.ActiveNodes=setdiff(obj.ActiveNodes,F);
            Res=obj.MsgTable(X,F);
            if(obj.animazione)
                    %Coloro l'arco attivo di verde
                    set(obj.MyHarchi(X,F), 'Color','g');
                    
                    set(obj.MyHnodi(F,2), 'FaceColor','g');
                    
                    pause;
                end
            
        else
            %Recupero tutti i MSG ENTRANTI IN X dai Nodi Variabili escluso X: Ne(F)\X=
            Ne=obj.Neighborhood(F);
            Ne=setdiff(Ne,X);
            Msg='';
            Res=['Sum( F' num2str(F) '(X' num2str(X) ','];
            Xs='';
            i=1;
            
            if(obj.animazione)
                %Coloro l'arco attivo di rosso
                set(obj.MyHarchi(X,F), 'Color','r');       
                set(obj.MyHnodi(F,2), 'FaceColor','r');
                pause;
                
            end
            
            %Nodi Variabile Vicini di F
            for Xm=Ne
             
                if(obj.animazione)
                    %Coloro l'arco attivo di rosso
                    set(obj.MyHarchi(Xm,F), 'Color','r');
                    pause;
                    set(obj.MyHnodi(Xm,2), 'FaceColor','r');
                    pause;
                end
                
                [R obj]=Ux_f(obj,Xm,F);
                if(i==1)
                    Msg=[Msg ' ' R];
                    %Marginal=[Marginal num2str(Xm)];
                    
                    Xs=[Xs 'X' num2str(Xm) ];
                else
                     Msg=[Msg '*' R];
                     %Marginal=[Marginal ',' num2str(Xm)];
                     
                     Xs=[Xs ',X' num2str(Xm)];
                end
                
                if(obj.animazione)
                    set(obj.MyHnodi(F,2), 'FaceColor','g');                    
                    pause;
                    %Coloro l'arco attivo di verde
                    set(obj.MyHarchi(Xm,F), 'Color','g');
                    pause;
                    
                end
                
                
                i=i+1;
            end 
            Res=[Res Xs ')*' Msg ',' Xs ')' ];
            obj.MsgTable(Xm,F)={Res};
            
            if(obj.animazione)
                %Coloro l'arco attivo di rosso
                set(obj.MyHarchi(Xm,F), 'Color','g');       
                set(obj.MyHnodi(Xm,2), 'FaceColor','g');
                pause;                
            end
            
        end
      end
      
   %Fattorizzo risolvendo le Dipendenze dal nodo Radice
   function [Res obj]=ForwardMsg(obj,X)
       %Cerco tutti i nodi Foglia eccetto X
       obj.ActiveNodes=1:obj.Nodes;
       obj.ActiveNodes=setdiff(obj.ActiveNodes,X);
       
       Ne=obj.Neighborhood(X);
       Res='';
       
       if(obj.animazione)
               %Coloro il nodo attivo di verde
               set(obj.MyHnodi(X,2), 'FaceColor','r');
               pause;
       end
       
       for fl=Ne
        if(obj.animazione)
                %Coloro l'arco attivo di rosso
                set(obj.MyHarchi(X,fl), 'Color','r');       
                set(obj.MyHnodi(fl,2), 'FaceColor','r');
                pause;
                
            end  
           
           [R obj]=Uf_x(obj,X,fl);           
           Res=[Res '*' R];
           if(obj.animazione)
                    %Coloro l'arco attivo di verde
                    set(obj.MyHarchi(X,fl), 'Color','g');
                    pause;
                    set(obj.MyHnodi(fl,2), 'FaceColor','g');                    
                    pause;
           end
           
       end
       
       if(obj.animazione)
               %Coloro il nodo attivo di verde
               set(obj.MyHnodi(X,2), 'FaceColor','g');
               pause;
       end
       
       
   end
      
   end% methods
end% classdef

