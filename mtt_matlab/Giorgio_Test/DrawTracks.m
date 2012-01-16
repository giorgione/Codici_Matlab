function [Track]=DrawTracks(Track,Myfig)
    
figure(Myfig); hold on
    for j = 1:length(Track)
        Tempo=Track(j).det;
        if isempty(Track(j).colore)
            Track.colore=rand(3,1);
        end
       
        plot(Tempo,Track(j).tid*ones(size(Tempo)),Track(j).colore)
       
    end
end