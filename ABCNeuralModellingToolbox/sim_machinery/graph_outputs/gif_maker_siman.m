function [] = gif_maker_siman(R,d)
if nargin<2
    d = R.d;
    d = sprintf('%d',[d(1:3)]);
end
% d = '2017612';
% This program creates a movie/slideshow from a set of images, and save it as an animated GIF file.
% Notice that the quality an image may decrease due to the GIF format.
%
% Written by Moshe Lindner , Bar-Ilan University, Israel.
% September 2010 (C)
for K = 1:3
    if K == 1
        file_path = [R.rootn 'outputs\' R.out.tag '\' R.out.dag '\dist_track\'];
    elseif K == 2
        file_path = [R.rootn 'outputs\' R.out.tag '\' R.out.dag '\r2track\'];
    elseif K == 3
        file_path = [R.rootn 'outputs\' R.out.tag '\' R.out.dag '\feattrack\'];
    end
    a = dir([file_path,'\*.jpg']);
    a = {a.name};
    file_name=sort_nat(a);
    
    if K== 1
        file_path2 = [R.rootn 'outputs\' R.out.tag '\' R.out.dag '\gifs\'];
        file_name2 = [R.out.tag '_dist_track_gif.gif'];
    elseif K == 2
        file_path2 = [R.rootn 'outputs\' R.out.tag '\' R.out.dag '\gifs\'];
        file_name2 = [R.out.tag '_r2track_gif.gif'];
    elseif K == 3
        file_path2 = [R.rootn 'outputs\' R.out.tag '\' R.out.dag '\gifs\'];
        file_name2 = [R.out.tag '_feat_track_gif.gif'];
    end
    mkdir(file_path2)
    loops = R.plot.gif.loops;
    delay=R.plot.gif.delay;
    
    delay1=R.plot.gif.start_t;
    
    delay2=R.plot.gif.end_t;
    
    h = waitbar(0,['0% done'],'name','Progress') ;
    for i=1:length(file_name)
        if strcmpi('gif',file_name{i}(end-2:end))
            [M  c_map]=imread([file_path,file_name{i}]);
        else
            a=imread([file_path,file_name{i}]);
            [M  c_map]= rgb2ind(a,256);
        end
        if i==1
            imwrite(M,c_map,[file_path2,file_name2],'gif','LoopCount',loops,'DelayTime',delay1)
        elseif i==length(file_name)
            imwrite(M,c_map,[file_path2,file_name2],'gif','WriteMode','append','DelayTime',delay2)
        else
            imwrite(M,c_map,[file_path2,file_name2],'gif','WriteMode','append','DelayTime',delay)
        end
        waitbar(i/length(file_name),h,[num2str(round(100*i/length(file_name))),'% done']) ;
    end
    close(h);
    %     msgbox('Finished Successfully!')
end