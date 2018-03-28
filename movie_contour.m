  runName='RunC';
  files=dir('*.mat');

valVec = [];
tVec = [];

  %Processor geometry
  dim_x=30;
  dim_y=20;

  % Total number of grid points
  TotalNx=600;
  TotalNy=600;

  %Print interval (seconds)
  dt_print=0.01d0;

  % The domain in x and y directions
  x_start=-6000.0d0;
  x_end=6000.0d0;
  y_start=-6000.0d0;
  y_end=6000.0d0;

  % Calculate grid sizes and local nr of grid points
  dx=(x_end-x_start)/TotalNx;
  Nx=TotalNx/dim_x;
  dy=(y_end-y_start)/TotalNy;
  Ny=TotalNy/dim_y;


fid = figure;
set(fid,'Position', [1 1  400  550])
set(fid,'PaperPosition',[5.5    7.0    10.0   15.66])

%% Set up the movie.
writerObj = VideoWriter(strcat(runName,'.avi')); % Name it.
writerObj.FrameRate = 10; % How many frames per second.
open(writerObj); 
%%


for file = files'
  fileName=file.name;
  load(fileName);
  
  xRed=x;
  yRed=y;
  ni2Red=ni2;
  ni1Red=ni1;

  
  f1=subplot(3,1,1);
%  ni2=ifft(exp(1i*(-x0-v0*t)*repmat(k.',1,Ny)).*fft(ni2,[],1),[],1);
  pcolor(xRed/1000,yRed/1000,ni2Red'/1e6);
  shading flat;
  f2=colorbar;
  axis equal;axis tight
%  xlabel('x (km)')
  ylabel('y (km)')
  title(['n_{i2} (cm^{-3})  t=' num2str(time,'%3.2f') ' s'],'FontWeight','normal')
  set(f1,'Position',[0.060 0.11+0.62 0.671 0.20])
  set(f2,'Position',[0.75 0.11+0.62 0.040 0.20])
  
  f1=subplot(3,1,2);
  
  pcolor(xRed/1000,yRed/1000,log10(abs(ni2Red'/1e6/1e5)+1e-5));
  shading flat;
  caxis([-2 2])
  f2=colorbar;
  hold on
  contour(xRed/1000,yRed/1000,ni2Red',[1e11 1e11],'r','Linewidth',1);

%  contourf(x/1000,y/1000,log(ni2'),log([1e12 1e12]),'b');
%  hold on
%  contourf(x/1000,y/1000,log(ni2'),log([1e13 1e13]),'k');
%  contourf(x/1000,y/1000,log(ni2'),log([1e14 1e14]),'r');
  
%  shading flat;
  axis equal;axis tight
%  xlabel('x (km)')
  ylabel('y (km)')
  title(['log_{10}(n_{i2}/10^{5})'],'FontWeight','normal')
  set(f1,'Position',[0.060 0.11+0.31 0.671 0.20])
  set(f2,'Position',[0.75 0.11+0.31 0.040 0.20])
  
  f1=subplot(3,1,3);
%  ni1=ifft(exp(1i*(-x0-v0*t)*repmat(k.',1,Ny)).*fft(ni1,[],1),[],1);
%  pcolor(xRed/1000,yRed/1000,ni1Red'/1e6);
		line(xRed/1000,ni2(TotalNy/2,:));
%  shading flat;
%  f2=colorbar;
%  axis equal;axis tight
  xlabel('x (km)')
  ylabel('n_{i2} (cm^{-3})')
  title(['n_{i2} Line Plot'],'FontWeight','normal')
  set(f1,'Position',[0.060 0.11 0.671 0.20])
  set(f2,'Position',[0.75 0.11 0.040 0.20])

		ni2Cut=ni2(TotalNy/2,TotalNx/2:TotalNx);
		[val,idx]=max(ni2Cut);
		valVec = [valVec xRed(idx+TotalNx/2)];
		tVec = [tVec time];
  pause(0.01);
 
  frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
  writeVideo(writerObj, frame);
%  return
end

close(writerObj); % Saves the movie.

		dlmwrite(strcat(runName,'Bounce.txt'),[transpose(tVec) transpose(valVec)],'delimiter','\t','precision',4);
		
		theFigure=figure;
		set(theFigure,'visible','off');
		
		line(tVec,valVec);
		title(strcat(runName,' Bounce Frequency'));
		xlabel('Time (s)');
		ylabel('Radius (m)');
		
		saveas(theFigure,strcat(runName,'Bounce.jpg'));
		saveas(theFigure,strcat(runName,'Bounce.fig'));
		clf(theFigure);

		dlmwrite(strcat(runName,'Bounce.txt'),[transpose(tVec) transpose(valVec)],'delimiter','\t','precision',4);
