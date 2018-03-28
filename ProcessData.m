
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

  RunID='RunC';

ni1Files=dir('ni1_*.dat');
ni2Files=dir('ni2_*.dat');
BzFiles=dir('Bz_*.dat');
vi1xFiles=dir('vi1x_*.dat');
vi1yFiles=dir('vi1y_*.dat');
vi2xFiles=dir('vi2x_*.dat');
vi2yFiles=dir('vi2y_*.dat');
x_vxFiles=dir('x_vx_*.dat');

nFiles=size(ni1Files,1);

for i=1:nFiles-1
    disp(strcat('Processing file ',num2str(i),'...'))
    ni1_local=load(ni1Files(i).name);
    ni2_local=load(ni2Files(i).name);
    Bz_local=load(BzFiles(i).name);
    vi1x_local=load(vi1xFiles(i).name);
    vi1y_local=load(vi1yFiles(i).name);
    vi2x_local=load(vi2xFiles(i).name);
    vi2y_local=load(vi2yFiles(i).name);
    
    x_vx=load(x_vxFiles(i).name);
    x_frame=x_vx(1);
    x0=0;
    sh=round((x0+x_frame)/dx);
  
    x=((1:TotalNx)-TotalNx/2)*dx-x0;
    y=((1:TotalNy)-TotalNy/2)*dy;
    
    ni1_local=reshape(ni1_local,Nx,Ny,dim_y,dim_x);
    ni2_local=reshape(ni2_local,Nx,Ny,dim_y,dim_x);
    Bz_local=reshape(Bz_local,Nx,Ny,dim_y,dim_x);
    vi1x_local=reshape(vi1x_local,Nx,Ny,dim_y,dim_x);
    vi1y_local=reshape(vi1y_local,Nx,Ny,dim_y,dim_x);
    vi2x_local=reshape(vi2x_local,Nx,Ny,dim_y,dim_x);
    vi2y_local=reshape(vi2y_local,Nx,Ny,dim_y,dim_x);
    
    ni1=zeros(TotalNx,TotalNy);
    ni2=zeros(TotalNx,TotalNy);
    Bz=zeros(TotalNx,TotalNy);
    vi1x=zeros(TotalNx,TotalNy);
    vi1y=zeros(TotalNx,TotalNy);
    vi2x=zeros(TotalNx,TotalNy);
    vi2y=zeros(TotalNx,TotalNy);
    
    for i2=1:dim_y
        for i1=1:dim_x
            ni1((i1-1)*Nx+1:i1*Nx,(i2-1)*Ny+1:i2*Ny)=ni1_local(:,:,i2,i1);
            ni2((i1-1)*Nx+1:i1*Nx,(i2-1)*Ny+1:i2*Ny)=ni2_local(:,:,i2,i1);
            Bz((i1-1)*Nx+1:i1*Nx,(i2-1)*Ny+1:i2*Ny)=Bz_local(:,:,i2,i1);
            vi1x((i1-1)*Nx+1:i1*Nx,(i2-1)*Ny+1:i2*Ny)=vi1x_local(:,:,i2,i1);
            vi1y((i1-1)*Nx+1:i1*Nx,(i2-1)*Ny+1:i2*Ny)=vi1y_local(:,:,i2,i1);
            vi2x((i1-1)*Nx+1:i1*Nx,(i2-1)*Ny+1:i2*Ny)=vi2x_local(:,:,i2,i1);
            vi2y((i1-1)*Nx+1:i1*Nx,(i2-1)*Ny+1:i2*Ny)=vi2y_local(:,:,i2,i1);
        end
    end
    
    ni1=circshift(ni1,[sh 0]);
    ni2=circshift(ni2,[sh 0]);
    Bz=circshift(Bz,[sh 0]);
    vi1x=circshift(vi1x,[sh 0]);
    vi1y=circshift(vi1y,[sh 0]);
    vi2x=circshift(vi2x,[sh 0]);
    vi2y=circshift(vi2y,[sh 0]);
    
    w1=(circshift(vi1y,[1 0])-circshift(vi1y,[-1 0]))/(2*dx)...
       -(circshift(vi1x,[0 1])-circshift(vi1x,[0 -1]))/(2*dy);
    w2=(circshift(vi2y,[1 0])-circshift(vi2y,[-1 0]))/(2*dx)...
       -(circshift(vi2x,[0 1])-circshift(vi2x,[0 -1]))/(2*dy);
   
    div1=(circshift(vi1y,[0 1])-circshift(vi1y,[0 -1]))/(2*dy)...
       +(circshift(vi1x,[1 0])-circshift(vi1x,[-1 0]))/(2*dx);
   
    div2=(circshift(vi2y,[0 1])-circshift(vi2y,[0 -1]))/(2*dy)...
       +(circshift(vi2x,[1 0])-circshift(vi2x,[-1 0]))/(2*dx);

    n=ni1+ni2;
gradn_x=(circshift(n,[1 0])-circshift(n,[-1 0]))/(2*dx);
gradn_y=(circshift(n,[0 1])-circshift(n,[0 -1]))/(2*dy);

    mu0=4*pi*1e-7;
    e=1.602e-19;
    BSq=Bz.^2/(2*mu0);
    gradBSq_x=(circshift(BSq,[1 0])-circshift(BSq,[-1 0]))/(2*dx);
    gradBSq_y=(circshift(BSq,[0 1])-circshift(BSq,[0 -1]))/(2*dy);
   
    poisson=-1./(n.^2*e).*(gradn_x.*gradBSq_y-gradn_y.*gradBSq_x);
    curBx=Bz.*(ni1.*vi1x+ni2.*vi2x)./(ni1+ni2);
    curBy=Bz.*(ni1.*vi1y+ni2.*vi2y)./(ni1+ni2);
    divCurB=(circshift(curBx,[1 0])-circshift(curBx,[-1 0]))/(2*dx)...
        +(circshift(curBy,[0 1])-circshift(curBy,[0 -1]))/(2*dy);
    
    curx=(ni1.*vi1x+ni2.*vi2x)./(ni1+ni2);
    cury=(ni1.*vi1y+ni2.*vi2y)./(ni1+ni2);
    divCur=(circshift(curx,[1 0])-circshift(curx,[-1 0]))/(2*dx)...
        +(circshift(cury,[0 1])-circshift(cury,[0 -1]))/(2*dy);
    
    frameID=i;
    time=i*dt_print;
    
    fileName=strcat(RunID,'_',num2str(i,'%04i'),'.mat');
    
    save(fileName,'RunID','x','y','frameID','time','ni1','ni2','Bz','vi1x','vi1y','vi2x','vi2y','w1','w2','div1','div2',...
	 'BSq','gradBSq_x','gradBSq_y','poisson','curBx','curBy','divCurB','curx','cury','divCur','gradn_x','gradn_y');
    
    
end


