%% Make domain
Xloc = [xmin xmax]; %lon
Yloc = [ymin ymax]; %lat



%% Initialise data matrix for ECCO2 data.
theta = nan(12*(MaxYear-MinYear+1),50,xmax-xmin+1,ymax-ymin+1);
    
%Cycle through correct naming convention for cube92
year_ind = MinYear:MaxYear;
mon_ind = 1:12;
loop_ind=1;
clear monstr
for yy=MinYear:MaxYear; %number of years
 for mm=1:12 %number of months
  if mm<10
  monstr(loop_ind)=str2num([num2str(yy),'0',num2str(mm)]);
  elseif mm>=10
  monstr(loop_ind)=str2num([num2str(yy),num2str(mm)]);
  end
  loop_ind=loop_ind+1;
  end
end


for TimeInd = 1:12*(MaxYear-MinYear+1)

theta(TimeInd,:,:,:) = permute(ncread([external_dir,'/ecco2/THETA.nc/THETA.1440x720x50.' num2str(monstr(TimeInd)) '.nc'],'THETA',[Xloc(1) Yloc(1) 1 1],[Xloc(2)-Xloc(1)+1 Yloc(2)-Yloc(1)+1 Inf Inf]),[3 1 2]);
disp([num2str(TimeInd) ' of ' num2str(12*(MaxYear-MinYear+1)) ' month done.'])
end

save(['cube92_iaf_theta_',RunName,'.mat'],'theta','-v7.3')

