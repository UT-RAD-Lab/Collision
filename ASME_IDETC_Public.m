%2020.09.04 - This code is written to detect the number of collisions that
%occur as robots move around a domain. This code accompanies the
%article: Evaluation and Modification of Kinetic Gas Collision Theory as
%Applied to Encounter Rate Dynamics for Multi-Robot Groups and Robot Swarms
%by Adam Schroeder (adam.schroeder@utoledo.edu), Mark Rooney, and Glenn
%Lipscomb

tic
clc
clear
clear plot

%general settings
travel=0.22; %travel speed = units distance/second
diameter=.16; %diameter of robots, units distance
levyalpha=1.0; %levy alpha parameter for changing the distribution that is used to sample a path length
time=200; %simulation time, units second
smax=100; %number of statistical runs

Robots = [2 3 4 5 6 7 8 9 10 11 12 13 14 15]; %by default, the number of robots in this vector are few so it will run quickly
%Robots = [2 3 4 5 6 7 8 9 10 11 12 13 14 15 25 50 75 100]; %these where the robot quantities used in the paper
deltat=0.05; %the domain is checked for collisions every time step, units second

%search area boundaries
xmin=0; %units distance
xmax=3;
ymin=0;
ymax=3;
Asys= (xmax-xmin)*(ymax-ymin); %system area

r=1; %multiple run index
rmax=numel(Robots); %maximum number of runs

while r<=rmax
    travel_effective=travel*deltat; %adjusted to time step
    Collisions=zeros(1,time/deltat); %will track collisions
    
    A=1:Robots(r);
    
    %% Calculates all potential collision pairs
    Co = nchoosek(A,2);
    Collisions_previous=zeros(1,length(Co)); %start with no collisions previously
    
    s=1;
    
    
    %% Begins simulation with collision tracking
    while s<=smax;
        %initiate robots
        Nx=(xmax-diameter).*rand(1,Robots(r))+diameter/2; %initialize population randomly
        Ny=(xmax-diameter).*rand(1,Robots(r))+diameter/2; %initialize population randomly
        
        pathremaining=((time*travel)+1).*ones(1,Robots(r)); %will give remaining path length
        robotheading=zeros(1,Robots(r)); %will track robot headings
        for i=1:length(robotheading)
            robotheading(i)=2*pi*rand;
        end
        
        t=0; %time index for simulation time
        loops=0; %loop index for number of times through while loop
        
        while t<=time
            Nx_prev=Nx; %will store previous positions
            Ny_prev=Ny; %will store previous positions
            for k=1:numel(Nx); %now update robot positions
                if pathremaining(k)<=0 %if a new path needs calculated
                    angle_rand=2*pi*rand;
                    x_rand=cos(angle_rand);
                    y_rand=sin(angle_rand);
                    robotheading(k)=atan2(y_rand,x_rand); %heading is summation of noise and reverse gradient
                    
                    future_Nx=Nx(1,k)+travel_effective*cos(robotheading(k));
                    future_Ny=Ny(1,k)+travel_effective*sin(robotheading(k));
                    
                    Nx(1,k)=future_Nx;
                    Ny(1,k)=future_Ny;
                    pathremaining(k)=rand^(-1/levyalpha)-travel_effective;
                else %if a new path doesn't need calculated
                    future_Nx=Nx(1,k)+travel_effective*cos(robotheading(k));
                    future_Ny=Ny(1,k)+travel_effective*sin(robotheading(k));
                    Nx(1,k)=future_Nx;
                    Ny(1,k)=future_Ny;
                    pathremaining(k)=pathremaining(k)-travel_effective; %decrease length of path remaining
                end %if a new path needs created
            end %end for robot loop
            
            %% Enforces periodic boundaries
            for k=1:numel(Nx)
                if Nx(1,k)>xmax
                    Nx(1,k)=Nx(1,k)-xmax;
                else if Nx(1,k)<0
                        Nx(1,k)=Nx(1,k)+xmax;
                    end
                end
                if Ny(1,k)>ymax
                    Ny(1,k)=Ny(1,k)-ymax;
                else if Ny(1,k)<0
                        Ny(1,k)=Ny(1,k)+ymax;
                    end
                end
                %
            end %end for loop
            
            t=t+deltat; %moved from below
            loops=loops+1;
            
            %% Determines if and between whom collisions have occured
            Nx_Diff=abs(Nx(:,Co(:,1))-Nx(:,Co(:,2)));
            Ny_Diff=abs(Ny(:,Co(:,1))-Ny(:,Co(:,2)));
            Nx_Diff(Nx_Diff>(xmax-diameter))=Nx_Diff(Nx_Diff>(xmax-diameter))-xmax;
            Ny_Diff(Ny_Diff>(ymax-diameter))=Ny_Diff(Ny_Diff>(ymax-diameter))-ymax;
            N_Diff=(Nx_Diff.^2+Ny_Diff.^2).^0.5;
            
            Collisions_locations=N_Diff<diameter;
            Collisions_locations_new=Collisions_locations&Collisions_previous~=1;
            Collisions(loops)=sum((Collisions_locations_new),2);
            Collisions_previous=Collisions_locations;
            
            Collisions_robots=[Collisions_locations_new.*Co(:,1)' ;Collisions_locations_new.*Co(:,2)'];
            Collisions_robots_reduced=[nonzeros(Collisions_robots(1,:))';nonzeros(Collisions_robots(2,:))'];
        end
        
        Collisions_total=sum(Collisions);
        collisions_intermediate(r,s)=Collisions_total;
        s=s+1;
    end %while statistical run
    
    Collisions_averaged(1,r)=mean(collisions_intermediate(r,:))./time;
    Collisions_deviation(1,r)=std(collisions_intermediate(r,:))./time;
    Collisions_confidence(1,r)=1.96*Collisions_deviation(r)/smax^0.5;
    
    r=r+1;
end %while multiple run

figure(1)
hold on
xtheory=0:Robots(end);

ytheory=(diameter*4/pi)*travel*(xtheory-1).*xtheory/(Asys);
ytheory_alt_crel=(diameter*2^0.5)*travel*(xtheory-1).*xtheory/(Asys);
ytheory_alt_low=(diameter*4/pi)*travel*(xtheory).*xtheory/(Asys);
plot(xtheory, ytheory)
hold on
plot(xtheory, ytheory_alt_crel)
plot(xtheory, ytheory_alt_low)
errorbar(Robots,Collisions_averaged(1,:),Collisions_deviation(1,:),'o')
xlabel('Robots','FontSize',14)
ylabel('Group Collisions per Second','FontSize',14)
ylim([0 1.5]);
xlim([0 15]);
lgd2=legend('Analytical Constant Velocity','Analytical Boltzmann','Analytical without N-1 Adjustment','Numerical')
lgd2.FontSize = 14;

Collisions_per_Robot(1,:)=Collisions_averaged(1,:)./Robots;

toc
