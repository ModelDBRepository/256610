
%%%%% Model for the Delayed Recognition Span Task in the spatial condition
%%%%% Based on the Delayed Response Task setup + short-term (calcium-mediated) synaptic facilitation

%%%%% The output shows examples of young, middle-aged and aged simulations (vce=5,7,9 respectively).

%%% Created by Sara Ibañez based on:
%%%  Wimmer K, Nykamp DQ, Constantinidis C, Compte A. Nat Neurosci. 17(3):431-439 (2014).

clear 

%%%%% PARAMETERS
Ne=640;                 % number of "excitatory cells" in the rate model
Ni=160;                 % number of "inhibitory cells" in the rate model
	
dt=2;                   % integration step in ms
tauE=20;                % time constant of rate equation for excitatory cells in ms
tauI=10;                % time constant of rate equation for inhibitory cells in ms
kappa=20;               % parameter defining concentration of e-to-e connectivity

GEEa=22.963;            % strength of excitation to excitatory neurons (AMPA receptors)
GEEn=55.297;            % strength of excitation to excitatory neurons (NMDA receptors)
GIE=25.19;              % strength of inhibition to excitatory neurons (GABAa receptors)
GEIa=118;               % strength of excitation to inhibitory neurons (AMPA receptors)
GEIn=141.5;             % strength of excitation to inhibitory neurons (NMDA receptors)
GII=205.54;             % strength of inhibition to inhibitory neurons (GABAa receptors)

sigE=1;                 % standard deviation of additive noise in rate equation of excitatory cells
sigI=3;                 % standard deviation of additive noise in rate equation of inhibitory cells
	
stim = 40000;           % strength of external stimulus

stimon = 1000;          % time when the first stimulus is applied (at 0 degrees) in ms
stimoff = 1500;         % time when the first stimulus ceases in ms
delayend = 3500;        % time when the first delay ends in ms
stimon2 = 3501;         % time when the first & second stimuli are applied (at 0 & 90 degrees) in ms
stimoff2 = 4001;        % time when the first & second stimuli cease in ms  	
delayend2 = 6001;       % time when the second delay ends in ms	
stimon3 = 6002;         % time when the first, second & third stimuli are applied (at 0, 90 & -90 degrees) in ms
stimoff3 = 6502;	  	% time when the first, second & third stimuli cease in ms 
delayend3 = 8502;       % time when the third delay ends in ms	
stimon4 = 8503;         % time when the first, second, third & fourth stimuli are applied (at 0, 90, -90 & 180 degrees) in ms
stimoff4 = 9003;        % time when the first, second, third & fourth stimuli cease in ms 
totalTime=stimoff4;     % total time of the simulation in ms



%%%%% PRELIMINARY CALCULATIONS
rE=zeros(Ne,1);
rI=zeros(Ni,1);
nsteps=floor(totalTime/dt);
tt=dt*(1:nsteps);

% E-to-E connectivity
theta0 = (0:Ne-1)/Ne*2*pi;
v0 = exp(kappa*cos(theta0));
v0 = v0/sum(v0);
WE = gallery('circul',v0);

% stimulus parameters
theta=theta0-pi+0.8;

v1 = exp(kappa*cos(theta));
v1 = v1/sum(v1);
stimulus1 = stim*v1';

v2 = exp(kappa*sin(theta));
v2 = v2/sum(v2);
stimulus2 = stim*v2';

v3 = exp(kappa*sin(-theta));
v3 = v3/sum(v3);
stimulus3 = stim*v3';

v4 = exp(kappa*-cos(theta));
v4 = v4/sum(v4);
stimulus4 = stim*v4';

stimulus_sum1 = stimulus1 + stimulus2;
stimulus_sum2 = stimulus1 + stimulus2 + stimulus3;
stimulus_sum3 = stimulus1 + stimulus2 + stimulus3 + stimulus4;

stimon = floor(stimon/dt);
stimoff = floor(stimoff/dt);
delayend = floor(delayend/dt);
stimon2 = floor(stimon2/dt);
stimoff2 = floor(stimoff2/dt);
delayend2 = floor(delayend2/dt);
stimon3 = floor(stimon3/dt);
stimoff3 = floor(stimoff3/dt);
delayend3 = floor(delayend3/dt);
stimon4 = floor(stimon4/dt);
stimoff4 = floor(stimoff4/dt);

% external bias current
I0E=80;  %pA	  
I0I=15;  %pA

% input-output function for excitatory and inhibitory cells, as used previously (used previously in: Brunel, Cereb Cortex 13:1151 (2003))
fE = inline('vce/(98^2)*x.*x.*(x>0).*(x<98)+vce*sqrt((4/98)*x-3).*(x>=98)','x','vce');
fI = inline('50/(20^2)*x.*x.*(x>0).*(x<20)+50*sqrt((4/20)*x-3).*(x>=20)','x');

% population vector decoder given the rates r for neurons with selectivity th
decode = inline('atan2(sum(r.*sin(th)),sum(r.*cos(th)))','r','th');

% Short-term synaptic facilitation parameters 
UE=0.001;       % baseline value for the utilization parameter
tf=1500;        % decaying time to the baseline value of the utilization parameter in ms
    


%%%%%% SIMULATION LOOP
RE=zeros(Ne,nsteps);
centerBang=zeros(nsteps,1); 
u=zeros(Ne,nsteps);

vce=[5,7,9];

for j=1:3
    
for i=1:nsteps
    % additive noise for each population
    noiseE = sigE*randn(Ne,1);
    noiseI = sigI*randn(Ni,1);
        
    if i==1
        IEa=0;
        IEn=0;
        IEg=0;
        IIa=0;
        IIn=0;
        IIg=0;
    else           
        IEa=IEa+(-IEa+GEEa*WE*uE.*rE)*dt/2;
        IEn=IEn+(-IEn+GEEn*WE*uE.*rE)*dt/100;
        IEg=IEg+(-IEg-GIE*mean(rI)*ones(Ne,1))*dt/10;
        IIa=IIa+(-IIa+GEIa*mean(rE)*ones(Ni,1))*dt/2;
        IIn=IIn+(-IIn+GEIn*mean(rE)*ones(Ni,1))*dt/100;
        IIg=IIg+(-IIg-GII*mean(rI)*ones(Ni,1))*dt/10;  
    end
    IE=IEa+IEn+IEg+I0E*ones(Ne,1);
    II=IIa+IIn+IIg+I0I*ones(Ni,1);

    if i>stimon && i<stimoff
        IE=IE+stimulus1;            % first stimulus presentation 
    end      
    if i>stimon2 && i<stimoff2
        IE=IE+stimulus_sum1;        % first & second stimuli presentation
    end
    if i>stimon3 && i<stimoff3
        IE=IE+stimulus_sum2;        % first, second & third stimuli presentation
    end
    if i>stimon4 && i<stimoff4
        IE=IE+stimulus_sum3;        % first, second, third & fourth stimuli presentation
    end
        
    % integration with time-step dt: Newton method
    rE = rE + (fE(IE,vce(j)) - rE + noiseE)*dt/tauE;        
    rI = rI + (fI(II) - rI + noiseI)*dt/tauI;
        
    % utilization parameter (fraction of resources used by each spike)
    if 1 <= i && i < 500
        uE=UE*ones(Ne,1);
    else
        uE = uE + (-uE + UE + tf*UE*rE - tf*UE*rE.*uE)*dt/tf;          
    end
        
    % get decoded angle from network activity  
    ang=decode(rE,theta');

    % Output variables  
    centerBang(i)=ang;
    RE(:,i)=rE;
    u(:,i)=uE;
end



%%%%% Plot excitatory cells firing rate
if vce(j)==5    
    figure
    mesh(tt,theta,RE)
    view(0,90)
    xlabel('time (s)')
    ylabel('excitatory cells preferred direction (\circ)')
    xlim([0 9000])
    ylim([-2.3 3.8])
    set(gca,'xtick',[0 1000 2000 3000 4000 5000 6000 7000 8000 9000])
    set(gca,'xticklabel',{'0','1','2','3','4','5','6','7','8','9'})
    set(gca,'ytick',[-3.14 -1.57 0 1.57 3.14])
    set(gca,'yticklabel',{'-180','-90', '0','90','180'})
    colorbar
    title('Excitatory cells firing rate for young case (vce= 5 Hz)')
elseif vce(j)==7 
    figure
    mesh(tt,theta,RE)
    view(0,90)
    xlabel('time (s)')
    ylabel('excitatory cells preferred direction (\circ)')
    xlim([0 9000])
    ylim([-2.3 3.8])
    set(gca,'xtick',[0 1000 2000 3000 4000 5000 6000 7000 8000 9000])
    set(gca,'xticklabel',{'0','1','2','3','4','5','6','7','8','9'})
    set(gca,'ytick',[-3.14 -1.57 0 1.57 3.14])
    set(gca,'yticklabel',{'-180','-90', '0','90','180'})
    colorbar
    title('Excitatory cells firing rate for middle-aged case (vce= 7 Hz)')
else
    figure
    mesh(tt,theta,RE)
    view(0,90)
    xlabel('time (s)')
    ylabel('excitatory cells preferred direction (\circ)')
    xlim([0 9000])
    ylim([-2.3 3.8])
    set(gca,'xtick',[0 1000 2000 3000 4000 5000 6000 7000 8000 9000])
    set(gca,'xticklabel',{'0','1','2','3','4','5','6','7','8','9'})
    set(gca,'ytick',[-3.14 -1.57 0 1.57 3.14])
    set(gca,'yticklabel',{'-180','-90', '0','90','180'})
    colorbar
    title('Excitatory cells firing rate for aged case (vce= 9 Hz)')
end
end











