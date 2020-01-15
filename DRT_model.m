
%%%%% Model for the Delayed Response Task (Oculomotor task)
%%%%% FIRING RATE NETWORK MODEL FOR THE BUMP ATTRACTOR (with continuous ring connectivity)

%%%%% The output shows examples of an under-excited network (persistent activity ends before the end of the delay),
%%%%% a network maintaining persistent activity tuned to the stimulus location (TPA-S) until the end of the delay, 
%%%%% and a partially over-excited network (all excitatory cells start firing at the same rate), 
%%%%% for vce = 3, 5, and 9 respectively.

%%%  Modified by Sara Ibañez from:
%%%  Wimmer K, Nykamp DQ, Constantinidis C, Compte A. Nat Neurosci. 17(3):431-439 (2014).


clear

%%%%% PARAMETERS
Ne=640;                 % number of "excitatory cells" in the rate model
Ni=160;                 % number of "inhibitory cells" in the rate model

totalTime=4200;         % total time of the simulation in ms
dt=2;                   % integration step in ms

tauE=20;                % time constant of rate equation for excitatory cells in ms
tauI=10;                % time constant of rate equation for inhibitory cells in ms

kappa=1.5;              % parameter defining concentration of e-to-e connectivity

 
GEEa=38.861;            % strength of excitation to excitatory neurons (AMPA receptors)
GEEn=69.312;            % strength of excitation to excitatory neurons (NMDA receptors)
GIE=31.318;             % strength of inhibition to excitatory neurons (GABAa receptors)
GEIa=195.77;            % strength of excitation to inhibitory neurons (AMPA receptors)
GEIn=215.14;            % strength of excitation to inhibitory neurons (NMDA receptors)
GII=220.94; 

sigE=1;                 % standard deviation of additive noise in rate equation of excitatory cells
sigI=3;                 % standard deviation of additive noise in rate equation of inhibitory cells

stimon = 1000;          % time when external stimulus is applied in ms
stimoff = 1500;         % time when external stimulus ceases in ms
stim =40000;            % strength of external stimulus
delayend=3500;          % time when delay ends in ms, and external input is applied to erase memory



%%%%% PRELIMINARY CALCULATIONS
rE=zeros(Ne,1);
rI=zeros(Ni,1);
nsteps=floor(totalTime/dt);
tt=dt*(1:nsteps);

% E-to-E connectivity
theta = (0:Ne-1)/Ne*2*pi;
v = exp(kappa*cos(theta));
v = v/sum(v);
WE = gallery('circul',v);

% stimulus parameters
theta=theta-pi;
v = exp(kappa*cos(theta));
v = v/sum(v);
stimulus = stim*v';
stimon = floor(stimon/dt);
stimoff = floor(stimoff/dt);
delayend = floor(delayend/dt);

% external bias current
I0E=80;  %pA	  
I0I=15;  %pA

% input-output function for excitatory & inhibitory cells (used previously in: Brunel, Cereb Cortex 13:1151 (2003))
fE = inline('vce/(98^2)*x.*x.*(x>0).*(x<98)+vce*sqrt((4/98)*x-3).*(x>=98)','x','vce');
fI = inline('50/(20^2)*x.*x.*(x>0).*(x<20)+50*sqrt((4/20)*x-3).*(x>=20)','x');

% population vector decoder given the rates r for neurons with selectivity th
decode = inline('atan2(sum(r.*sin(th)),sum(r.*cos(th)))','r','th');



%%%% SIMULATION LOOP
RE=zeros(Ne,nsteps);
centerBang=zeros(nsteps,1);

vce=[3,5,9];

for j=1:3

for i=1:nsteps
  % additive noise for each population
  noiseE = sigE*randn(Ne,1);
  noiseI = sigI*randn(Ni,1);
  
  % current input to each population
  if i==1
      IEa=0;
      IEn=0;
      IEg=0;
      IIa=0;
      IIn=0;
      IIg=0;
  else
      IEa=IEa+(-IEa+GEEa*WE*rE)*dt/2;
      IEn=IEn+(-IEn+GEEn*WE*rE)*dt/100;
      IEg=IEg+(-IEg-GIE*mean(rI)*ones(Ne,1))*dt/10;
      IIa=IIa+(-IIa+GEIa*mean(rE)*ones(Ni,1))*dt/2;
      IIn=IIn+(-IIn+GEIn*mean(rE)*ones(Ni,1))*dt/100;
      IIg=IIg+(-IIg-GII*mean(rI)*ones(Ni,1))*dt/10;
  end 
  IE=IEa+IEn+IEg+I0E*ones(Ne,1);
  II=IIa+IIn+IIg+I0I*ones(Ni,1);
 
  % external task-dependent inputs
  if i>stimon && i<stimoff
      IE=IE+stimulus;         % cue stimulus before delay
  end
  if i>delayend && i<delayend+(stimoff-stimon)      
      IE=IE-stim;             % erasing global input after delay
  end

  % integration with time-step dt: Newton method
  rE = rE + (fE(IE,vce(j)) - rE + noiseE)*dt/tauE;
  rI = rI + (fI(II) - rI + noiseI)*dt/tauI;
  
  % get decoded angle from network activity  
  ang=decode(rE,theta');
  
  % Excitatory cells firing rate  
  RE(:,i)=rE;
  centerBang(i)=ang;
end



%%%%% Plot excitatory cells firing rate
if vce(j)==3 
    figure
    mesh(tt,theta,RE)
    view(0,90)
    xlabel('time (s)')
    ylabel('excitatory cells preferred direction (\circ)')
    xlim([0 4200])
    ylim([-3.14 3.14])
    set(gca,'xtick',[0 1000 2000 3000 4000]);
    set(gca,'xticklabel',{'0','1','2','3','4'});
    set(gca,'ytick',[-3.14 -1.57 0 1.57 3.14]);
    set(gca,'yticklabel',{'-180','-90', '0','90','180'});
    colorbar
    title('Excitatory cells firing rate. Under-excited network (vce= 3 Hz)')
elseif vce(j)==5
    figure
    mesh(tt,theta,RE)
    view(0,90)
    xlabel('time (s)')
    ylabel('excitatory cells preferred direction (\circ)')
    xlim([0 4200])
    ylim([-3.14 3.14])
    set(gca,'xtick',[0 1000 2000 3000 4000]);
    set(gca,'xticklabel',{'0','1','2','3','4'});
    set(gca,'ytick',[-3.14 -1.57 0 1.57 3.14]);
    set(gca,'yticklabel',{'-180','-90', '0','90','180'});
    colorbar
    title('Excitatory cells firing rate. Network maintaining TPA-S (vce= 5 Hz)')
elseif vce(j)==9
    figure
    mesh(tt,theta,RE)
    view(0,90)
    xlabel('time (s)')
    ylabel('excitatory cells preferred direction (\circ)')
    xlim([0 4200])
    ylim([-3.14 3.14])
    set(gca,'xtick',[0 1000 2000 3000 4000]);
    set(gca,'xticklabel',{'0','1','2','3','4'})
    set(gca,'ytick',[-3.14 -1.57 0 1.57 3.14]);
    set(gca,'yticklabel',{'-180','-90', '0','90','180'});
    colorbar
    title('Excitatory cells firing rate. Partially over-excited network (vce= 9 Hz)')
end
end
    
    
    
    
    