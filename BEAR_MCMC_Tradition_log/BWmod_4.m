function  [Ctssm,Mwashoff,Mres] = Bwmod_3(c3,c4,c2,c1,Q)

%modele du Genie Rural a 4 parametres Journalier
% c3: coeffiicient
% c4: Exponent
% c2: build up rate constant (1/day)
% c1: Maximum build-up possible(kg)

% Mwashoff(k): eroded pollutant mass at k (kg)
% Mres(k): available pollutant mass at k (kg)
% Mbuildup: the mass of build up at te beginning of the rainfall event(kg)

% Q: input-steamflow (m3/S)
% Ctssm: Modelled concentration of Ts (mg/l)
% rain_event[a,b,c]: a-event number; b-bth day in event a; c-1 for wet,0
% for dry.
% build-up:exponential function; wash-off: exponential function


    nDATA = length(Q);

    %% initial simulation
    k=1; Mres(1)=c1*0.9;
            temp1 = Q(k); %% effective surface flow
            Mwashoff(k) = c3*temp1^c4*Mres(k);  %%wash-off 
            
            Ctssm(k) = Mwashoff(k)/Q(k)/86.4;                        
            Mres(k+1)= max(0,Mres(k)-Mwashoff(k));
        
    %% continuous simulation   
	for k = 2:nDATA
              
            t = -log(1-Mres(k)/c1)/c2 + 1;
            Mbuildup(k)= c1*(1-exp(-1*c2*t));%%build-up
            Mres(k)= Mbuildup(k);  
            
            temp1 = Q(k);
            Mwashoff(k) = c3*temp1^c4*Mres(k);
            Ctssm(k) = Mwashoff(k)/Q(k)/86.4;
            
    
%             V(k) = Q(k) *60*60*24 ;
                                   
            Mres(k+1)=max(0,Mres(k)-Mwashoff(k));
    end
    
end