function new_states = QAR_step(old_states, graph,baseDrive)

%Take an existing State (QAR_step) and step model forward


%create new empty state filled with Q as default
new_states = zeros(size(old_states));

%The model can run on weighted or directed graphs

% -2 and −1 for the refractory state
% 1 corresponds to the active state


%loop through each neuron

for ii = 1:length(old_states)
    
    %update Activt to Refrac
    if old_states(ii) > 0 
      
        new_states(ii) = -2; %turn into refractory stage
        
   %Update each neighbour
   activeID = rand(size(new_states))< graph(:,ii);
   
   %only those currently in Q state can be infected
   activeID = activeID & old_states == 0;
   
    %then activate any connected neighbours 
    new_states(activeID) = 1;
        
       
    %if in refractory increase one step closer to Q 
    elseif  old_states(ii) < 0 
        new_states(ii) = old_states(ii) +1;
        
        
    %Q can be randomly turned to active with rate controlled by baseDrive
    elseif  old_states(ii) == 0 && rand < baseDrive
        new_states(ii) = 1;
             
    end
        
        
end
