function history = QAR(initstates,graph,T,baseDrive)

%large matrix to store all states for simulation
history = zeros(T,numel(initstates));

%grab the initial state
states = initstates;

%loop through each time
for tt = 1:T
    
    if ~mod(tt,100)
        disp(100*tt/T)
    end
    
    %step through QAR
    states = QAR_step(states,graph,baseDrive); 
    
    %update history
    history(tt,:) = states;
end

end