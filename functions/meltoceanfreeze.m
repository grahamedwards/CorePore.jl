function mof = meltoceanfreeze(Benthic,ocean2freeze,freeze2melt)
% MELTOCEANFREEZE calculates climate conditions for an input of the current 
% benthic value and the threshold variables.
% 
% Inputs: Value of benthic record (benthic) at timestep under investigation, 
% boundary value between ocean and freezing conditions, boundary value 
% between freezing and melting conditions.
% 
% Outputs: -1,0,1 --> -1 = freezing, 0 = ocean, 1 = melting

% verify that ocean2freeze is smaller than freeze2melt
if ocean2freeze >= freeze2melt
    error('ocean2freeze must be less than freeze2melt')
end


if Benthic < ocean2freeze
    mof=0; % ocean conditions

elseif Benthic < freeze2melt && Benthic >= ocean2freeze
    mof=-1; % freezing conditions

elseif Benthic >= freeze2melt
    mof=1; % melting conditions

end

end




