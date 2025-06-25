function heat = Q_to_heat_matrix(Q, m, c)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENGR 130 - Fall 2024
%
% Function Purpose
% Calculate the temperature of a material given its thermal energy
%
% Function Call
% heat = Q_to_heat(Q, m, c)
%
% Input Arguments
% Q: Thermal energy of object (J)
% m: Mass of object (kg)
% c: Specific heat capacity of object (J / (kg * K))
%
% Output Arguments
% heat: The temperature of the object (K)
%
% Assignment Information
%   Assignment:     Project 4
%   Author:         Tucker Bremer, bremer9@purdue.edu
%   Team ID:        001-15
%  	Contributor:    Name, login@purdue [repeat for each]
%   My contributor(s) helped me:	
%     [ ] understand the assignment expectations without
%         telling me how they will approach it.
%     [ ] understand different ways to think about a solution
%         without helping me plan my solution.
%     [ ] think through the meaning of a specific error or
%         bug present in my code without looking at my code.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    heat = Q ./ (m .* c);
end

%% ACADEMIC INTEGRITY STATEMENT
% I have not used source code obtained from any other unauthorized
% source, either modified or unmodified.  Neither have I provided
% access to my code to another. The project I am submitting
% is my own original work.