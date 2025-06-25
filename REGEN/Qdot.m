function rate = Qdot(K, A, t_hot, t_cold, d)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENGR 130 - Fall 2024
%
% Function Purpose
% Calculate heat transfer rate
%
% Function Call
% rate = Qdot(K, A, t_hot, t_cold, d)
%
% Input Arguments
% K: Conductive heat transfer coeffecient of the material (W / (m^2 * K))
% A: Surface area of contact between two objects (m^2)
% t_hot: Temperature of hot object (K)
% t_cold: Temperature of cold object (K)
%
% Output Arguments
% rate: Heat transfer rate (J/s)
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
rate = (K * A * (t_hot - t_cold)) / d;
end

%% ACADEMIC INTEGRITY STATEMENT
% I have not used source code obtained from any other unauthorized
% source, either modified or unmodified.  Neither have I provided
% access to my code to another. The project I am submitting
% is my own original work.