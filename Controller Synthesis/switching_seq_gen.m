function seq = switching_seq_gen(start, L, possible_switch)
%This function generates all admissible sequences of length L+1 based on
%possible starting modes and admissible transitions.

if ~ischar(start)
    start = num2str(start);
end
if ~ischar(possible_switch)
    possible_switch = num2str(possible_switch);
end

if L == 1
    seq = start;
else
    
    j = 1;
    steps = L;
    n = size(start,1);
    prev_seq = start;
    for i = 1:n
        pos = find(possible_switch(:,1)==prev_seq(i,end));
        for k = 1:length(pos)
            next_seq(j,:) = [prev_seq(i,:), possible_switch(pos(k),end)];
            j = j+1;
        end
    end
    
    seq = switching_seq_gen(next_seq,steps-1,possible_switch);
end