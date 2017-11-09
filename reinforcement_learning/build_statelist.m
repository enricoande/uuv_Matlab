function statelist = build_statelist(s1,s2)
% build_statelist.m      e.anderlini@ucl.ac.uk     23/10/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function returns a statelist (it could be either a list of discrete
% states or the position of the centres of radial basis functions).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

l1 = length(s1);
l2 = length(s2);

statelist = zeros(l1*l2,2);
k = 0;
for i=1:l1
    for j=1:l2
        k = k + 1;
        statelist(k,:) = [s1(i);s2(j)];
    end
end

end