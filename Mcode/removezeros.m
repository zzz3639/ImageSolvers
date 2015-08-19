function csrres=removezeros(csmres,fold)
%   Modified in 2015.06.20 by ZHANG Haowen
    Int=csmres(:,3);
    Int=Int/sum(sum(Int));
    [u,v]=sort(Int);
    csmres=csmres(v(end:-1:1),:);
    Int=Int(v(end:-1:1),:);
    C=cumsum(Int);
    if size(csmres,1)==0
        csrres=zeros(0,3);
    else
        csrres=csmres(1:sum(C<fold)+1,:);
end