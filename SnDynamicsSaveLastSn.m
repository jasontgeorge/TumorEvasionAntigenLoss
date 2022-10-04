function [snPlus1,SEnd]=SnDynamicsSaveLastSn(s,q,p,f)
    %Description: This function takes as inputs s,q,p,f and outputs a
    %simulated next antigen, provided a tie occurs.
    
    %If the Recognizer wins, s_{n+1}=Inf <- motivated by fact that recognizer generally does better for higher s
    %If the Recognizer loses, s_{n+1}=-Inf;
    %If the Recognizer ties, s_{n+1}>0
    
    %I. Number of recognized antigens
    r=binornd(s,q);
    if r==0 %Recognizer loses / Evader wins
        SEnd=s;
        snPlus1=-Inf;
    else %Recognizer does not lose
        l=binornd(r,p);
    if l==r %Tie
        SEnd=s-l+f;
        snPlus1=s-l+f; %For now, just use mean f
    elseif l<r %Evader loss / Recognizer win
        SEnd=s;
        %SEnd=s-l;
        %SEnd=s-l+f;
        snPlus1=Inf;
    end
    %If snPlus1==0, then in the next round, r=0 is guaranteed.
    end
end