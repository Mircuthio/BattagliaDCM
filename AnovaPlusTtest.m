function [s_all,hfig,hfig2]=AnovaPlusTtest(Xt,groups,infvec,info_string,normalization)
% Xt=globTotLen{end};
% fprintf('\nTest2 Gaussian 1 versus Gaussian 1\n');
% [h,p,ci,stats]=ttest(Xt1);
try
    infvec;
catch
    infvec={'Low Cooperation: No Recognition and No Signaling .','Agent B Recognition, No Signaling','Agent B Recognition & Signaling','Agent A,B Recognition, No Signaling','Agent A,B Recognition & Signaling'};
    conds=unique(groups);
    infvec=cell(length(conds),1);
    for i=1:length(infvec)
        infvec{i}=['Condition ' num2str(groups(i))];
    end 
end
try
    s=info_string;
catch
    s='';
end
try
    nn=normalization;
catch
    nn=1;
end
conditions=unique(groups);
Xt1=Xt / nn;
cend=conditions(end);
s_all=cell(0,0);
for ind=1:length(conditions)-1
    cind=conditions(ind);
    [h,p,ci,t_stats]=ttest2(Xt1(groups==cind),Xt1(groups==cend));
    sn=sprintf('%s tTest %s versus %s\n',s,infvec{end},infvec{ind});
    fprintf('\n%s\n',sn);
    
    st=printTTest(h,p,ci,t_stats);
    s_all=[s_all {sn} st];
end
% [p,table,a_stats]=anova1(Xt1);
[~,table,a_stats]=anova1(Xt1,groups);
hfig=gcf;
hfig2=figure;
[c,m,hfig]=multcompare(a_stats);

% etasquared=table{2,2}/table{4,2};
% effectsize=sqrt(etasquared/(1-etasquared));

result=sprintf('%s Anova F[%g,%g]=%g, p=%g',s,table{2,3},table{3,3},table{2,5},table{2,6});
s_all =[result s_all];
fprintf('\n%s\n',result);
