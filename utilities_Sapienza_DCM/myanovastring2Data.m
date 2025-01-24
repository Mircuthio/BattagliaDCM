    %% Anova Test
function string = myanovastring2Data(Y)
    
s=[];
% [p,table,a_stats]=anova1(Xt1);
[~,table,a_stats]=anova1(Y);
    % figure;
    % [c,m,hfig]=multcompare(a_stats);
    
%     fprintf('\n%s F[%g,%g]=%g, p=%g\n',s,table{2,3},table{3,3},table{2,5},table{2,6});
    alpa = [0.01,0.05];
    [h,p,ci,stats] = ttest2(Y(:,1),Y(:,2),'Alpha',alpa(1));
    [h1,p1,ci1,stats1] = ttest2(Y(:,1),Y(:,2),'Alpha',alpa(2));
    % disp(h), disp(stats.tstat), disp(p), disp(stats.df)
    % disp(h1), disp(stats1.tstat), disp(p1), disp(stats1.df)

    string=sprintf('\n%s F[%g,%g]=%g, p=%g\n',s,table{2,3},table{3,3},table{2,5},table{2,6});
        fprintf('\n%s t(%g)=%g, p=%g, alpha=%g\n',s,stats1.df,stats.tstat,p,alpa(1));
    fprintf('\n%s t(%g)=%g, p=%g, alpha=%g\n',s,stats1.df,stats.tstat,p,alpa(2));
    fprintf('\n%s F[%g,%g]=%g, p=%g\n',s,table{2,3},table{3,3},table{2,5},table{2,6});