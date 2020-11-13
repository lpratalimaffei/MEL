REP = importdata('Tabular_data.out');


REP_DATA = REP.data(1:4460,:);

check_REP = REP_DATA(REP_DATA(:,end)==1,:);


pen_count = 0;
save_REP  = [];
obj_start = REP_DATA(1,end);
save_REP(1,:) = [1 obj_start pen_count];

count = 2;


for i = 2:size(REP_DATA,1)
    
    if REP_DATA(i,end) == 1
        pen_count = pen_count +1 ;
    end
    
    if REP_DATA(i,end) < obj_start
        
        obj_start = REP_DATA(i,end);
        save_REP(count,:) = [i obj_start pen_count];
        count = count +1 ;
    end
    
end



hold on; 
%xlim([0 4500])
plot(save_REP(:,1),save_REP(:,2), '-b')

set(gca,'FontSize',20)
set(gca, 'Layer', 'top')

xlabel('Evaluation number','fontweight','bold' )
ylabel('1-CM','fontweight','bold' )

yyaxis right
plot(save_REP(:,1),save_REP(:,3), '--b')
ylabel('Number of Penalties','fontweight','bold' )
ylim([0 2500])



lgd = legend({'Objective Function' 'Penalties'}, 'FontSize', 18, 'Location', 'northeast');
lgd.FontSize = 18;

