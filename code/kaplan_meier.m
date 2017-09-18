function D = kaplan_meier(filename)


D = dlmread(filename,'\t',1,0);

t = D(:,1)/30.0; % Survival times
censored = 1-D(:,2); % Event
age = abs(D(:,3))/365.0; %Age
group = D(:,4); % Cluster

% Possible line colors
c = {'b','r','k','g','m','y','b','g'};
c = {[0 0 1], [0 0 0], [0 0.4 0], [1 0 0.0], [1 1 0], [0.6 0 0.4], [0 1 0.0], [0 1 1], [1 0 1], [1 0.2 0.2]};


% Possible legend titles
for i=1:max(group)
    clust{i} = ['Clust ' num2str(i) ' (' num2str(length(find(group==i))) ',' num2str(length(find(group==i & censored==0))) ')'];
    clus{i} = ['Clust ' num2str(i)]; 
    ages{i} = age(find(group==i & age ~= 0));
    fprintf('Average age for cluster %d is: %f and median is: %f \n',i,mean(ages{i}),median(ages{i}));
end;

% Survival curves
figure(1)
for i=min(group):1:max(group)
    % Calculate and plot the empirical cdf and confidence bounds
    % ecdf is the Empirical (Kaplan-Meier) cumulative distribution function
    [f,x,flo,fup] = ecdf(t(find(group==i)),'censoring',censored(find(group==i))); 
    stairs(x,1-f,'Color',c{i},'LineWidth',2);
    hold on
end;
xp = xlabel('Time (Months)');
yp = ylabel('Survival probability');
l = legend(clust{1:max(group)});
set(xp,'FontSize',14);
set(yp,'FontSize',14);
set(l,'FontSize',14);
set(gca,'FontSize',14);

% Box plot average age
figure(2)
boxplot(age(find(age ~= 0)),group(find(age ~= 0)),'labels',clus);
h = findobj(gca,'tag','Box');
med = findobj(gca,'tag','Median');
uw = findobj(gca,'tag','Upper Whisker');
uav = findobj(gca,'tag','Upper Adjacent Value');
lw = findobj(gca,'tag','Lower Whisker');
lav = findobj(gca,'tag','Lower Adjacent Value');

for jj=1:length(h)
    set(h(length(h)-jj+1),'LineWidth',2,'Color',c{jj});
    set(med(length(med)-jj+1),'Color',c{jj},'LineWidth',2);
    set(uw(length(uw)-jj+1),'Color',c{jj},'LineWidth',2);
    set(uav(length(uav)-jj+1),'Color',c{jj},'LineWidth',2);
    set(lw(length(lw)-jj+1),'Color',c{jj},'LineWidth',2);
    set(lav(length(lav)-jj+1),'Color',c{jj},'LineWidth',2);
end;
set(findobj(gca,'Type','text'),'FontSize',14);
yp = ylabel('Age (years)');
set(yp,'FontSize',14);
set(gca,'FontSize',14);
