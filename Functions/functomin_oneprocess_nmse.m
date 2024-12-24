function cost=functomin_oneprocess_nmse(rf,rferr,inf,inferr,funcorr,funcorrerr,Fiber_length,var2,interval)


%no automatic integer setting
vel=var2(1);
ceff=var2(3);
a=var2(2);

t=(-(log(1-rf)*(ceff+1)*(ceff+2))/(2*a*vel)).^(1/(ceff+2));
infth=a.*t.^ceff;
valuesinf=nansum((infth-inf).^2)./nansum((infth-mean(inf)).^2); %nmse

valuescor=[];
for i=1:length(funcorr(:,1))
    
max_length1=2*vel*t(i);
if max_length1>Fiber_length-1
    max_length1=Fiber_length-1;
end
x=0:max_length1;
y=funcorr(i,x+Fiber_length);
phi=exp(-2*vel*a*t(i).^(ceff+2)/((ceff+1)*(ceff+2)));
comp=t(i).*(1 - x./max_length1);
xfit_unbiased=1-2*phi+(phi^2).*exp(2*vel*a.*comp.^(ceff+2)/((ceff+1)*(ceff+2)));

valuescor=[valuescor,nansum((xfit_unbiased-y).^2)./nansum((xfit_unbiased-mean(y)).^2)];

end

cost=(nansum(valuesinf)+nansum(valuescor));
                                
end