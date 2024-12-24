function cost=functomin_twoprocesses_nmse(rf,rferr,inf,inferr,funcorr,funcorrerr,Fiber_length,var2,interval)


%no automatic integer setting
vel=var2(1);
ceff=var2(3);
a=var2(2);

ceff2=var2(6);
a2=var2(5);
vel2=var2(4);
teta=var2(7:end);

t=(-(log(1-rf)*(ceff+1)*(ceff+2))/(2*a*vel)).^(1/(ceff+2));
t2=(-(log(1-rf)*(ceff2+1)*(ceff2+2))/(2*a2*vel2)).^(1/(ceff2+2));
infth=teta.*a.*t.^ceff+(1-teta).*a2.*t2.^ceff2;
valuesinf=nansum((infth-inf).^2)./nansum((infth-mean(inf)).^2); %nmse

valuescor=[];
for i=1:length(funcorr(:,1))
    
max_length1=2*vel*t(i);
if max_length1>Fiber_length-1
    max_length1=Fiber_length-1;
end
x1=0:max_length1;
y1=funcorr(i,x1+Fiber_length);
phi1=exp(-2*vel*a*t(i).^(ceff+2)/((ceff+1)*(ceff+2)));
comp1=t(i).*(1 - x1./max_length1);
xfit_unbiased1=1-2*phi1+(phi1^2).*exp(2*vel*a.*comp1.^(ceff+2)/((ceff+1)*(ceff+2)));

max_length2=2*vel2*t2(i);
if max_length2>Fiber_length-1
    max_length2=Fiber_length-1;
end
x2=0:max_length2;
y2=funcorr(i,x2+Fiber_length);
phi2=exp(-2*vel2*a2*t2(i).^(ceff2+2)/((ceff2+1)*(ceff2+2)));
comp2=t2(i).*(1 - x2./max_length2);
xfit_unbiased2=1-2*phi2+(phi2^2).*exp(2*vel2*a2.*comp2.^(ceff2+2)/((ceff2+1)*(ceff2+2)));

if max_length1>max_length2
    if max_length1>Fiber_length-1
    max_length1=Fiber_length-1;
    end
    x=x1;
    y=y1;
    xfit_unbiased2(end:floor(max_length1)+1)=xfit_unbiased2(end);
    xfit_unbiased=teta(i)*xfit_unbiased1+(1-teta(i))*xfit_unbiased2;
else
    if max_length2>Fiber_length-1
    max_length2=Fiber_length-1;
    end
    x=x2;
    y=y2;
    xfit_unbiased1(end:floor(max_length2)+1)=xfit_unbiased1(end);
    xfit_unbiased=teta(i)*xfit_unbiased1+(1-teta(i))*xfit_unbiased2;
end

valuescor=[valuescor,nansum((xfit_unbiased-y).^2)./nansum((xfit_unbiased-mean(y)).^2)];
end

cost=(nansum(valuesinf)+nansum(valuescor));
                        


end