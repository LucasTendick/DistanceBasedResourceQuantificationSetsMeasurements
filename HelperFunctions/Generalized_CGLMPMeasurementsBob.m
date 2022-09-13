function Mby = Generalized_CGLMPMeasurementsBob(d,m)

% See Generalized_CGLMPMeasurementsAlice for a description

Id = eye(d);

for l = 1:d
basis_vec(:,l) =  Id(:,l);
end


Mby = zeros(d,d,d,m);
Mby_vec = zeros(d,d,m);
for y = 1:m
 for b = 1:d   
     bet(y) = y/m;
     for k = 0:d-1
            Mby_vec(:,b,y) = Mby_vec(:,b,y)+exp(-2*pi*1i/d *k*((b-1)-bet(y)))*basis_vec(:,k+1);
     end
     Mby_vec(:,b,y) = (1/sqrt(d))*Mby_vec(:,b,y);
     Mby(:,:,b,y) = Mby_vec(:,b,y)*Mby_vec(:,b,y)';
 end   
end



end