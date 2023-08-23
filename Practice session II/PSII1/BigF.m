function F = BigF(Ia,Iss,I1,I2,b0,b1,b2);


a1 = -b1*Iss;
a0 = a1 + (b1-b0)*I1;
a2 = a1 + (b1-b2)*I2;


if Ia<=I1
    F=a0 + b0*Ia;
elseif Ia<=I2
    F=a1 + b1*Ia;
else
    F=a2 + b2*Ia;
end


