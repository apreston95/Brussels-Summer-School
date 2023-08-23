function F = PSIIbx(Ia,cti,I1,I2,Iss,b0,b1,b2);


 Ii=BigF(Ia,Iss,I1,I2,b0,b1,b2) +cti;

F = Ii - Ia;