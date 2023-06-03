Kp=0.48*0.6
Ti = 
Td =0
PIDC = pid(Kp,1/Ti,Td)
T = feedback(PIDC*Gs,1)
t = 0:Tp:1000;
step(T,t)