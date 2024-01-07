% puntos iniciales

x1 = 30;
y1 = 220;

x2 = 50;
y2 = 100;

x3 = 100;
y3 = -30;

x4 = 260;
y4 = 10;

% creacion de variables

% se crean las matrices para los valores en x e y
% en y simplemente se guardaran los resultados de las ecuaciones
% para x se calculara una ecuacion cubica en cada punto


Y = [y1;y2;y3;y4];
X = [x1^3 x1^2 x1 1; x2^3 x2^2 x2 1; x3^3 x3^2 x3 1; x4^3 x4^2 x4 1];

% 

% a traves de la funcion inv que hace una matriz inversa con los 
% arreglos, se calculan los coeficientes y se guardan en a
A = inv(X) * Y;

% se crean valores 201 valores del x1 al x4 para poder realizar calculos en
% ellos y finalmente hacer 
x = linspace(x1, x4,201);

% 
fcurva = A(1)*x.^3+A(2)*x.^2 + A(3)*x +A(4);
 
xx=linspace(x1,x4,201);
yy=(A(1)*x.^3+A(2)*x.^2+A(3)*x+A(4));
dx=[];
dy=[];
ds=[];

% FOR LOOP
% se calcula la longitud de arco al tomar las derivadas como el cambio 
% de los valores de x e y siguientes con el actual
% despues se meten en la formula de longitud 
% de arco


for i=1:length(x)-1

    dx=[dx,xx(i+1)-xx(i)];
    dy=[dy,yy(i+1)-yy(i)];
    ds=[ds,sqrt(dx(i)^2+dy(i)^2)];
    der1x = dx./ds;
    der1y = dy./ds;

% se consigue ds
% se consigue der1x y der1y
%

end

longituddd=sum(ds);
disp(longituddd)


ddx = [];
ddy = [];
dr = [];

for i=1:length(der1x)-1

    ddx=[ddx,(der1x(i+1)-der1x(i))/ds(i)];
    ddy=[ddy,(der1y(i+1)-der1y(i))/ds(i)];
    dr = [dr,1/sqrt(ddx(i).^2+ddy(i).^2)];

end

% radio = min(dr);
% disp(radio)
c = [];

% Primera coordenada
for j = 2:length(der1y)-1
    if (der1y(j)>0 && der1y(j-1)<0) || (der1y(j)>0 && der1y(j+1)<0)
        c=[c,j];
    end
end
% disp(c)


indice1 = find(dr==min(dr(c(1) - 5: c(1)+5)));
indice2 = find(dr==min(dr(c(2) - 5: c(2)+5)));
R1 = dr(indice1);
R2 = dr(indice2);

t1 = find(x>(x(indice1)-R1-0.6) & x<(x(indice1)-R1+0.6));
t2 = find(x>(x(indice2)-R2-0.6) & x<(x(indice2)-R2+0.6));

longitud1 = sum(ds(1:t1));
disp(longitud1)
longitud2 = sum(ds(1:t2));
disp(longitud2);

% Recta 1

m = der1y(t1)/der1x(t1);
b1 = fcurva(t1) - m*x(t1); 
coordx = [x(t1-35):x(t1+35)];
tangente1 = m*coordx + b1;

perp = -(1/m)*coordx + (fcurva(t1)+x(t1)/m);

parallel = tangente1-20;
xintercept = ((fcurva(t1)+x(t1)/m)-(b1-20))/(m+(1/m));
yintercept = m*xintercept + b1-20;

parallel2 = tangente1-30;
xintercept2 = ((fcurva(t1)+x(t1)/m)-(b1-30))/(m+(1/m));
yintercept = m*xintercept2 + b1-30;

distance = sqrt((coordx-xintercept).^2+(parallel-yintercept).^2);
distance2 = sqrt((coordx-xintercept).^2+(parallel2-yintercept).^2);
long1 = find(distance<41 & distance>39);
long2 = find(distance2<41 & distance2>39);

gradax1 = coordx(long1(1));
graday1 = parallel(long1(1));
gradax2 = coordx(long1(2));
graday2 = parallel(long1(2));

gradax3 = coordx(long2(1));
graday3 = parallel2(long2(1));
gradax4 = coordx(long2(2));
graday4 = parallel2(long2(2));

% Recta 2

m = der1y(t2)/der1x(t2);
b1 = fcurva(t2) - m*x(t2); 
coordx = [x(t2-35):x(t2+35)];
tangente1 = m*coordx + b1;

perp = -(1/m)*coordx + (fcurva(t2)+x(t2)/m);

parallel = tangente1+20;
xintercept = ((fcurva(t2)+x(t2)/m)-(b1+20))/(m+(1/m));
yintercept = m*xintercept + b1+20;

parallel2 = tangente1+30;
xintercept2 = ((fcurva(t2)+x(t2)/m)-(b1+30))/(m+(1/m));
yintercept = m*xintercept2 + b1+30;

distance = sqrt((coordx-xintercept).^2+(parallel-yintercept).^2);
distance2 = sqrt((coordx-xintercept).^2+(parallel2-yintercept).^2);
long1 = find(distance<41 & distance>39);
long2 = find(distance2<42 & distance2>39);

gradax21 = coordx(long1(1));
graday21 = parallel(long1(1));
gradax22 = coordx(long1(2));
graday22 = parallel(long1(2));

gradax23 = coordx(long1(1));
graday23 = parallel2(long1(1));
gradax24 = coordx(long2(2));
graday24 = parallel2(long2(2));

theta = 9.9;
mu = 0.9;
g = 9.8;
R = dr;
masa = 800;

vnper1 = sqrt(mu*g*dr(t1));
vnper2 = sqrt(mu*g*dr(t2));
vper1 = sqrt(dr(t1)*g * (sind(theta)+mu*cosd(theta))/(cosd(theta)+mu*cosd(theta)-mu*sind(theta)));
vper2 = sqrt(R(t2)*g *(sin(theta)+mu*cos(theta))/(cos(theta)+mu*cos(theta)-mu*sin(theta)));

peralte = "Para una pista son peralte escriba 1. Para un pista con peralte escriba 2: ";
opt1 = input(peralte);

[vnper1, vnper2, vper1, vper2];

vel0 = "Ingrese un valor de velocidad (en m/s): ";
opt2 = input(vel0);

parameter = 0;
d = 0;
calor = 0;
distancia = 0;

if opt1 == 1
    if (opt2 < vnper1) && (opt2 < vnper2)
        parameter = 1;
        d = longituddd;
        calor = masa *mu*g*d;
        distancia =  489.7897;      
    else if (opt2 >= vnper1)
        parameter = 2;
        d = opt2^2/(2*mu*g);
        calor = masa*mu*g*d;
        distancia = 250.2758;
    else if (opt2 < vnper1) && (opt2 >= vnper2)
        parameter = 3;
        d = opt2^2/(2*mu*g);
        calor = masa*mu*g*d;
        distancia = 404.0603;
    end
    end
    end

else if opt1 == 2
        if (opt2 < vper1) && (opt2 < vper2)
            parameter = 1;
            d = longituddd;
            calor = masa *mu*g*d;
            distancia =  489.7897;
        else if (opt2 >= vper1)
            parameter = 2;
            d = opt2^2/(2*mu*g);
            calor = masa*mu*g*d;
            distancia = 250.2758;
        else if (opt2 < vper1) && (opt2 >= vper2)
            parameter = 3;
            d = opt2^2/(2*mu*g);
            calor = masa*mu*g*d;
            distancia = 404.0603;

        end
        end
        end
end

end

pause_time = 0.07;
dim = length(x);
y = fcurva;

figure(2),clf
hold on
plot(x,fcurva,'LineWidth',20, 'Color', 'black')
plot(x,fcurva, 'LineWidth',7,'Color','white')

patch('Faces',[1,2,3,4],'Vertices', [gradax3 graday3; gradax4 graday4; gradax2 graday2; gradax1 graday1])% la de 1 y 2 no se si estenb bien

patch('Faces',[1,2,3,4],'Vertices',[gradax23 graday23; gradax24 graday24; gradax22 graday22; gradax21 graday21])%la de 22 y 21 no se si esten bien
xlabel('x(m)')
ylabel('y(m)')

txt1 = 'Sobrepasó el límite de velocidad';
txt2 = 'Calor disipado: ';
txt3 = 'KJ';
txt4 = 'La distancia recorrida fue: ';
txt5 = 'metros';

for step = 1:dim
    if parameter == 1
        point = scatter(x(step),y(step),40,[1 0 0],'filled');
        pause(pause_time)
        delete(point)
     
    else if parameter == 2

            if x(step)>= x(t1)
                text(80,100,txt1)
                %text(170,80,[txt2 num2str(calor/1000, '%.5g') txt3])
                point = scatter(x(step),m*x(step) + b1,40,[1 0 0], 'filled');
                pause(pause_time)

                if sqrt((x(step)-x(t1))^2 + (m*(x(step)-x(t1)))^2) >= ((y(step)-y(t1))^2 + (m*(y(step)-y(t1)))^2) 
                    break
                else
                    
                    delete(point)
                end

            else 
                point = scatter(x(step),y(step),40,[1 0 0],'filled');
                pause(pause_time)
                delete(point)
            end


        else if parameter == 3
             if x(step)>= x(t2)
                 text(80,100,txt1)
                 text(170,80,[txt2 num2str(calor/1000, '%0.5g') txt3])
                 point = scatter(x(step),m*x(step) + b1,40,[1 0 0], 'filled'); 
                 delete(point)
                 if sqrt((x(step)-x(t2))^2 + (m*(x(step)-x(t2)))^2)>= ((y*(step)-y(t2))^2 + (m*(y(step)-y(t2)))^2)
                    break
                 else
%                      continue
                     delete(point)
                 end


             else
                 text(80,100,txt1)
                 %text(170,80,[txt2 num2str(calor/1000, '%0.5g') txt3])
                 
                 point = scatter(x(step),y(step),40,[1 0 0],'filled');
                 pause(pause_time)
                 delete(point)
             end
       end
    end
    end

end

text(80,120,[txt2 num2str(calor/1000, '%0.5g') txt3])
text(80,140,[txt4 num2str(distancia, '%0.5g') txt5])

                 
hold off
