function [z_dot] = multipleDOF_Car_Suspension_Function(t,z,y,time,m,c,k)


%t isn't actually an array. The function gets passed one individual value
 %of t at a time. check it out with
t;
length(t);

% y = 0;
% if t >= d/v && t <= (d+L)/v   
%    y = h*sin( (pi*v/L) * (t - d/v));
% end

y = interp1(time,y,t); %values of y at whatever t must be interpolated,
%since there might not be a given y value for that value of t

z_dot(1) = z(3);
z_dot(2) = z(4);
z_dot(3) = (1/m(1)) * ( c(1)*z(4) + k(1)*z(2) - c(1)*z(3) - k(1)*z(1) );
z_dot(4) = (1/m(2)) * ( k(2)*y + c(1)*z(3) + k(1)*z(1) - c(1)*z(4) - (k(1)+k(2))*z(2) );

z_dot = z_dot';
end