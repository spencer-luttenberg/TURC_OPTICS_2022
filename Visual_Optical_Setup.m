clc;
cnt = 0;
num_hit_rays = 0;

x_pts = zeros(1, 7);
y_pts = zeros(1, 7);
i = 1;
% for y = 4:0.025:5
%     num_hit_rays = 0;
%     cnt = 0;

hold on
for x = -37.5:0.1:37.5
    x_pts(1) = 0;
    y_pts(1) = x;
    ray = input_lens(x, 77.52, 11.45, 1.78, 1.50091);
    x_pts(2) = ray(3,1);
    y_pts(2) = ray(1,1);
    ray = traverse(ray, 121.6);
    x_pts(3) = ray(3,1);
    y_pts(3) = ray(1,1);
    ray = concave(ray, 1.74377, 9.42, 12,4.12, 2.20, 0);
    x_pts(4) = ray(3,1);
    y_pts(4) = ray(1,1);
    ray = traverse(ray, 4.825);
    x_pts(5) = ray(3,1);
    y_pts(5) = ray(1,1);    
    ray = convex(ray, 1.74377, 7.85, 10, 1.45, 3.25, 4.825);
    x_pts(6) = ray(3,1);
    y_pts(6) = ray(1,1);
    ray = traverse(ray, 150 - ray(3,1));
    x_pts(7) = ray(3,1);
    y_pts(7) = ray(1,1);
    if (ray(1,1)) >= -0.75 && (ray(1,1)) <= 0.75
        if ray(2,1) >= -12 && ray(2,1) <= 12
        num_hit_rays = num_hit_rays + 1;
        end
    end
    cnt = cnt + 1;
    plot(x_pts, y_pts,'Color',[1;1;0]);
end

%    array_of_distances(i) = y;
%    array_of_percent(i) = num_hit_rays/cnt*100;
%   
%    i = i + 1;
%   
% end
disp((num_hit_rays/cnt)*100);

%plot(array_of_distances, array_of_percent);








% ray = input_lens(-10, 77.52, 11.45, 1.78, 1.50091)
% r_new = vpa(concave(ray, 1.74377, 7.07, 9, 3.63, 2.25, 130))






function ray_new = traverse(ray, distance)
ray_new = [ray(1,1) + distance*tand(ray(2,1)); ray(2,1); ray(3,1) + distance];

end


function ray = input_lens(input_h, radius, tc, te, ng)
theta_r = asind(input_h/radius);
a = radius*cosd(theta_r) - radius + tc - te;
tx = a + te;          
theta_1 = theta_r;
theta_2 = asind((1/ng)*sind(theta_1));
theta_in = theta_r - theta_2;
ray(1,1) = input_h - tand(theta_in)*tx;
ray(2,1) = -(asind(ng*sind(theta_in)));
ray(3,1) = tc;
end


function ray_output = concave(ray, ng, r, lens_h, te, tc, d_from_to)
ax = te - tc;
dc = d_from_to - r + ax;
d1 = r - ax;
%theta__r = asind(h/r);
ray_new = ray;
%ray_new = traverse(ray, d_from_to);

z = tand(ray_new(2,1));
h_in = ray_new(1,1);
a1 = (-z*h_in + sqrt(r^2 + z^2*r^2 - h_in^2))/(1+z^2);
%see if it intersects the circle of concave lens
if isreal(a1) && a1 > d1
    x = a1;
    y = x*tand(ray_new(2,1)) + ray_new(1,1);
    theta_r = asind(y/r);
    theta_i = ray(2,1);
    tx = r - x + tc;
    
    if theta_i < 0
       if theta_r > 0
           theta_1 = theta_r + abs(theta_i);
           theta_2 = asind(sind(theta_1)*(1/ng)) - theta_r;
           ray_output(2,1) = (asind(ng*sind(-theta_2)));
%            vpa(ray_output(2,1))
           ray_output(1,1) = y - tx*tand(theta_2);
       end
       if theta_r < 0
           if abs(theta_r) > abs(theta_i)
               theta_1 = abs(theta_r) - abs(theta_i);
           else     
               theta_1 = abs(theta_i) - abs(theta_r);
           end
               theta_2 = asind(sin(theta_1)*(1/ng)) + abs(theta_r);
               ray_output(1,1) = y - tx*tand(theta_2);
               ray_output(2,1) = -(asind(ng*sind(theta_2)));     
       end
    else
        if theta_r > 0
            if theta_i > theta_r
                theta_1 = theta_i - theta_r;
            end
            if theta_r > theta_i
                theta_1 = theta_r - theta_i;
            end
            theta_2 = asind(sind(theta_1)*(1/ng)) + theta_r;
            ray_output(1,1) = y + tx*tand(theta_2);
            ray_output(2,1) = (asind(ng*sind(theta_2))); 
            
        else
            theta_1 = abs(theta_r) + theta_i;
            theta_2 = asind(sind(theta_1)*(1/ng)) - abs(theta_r);
            ray_output(1,1) = y + tx*tand(theta_2);
            ray_output(2,1) = (asind(ng*sind(theta_2))); 
        end  
    end 
    ray_output(3,1) = dc + r + tc + ray(3,1);
    
else
    ray_output = traverse(ray_new, r+tc);
end
end

function r_out = convex(ray, ng, r, lens_h, te, tc, d_from_to)
%traverse the distance
%ray = traverse(ray, d_from_to);
%create a line for intersection -> (0,0) is the center of the lens
 b = ray(1,1);
 m = tand(ray(2,1));
 ax = tc - te;
 p_intersect = ((-b*m)+r-sqrt(r^2-2*b*m*r-b^2))/(m^2+1);
 if isreal(p_intersect) && p_intersect <= ax
     %x and y of positions with reference to front center of lens
     x = p_intersect;
     y = m*x + b;
     theta_r = asind(y/r);
     theta_i = ray(2,1);
     %tx = distance from enter to end of glass
     tx = tc - x;
     if theta_r > 0
         if theta_i >= 0
             theta_1 = theta_i + theta_r;
             theta_2 = asind((1/ng)*sind(theta_1)) - theta_r;
             %add to h 
             r_out(1,1) = ray(1,1) + tand(theta_2)*tx;
             %we can use sign conventions and keep the same angle sign
             r_out(2,1) = asind(ng*sind(theta_2));
         else
             if abs(theta_i) > theta_r
                 %get 'positive' angle from axis inside glass
                 theta_1 = abs(theta_i) - theta_r;
                 theta_2 = asind((1/ng)*sind(theta_1)) + theta_r;
             end
             if abs(theta_i) < theta_r
                 %get 'positive' angle from axis inside glass
                 theta_1 = theta_r - abs(theta_i);
                 theta_2 = theta_r - asind((1/ng)*sind(theta_1));
             end
             if abs(theta_i) == theta_r
                 theta_2 = theta_r;
             end
             %output new height
             r_out(1,1) = ray(1,1) - tand(theta_2)*tx;
             %set new angle
             r_out(2,1) = -asind(ng*sind(theta_2));        
         end
         
         r_out(3,1) = ray(3,1) + tc;
     end
     if theta_r < 0 
         if theta_i >= 0
             if theta_i > abs(theta_r)
                 theta_1 = theta_i - abs(theta_r);
                 theta_2 = asind((1/ng)*sind(theta_1)) + abs(theta_r);
             end
             if theta_i < abs(theta_r)
                 theta_1 = abs(theta_r) - theta_i;
                 theta_2 = abs(theta_r) - asind((1/ng)*sind(theta_1));
             end        
             if theta_i == theta_r
                 theta_2 = theta_i;
             end
             r_out(1,1) = ray(1,1) + tand(theta_2)*tx;
             r_out(2,1) = asind(ng*sind(theta_2));
         end
         
         if theta_i < 0  
             theta_1 = abs(theta_i) + abs(theta_r);
             theta_2 = asind((1/ng)*sind(theta_1)) - abs(theta_r);
             %'negative' angle
             r_out(1,1) = ray(1,1) - tx*tand(theta_2);
             r_out(2,1) = -asind(ng*sind(theta_2));
             
         end
         r_out(3,1) = ray(3,1) + tc;
     end
     if theta_r == 0
         if theta_i == 0 
             r_out = traverse(ray, ray(3,1) + tc);
         else
             theta_1 = theta_i;
             theta_2 = asind((1/ng)*sind(theta_1));
             r_out(1,1) = ray(1,1) + tx*tand(theta_2);
             r_out(2,1) = asind(ng*sind(theta_2));
         end
     end     
     else
     r_out = traverse(ray, ray(3,1) + tc);
 end
end


