function is_cw = is_clockwise(polygon)
% Checks if polygon is clockwise, by evaluating (x2 ? x1)(y2 + y1) at each
% point (positive if clockwise, negative if counterclockwise), and summing
% up these values across the entire polygon.
% e.g. (from Roberto Bonvallet on StackOverflow)
% point[0] = (5,0)   edge[0]: (6-5)(4+0) =   4
% point[1] = (6,4)   edge[1]: (4-6)(5+4) = -18
% point[2] = (4,5)   edge[2]: (1-4)(5+5) = -30
% point[3] = (1,5)   edge[3]: (1-1)(0+5) =   0
% point[4] = (1,0)   edge[4]: (5-1)(0+0) =   0
%                                          ---
%                                          -44  counter-clockwise (if
%                                          positive y goes upwards)
    poly1 = polygon;                % Polygon
    poly2 = circshift(poly1,-1,1);  % "Next corners" of polygon
    is_cw = sum( (poly2(:,1) - poly1(:,1)) .* (poly2(:,2) + poly1(:,2)) ) < 0;  % Sign is flipped since positive y goes downwards in MATLAB
end