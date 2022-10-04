function EmbededLegend(x, y, xPos, label, labelColor)
%Description: Finds and plots the specified label tangent to a specified 
%plot (x,y) at a specified x-location    

%Inputs
%x      - vector of x used in plot
%y      - vecotr of y used in plot
%xPos   - scalar x position for label
%label  - string containing label information

% Obtain local slope
d = (y(xPos+1)-y(xPos))/(x(xPos+1)-x(xPos));
X = diff(get(gca, 'xlim'));
Y = diff(get(gca, 'ylim'));
p = pbaspect;
a = atan(d*p(2)*X/p(1)/Y)*180/pi;

% Display the text
text(x(xPos), y(xPos), label ,'Color',labelColor, 'BackgroundColor', 'w', 'rotation', a);

end