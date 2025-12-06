diameter = 1.5/39.37;
length = 0.5/39.37;
convergelength = 0.625/39.37;
angle = 45;
throatDiameter = 0.213/39.37;
throatLength = 0.25/39.37;
overallLength = 1.376/39.37;

x = linspace(0,overallLength,1000)';
y = ([linspace(diameter,diameter,363),linspace(diameter,throatDiameter,454),linspace(throatDiameter,throatDiameter,183)]/2)';
throatArea = pi*(throatDiameter/2)^2;
areaRatio = (pi*(y.^2))/throatArea;

output = [areaRatio,x,y];