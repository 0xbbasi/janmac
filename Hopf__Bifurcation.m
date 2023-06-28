
clc;
clear all;
close all;function Hopf_Bifurcation()

    % parameters
    a = 1;
    b = 3;
    c = 1;
    d = 5;
    r = 0.006;
    s = 4;
    I = 0.8;
    xr = (c/d)^(1/2); % calculate xr
   
    % create meshgrid for phase plane plot
    x1 = linspace(-4, 4, 20);
    x2 = linspace(-4, 4, 20);
    [X1, X2] = meshgrid(x1, x2);
       
    % evaluate the system's derivatives at each point on the grid
    x1dot = X2 - a*X1.^3 + b*X1.^2 + I - X2;
    x2dot = c - d*X1.^2 - X2;
    x3dot = r*(s*(X1 - xr) - X3);
    
    % plot the phase plane
    quiver(X1, X2, x1dot, x2dot);
    xlabel('x');
    ylabel('y');
    axis([-4 4 -4 4]);
    hold on;
    
    % find the Jacobian matrix
    syms x1 x2 x3
    Jr = jacobian([x2 - a*x1^3 + b*x1^2 + I - x3; c - d*x1^2 - x2;
                   r*(s*(x1 - xr) - x3)], [x1, x2, x3]);
    
    % find the equilibrium points
    eqPoints = vpasolve([x2 - a*x1^3 + b*x1^2 + I - x3 == 0, c - d*x1^2 - x2 == 0,
                        r*(s*(x1 - xr) - x3) == 0], [x1, x2, x3]);
    
    % calculate the eigenvalues at each equilibrium point
    for i = 1:length(eqPoints.x1)
        JrNum = double(subs(Jr, [x1, x2, x3], [eqPoints.x1(i), eqPoints.x2(i), eqPoints.x3(i)]));
        eigs_Jr = eig(JrNum);
        real_eigs(:, i) = real(eigs_Jr); % store the real parts of the eigenvalues
    end
    
    % check for Hopf bifurcation
    for i = 1:length(eqPoints.x1)
        if eqPoints.x3(i) == 0 && all(real_eigs(:, i) < 0)
            plot(eqPoints.x1(i), eqPoints.x2(i), 'ro');
            disp('Hopf bifurcation detected');
            disp(['Equilibrium point: (' num2str(eqPoints.x1(i)) ', ' num2str(eqPoints.x2(i)) ', ' num2str(eqPoints.x3(i)) ')']);
            break;
        end
    end
    
end
