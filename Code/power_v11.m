% méthode de la puissance itérée avec déflation

% Données
% A          : matrice dont on cherche des couples propres
% m          : nombre maximum de valeurs propres que l'on veut calculer
% percentage : pourcentage recherché de la trace
% eps        : seuil pour déterminer si un couple propre a convergé (méthode de la puissance itérée)
% maxit      : nombre maximum d'itérations pour calculer une valeur propre (méthode de la puissance itérée)

% Résultats
% V : matrice des vecteurs propres
% D : matrice diagonale contenant les valeurs propres (ordre décroissant)
% n_ev : nombre de couples propres calculés
% itv : nombre d'itérations pour chaque couple propre
% flag : indicateur sur la terminaison de l'algorithme
%  flag = 0  : on a convergé (on a calculé le pourcentage voulu de la trace)
%  flag = 1  : on a atteint le nombre maximum de valeurs propres sans avoir atteint le pourcentage
%  flag = -3 : on n'a pas convergé en maxit itérations pour calculer une valeur propre
function [ V, D, n_ev, itv, flag ] = power_v11( A, m, percentage, eps, maxit )
    
    n = size(A,1);

    % initialisation des résultats
    W = [];
    V = [];
    itv = [];
    n_ev = 0;

    % trace de A
    tA = trace(A);
 
    % somme des valeurs propres 
    eig_sum = 0.0;

    % indicateur de la convergence (si le pourcentage est atteint)
    convg = 0;
    
    % numéro du couple propre de l'itération courante
    k = 0;

    while (~convg && k < m)
        
        k = k + 1;

        % Méthode de la puissance itérée
        v = randn(n,1);
        z = A*v;
        beta = v'*z;

        % conv = || beta * v - A*v||/|beta| < eps
        norme = norm(beta*v - z, 2)/norm(beta,2);
        nb_it = 1;
        
        while(norme > eps && nb_it < maxit)
          y = A*v;
          v = y / norm(y,2);
          z = A*v;
          beta = v'*z;
          norme = norm(beta*v - z, 2)/norm(beta,2);
          nb_it = nb_it + 1;
        end
        
        % Si le calcul de ce couple propre a échoué, l'algorithme a échoué
        if(nb_it == maxit)
          flag = -3;
          D = 0;
          return;
        end
        
        % On sauvegarde le couple propre
        W(k) = beta;
        V(:,k) = v;
        itv(k) = nb_it;
        eig_sum = eig_sum + beta;

        % On applique la déflation
        A = A - beta * (v*v');

        % On vérifie la condition de convergence en regardant si on atteint le pourcentage
        convg = eig_sum/tA > percentage;
        
    end

    % Si on a atteint le pourcentage
    if (convg)
      n_ev = k;
      flag = 0;
      W = W';
      D = diag(W);
    else
      % Si on a atteint pas le pourcentage
      flag = 1;
      D = 0;
    end
    
end
