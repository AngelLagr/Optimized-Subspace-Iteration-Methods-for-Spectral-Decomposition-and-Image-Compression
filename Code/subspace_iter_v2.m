% version améliorée de la méthode de l'espace invariant (v2)
% avec utilisation de la projection de Raleigh-Ritz 
% avec une accélération bloc

% Données
% A          : matrice dont on cherche des couples propres
% m          : taille maximale de l'espace invariant que l'on va utiliser
% percentage : pourcentage recherché de la trace 
% p          : nombre de produits Y = A^p . V
% eps        : seuil pour déterminer si un vecteur de l'espace invariant a convergé
% maxit      : nombre maximum d'itérations de la méthode

% Résultats
% V : matrice des vecteurs propres
% D : matrice diagonale contenant les valeurs propres (ordre décroissant)
% n_ev : nombre de couples propres calculées
% it : nombre d'itérations de la méthode
% itv  : nombre d'itérations pour chaque couple propre
% flag : indicateur sur la terminaison de l'algorithme
    %  flag = 0  : on converge en ayant atteint le pourcentage de la trace recherché
    %  flag = 1  : on converge en ayant atteint la taille maximale de l'espace
    %  flag = -3 : on n'a pas convergé en maxit itérations

function [ V, D, n_ev, it, itv, flag ] = subspace_iter_v2( A, m, percentage, p, eps, maxit )

    % calcul de la norme de A (pour le critère de convergence d'un vecteur (gamma))
    normA = norm(A,'fro');

    % trace de A
    traceA = trace(A);

    % valeur correspondant au pourcentage spectral à atteindre
    vtrace = percentage*traceA;

    n = size(A,1);
    W = zeros(m,1);
    itv = zeros(m,1);

    % numéro de l'itération courante
    k = 0;
    % somme des valeurs propres à l'itération courante
    eigsum = 0.0;
    % nombre de vecteurs ayant convergé
    nb_c = 0;

    % variable indicant la convergence la convergence
    conv = 0;

    % On génère aléatoirement un ensemble initial de m vecteurs orthogonaux
    Vr = randn(n, m);
    Vr = mgs(Vr);

    % rappel : conv = (eigsum >= trace) | (nb_c == m)
    while (~conv && k < maxit)
        
        k = k+1;
        %% Y <- A^p*V
        for i=1:p
            Vr = A*Vr;
        end
        Y = Vr;
        %% orthogonalisation
        Vr = mgs(Y);
        
        %% Projection de Rayleigh-Ritz
        [Wr, Vr] = rayleigh_ritz_projection(A, Vr);
        
        %% Quels vecteurs ont convergé à cette itération
        analyse_cvg_finie = 0;
        % nombre de vecteurs ayant convergé à cette itération
        nbc_k = 0;
        % nb_c est le dernier vecteur à avoir convergé à l'itération précédente
        i = nb_c + 1;
        
        while(~analyse_cvg_finie)
            % Tous les vecteurs du sous espace ont convergé
            % L'algorithme se termine sans avoir obtenu le pourcentage
            % spectral souhaité
            if(i > m)
                analyse_cvg_finie = 1;
            else
                % On regarde si le vecteur i a convergé
                
                % Calcul de la norme du résidu
                aux = A*Vr(:,i) - Wr(i)*Vr(:,i);
                res = sqrt(aux'*aux);
                
                if(res >= eps*normA)
                    % Si le vecteur i n'a pas convergé les vecteurs suivant
                    % non plus donc l'itération est terminée
                    analyse_cvg_finie = 1;
                else
                    % Si le vetceur i a convergé on continue avec un de
                    % plus
                    nbc_k = nbc_k + 1;
                    % On le stock et on garde aussi sa valeur propre
                    W(i) = Wr(i);
                    
                    itv(i) = k;
                    
                    % On met à jour la somme des valeurs propres
                    eigsum = eigsum + W(i);
                    
                    % si cette somme actualisée permet d'ateindre le
                    % pourcentage spectral souhaité, on s'arrête
                    if(eigsum >= vtrace)
                        analyse_cvg_finie = 1;
                    else
                        % Sinon, on passe au vecteur suivant
                        i = i + 1;
                    end
                end
            end
        end
        
        % On met à jour le nombre de vecteurs ayant convergé
        nb_c = nb_c + nbc_k;
        
        % Condition de convergence ( deux cas )
        conv = (nb_c == m) | (eigsum >= vtrace);
        
    end
    
    if(conv)
        % Mise à jour des résultats
        n_ev = nb_c;
        V = Vr(:, 1:n_ev);
        W = W(1:n_ev);
        D = diag(W);
        it = k;
    else
        % L'algorithme n'a pas convergé
        D = zeros(1,1);
        V = zeros(1,1);
        n_ev = 0;
        it = k;
    end

    % On indique comment le programme a terminé
    if(eigsum >= vtrace)
        flag = 0;
    else if (n_ev == m)
            flag = 1;
        else
            flag = -3;
        end
    end
end
