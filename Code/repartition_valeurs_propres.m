close all;
% Génération des matrices pour chaque type
n = 300; % taille de la matrice

imat =1;
[type1_matrix,~,~] = matgen_csad(imat,n);
imat =2;
[type2_matrix,~,~] = matgen_csad(imat,n);
imat =3;
[type3_matrix,~,~] = matgen_csad(imat,n);
imat =4;
[type4_matrix,~,~] = matgen_csad(imat,n);

% Calcul des valeurs propres pour chaque type de matricebar
type1_eigenvalues = eig(type1_matrix);
type2_eigenvalues = eig(type2_matrix);
type3_eigenvalues = eig(type3_matrix);
type4_eigenvalues = eig(type4_matrix);

% Tracé des spectres des valeurs propres sous forme d'histogrammes
edges = 0:1:200;
figure;
histogram(type1_eigenvalues,edges, Normalization='count');
ytickformat(gca,'percentage');
title('Spectre de la matrice de type 1');

edges = 0:0.1:1;
figure;
histogram(type2_eigenvalues,edges, Normalization='count');
ytickformat(gca,'percentage');
title('Spectre de la matrice de type 2');

edges = 0:0.1:1;
figure;
histogram(type3_eigenvalues,edges, Normalization='count');
ytickformat(gca,'percentage');
title('Spectre de la matrice de type 3');

edges = 0:1/200:1;
figure;
histogram(type4_eigenvalues,edges, Normalization='count');
ytickformat(gca,'percentage');
title('Spectre de la matrice de type 4');

%Tracés des droites de répartition

% Matrice 1
figure;
plot(type1_eigenvalues, 'b.-');
xlabel('Indice');
ylabel('Valeur propre');
title('Distribution des valeurs propres de la matrice 1');
grid on;
saveas(gcf, 'répartition_matrice1.png');

% Matrice 2
figure;
plot(type2_eigenvalues, 'b.');
xlabel('Indice');
ylabel('Valeur propre');
title('Distribution des valeurs propres de la matrice 2');
grid on;
saveas(gcf, 'répartition_matrice2.png');

% Matrice 3
figure;
plot(type3_eigenvalues, 'b.-');
xlabel('Indice');
ylabel('Valeur propre');
title('Distribution des valeurs propres de la matrice 3');
grid on;
saveas(gcf, 'répartition_matrice3.png');

% Matrice 4
figure;
plot(type4_eigenvalues, 'b.-');
xlabel('Indice');
ylabel('Valeur propre');
title('Distribution des valeurs propres de la matrice 4');
grid on;
saveas(gcf, 'répartition_matrice4.png');
