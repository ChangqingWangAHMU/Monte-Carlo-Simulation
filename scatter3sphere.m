function scatter3sphere(r,size_sphere)
figure;axis([-40 40 -40 40 -40 40]);view([1 1 1]);
markerwidth=5.425;
scatter3(r(:,1),r(:,2),r(:,3),(size_sphere'/2*markerwidth).^2,'filled');
axis([-40 40 -40 40 -40 40]);grid on;