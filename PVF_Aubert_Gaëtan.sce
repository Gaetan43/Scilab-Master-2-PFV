clear 
close
clc

//on appelle toutes les fonctions (sci) contenue sur l'ordi, plutot dans mon dossier du CC que l'on a créé !
getd("C:\Users\gaetan\Desktop\TP scilab Master 2\Contrôle Continue\données CC Gaetan");

// on crée la boucle qui va chercher les fichier à traiter 1 par 1 dans C:\Users\gaetan\Desktop\TP scilab Master 2\Contrôle Continue\données CC\SPRINT_RADAR

Files = findfiles('C:\Users\gaetan\Desktop\TP scilab Master 2\Contrôle Continue\données CC Gaetan\SPRINT_RADAR', "*.txt")
B = size(Files)

//on mets en place la récursivité sur tous les fichier txt du dossier 
for n =1:B(:,1)

//on va cherche toutes nos donnés a traiter, soit les vitesses et temps enregistrés du fichier n
data = csvRead (Files(n,1),"	  ",".","double")

//on va cherche toutes les donées anthropométriques des sujets contenues dans le fichier "masse_taille.xls"
sheets= readxls("C:\Users\gaetan\Desktop\TP scilab Master 2\Contrôle Continue\données CC Gaetan\masse_taille.xls");

S1=sheets(1);// pointe sur la 1er page du fichier xls
masse_taille= S1.value; // récupère les valeurs numériques 
// donc les valeurs importés du fichier xls sont contenues dans la matrice de la variables "masse"
Nom_Prénom= S1.text; // extrait le nom des mecs 

// on ressort nom, masse, taille du sujet que l'on traite ... a mettre en boucle de sorte que sa prennent tout seul pour la suite...!!!!!!!!
NomSujet = Nom_Prénom (n,1);
MasseSujet = masse_taille (n,2);
TailleSujet = masse_taille (n,3);

t = data(2:length(data(:,1)),1); // temps en seconde 
vitessekm_h = data(2:length(data(:,1)),2); // vitesse en km/h 
v =vitessekm_h/3.6;// on met vitesse en m/s

// les valeurs vers lesquelles nous souhaitons que les coeficients se rapprochent.
c0 = [8;0.1;1.3]';

// On lance la fonction de l'optimisation 
[f, coef] = leastsq(list(minimisation,t,v),c0);
// fonction pour minimiser la fonction 1-Exp par rapport au data enregistré vitesse et temps. Ou f =la valeur de la somme des carrés des écarts  et coef = le vecteur contenant la valeur « optimisée » des paramètres du modèle.

VitesseModélisé = yExp(t,coef);// les valeurs de de la fonction qui a été optimisé, soit on obtient les valeurs de vitesse modélisés.
subplot(1,1,1)
plot(t,VitesseModélisé,'b',t,v,'r');// courbe des valeurs de vitesses modélisées avec les coefficients trouvés via leastsq..

// on mets dans délai, Tau et Vmax leur valeur trouvé lors de l'optimisation.
délai = coef (1,2);
Tau = coef(1,3);
Vmax = coef (1,1);
Time = t - délai; 

distance = Vmax*((Time)+Tau*exp(-(Time/Tau)))-Vmax*Tau ; 
//plot(t,distance)

SampleF = 46.875 ;

// on ajoute nos valeur pour tous ceux qui est Force de resistance à l'air 

Pb = 760;
T° = 20;
Cd = 0.9;
Af = (0.2025*TailleSujet^0.725*MasseSujet^0.425)*0.266;
rho = 1.293*Pb/760*273/(273+T°);
K = 0.5*rho*Af*Cd ;
Vent = 0 ;

// calcul accélération réel 
for i = 1:length(v)-1;
    Acc_H (i) = ((v(i+1,1)-v(i,1))/(t(i+1,1)-t(i,1)));
end

// calcul accélération modélisé
Acc_H_Mod = (Vmax/Tau)*exp(-((t-délai)/Tau));

// calcule force horizontale à partir Acc_H_Mod
F_Hor = MasseSujet*Acc_H_Mod+K*(VitesseModélisé-Vent)^2;

// force normalisé au poids du sujet.
F_Hor_Poids = F_Hor/MasseSujet;

// Calcule de la puissance avec P = V_Mod*F
Puissance = VitesseModélisé .* F_Hor_Poids;

//calcul V0 (m/s), F0 (N), Sfv, R² 
[polyFunction] = polyfitT(VitesseModélisé,F_Hor,1);

subplot(2,1,2)
plot(VitesseModélisé,F_Hor,'r-');

F0 = polyFunction(1,1);
Sfv = polyFunction(2,1);
V0 = -F0/polyFunction(2,1);

// calcul de la Force total en N
Ftot = sqrt((F_Hor)^2+(MasseSujet*9.81)^2);

// calcul de RF, et cela lorsque le temps (t)>= 0.3s
Find_t = find(t>=0.3);
Début_RF = Find_t (1,1);
RF = (F_Hor(Début_RF:length(Ftot),:)./Ftot(Début_RF:length(Ftot),:));

// calcul RF sur 10m et moyenne 
Find_10m = find (distance <= 10);
RF_10m = RF(1:length(Find_10m),1)
Mean_RF_10m = mean(RF_10m)

//calcul de la relation puissance-vitesse
Relation_Pui_V = polyfitT(VitesseModélisé,Puissance,2)

subplot(3,1,2)
plot(VitesseModélisé,Puissance,'b')

//calcul Pmax_Théorique

Vopt = -Relation_Pui_V(2,1)/(2*Relation_Pui_V(3,1))
Pmax_Théorique = Relation_Pui_V(3,1)* (Vopt)^2+Relation_Pui_V(2,1)*Vopt+Relation_Pui_V(1,1)

//calcul Puissance_Modélisé

Puissance_Modélisé = Relation_Pui_V(3,1)*(VitesseModélisé)^2+Relation_Pui_V(2,1)*VitesseModélisé + Relation_Pui_V(1,1)

//calcul DRF:
[DRFpente] = polyfitT(VitesseModélisé(Début_RF:length(Ftot)),RF,1);
DRF = DRFpente (2,1)*100;

// temps à 5m 
Find_Time_5m = find(distance<=5);
Time_5m = t(length (Find_Time_5m),1);


// temps à 10m 
Find_Time_10m = find(distance<=10);
Time_10m = t(length (Find_Time_10m),1);

// temps à 20m 
Find_Time_20m = find(distance<=20);
Time_20m = t(length (Find_Time_20m),1);

// temps à 30m 
Find_Time_30m = find(distance<=30);
Time_30m = t(length (Find_Time_30m),1);

// temps à 40m 
Find_Time_40m = find(distance<=40);
Time_40m = t(length (Find_Time_40m),1);

//distance en 2s
Find_Distance_2s = find(t<=2);
Distance_2s = distance(length(Find_Distance_2s),1);

//distance en 4s
Find_Distance_4s = find(t<=4);
Distance_4s = distance(length(Find_Distance_4s),1);

//top speed 
VitesseMax = max(VitesseModélisé)

//Acceleration time constant

Acc_Constant = Tau

// calcule accéleration moyenne, et tous les valeur que l'on veut pour feuille bilan  
acc_M = mean(Acc_H_Mod);
F0_N_kg = F0/MasseSujet ;
Pmax_Hor = max(Puissance) ;
Pmax = Pmax_Hor*MasseSujet ;
RF_Peak = max(RF)


//data ou on aura toute les valeurs pour recréer le fichier txt.
DATA_FINAL(n,:) = [MasseSujet, V0, F0, F0_N_kg, Pmax, Pmax_Hor, Sfv, Mean_RF_10m, RF_Peak, DRF, Time_5m, Time_10m, Time_20m, Time_30m, Time_40m, Distance_2s, Distance_4s, VitesseMax, Acc_Constant]

end

csvWrite(DATA_FINAL,"C:\Users\gaetan\Desktop\DATAFINALCC.txt","  ",".")

