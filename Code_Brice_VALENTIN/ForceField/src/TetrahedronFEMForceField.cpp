#include <ForceField_kernels.h>

using namespace sofa::defaulttype;

void TetrahedronFEMForceField3f_initialize( TVecVec3 & positions0, TVecTetra & tetrahedra, TReal youngModulus, TReal poissonRatio, std::vector<Element> & elems)
{
//initialisation des matrices de rigidite, locales a chaque element

    //copie des indices de chaque tetrahedre dans le tableau "vIx" de "elems"
    for (unsigned int i=0; i<tetrahedra.size(); ++i)
    {
        for (int j=0; j<4; ++j)
        {
            elems[i].vIx[j] = tetrahedra[i][j];
        }
    }

    //Conversion des parametres du modele FE avec les coefficients de LAME
    TReal lambda, mu;
    lambda = (youngModulus*poissonRatio)/((1+poissonRatio)*(1-2*poissonRatio));
    mu = youngModulus/(2*(1+poissonRatio));

    //Creation de la matrice de rigidite interne De
    TMat6x6 De;
    for (int i=0; i<3; ++i)
    {
        De[i][i] = lambda+2*mu;
    }
    for (int i=3; i<6; ++i)
    {
        De[i][i] = mu;
    }
    De[0][1] = lambda;
    De[0][2] = lambda;
    De[1][0] = lambda;
    De[1][2] = lambda;
    De[2][0] = lambda;
    De[2][1] = lambda;

    //Declaration des matrices Ve, Ve_inv, Be
    TMat4x4 Ve, Ve_inv;
    TMat6x12 Be;

    //Boucle sur tous les elements pour pouvoir calculer Ke : Ke = |Ve|*Be'*De*Be/6
    for (unsigned int j=0; j<elems.size(); ++j)
    {
    //Construction de la matrice Be par etapes **********************************************
    //1.Construction de la matrice Ve
        //Boucle sur les lignes de Ve, soit sur les 4 points du tetra
        for (int i=0; i<4; ++i)
        {
            Ve[i][0] = 1;
            Ve[i][1] = positions0[elems[j].vIx[i]][0]; //position initiale x
            Ve[i][2] = positions0[elems[j].vIx[i]][1]; //position initiale y
            Ve[i][3] = positions0[elems[j].vIx[i]][2]; //position initiale z
        }

    //2.Construction de la matrice Ve_inv, inverse de Ve
        TReal detVe = invertMatrix(Ve_inv, Ve); //recuperation du determinant de Ve "detVe" lors de l'inversion de Ve
        if (detVe<0)
        {
            //warning pour prevenir l'utilisateur que les tetrahedres sont mal orientes
            std::cout<<"Les tetrahedres sont mal orientes, il faut changer la valeur de flip_tetra"<<endl;
        }

    //3.Construction de la matrice Be a la main
        Be[0][0] = Ve_inv[1][0];
        Be[0][3] = Ve_inv[1][1];
        Be[0][6] = Ve_inv[1][2];
        Be[0][9] = Ve_inv[1][3];

        Be[1][1] = Ve_inv[2][0];
        Be[1][4] = Ve_inv[2][1];
        Be[1][7] = Ve_inv[2][2];
        Be[1][10] = Ve_inv[2][3];

        Be[2][2] = Ve_inv[3][0];
        Be[2][5] = Ve_inv[3][1];
        Be[2][8] = Ve_inv[3][2];
        Be[2][11] = Ve_inv[3][3];

        Be[3][0] = Ve_inv[2][0];
        Be[3][1] = Ve_inv[1][0];
        Be[3][3] = Ve_inv[2][1];
        Be[3][4] = Ve_inv[1][1];
        Be[3][6] = Ve_inv[2][2];
        Be[3][7] = Ve_inv[1][2];
        Be[3][9] = Ve_inv[2][3];
        Be[3][10] = Ve_inv[1][3];

        Be[4][0] = Ve_inv[3][0];
        Be[4][2] = Ve_inv[1][0];
        Be[4][3] = Ve_inv[3][1];
        Be[4][5] = Ve_inv[1][1];
        Be[4][6] = Ve_inv[3][2];
        Be[4][8] = Ve_inv[1][2];
        Be[4][9] = Ve_inv[3][3];
        Be[4][11] = Ve_inv[1][3];

        Be[5][1] = Ve_inv[3][0];
        Be[5][2] = Ve_inv[2][0];
        Be[5][4] = Ve_inv[3][1];
        Be[5][5] = Ve_inv[2][1];
        Be[5][7] = Ve_inv[3][2];
        Be[5][8] = Ve_inv[2][2];
        Be[5][10] = Ve_inv[3][3];
        Be[5][11] = Ve_inv[2][3];
    //Fin de la construction de la matrice Be **********************************************

    //Remplissage de la matrice Ke pour chaque element : Ke = |Ve|*Be'*De*Be/6
        elems[j].Ke = detVe*Be.transposed()*De*Be/6;
    }
}


void TetrahedronFEMForceField3f_addForce(TVecVec3 & positions0, std::vector<Element> & elems, TVecVec3 & f, const TVecVec3 & positions)
{
// accumulation des forces internes a l'objet : f = f - K*u
// avec u = deplacement depuis la positions initiale
// et K = la matrice de rigidite (Stiffness)

// 2 methodes sont presentes ici :
// 1ere methode : on alloue la matrice K
// 2e methode : on n'alloue plus la matrice K car elle est creuse, donc on calcule directement f=f-(G*Ke*G')*u

// pour utiliser une methode : decommenter la methode souhaitee et commenter la methode non utilise

/*
// 1ere methode **********************************************************************************************************************

    //Allocation dynamique de la memoire pour la matrice K de taille (3*(nombre de points),3*(nombre de points))
    TReal **K = new TReal* [3*positions0.size()];
    for (int i=0; i<3*positions0.size(); ++i)
    {
        K[i] = new TReal[3*positions0.size()];
    }

    //Initialisation de chaque élément de K a 0
    for (int i=0; i<3*positions0.size(); ++i)
    {
        for (int j=0; j<3*positions0.size(); ++j)
        {
            K[i][j] = 0;
        }
    }

    //Boucle sur tous les elements
    for (int a=0; a<elems.size(); ++a)
    {
        //Definitions de variables pour simplifier la visibilite du code
        Element & e = elems[a]; //element courant
        int & i = e.vIx[0]; //indice i du tetrahedre defini par les indices [i,j,k,l]
        int & j = e.vIx[1]; //indice j du tetrahedre defini par les indices [i,j,k,l]
        int & k = e.vIx[2]; //indice k du tetrahedre defini par les indices [i,j,k,l]
        int & l = e.vIx[3]; //indice l du tetrahedre defini par les indices [i,j,k,l]

        //double boucles pour parcourir les matrices A,B,...P de Ke qui sont de taille 3*3 avec Ke=[A B C D; E F G H; I J K L; M N O P]
        for (int x=0; x<3; ++x)
        {
            for (int y=0; y<3; ++y)
            {
                //Premiere ligne de Ke
                K[3*i+x][3*i+y] += e.Ke[3*0+x][3*0+y]; //a l'endroit (3*i,3*i) de K
                K[3*i+x][3*j+y] += e.Ke[3*0+x][3*1+y]; //a l'endroit (3*i,3*j) de K
                K[3*i+x][3*k+y] += e.Ke[3*0+x][3*2+y]; //a l'endroit (3*i,3*k) de K
                K[3*i+x][3*l+y] += e.Ke[3*0+x][3*3+y]; //a l'endroit (3*i,3*l) de K

                //Deuxieme ligne de Ke
                K[3*j+x][3*i+y] += e.Ke[3*1+x][3*0+y]; //a l'endroit (3*j,3*i) de K
                K[3*j+x][3*j+y] += e.Ke[3*1+x][3*1+y]; //a l'endroit (3*j,3*j) de K
                K[3*j+x][3*k+y] += e.Ke[3*1+x][3*2+y]; //a l'endroit (3*j,3*k) de K
                K[3*j+x][3*l+y] += e.Ke[3*1+x][3*3+y]; //a l'endroit (3*j,3*l) de K

                //Troisieme ligne de Ke
                K[3*k+x][3*i+y] += e.Ke[3*2+x][3*0+y]; //a l'endroit (3*k,3*i) de K
                K[3*k+x][3*j+y] += e.Ke[3*2+x][3*1+y]; //a l'endroit (3*k,3*j) de K
                K[3*k+x][3*k+y] += e.Ke[3*2+x][3*2+y]; //a l'endroit (3*k,3*k) de K
                K[3*k+x][3*l+y] += e.Ke[3*2+x][3*3+y]; //a l'endroit (3*k,3*l) de K

                //Quatrieme ligne de Ke
                K[3*l+x][3*i+y] += e.Ke[3*3+x][3*0+y]; //a l'endroit (3*l,3*i) de K
                K[3*l+x][3*j+y] += e.Ke[3*3+x][3*1+y]; //a l'endroit (3*l,3*j) de K
                K[3*l+x][3*k+y] += e.Ke[3*3+x][3*2+y]; //a l'endroit (3*l,3*k) de K
                K[3*l+x][3*l+y] += e.Ke[3*3+x][3*3+y]; //a l'endroit (3*l,3*l) de K
            }
        }
    }

    //Construction de u de type TVecVec3 comme positions0 et positions
    TVecVec3 u;
    u.resize(positions0.size());

    //Calcul de u : difference des positions courantes avec les positions initiales
    for (int a=0; a<positions0.size(); ++a)
    {
        u[a] = positions[a]-positions0[a];
    }

    //Calcul de f : forces internes
    for (int a=0; a<positions0.size(); ++a)
    {
        for (int b=0; b<3; ++b)
        {
            for (int c=0; c<positions0.size(); ++c)
            {
                for (int d=0; d<3; ++d)
                {
                    f[a][b] -= K[3*a+b][3*c+d]*u[c][d];
                }
            }
        }
    }

    //Destruction de la matrice K
    for (int a=0; a<3*positions0.size(); ++a)
        {
        delete [] K[a];
        }
    delete [] K;

// FIN 1ere methode **********************************************************************************************************************
*/

// 2e methode ***************************************************************************************************************************

    //Construction de u de type TVecVec3 comme positions0 et positions
    TVecVec3 u;
    u.resize(positions0.size());

    //Calcul de u : difference des positions courantes avec les positions initiales
    for (unsigned int a=0; a<positions0.size(); ++a)
    {
        u[a] = positions[a]-positions0[a];
    }

    //Boucle sur tous les elements
    for (unsigned int a=0; a<elems.size(); ++a)
    {
        //Definitions de variables pour simplifier la visibilite du code
        Element & e = elems[a]; //element courant
        int & i = e.vIx[0]; //indice i du tetrahedre defini par les indices [i,j,k,l]
        int & j = e.vIx[1]; //indice j du tetrahedre defini par les indices [i,j,k,l]
        int & k = e.vIx[2]; //indice k du tetrahedre defini par les indices [i,j,k,l]
        int & l = e.vIx[3]; //indice l du tetrahedre defini par les indices [i,j,k,l]

        //double boucles pour parcourir les matrices A,B,...P de Ke qui sont de taille 3*3 avec Ke=[A B C D; E F G H; I J K L; M N O P]
        for (int x=0; x<3; ++x)
        {
            for (int y=0; y<3; ++y)
            {
                //Premiere ligne de Ke se retrouve a la ligne i de f
                f[i][x] -= e.Ke[3*0+x][3*0+y]*u[i][y]; //premiere colonne de Ke avec ligne i de u
                f[i][x] -= e.Ke[3*0+x][3*1+y]*u[j][y]; //deuxieme colonne de Ke avec ligne j de u
                f[i][x] -= e.Ke[3*0+x][3*2+y]*u[k][y]; //troisieme colonne de Ke avec ligne k de u
                f[i][x] -= e.Ke[3*0+x][3*3+y]*u[l][y]; //quatrieme colonne de Ke avec ligne l de u

                //Deuxieme ligne de Ke se retrouve a la ligne j de f
                f[j][x] -= e.Ke[3*1+x][3*0+y]*u[i][y]; //premiere colonne de Ke avec ligne i de u
                f[j][x] -= e.Ke[3*1+x][3*1+y]*u[j][y]; //deuxieme colonne de Ke avec ligne j de u
                f[j][x] -= e.Ke[3*1+x][3*2+y]*u[k][y]; //troisieme colonne de Ke avec ligne k de u
                f[j][x] -= e.Ke[3*1+x][3*3+y]*u[l][y]; //quatrieme colonne de Ke avec ligne l de u

                //Troisieme ligne de Ke se retrouve a la ligne k de f
                f[k][x] -= e.Ke[3*2+x][3*0+y]*u[i][y]; //premiere colonne de Ke avec ligne i de u
                f[k][x] -= e.Ke[3*2+x][3*1+y]*u[j][y]; //deuxieme colonne de Ke avec ligne j de u
                f[k][x] -= e.Ke[3*2+x][3*2+y]*u[k][y]; //troisieme colonne de Ke avec ligne k de u
                f[k][x] -= e.Ke[3*2+x][3*3+y]*u[l][y]; //quatrieme colonne de Ke avec ligne l de u

                //Quatrieme ligne de Ke se retrouve a la ligne l de f
                f[l][x] -= e.Ke[3*3+x][3*0+y]*u[i][y]; //premiere colonne de Ke avec ligne i de u
                f[l][x] -= e.Ke[3*3+x][3*1+y]*u[j][y]; //deuxieme colonne de Ke avec ligne j de u
                f[l][x] -= e.Ke[3*3+x][3*2+y]*u[k][y]; //troisieme colonne de Ke avec ligne k de u
                f[l][x] -= e.Ke[3*3+x][3*3+y]*u[l][y]; //quatrieme colonne de Ke avec ligne l de u
            }
        }
    }

// FIN 2e methode **********************************************************************************************************************

}

void TetrahedronFEMForceField3f_addDForce( double factor, std::vector<Element> & elems, TVecVec3 & df,const TVecVec3 & dx)
{
    //df = df + factor*(dH(u)/du)*dx

    //Boucle sur tous les elements
    for (unsigned int a=0; a<elems.size(); ++a)
    {
        //Definitions de variables pour simplifier la visibilite du code
        Element & e = elems[a]; //element courant
        int & i = e.vIx[0]; //indice i du tetrahedre defini par les indices [i,j,k,l]
        int & j = e.vIx[1]; //indice j du tetrahedre defini par les indices [i,j,k,l]
        int & k = e.vIx[2]; //indice k du tetrahedre defini par les indices [i,j,k,l]
        int & l = e.vIx[3]; //indice l du tetrahedre defini par les indices [i,j,k,l]

        //double boucles pour parcourir les matrices A,B,...P de Ke qui sont de taille 3*3 avec Ke=[A B C D; E F G H; I J K L; M N O P]
        for (int x=0; x<3; ++x)
        {
            for (int y=0; y<3; ++y)
            {
                //Premiere ligne de Ke se retrouve a la ligne i de df
                df[i][x] += factor*e.Ke[3*0+x][3*0+y]*dx[i][y]; //premiere colonne de Ke avec ligne i de dx
                df[i][x] += factor*e.Ke[3*0+x][3*1+y]*dx[j][y]; //deuxieme colonne de Ke avec ligne j de dx
                df[i][x] += factor*e.Ke[3*0+x][3*2+y]*dx[k][y]; //troisieme colonne de Ke avec ligne k de dx
                df[i][x] += factor*e.Ke[3*0+x][3*3+y]*dx[l][y]; //quatrieme colonne de Ke avec ligne l de dx

                //Deuxieme ligne de Ke se retrouve a la ligne j de df
                df[j][x] += factor*e.Ke[3*1+x][3*0+y]*dx[i][y]; //premiere colonne de Ke avec ligne i de dx
                df[j][x] += factor*e.Ke[3*1+x][3*1+y]*dx[j][y]; //deuxieme colonne de Ke avec ligne j de dx
                df[j][x] += factor*e.Ke[3*1+x][3*2+y]*dx[k][y]; //troisieme colonne de Ke avec ligne k de dx
                df[j][x] += factor*e.Ke[3*1+x][3*3+y]*dx[l][y]; //quatrieme colonne de Ke avec ligne l de dx

                //Troisieme ligne de Ke se retrouve a la ligne k de df
                df[k][x] += factor*e.Ke[3*2+x][3*0+y]*dx[i][y]; //premiere colonne de Ke avec ligne i de dx
                df[k][x] += factor*e.Ke[3*2+x][3*1+y]*dx[j][y]; //deuxieme colonne de Ke avec ligne j de dx
                df[k][x] += factor*e.Ke[3*2+x][3*2+y]*dx[k][y]; //troisieme colonne de Ke avec ligne k de dx
                df[k][x] += factor*e.Ke[3*2+x][3*3+y]*dx[l][y]; //quatrieme colonne de Ke avec ligne l de dx

                //Quatrieme ligne de Ke se retrouve a la ligne l de df
                df[l][x] += factor*e.Ke[3*3+x][3*0+y]*dx[i][y]; //premiere colonne de Ke avec ligne i de dx
                df[l][x] += factor*e.Ke[3*3+x][3*1+y]*dx[j][y]; //deuxieme colonne de Ke avec ligne j de dx
                df[l][x] += factor*e.Ke[3*3+x][3*2+y]*dx[k][y]; //troisieme colonne de Ke avec ligne k de dx
                df[l][x] += factor*e.Ke[3*3+x][3*3+y]*dx[l][y]; //quatrieme colonne de Ke avec ligne l de dx
            }
        }
    }
}
