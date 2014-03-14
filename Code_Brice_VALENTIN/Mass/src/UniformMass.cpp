#include <Mass_kernels.h>
#include <stdio.h>

void UniformMass3f_addForce( const TVec3 & gavity, TReal massDensity, TVecVec3 & f )
{
//application des forces externes sur le modele element fini
//ici : forces externes = gravite : f = m*g

    for (unsigned int i=0; i<f.size(); ++i)
    {
        //l'operateur "*" est defini entre un Vec3 et un scalaire
        //l'operateur "+=" est defini entre deux Vec3
        f[i] += gavity*massDensity;
    }
}

void UniformMass3f_accFromF( TReal mass, TVecVec3 & a, const TVecVec3 & f )
{
//resolution du systeme lineaire M*a = f
//M est une matrice diagonale dont chaque valeur de la diagonale est egale a "mass"

    for (unsigned int i=0; i<f.size(); ++i)
    {
        //l'operateur "/" est defini entre un Vec3 et un scalaire
        a[i] = f[i]/mass;
    }
}

void UniformMass3f_addMDx( TReal mass, TVecVec3 & res, const TVecVec3 & dx )
{
//calcul du produit de la matrice de masse (diagonale avec elements egaux a "mass") et du vecteur "dx"
    for (unsigned int i=0; i<res.size(); ++i)
    {
        //l'operateur "*" est defini entre un Vec3 et un scalaire
        res[i] = dx[i]*mass;
    }
}
