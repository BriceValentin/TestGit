#include <State_kernels.h>
#include <string.h>

void MechanicalObject3f_vPEq1( TVecVec3 & res, int index, const TVec3 val )
{
//accumule la valeur "val" sur l'element "index" du vecteur "res"

    res[index] += val;
}

void MechanicalObject3f_vClear( TVecVec3 & res )
{
//assigne le vecteur "null" dans toutes les cases du vecteur "res"

    for (unsigned int i=0; i<res.size(); ++i)
    {
        //utilisation de la fonction "clear" qui met à 0 les cases de res[i] de type Vec3
        res[i].clear();
    }
}

void MechanicalObject3f_vEqBF( TVecVec3 & res, const TVecVec3 & b, float f )
{
//calcul de res = b*f

    for (unsigned int i=0; i<res.size(); ++i)
    {
        //l'operateur "*" est defini entre un Vec3 et un scalaire (pas besoin d'une double boucle)
        res[i] = b[i]*f;
    }
}

void MechanicalObject3f_vPEqBF( TVecVec3 & res, const TVecVec3 & b, float f )
{
//calcul de res = res + b*f

    for (unsigned int i=0; i<res.size(); ++i)
    {
        //l'operateur "*" est defini entre un Vec3 et un scalaire
        //et l'operateur "+" est defini entre deux Vec3
        //(pas besoin d'une double boucle)
        res[i] = res[i] + b[i]*f;
    }
}

void MechanicalObject3f_vOp( TVecVec3 & res, const TVecVec3 & a, const TVecVec3 & b, float f )
{
//calcul de res = a + b*f

    for (unsigned int i=0; i<res.size(); ++i)
    {
        //l'operateur "*" est defini entre un Vec3 et un scalaire
        //et l'operateur "+" est defini entre deux Vec3
        //(pas besoin d'une double boucle)
        res[i] = a[i] + b[i]*f;
    }
}

TReal MechanicalObject3f_vDot( const TVecVec3 & a, const TVecVec3 & b)
{
//calcul du produit scalaire de a et b

    TReal res=0; //initialisation de res a 0
    for (unsigned int i=0; i<a.size(); ++i)
    {
        //l'operateur "*" est defini entre deux Vec3
        res += a[i]*b[i]; // accumulation dans res
    }
    return res;

}
