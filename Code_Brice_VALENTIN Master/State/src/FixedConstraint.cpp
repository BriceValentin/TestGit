#include <State_kernels.h>

void FixedConstraint3f_projectResponseIndexed( const TVecInt & indices, TVecVec3 & dx )
{
//assigne le vecteur "null" sur tous les points fixes (stockes dans le vecteur "indices")

    for (unsigned int i=0; i<indices.size(); ++i)
    {
        //utilisation de la fonction "clear" qui met a 0 les cases de dx[indices[i]] de type Vec3
        dx[indices[i]].clear();
    }
}
