#include <State_kernels.h>
#include <ForceField_kernels.h>
#include <Mass_kernels.h>
#include <stdio.h>

#include <Sofa/Vec.h>
#include <Sofa/Mat.h>
#include <DataStructures.h>


void computeForce(const TVec3 & gravity, FEMMesh* mesh, TVecVec3 &f, double rayleighMass)
{
    //question 4.1.1 : M*a = p - (alpha*M)*v - K*u , avec p forces externes

    //redimensionnement de f, de la meme taille que le vecteur positions
    f.resize(mesh->positions.size());

    //on efface les forces du pas de temps precedent
    MechanicalObject3f_vClear(f);

    //accumulation des forces externes p
        //accumulation de la gravite
        UniformMass3f_addForce(gravity, mesh->massDensity, f);
        //accumulation des forces d'interaction avec la souris si l'index est compris entre 0 et la taille de f
        if ( (mesh->externalForce.index >= 0) && ((unsigned)mesh->externalForce.index < f.size() )  )
        {
            f[mesh->externalForce.index] += mesh->externalForce.value;
        }

    //accumulation des forces internes du modele element fini : f=f-K*u
    TetrahedronFEMForceField3f_addForce(mesh->positions0, mesh->femElem, f, mesh->positions);

    //accumulation de l'amortissement en appliquant le facteur rayleighMass : f=f-alpha*M*v
    MechanicalObject3f_vPEqBF(f, mesh->velocity, -(mesh->massDensity)*rayleighMass);

    //application de l'operateur de projection pour appliquer les conditions limites, il annule toutes les forces appliquees sur les points fixes
    FixedConstraint3f_projectResponseIndexed(mesh->fixedParticles, f);
}

void timeIntegrator_EulerExplicit(const SimulationParameters* params, FEMMesh* mesh, double h)
{
    //appel de la fonction computeForce
    computeForce(params->gravity, mesh, mesh->f, params->rayleighMass);

    //redimensionnement de a, de la meme taille que le vecteur f
    mesh->a.resize(mesh->f.size());

    //calcul de a avec M*a=f
    MechanicalObject3f_vEqBF(mesh->a, mesh->f, 1/(mesh->massDensity));

    //integration des vitesses
    MechanicalObject3f_vOp(mesh->velocity, mesh->velocity, mesh->a, h);

    //integration des positions
    MechanicalObject3f_vOp(mesh->positions, mesh->positions, mesh->velocity, h);
}


