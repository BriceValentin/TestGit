#include <State_kernels.h>
#include <ForceField_kernels.h>
#include <Mass_kernels.h>
#include <stdio.h>

int simulation_cg_iter;

void computeForce(const TVec3 & gravity, FEMMesh* mesh, TVecVec3 &b, double rayleighMass, double rayleighStiffness, double h)
{
    //redimensionnement de b, de la meme taille que le vecteur positions
    b.resize(mesh->positions.size());

    //on efface les forces du pas de temps precedent
    MechanicalObject3f_vClear(b);

    //accumulation des forces externes p
        //accumulation de la gravite
        UniformMass3f_addForce(gravity, mesh->massDensity, b);
        //accumulation des forces d'interaction avec la souris si l'index est compris entre 0 et la taille de b
        if ( (mesh->externalForce.index >= 0) && ((unsigned)mesh->externalForce.index < b.size() )  )
        {
            b[mesh->externalForce.index] += mesh->externalForce.value;
        }

    //accumulation des forces internes du modele element fini : b=b-K*u
    TetrahedronFEMForceField3f_addForce(mesh->positions0, mesh->femElem, b, mesh->positions);

    //accumulation de l'amortissement en appliquant le facteur rayleighMass : b=b-alpha*M*v
    MechanicalObject3f_vPEqBF(b, mesh->velocity, -(mesh->massDensity)*rayleighMass);

    //accumulation du facteur rayleighStiffness : b=b-(dH(u)/du)*(h-beta)*v
    TetrahedronFEMForceField3f_addDForce(h-rayleighStiffness, mesh->femElem, b, mesh->velocity);

    //application de l'operateur de projection pour appliquer les conditions limites, il annule toutes les forces appliquees sur les points fixes
    FixedConstraint3f_projectResponseIndexed(mesh->fixedParticles, b);
}

void  mulMatrixVector(const  SimulationParameters* params,  FEMMesh*  mesh,  std::vector<TVec3>&  result,  const  std::vector<TVec3>&  input, double h)
{
    //multiplication de la matrice A et un vecteur x
    //avec A = (1+h*alpha)*M + (h^2+h*beta)*(dH(u)/du)
    //et x = input

    //Decomposition en 2 etapes :
    //result = [(1+h*alpha)*M]*input
    //result = result + [(h^2+h*beta)*(dH(u)/du)]*input

    //Etape 1 : result = [(1+h*alpha)*M]*input
    MechanicalObject3f_vEqBF(result, input, (1+h*(mesh->massDensity))*(params->rayleighMass));

    //Etape 2 : result = result + [(h^2+h*beta)*(dH(u)/du)]*input
    TetrahedronFEMForceField3f_addDForce(h*h+h*(params->rayleighStiffness), mesh->femElem, result, input);

    //application de l'operateur de projection pour appliquer les conditions limites, il annule toutes les forces appliquees sur les points fixes
    FixedConstraint3f_projectResponseIndexed(mesh->fixedParticles, result);
}

void linearSolver_ConjugateGradient(const  SimulationParameters*  params,  FEMMesh*  mesh, TVecVec3  &  x,  TVecVec3  &  b, double h)
{
    //Algorithme du Gradient Conjugué (CG)

    //initialisation de la variable itérative i
    int i = 0;

    //Construction de r de type TVecVec3 comme b et Ax
    //puisque r = b - A*x
    TVecVec3 r;
    r.resize(b.size());

    //Calcul de A*x avec la fonction mulMatrixVector définie au-dessus
    TVecVec3 Ax;
    Ax.resize(b.size());
    mulMatrixVector(params, mesh, Ax, x, h);

    //Calcul de r : r = b - A*x
    MechanicalObject3f_vOp( r, b, Ax, -1);

    //initialisation de d
    TVecVec3 d;
    d.resize(r.size());
    d = r;

    //delta_new = produit scalaire r'r
    TReal delta_new;
    delta_new = MechanicalObject3f_vDot(r, r);

    //initialisation de delta_0
    TReal delta_0;
    delta_0 = delta_new;

    //declaration des variables intermediaires utilisees dans la boucle while
    TVecVec3 q;
    q.resize(r.size());
    TReal alpha, delta_old, beta;

    //boucle while de l'algorithme du gradient conjugue
    while (i<params->maxIter && delta_new>((params->tolerance)*(params->tolerance)*delta_0))
    {
        //calcul de q : q = A*d
        mulMatrixVector(params, mesh, q, d, h);

        //calcul de alpha : alpha = delta_new/(d'q)
        alpha = delta_new/MechanicalObject3f_vDot(d, q);

        //mise a jour de x : x = x + alpha*d
        MechanicalObject3f_vPEqBF( x, d, alpha);

        //si i est divisible par 50 (utilisation d'un modulo)
        if ((i%50) == 0)
        {
            //mise a jour de Ax
            mulMatrixVector(params, mesh, Ax, x, h);

            //mise a jour de r : r = b - A*x
            MechanicalObject3f_vOp( r, b, Ax, -1);
        }
        else
        {
            //mise a jour de r : r = r - alpha*q
            MechanicalObject3f_vPEqBF( r, q, -alpha);
        }

        //variable tampon
        delta_old = delta_new;

        //mise a jour de delta_new = produit scalaire r'r
        delta_new = MechanicalObject3f_vDot(r, r);

        //calcul de beta : beta = delta_new/delta_old
        beta = delta_new/delta_old;

        //mise a jour de d : d = r + beta*d
        MechanicalObject3f_vOp( d, r, d, beta);

        //incrementation de i
        i = i + 1;
    }

    //nombre d’iterations necessaire pour atteindre convergence dans le CG
    simulation_cg_iter = i;
}


void timeIntegrator_EulerImplicit(const SimulationParameters* params, FEMMesh* mesh, double h)
{
    //derivee de H(u) par rapport a u => K car H(u)=K*u

    //appel de computeForce pour calculer le terme b de A*x=b
    computeForce(params->gravity, mesh, mesh->b, params->rayleighMass, params->rayleighStiffness, h);

    //redimensionnement de a, de la meme taille que le vecteur b
    mesh->a.resize(mesh->b.size());

    //resolution du systeme A*x=b par la methode du gradient conjugue definie dans la fonction linearSolver_ConjugateGradient
    linearSolver_ConjugateGradient(params, mesh, mesh->a, mesh->b, h);

    //integration des vitesses
    MechanicalObject3f_vOp(mesh->velocity, mesh->velocity, mesh->a, h);

    //integration des positions
    MechanicalObject3f_vOp(mesh->positions, mesh->positions, mesh->velocity, h);
}
