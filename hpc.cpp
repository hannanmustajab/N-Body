/*
    N-Body problem solution using openMP.
    Based on code provided by: Guido Giuntoli
    link: https://www.linkedin.com/pulse/2-optimizing-c-n-body-problem-openmp-guido-giuntoli/
*/

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <omp.h>

class Particle
{
public:
    double pos[3], vel[3];
};

class Problem
{
public:
    Problem(double mass, double dt, unsigned numParticles) : mMass(mass),
                                                             mInverseMass(1.0 / mass),
                                                             mDt(dt),
                                                             mNumParticles(numParticles),
                                                             mParticles(new Particle[numParticles]) {}

    ~Problem()
    {
        delete[] mParticles;
    }

    void integrate();

private:
    const double mG = 6.6743e-11;
    const double mMass;
    const double mInverseMass;
    const double mDt;
    const unsigned mNumParticles;
    Particle *const mParticles;
};

void Problem::integrate()
{
    const double Const = mG * mMass * mMass;

    for (int pi = 0; pi < mNumParticles; pi++)
    {

        double force[3] = {};

        // Calculate total force in body pi

        for (int pj = 0; pj < mNumParticles; pj++)
        {

            if (pj != pi)
            {

                const double dij[3] = {
                    mParticles[pj].pos[0] - mParticles[pi].pos[0],
                    mParticles[pj].pos[1] - mParticles[pi].pos[1],
                    mParticles[pj].pos[2] - mParticles[pi].pos[2]};

                const double dist2 = dij[0] * dij[0] +
                                     dij[1] * dij[1] +
                                     dij[2] * dij[2];

                const double ConstDist2 = Const / dist2;
                const double idist = 1 / sqrt(dist2);

                // F = C * m * m / ||x2 - x1||^2 * (x2 - x1) / ||x2 - x1||

                force[0] += ConstDist2 * dij[0] * idist;
                force[1] += ConstDist2 * dij[1] * idist;
                force[2] += ConstDist2 * dij[2] * idist;
            }
        }

        // dv / dt = a = F / m

        mParticles[pi].vel[0] += force[0] * mInverseMass * mDt;
        mParticles[pi].vel[1] += force[1] * mInverseMass * mDt;
        mParticles[pi].vel[2] += force[2] * mInverseMass * mDt;
    }

    // Update pos, should be done after forces have being computed

    for (int pi = 0; pi < mNumParticles; pi++)
    {

        // dx / dt = v
        mParticles[pi].pos[0] += mParticles[pi].vel[0] * mDt;
        mParticles[pi].pos[1] += mParticles[pi].vel[1] * mDt;
        mParticles[pi].pos[2] += mParticles[pi].vel[2] * mDt;
    }
}

int main()
{
    const int nTimeSteps = 100;
    const double Mass = 1e12;
    const double dt = 1e-6;
    const unsigned numParticles = 10000;
    Problem problem(Mass, dt, numParticles);

    double start_time = omp_get_wtime();
    for (int ts = 0; ts < nTimeSteps; ts++)
        problem.integrate();
    double time = omp_get_wtime() - start_time;
    printf("Time: \t %f \n", time);

    return 0;
}
