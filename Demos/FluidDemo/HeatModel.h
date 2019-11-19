#pragma once
#ifndef __HeatModel_h__
#define __HeatModel_h__

#include "Demos/Simulation/ParticleData.h"
#include <vector>
#include "Demos/Simulation/NeighborhoodSearchSpatialHashing.h"
#include "PositionBasedDynamics/SPHKernels.h"

#include "TimeStepFluidModel.h"
#include "Demos/Simulation/TimeManager.h"
#include "PositionBasedDynamics/PositionBasedFluids.h"
#include "PositionBasedDynamics/TimeIntegration.h"
#include "Demos/Utils/Timing.h"
#include "HeatModel.h"

namespace PBD
{
	class HeatModel
	{

	protected:
		static Real m_radius;
		static Real m_k;
		static Real m_l;
		static Real m_W_zero;
	public:

		static void HeatFlux(int Particle, ParticleData &pd, const unsigned int nParticles, unsigned int **neighbors, unsigned int *numNeighbors, const Real density[],  const Real boundaryMass[], const Vector3r boundaryPos[], std::vector<bool> &boundaryActive, Real *lastTemp ,Real difussion = 100, bool Border = false , Real BorderTemp = 25)
		{

			Vector3r particlePos = pd.getPosition(Particle);
			Real currentTemp = pd.getTemp(Particle);
			Real meantemp = 0.0;
			Real fluidNeighbors = 0;
			Real mass = pd.getMass(Particle);
			Real particleVol = mass / density[Particle];
			Real t = TimeManager::getCurrent()->getTimeStepSize();

			Real factor = 1;
			
			for (unsigned int j = 0; j < numNeighbors[Particle]; j++)
			{
				const unsigned int neighborIndex = neighbors[Particle][j];			
				Real heaviside = 0;
				
				if (neighborIndex < nParticles)		// Test if fluid particle
				{
					/* for condesation of the gas particles
					if (pd.getPosition(neighborIndex)[1] > 1.5)
						continue;
					//*/

					Vector3r Xij = particlePos - pd.getPosition(neighborIndex);
					Real xij = Xij.norm();
					Real nTemp = pd.getTemp(neighborIndex);
					Real neighborVol = pd.getMass(neighborIndex) / density[neighborIndex];
					/*if (nTemp < currentTemp)
						heaviside = 1;*/
					//meantemp += neighborVol*((2*(nTemp - currentTemp)*dirac + (nTemp - currentTemp)*(1-dirac))/xij)*CubicKernel::dW(Xij);
					//meantemp +=neighborVol *(factor*(nTemp - currentTemp)*heaviside + (nTemp - currentTemp)*(1.0 - heaviside))*CubicKernel::W(Xij);
					meantemp += neighborVol * (nTemp - currentTemp)*CubicKernel::W(Xij);
					fluidNeighbors++;

				}
				else if (Border && boundaryActive[neighborIndex - nParticles]) {
					//std::cout << pd.getTemp(neighborIndex)<<"\n";
					Vector3r Xij = particlePos - boundaryPos[neighborIndex - nParticles];
					Real xij = Xij.norm();				
					/*if (BorderTemp > currentTemp)
						heaviside = 1;*/
					meantemp += particleVol *(BorderTemp - currentTemp)*CubicKernel::W(Xij);
					fluidNeighbors++;

				}
			}

			

				//pd.getTemp(Particle) = currentTemp + t*(meantemp)/difussion;
				Real dt = t*(meantemp)*difussion;

				/*
				if ((pd.getState(Particle) == -2 || pd.getState(Particle) == 1) && pd.getPosition(Particle)[1] > 2.0) {
					float damping = 3.0;//3.5;
					if (pd.getTemp(Particle) >= 15)
						dt -= damping;
				}//*/

				pd.getLatent(Particle)[0] = t*dt*pd.getC(Particle);
				*lastTemp = currentTemp+dt;
	

		}

		static void Cleary(int Particle, ParticleData &pd, const unsigned int nParticles, unsigned int **neighbors, unsigned int *numNeighbors,const Real density[],const Real boundaryMass[], const Real boundaryVolume[], const Vector3r boundaryPos[], std::vector<bool> &boundaryActive, Real *lastTemp ,  bool Border = false, Real BorderTemp = 25)
		{
			/*Conduction Modelling Using SmoothedParticle Hydrodynamics
			Paul W. Cleary⁄ and Joseph J. Monaghany*/

			
			Real currentTemp = pd.getTemp(Particle);

	
			Real vgradient = 0;
			Real mass = pd.getMass(Particle);
			Real particleVol = mass / density[Particle];
			Real C = pd.getC(Particle); //coefficient
			Real phi = particleVol / (mass * C);
			Real ki = pd.getConductivity(Particle); 
			Real sum = 0;
			Vector3r particlePos = pd.getPosition(Particle);
			Real t = TimeManager::getCurrent()->getTimeStepSize();
			//std::cout << particleVol << "\n";
			
			for (unsigned int j = 0; j < numNeighbors[Particle]; j++)
			{
				const unsigned int neighborIndex = neighbors[Particle][j];

				if (neighborIndex < nParticles)		// Test if fluid particle
				{
					/*Real kj = pd.getConductivity(neighborIndex);
					Real deltaT = pd.getTemp(neighborIndex) - currentTemp;
					Real neighborVol = pd.getMass(neighborIndex) / density[neighborIndex];
					Vector3r Xij = particlePos - pd.getPosition(neighborIndex);
					Real xij = Xij.norm();
					Real Kparameter = (4 * ki*kj) / (ki + kj);
					sum += Kparameter*neighborVol*(deltaT / xij)*CubicKernel::W(Xij);*/
					//std::cout <<xij<< " -> "<< CubicKernel::dW(Xij) << "\n";

					Real kj = pd.getConductivity(neighborIndex);
					Real deltaT = pd.getTemp(neighborIndex) - currentTemp;
					Real neighborVol = pd.getMass(neighborIndex) / (density[neighborIndex]);
					Vector3r Xij = particlePos - pd.getPosition(neighborIndex);
					Real xij = Xij.norm();
					Real Kparameter = (4 * ki*kj) / (ki + kj);

					Real val = Kparameter*neighborVol*(deltaT)*CubicKernel::W(Xij);

					sum += val;
				}
				else if (Border && boundaryActive[neighborIndex - nParticles]) {
					Real kj = pd.getConductivity(Particle);
					Real deltaT = BorderTemp - currentTemp;
					Real neighborVol = boundaryVolume[neighborIndex - nParticles];
					Vector3r Xij = particlePos - boundaryPos[neighborIndex - nParticles];
					Real xij = Xij.norm();
					Real Kparameter = (4 * ki*kj) / (ki + kj);
					sum += Kparameter*neighborVol*(deltaT)*CubicKernel::W(Xij);
				}

				
			}
			Real dt = t*sum;// pd.getC(Particle);

			pd.getLatent(Particle)[0] = t*dt*pd.getC(Particle);
			*lastTemp = currentTemp + dt;
		}



		static void CalculateLatentHeat(ParticleData &pd ,Real melting, Real evaporation,Vector2r Thresholds, std::vector<Real> &lasttemp) {


			const unsigned int nParticles = pd.size();
		#pragma omp parallel default(shared)
					{
		#pragma omp for schedule(static)  
						for (int i = 0; i < (int)nParticles; i++)
						{
							switch ((int)pd.getState(i))
							{
								
							/*case 0:
								if (lasttemp[i] < melting && pd.getLatent(i)[1] > 0) {
									lasttemp[i] = melting;
									pd.getLatent(i)[1] += pd.getLatent(i)[0];
								}
								break;
							case 1:
								if (lasttemp[i] >melting && pd.getLatent(i)[1]<Thresholds[0]) {
									lasttemp[i] = melting;
									pd.getLatent(i)[1] += pd.getLatent(i)[0];
								}
								if (lasttemp[i] < evaporation && pd.getLatent(i)[2] > 0) {
									lasttemp[i] = evaporation;
									pd.getLatent(i)[2] += pd.getLatent(i)[0];
								}
								break;
							case 2:
								if (lasttemp[i]  > evaporation && pd.getLatent(i)[2] < Thresholds[1]) {
									lasttemp[i] = evaporation;
									pd.getLatent(i)[2] += pd.getLatent(i)[0];
								}*/
							case 0:
								if (lasttemp[i] >= melting ) {
									lasttemp[i] = melting;
								}
								break;
							case 1:
								if (lasttemp[i]<=melting ) {
									lasttemp[i] = melting;
									break;
								}
								else if (lasttemp[i] >= evaporation) {
									lasttemp[i] = evaporation;
								}
								break;
							case 2:
								if (lasttemp[i] <= evaporation) {
									lasttemp[i] = evaporation;
								}
								break;
							case -1:

								pd.getLatent(i)[1] += pd.getLatent(i)[0];
								
								if (pd.getLatent(i)[1] < 0) {
									pd.getLatent(i)[1] = 0;
								}
								else if(pd.getLatent(i)[1]  >= Thresholds[0] ) {
									pd.getLatent(i)[1] = Thresholds[0];
								}
								else
									lasttemp[i] = melting;
								break;
							case -2:

								pd.getLatent(i)[2] += pd.getLatent(i)[0];

								if (pd.getLatent(i)[2] < 0) {
									pd.getLatent(i)[2] = 0;
								}
								else if (pd.getLatent(i)[2]  >= Thresholds[1]) {
									pd.getLatent(i)[2] = Thresholds[1];
								}
								else
									lasttemp[i] = evaporation;
							default:
								break;
							}

						}

					}



		}

		static Real sigmaF(Real LatentHeat, Real gasThreshold, Real begin, Real offset) {
			Real Lf = gasThreshold;
			Real a = 1.0 + (12.0 / Lf);
			Real b = (0.7*Lf)- LatentHeat;
			Real denominator = 1 + pow(a, b);
			Real total = (-begin+offset) / denominator + begin;
			return total;
		}

	};



}

#endif

