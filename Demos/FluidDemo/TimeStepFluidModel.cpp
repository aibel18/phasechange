#include "TimeStepFluidModel.h"
#include "Demos/Simulation/TimeManager.h"
#include "PositionBasedDynamics/PositionBasedFluids.h"
#include "PositionBasedDynamics/TimeIntegration.h"
#include "PositionBasedDynamics/SPHKernels.h"
#include "Demos/Utils/Timing.h"
#include "HeatModel.h"
#include "Demos/Simulation/Constraints.h"

using namespace PBD;
using namespace std;

bool active = true;

TimeStepFluidModel::TimeStepFluidModel()
{
}

TimeStepFluidModel::~TimeStepFluidModel(void)
{
}

void TimeStepFluidModel::step(FluidModel &model)
{
	START_TIMING("simulation step");
	TimeManager *tm = TimeManager::getCurrent ();
	const Real h = tm->getTimeStepSize();
	ParticleData &pd = model.getParticles();

	START_TIMING("acceleration");
	clearAccelerations(model);
	STOP_TIMING_AVG;


	// Update time step size by CFL condition
	//updateTimeStepSizeCFL(model, 0.0001, 0.005);

	// Time integration
	for (unsigned int i = 0; i < pd.size(); i++)
	{ 
		//*
		if (!pd.getActive(i))
			continue;//*/

		model.getDeltaX(i).setZero();
		pd.getLastPosition(i) = pd.getOldPosition(i);
		pd.getOldPosition(i) = pd.getPosition(i);
		TimeIntegration::semiImplicitEuler(h, pd.getMass(i), pd.getPosition(i), pd.getVelocity(i), pd.getAcceleration(i));

	}



	// Perform neighborhood search
	START_TIMING("neighborhood search");
	model.getNeighborhoodSearch()->neighborhoodSearch(&model.getParticles().getPosition(0), model.numBoundaryParticles(), &model.getBoundaryX(0));
	STOP_TIMING_AVG;

	

	// Solve density constraint
	START_TIMING("constraint projection");
	constraintProjection(model);
	STOP_TIMING_AVG;


	// Update velocities	
	for (unsigned int i = 0; i < pd.size(); i++)
	{
		if (m_velocityUpdateMethod == 0)
			TimeIntegration::velocityUpdateFirstOrder(h, pd.getMass(i), pd.getPosition(i), pd.getOldPosition(i), pd.getVelocity(i));
		else
			TimeIntegration::velocityUpdateSecondOrder(h, pd.getMass(i), pd.getPosition(i), pd.getOldPosition(i), pd.getLastPosition(i), pd.getVelocity(i));
	}

	// Compute viscosity
	START_TIMING("viscosity");
	computeXSPHViscosity(model);
	STOP_TIMING_AVG;

	// Compute Heat Transfer
	if (active){
	START_TIMING("heat Transfer");
	HeatTransfer(model);
	STOP_TIMING_AVG;

	START_TIMING("manage Constraints");
	manageConstraints(model);
	STOP_TIMING_AVG;
	}

	// Compute new time	
	tm->setTime (tm->getTime () + h);
	model.getNeighborhoodSearch()->update();
	STOP_TIMING_AVG;
}



/** Clear accelerations and add gravitation.
 */
void TimeStepFluidModel::clearAccelerations(FluidModel &model)
{
	ParticleData &pd = model.getParticles();
	const unsigned int count = pd.size();
	const Vector3r evaporate(0.0, 4, 0.0);
	//Real diam = model.getParticleRadius();
	//Real density = model.getDensity0();

	

#pragma omp parallel default(shared)
	{
#pragma omp for schedule(static)  
		for (int i = 0; i < count; i++)
		{

			if (!pd.getActive(i))
				continue;//*/

			Vector3r grav;

			//*
			if (pd.getState(i) == -2|| pd.getState(i) == 2) {
				//Real val = sigmoide(pd.getLatent(i)[2], model.getLatent()[1], -9.8, -0.0);
				Real val = logarithmic(pd.getLatent(i)[2], model.getLatent()[1], -9.8, 1.0,0.015);
				/*if (i == 0)
					printf("%f val: %f\n", pd.getLatent(i)[2],val);*/
				grav = Vector3r(0.0, val , 0.0);
			}
			else //*/
				grav = Vector3r(0.0, -9.8, 0.0);

			// Clear accelerations of dynamic particles
			if ((pd.getMass(i) != 0.0))
			{
				Vector3r &a = pd.getAcceleration(i);
				a = grav;

			}
		}
	}
}

/** Determine densities of all fluid particles. 

void TimeStepFluidModel::computeDensities(FluidModel &model)
{
	ParticleData &pd = model.getParticles();
	const unsigned int numParticles = pd.size();
	unsigned int **neighbors = model.getNeighborhoodSearch()->getNeighbors();
	unsigned int *numNeighbors = model.getNeighborhoodSearch()->getNumNeighbors();

#pragma omp parallel default(shared)
	{
#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			Real &density = model.getDensity(i);
			Real density_err;
			PositionBasedFluids::computePBFDensity(i, numParticles, &pd.getPosition(0), &pd.getMass(0), &model.getBoundaryX(0), &model.getBoundaryPsi(0), numNeighbors[i], neighbors[i], model.getDensity0(), true, density_err, density);
		}
	}
}
*/
/** Update time step size by CFL condition.
*/
void TimeStepFluidModel::updateTimeStepSizeCFL(FluidModel &model, const Real minTimeStepSize, const Real maxTimeStepSize)
{
	const Real radius = model.getParticleRadius();
	const Real cflFactor = 1.0;
	Real h = TimeManager::getCurrent()->getTimeStepSize();

	// Approximate max. position change due to current velocities
	Real maxVel = 0.1;
	ParticleData &pd = model.getParticles();
	const unsigned int numParticles = pd.size();
	const Real diameter = 2.0*radius;
	for (unsigned int i = 0; i < numParticles; i++)
	{
		const Vector3r &vel = pd.getVelocity(i);
		const Vector3r &accel = pd.getAcceleration(i);
		const Real velMag = (vel + accel*h).squaredNorm();
		if (velMag > maxVel)
			maxVel = velMag;
	}

	// Approximate max. time step size 		
	h = cflFactor * .4 * (diameter / (sqrt(maxVel)));

	h = min(h, maxTimeStepSize);
	h = max(h, minTimeStepSize);

	TimeManager::getCurrent()->setTimeStepSize(h);
}

/** Compute viscosity accelerations.
*/
void TimeStepFluidModel::computeXSPHViscosity(FluidModel &model)
{
	ParticleData &pd = model.getParticles();
	const unsigned int numParticles = pd.size();	

	unsigned int **neighbors = model.getNeighborhoodSearch()->getNeighbors();
	unsigned int *numNeighbors = model.getNeighborhoodSearch()->getNumNeighbors();

	const Real viscosity = model.getViscosity();
	const Real h = TimeManager::getCurrent()->getTimeStepSize();

	const Real tempMean = (model.getEvaporation() - model.getMelting())*0.25;

	// Compute viscosity forces (XSPH)
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{

			if (!pd.getActive(i))
				continue;//*/

			//if (pd.getState(i) == 1 || pd.getState(i) == -1 || pd.getState(i) == -2 ){
			if (pd.getState(i) != 0 ) {
				const Vector3r &xi = pd.getPosition(i);

				Real viscosity2 = sigmoide(pd.getTemp(i) - model.getMelting(), tempMean, 1.0, model.getViscosity());

				Vector3r &vi = pd.getVelocity(i);
				const Real density_i = model.getDensity(i);
				for (unsigned int j = 0; j < numNeighbors[i]; j++)
				{
					const unsigned int neighborIndex = neighbors[i][j];					

					if (neighborIndex < numParticles)		// Test if fluid particle
					{
						if (!pd.getActive(neighborIndex))
							continue;//*/

						// Viscosity for neighbor
						if (pd.getState(neighborIndex) == 2)
							continue;

						const Vector3r &xj = pd.getPosition(neighborIndex);
						const Vector3r &vj = pd.getVelocity(neighborIndex);
						const Real density_j = model.getDensity(neighborIndex);
						vi -= viscosity2 * (pd.getMass(neighborIndex) / density_j) * (vi - vj) * CubicKernel::W(xi - xj);
						

						

					}
					// 				else 
					// 				{
					// 					const Vector3r &xj = model.getBoundaryX(neighborIndex - numParticles);
					// 					vi -= viscosity * (model.getBoundaryPsi(neighborIndex - numParticles) / density_i) * (vi)* CubicKernel::W(xi - xj);
					// 				}
				}
			}
		}
	}
}


void TimeStepFluidModel::HeatTransfer(FluidModel &model) {

	ParticleData &pd = model.getParticles();
	Real particleRadius = model.getParticleRadius();
	const unsigned int count = pd.size();
	const unsigned int nParticles = pd.size();
	unsigned int **neighbors = model.getNeighborhoodSearch()->getNeighbors();
	unsigned int *numNeighbors = model.getNeighborhoodSearch()->getNumNeighbors();
	Real C = 4185.5;
	std::vector<Real> last_temp;
	last_temp.resize(pd.size());

	/*
	Water (liquid): CP = 4185.5 J/(kg⋅K) (15 °C, 101.325 kPa)
	*/
	// Heat Transfers
#pragma omp parallel default(shared)
	{
#pragma omp for schedule(static)  
		for (int i = 0; i < (int)nParticles; i++) 
		{

			if (!pd.getActive(i))
				continue;//*/


			Real difusion = model.getDifusion();
			// the solid has major difusion
			if (pd.getState(i) == 0 || pd.getState(i) == -1)
				difusion = model.getDifusion()*2;
			else if (pd.getState(i) == 2 /*|| pd.getState(i) == -2*/)
				//difusion = model.getDifusion()*0.03125;
				//difusion = model.getDifusion()*0.0625;
				difusion = model.getDifusion()*0.5;


			switch ( model.heatTransferModel) {
			case 0:
				HeatModel::HeatFlux(i, pd, nParticles, neighbors, numNeighbors, &model.getDensity(0), &model.getBoundaryPsi(0), &model.getBoundaryX(0),model.getBoundaryActives(), &last_temp[i], difusion, model.getContact(), model.getTempContact());
				break;
			case 1:
				HeatModel::Cleary(i, pd, nParticles, neighbors, numNeighbors, &model.getDensity(0), &model.getBoundaryPsi(0), &model.getBoundaryVolume(0), &model.getBoundaryX(0), model.getBoundaryActives(), &last_temp[i], model.getContact(), model.getTempContact());
				break;
			}
		}
	}
	//Latent Heat Handle

	Real melting = model.getMelting();
	Real evaporation = model.getEvaporation();
	Vector2r Thresholds = model.getLatent();
	HeatModel::CalculateLatentHeat(pd, melting, evaporation, Thresholds, last_temp);

	//std::cout << " Gravity " << HeatModel::sigmaF(pd.getLatent(1)[2], model.getLatent()[1],-9.8, 4) << "\n" << " Latent " << pd.getLatent(1) << " Temp " << pd.getTemp(1) << "\n";


	//Heat Update
	#pragma omp parallel default(shared)
	{
	#pragma omp for schedule(static)  
		for (int i = 0; i < (int)nParticles; i++)
		{
			//if (i == 1)
				//std::cout << pd.getTemp(i) << "  -> " << last_temp[i] << std::endl;

			/* for condesation of the gas particles
			if (pd.getState(i) == 1 && pd.getPosition(i)[1] > 2.0) {
				/*last_temp[i] = 6.0;m
				pd.setTemp(i, 6.0);
				pd.setLatent(i, Vector3r(0, 0, 0));
				pd.evalState(i, model.getMelting(), model.getEvaporation());* /
				if (i == 1)printf("i%d: %f", i, pd.getVelocity(i)[1]);

				//continue;
			}
			//*/			

			if (!pd.getActive(i))
				continue;//*/

			pd.setTemp(i, last_temp[i]);
			pd.evalState(i, model.getMelting(), model.getEvaporation());			
			//if (i == 0)
				//printf("%f %f\n", pd.getTemp(i), pd.getPosition(i)[0]);
		}
	}
}

void TimeStepFluidModel::reset()
{

}


/** Solve density constraint.
*/
void TimeStepFluidModel::constraintProjection(FluidModel &model)
{
	const unsigned int maxIter = 8;
	unsigned int iter = 0;


	ParticleData &pd = model.getParticles();
	const unsigned int nParticles = pd.size();
	vector<newDistanceConstraint>  &constraints = model.getConstraints();
	const unsigned int nConstraints = constraints.size();
	unsigned int **neighbors = model.getNeighborhoodSearch()->getNeighbors();
	unsigned int *numNeighbors = model.getNeighborhoodSearch()->getNumNeighbors();

	while (iter < maxIter)
	{
		Real avg_density_err = 0.0;

		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static)  
			for (int i = 0; i < (int)nParticles; i++)
			{
				Real density = model.getDensity0();

				// the solid has menor density
				
				if (pd.getState(i) == 0)
					density = density * 0.9;
				else if ((pd.getState(i) == -2  || pd.getState(i) == 2)) {
					density = logarithmic2(pd.getLatent(i)[2], model.getLatent()[1], density, 700,120);
					//density = HeatModel::sigmaF(pd.getLatent(i)[2], model.getLatent()[1], density, 700);
					//if (i == 0)
						//printf("%f val: %f\n", pd.getLatent(i)[2], density);
					

					//std::cout << pd.getLatent(i)[2] << "   " << pd.getState(i) << "   " << factor << std::endl;
				}//*/

				Real density_err;
				PositionBasedFluids::computePBFDensity(i, nParticles, &pd.getPosition(0), &pd.getMass(0), pd.getState(i),pd.getTemp(i),model.getEvaporation(),&model.getBoundaryX(0), &model.getBoundaryPsi(0), numNeighbors[i], neighbors[i], density, true, density_err, model.getDensity(i));
				PositionBasedFluids::computePBFLagrangeMultiplier(i, nParticles, &pd.getPosition(0), &pd.getMass(0), &model.getBoundaryX(0), &model.getBoundaryPsi(0), model.getDensity(i), numNeighbors[i], neighbors[i], density, true, model.getLambda(i));
			}
		}
		
		#pragma omp parallel default(shared)
		{
		#pragma omp for schedule(static)  
			for (int i = 0; i < (int)nParticles; i++){

				Real density = model.getDensity0();

				// the solid has menor density
				//*
				if (pd.getState(i) == 0 )
					density = density * 0.9;
				if ((pd.getState(i) == -2  || pd.getState(i) == 2) ) {
					density = logarithmic2(pd.getLatent(i)[2], model.getLatent()[1], density, 700,120);
					//density = HeatModel::sigmaF(pd.getLatent(i)[2], model.getLatent()[1], density, 700);
					//std::cout << pd.getLatent(i)[2] << "   " << pd.getState(i) << "   " << factor << std::endl;
				}//*/
				
				Vector3r corr;
				PositionBasedFluids::solveDensityConstraint(i, nParticles, &pd.getPosition(0), &pd.getMass(0), pd.getState(i), pd.getLatent(i)[2], model.getLatent()[1],&model.getBoundaryX(0), &model.getBoundaryPsi(0), numNeighbors[i], neighbors[i], density, true, &model.getLambda(0), corr);
				//pd.getPosition(i) += corr;
				model.getDeltaX(i) = corr;
			}
		}

		//*
		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static)  
			for (int i = 0; i < (int)nParticles; i++)
			{

				if (!pd.getActive(i))
					continue;//*/

				pd.getPosition(i) += model.getDeltaX(i);

				/// collition with floor
				if (pd.getPosition(i)[1] < model.getParticleRadius())
					pd.getPosition(i)[1] = model.getParticleRadius();

			}
		
		}//*/

		 /// distance constraint 
		#pragma omp parallel default(shared)
			{
			#pragma omp for schedule(static) 
				for (int i = 0; i < (int)nConstraints; i++) {
					//constraints[i].updateConstraint(model);
					constraints[i].solvePositionConstraint(model);
				}
			}

		iter++;
	}
}

void TimeStepFluidModel::manageConstraints(FluidModel &model) {

	ParticleData &pd = model.getParticles();
	std::vector<newDistanceConstraint> &c = model.getConstraints();
	const unsigned int nc = model.getConstraints().size();
	const unsigned int nParticles = pd.size();
	unsigned int **neighbors = model.getNeighborhoodSearch()->getNeighbors();
	unsigned int *numNeighbors = model.getNeighborhoodSearch()->getNumNeighbors();
	Vector2r Thresholds = model.getLatent();

	
	/* Delete Constraints*/
	std::vector<newDistanceConstraint> newConstraint;

	//#pragma omp parallel default(shared)
	{
		//#pragma omp for schedule(static) 
		for (int i = 0; i < (int)nc; i++)
		{
			if (c[i].evalConstraint(model))
				newConstraint.push_back(c[i]);
		}
	}
	int sizeConstraint = newConstraint.size();
	c.resize(sizeConstraint);

	#pragma omp parallel default(shared)
	{
	#pragma omp for schedule(static) 
		for (int i = 0; i < sizeConstraint; i++)
		{
			c[i] = (newConstraint[i]);
		}
	}

	/* Create Constraints*/
	//#pragma omp parallel default(shared)
	//{
		//#pragma omp for schedule(static)  
		for (int i = 0; i < (int)nParticles; i++)
		{	

			if (!pd.getActive(i))
				continue;//*/

			if(pd.getState(i) != -1 || pd.getLatent(i)[1] >= 20 || pd.getLatent(i)[0] >= 0)
			//if (pd.getState(i) != -1 || pd.getState0(i) != 1 )
			//if (pd.getState(i) != 0 && pd.getState(i) != -1)
				continue;
			
			/*if (i == 1)
				std::cout << pd.getState(i) << "  " << pd.getState0(i) << endl;*/
			const Real numMaxConstraints = model.getNumMaxConstraints();
			const Real radioSolidification = model.getRadioSolidification();
			const Vector3r &xi = pd.getPosition(i);
			int numMax = 0;

			for (unsigned int j = 0; j < numNeighbors[i]; j++)
			{
				const unsigned int neighborIndex = neighbors[i][j];
				if (neighborIndex >= (int)nParticles)
					continue;

				if (!pd.getActive(neighborIndex))
					continue;//*/

				//if (i == 1)
					//std::cout << "->  " << pd.getState(neighborIndex) << "  " << pd.getState0(neighborIndex)<<  << endl;
				if ( pd.getState(neighborIndex) == 0 /*|| pd.getState(neighborIndex) == -1*/ )
				//if (neighborIndex < nParticles && pd.getLatent(i)[1] >= 0 && pd.getLatent(i)[1]  < Thresholds[1] == 0)
				{
					const Vector3r &xj = pd.getPosition(neighborIndex);
					//if (i == 1)
						//std::cout << "  " << pd.getState(neighborIndex) << "  " << pd.getState0(neighborIndex) <<"  "<< (xi - xj).norm() << endl;

					
					if ((xi - xj).norm() < radioSolidification && numMax < numMaxConstraints) {
						model.addDistanceConstraint(i, neighborIndex);
						numMax++;
					}
				}
			}
		}
	//}

}

