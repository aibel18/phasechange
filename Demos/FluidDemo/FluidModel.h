#ifndef __FluidModel_h__
#define __FluidModel_h__

#include "Demos/Simulation/ParticleData.h"
#include <vector>
#include "Demos/Simulation/NeighborhoodSearchSpatialHashing.h"
#include "Demos/Simulation/TetModel.h"
#include "Demos/Simulation/TriangleModel.h"



namespace PBD
{
	class newDistanceConstraint;

	class FluidModel 
	{

		public:

			FluidModel();
			virtual ~FluidModel();

			int heatTransferModel = 0;

		protected:	
			Real viscosity;
			Real m_density0;
			Real m_particleRadius;
			Real m_supportRadius;
			Real m_cloth_stiffness;
			Real m_dilation;

			Real t_difusion;
			Real t_melting;
			Real t_evaporation;

			Vector2r LatentThreshhold;

			bool XPBD;
			bool t_contact;
			Real t_tempContact;

			ParticleData m_particles;
			std::vector<bool> m_boundaryActive;
			std::vector<Real> m_boundaryVolume;
			std::vector<Vector3r> m_boundaryX;
			std::vector<Real> m_boundaryPsi;
			std::vector<Real> m_density;
			
			std::vector<Real> m_lambda;		
			std::vector<Vector3r> m_deltaX;
			NeighborhoodSearchSpatialHashing *m_neighborhoodSearch;		
			std::vector<newDistanceConstraint>  m_constraints;
			std::vector<int>  m_constraintsIndex;

			Real numMaxConstraints;
			Real radioSolidification;

			void initMasses();

			void resizeFluidParticles(const unsigned int newSize);
			void releaseFluidParticles();


		public:
			std::vector<int> particlesFixed;
			void cleanupModel();
			virtual void reset();
			Real getClothStiffness() const { return m_cloth_stiffness; }
			void setClothStiffness(Real val) { m_cloth_stiffness = val; }
			Real getDilation() const { return m_dilation; }
			void setDilation(Real val) { m_dilation = val; }
			Real getMelting() const { return t_melting; }
			void setMelting(Real val) { t_melting = val; }

			Vector2r getLatent() const { return LatentThreshhold; }
			void setLatent(Vector2r val) { LatentThreshhold = val; }

			Real getDifusion() const { return t_difusion; }
			void setDifusion(Real val) { t_difusion = val; }
			Real getEvaporation() const { return t_evaporation; }
			void setEvaporation(Real val) { t_evaporation = val; }


			bool getXPBD() const { return XPBD; }
			void setXPBD(bool val) { XPBD = val; }

			bool getContact() const { return t_contact; }
			void setContact(bool val) { t_contact = val; }
			Real getTempContact() const { return t_tempContact; }
			void setTempContact(Real val) { t_tempContact = val; }
			

			ParticleData &getParticles();
			std::vector<newDistanceConstraint>  &getConstraints();

			void importMyModel(const char* path, std::vector<Vector3r>* vertices, std::vector<unsigned int>* edges, float scale,float fac);
			void addModel(const std::string &nameFile, Real, Real, Real, Vector3r position,bool fixed,float fac);
			void addObject(const std::string &nameFile,Real, Real, Vector3r size, Vector3r position);
			void addSolid(Real,Real,Vector3r size, Vector3r position);
			void addFluid(Real, Real, Vector3r size, Vector3r position);
			void addGas(Real, Real, Vector3r size, Vector3r position);
			void addAir(Real, Real, Vector3r size, Vector3r position, int density);


			void initModel(const unsigned int nFluidParticles, Vector3r* fluidParticles, const unsigned int nBoundaryParticles, Vector3r* boundaryParticles);
			void initModel(const unsigned int nBoundaryParticles, Vector3r* boundaryParticles, std::vector<bool>* boundaryActives);

			const unsigned int numBoundaryParticles() const { return (unsigned int)m_boundaryX.size(); }
			Real getDensity0() const { return m_density0; }
			Real getSupportRadius() const { return m_supportRadius; }
			Real getParticleRadius() const { return m_particleRadius; }
			void setParticleRadius(Real val) { m_particleRadius = val; m_supportRadius = 4*m_particleRadius; }
			NeighborhoodSearchSpatialHashing* getNeighborhoodSearch() { return m_neighborhoodSearch; }	

			bool addDistanceConstraint(const unsigned int particle1, const unsigned int particle2);

			Real getViscosity() const { return viscosity; }
			void setViscosity(Real val) { viscosity = val; }

			Real getNumMaxConstraints() const { return numMaxConstraints; }
			void setNumMaxConstraints(Real val) { numMaxConstraints = val; }

			Real getRadioSolidification() const { return radioSolidification; }
			void setRadioSolidification(Real val) { radioSolidification = val; }

			FORCE_INLINE std::vector<bool>& getBoundaryActives()
			{
				return m_boundaryActive;
			}
			FORCE_INLINE bool getBoundaryActive(const unsigned int i)
			{
				return m_boundaryActive[i];
			}

			FORCE_INLINE void setBoundaryActive(const unsigned int i, bool val)
			{
				m_boundaryActive[i] = val;
			}
			FORCE_INLINE const Vector3r& getBoundaryX(const unsigned int i) const
			{
				return m_boundaryX[i];
			}

			FORCE_INLINE Vector3r& getBoundaryX(const unsigned int i)
			{
				return m_boundaryX[i];
			}

			FORCE_INLINE void setBoundaryX(const unsigned int i, const Vector3r &val)
			{
				m_boundaryX[i] = val;
			}
			FORCE_INLINE const Real& getBoundaryVolume(const unsigned int i) const
			{
				return m_boundaryVolume[i];
			}

			FORCE_INLINE Real& getBoundaryVolume(const unsigned int i)
			{
				return m_boundaryVolume[i];
			}

			FORCE_INLINE void setBoundaryVolume(const unsigned int i, const Real &val)
			{
				m_boundaryVolume[i] = val;
			}

			FORCE_INLINE const Real& getBoundaryPsi(const unsigned int i) const
			{
				return m_boundaryPsi[i];
			}

			FORCE_INLINE Real& getBoundaryPsi(const unsigned int i)
			{
				return m_boundaryPsi[i];
			}

			FORCE_INLINE void setBoundaryPsi(const unsigned int i, const Real &val)
			{
				m_boundaryPsi[i] = val;
			}

			FORCE_INLINE const Real& getLambda(const unsigned int i) const
			{
				return m_lambda[i];
			}

			FORCE_INLINE Real& getLambda(const unsigned int i)
			{
				return m_lambda[i];
			}

			FORCE_INLINE void setLambda(const unsigned int i, const Real &val)
			{
				m_lambda[i] = val;
			}

			FORCE_INLINE const Real& getDensity(const unsigned int i) const
			{
				return m_density[i];
			}

			FORCE_INLINE Real& getDensity(const unsigned int i)
			{
				return m_density[i];
			}

			FORCE_INLINE void setDensity(const unsigned int i, const Real &val)
			{
				m_density[i] = val;
			}

			FORCE_INLINE Vector3r &getDeltaX(const unsigned int i)
			{
				return m_deltaX[i];
			}

			FORCE_INLINE const Vector3r &getDeltaX(const unsigned int i) const
			{
				return m_deltaX[i];
			}

			FORCE_INLINE void setDeltaX(const unsigned int i, const Vector3r &val)
			{
				m_deltaX[i] = val;
			}

	};
}

#endif