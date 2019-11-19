#include "FluidModel.h"
#include "PositionBasedDynamics/PositionBasedDynamics.h"
#include "PositionBasedDynamics/SPHKernels.h"
#include "Demos/Simulation/Constraints.h"
#include <iostream>
#include <fstream>
#include "Demos/Utils/TetGenLoader.h"
#include "Demos/Utils/OBJLoader.h"
#include "Demos/Utils/Logger.h"

using namespace PBD;

FluidModel::FluidModel() :
	m_particles()
{	
	m_cloth_stiffness = 1.0;
	m_dilation = 0.001;
	m_density0 = 1000.0;
	m_particleRadius = 0.025;
	viscosity = 0.025;
	m_neighborhoodSearch = NULL;

	t_melting = 0;
	t_evaporation = 100.0;
	t_difusion = 100;
	LatentThreshhold = Vector2r(333, 2257);

	radioSolidification = 0.055;
	numMaxConstraints = 4;

	t_contact = false;
	t_tempContact = 0;

}

FluidModel::~FluidModel(void)
{
	cleanupModel();
}

void FluidModel::cleanupModel()
{
	m_particles.release();
	m_lambda.clear();
	m_density.clear();
	m_deltaX.clear();
	/*
	for (unsigned int i = 0; i < m_constraints.size(); i++)
		delete m_constraints[i];*/
	m_constraints.clear();
	delete m_neighborhoodSearch;
}

void FluidModel::reset()
{
	const unsigned int nPoints = m_particles.size();
	
	for(unsigned int i=0; i < nPoints; i++)
	{
		const Vector3r& x0 = m_particles.getPosition0(i);
		m_particles.getPosition(i) = x0;
		m_particles.setTemp(i, m_particles.getTemp0(i));
		m_particles.getLastPosition(i) = m_particles.getPosition(i);
		m_particles.getOldPosition(i) = m_particles.getPosition(i);
		m_particles.getVelocity(i).setZero();
		m_particles.getAcceleration(i).setZero();
		m_deltaX[i].setZero();
		m_lambda[i] = 0.0;
		m_density[i] = 0.0;
		m_particles.getC(i) = 1000;
		m_particles.getConductivity(i) = 50.0;	
		m_particles.getLatent(i).setZero();

		m_particles.evalState(i,t_melting,t_evaporation);

		switch ((int)m_particles.getState(i)) {

		case 0:
			//, Vector3r(0.0f, this->getLatent()[0], this->getLatent()[1]));
			break;
		case 1:
			this->m_particles.getLatent(i)[1] = this->getLatent()[0];
			break;
		case 2:
			this->m_particles.getLatent(i)[1] = this->getLatent()[0];
			this->m_particles.getLatent(i)[2] = this->getLatent()[1];
			this->m_particles.getVelocity(i)[1] = ((rand() % 900+150.0)*0.0002);

			this->m_particles.getVelocity(i)[1] = this->m_particles.getVelocity(i)[1] + 0.00390625;
			this->m_particles.getPosition(i)[1] = this->m_particles.getPosition(i)[1] + this->m_particles.getVelocity(i)[1] * 0.00390625;
			break;
		case -1:
			this->m_particles.getLatent(i)[1] = this->getLatent()[0]/2;
			break;
		case -2:
			this->m_particles.getLatent(i)[1] = this->getLatent()[0];
			this->m_particles.getLatent(i)[2] = this->getLatent()[1] / 2;
			break;

		}

		/*
			if (i<nPoints / 2) {
				m_particles.getTemp(i) = 80;
			}
			else {
				m_particles.getTemp(i) = 0;
			}
			m_particles.setState(i, 1);
		
		//*/
	}

	/// if exist the model 3d
	this->m_constraints.clear();
	if (m_constraintsIndex.size()>0) {	
		
		int numConstrainst = m_constraintsIndex.size();		
		for (int i = 0; i < numConstrainst;i+=2) {
			this->addDistanceConstraint(m_constraintsIndex[i], m_constraintsIndex[i+1]);
		}
	}

}

ParticleData & PBD::FluidModel::getParticles()
{
	return m_particles;
}
std::vector<newDistanceConstraint> & PBD::FluidModel::getConstraints()
{
	return m_constraints;
}

void FluidModel::initMasses()
{
	const int nParticles = (int) m_particles.size();
	const Real diam = 2.0*m_particleRadius;

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < nParticles; i++)
		{

			//if(m_particles.getActive(i)){

				m_particles.setMass(i, 0.8 * diam*diam*diam * m_density0);		// each particle represents a cube with a side length of r		
																				// mass is slightly reduced to prevent pressure at the beginning of the simulation
			/*}
			else {
				m_particles.setMass(i,0);
			}*/

		}
	}
}


/** Resize the arrays containing the particle data.
*/
void FluidModel::resizeFluidParticles(const unsigned int newSize)
{
	m_particles.resize(newSize);
	m_lambda.resize(newSize);
	m_density.resize(newSize);

	m_deltaX.resize(newSize);
}


/** Release the arrays containing the particle data.
*/
void FluidModel::releaseFluidParticles()
{
	m_particles.release();
	m_lambda.clear();
	m_density.clear();

	m_deltaX.clear();
}

void FluidModel::initModel(const unsigned int nFluidParticles, Vector3r* fluidParticles, const unsigned int nBoundaryParticles, Vector3r* boundaryParticles)
{
	releaseFluidParticles();
	resizeFluidParticles(nFluidParticles);

	// init kernel
	CubicKernel::setRadius(m_supportRadius);

	// copy fluid positions
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)nFluidParticles; i++)
		{
			m_particles.getPosition0(i) = fluidParticles[i];
		}
	}

	m_boundaryX.resize(nBoundaryParticles);
	m_boundaryActive.resize(nBoundaryParticles);
	m_boundaryPsi.resize(nBoundaryParticles);
	m_boundaryVolume.resize(nBoundaryParticles);

	// copy boundary positions
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)nBoundaryParticles; i++)
		{
			m_boundaryX[i] = boundaryParticles[i];
		}
	}

	// initialize masses
	initMasses();

	//////////////////////////////////////////////////////////////////////////
	// Compute value psi for boundary particles (boundary handling)
	// (see Akinci et al. "Versatile rigid - fluid coupling for incompressible SPH", Siggraph 2012
	//////////////////////////////////////////////////////////////////////////

	// Search boundary neighborhood
	NeighborhoodSearchSpatialHashing neighborhoodSearchSH(nBoundaryParticles, m_supportRadius);
	neighborhoodSearchSH.neighborhoodSearch(&m_boundaryX[0]);
	 
	unsigned int **neighbors = neighborhoodSearchSH.getNeighbors();
	unsigned int *numNeighbors = neighborhoodSearchSH.getNumNeighbors();

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int) nBoundaryParticles; i++)
		{
			Real delta = CubicKernel::W_zero();
			for (unsigned int j = 0; j < numNeighbors[i]; j++)
			{
				const unsigned int neighborIndex = neighbors[i][j];
				delta += CubicKernel::W(m_boundaryX[i] - m_boundaryX[neighborIndex]);
			}
			const Real volume = 1.0 / delta;
			m_boundaryPsi[i] = m_density0 * volume;
			m_boundaryVolume[i] = volume;
		}
	}


	// Initialize neighborhood search
	if (m_neighborhoodSearch == NULL)
		m_neighborhoodSearch = new NeighborhoodSearchSpatialHashing(m_particles.size(), m_supportRadius);
	m_neighborhoodSearch->setRadius(m_supportRadius);

	reset();
}

void FluidModel::initModel( const unsigned int nBoundaryParticles, Vector3r* boundaryParticles, std::vector<bool>* boundaryActives)
{
	// init kernel
	CubicKernel::setRadius(m_supportRadius);

	m_boundaryX.resize(nBoundaryParticles);
	m_boundaryActive.resize(nBoundaryParticles);
	m_boundaryVolume.resize(nBoundaryParticles);
	m_boundaryPsi.resize(nBoundaryParticles);

	// copy boundary positions
   #pragma omp parallel default(shared)
	{
  #pragma omp for schedule(static)  
		for (int i = 0; i < (int)nBoundaryParticles; i++)
		{
			m_boundaryX[i] = boundaryParticles[i];
			m_boundaryActive[i] = (*boundaryActives)[i];
		}
	}

	// initialize masses
	initMasses();

	//////////////////////////////////////////////////////////////////////////
	// Compute value psi for boundary particles (boundary handling)
	// (see Akinci et al. "Versatile rigid - fluid coupling for incompressible SPH", Siggraph 2012
	//////////////////////////////////////////////////////////////////////////

	// Search boundary neighborhood
	NeighborhoodSearchSpatialHashing neighborhoodSearchSH(nBoundaryParticles, m_supportRadius);
	neighborhoodSearchSH.neighborhoodSearch(&m_boundaryX[0]);

	unsigned int **neighbors = neighborhoodSearchSH.getNeighbors();
	unsigned int *numNeighbors = neighborhoodSearchSH.getNumNeighbors();

#pragma omp parallel default(shared)
	{
#pragma omp for schedule(static)  
		for (int i = 0; i < (int)nBoundaryParticles; i++)
		{
			Real delta = CubicKernel::W_zero();
			for (unsigned int j = 0; j < numNeighbors[i]; j++)
			{
				const unsigned int neighborIndex = neighbors[i][j];
				delta += CubicKernel::W(m_boundaryX[i] - m_boundaryX[neighborIndex]);
			}
			const Real volume = 1.0 / delta;
			m_boundaryPsi[i] = m_density0 * volume;
			m_boundaryVolume[i] = volume;
		}
	}
	// Initialize neighborhood search
	if (m_neighborhoodSearch == NULL)
		m_neighborhoodSearch = new NeighborhoodSearchSpatialHashing(m_particles.size(), m_supportRadius);
	m_neighborhoodSearch->setRadius(m_supportRadius);

	reset();
}
void FluidModel::importMyModel(const char* path, std::vector<Vector3r>* vertices, std::vector<unsigned int>* edges, float scale, float fac) {

	std::ifstream file(path);
	int nVertices = 0;
	int nEdges = 0;

	if (!file)
		return;

	file >> nVertices >> nEdges;

	float x, y, z,minY=0.1f;
	for (int i = 0; i < nVertices; i++) {

		file >> x >> y >> z;
		if (y < minY) {
			minY = y;
		}
		//if( y < fac)
			vertices->push_back(Vector3r(x, y, z));
	}
	for (int i = 0; i < vertices->size(); i++)
	{
		if ((*vertices)[i][1] < fac) {
			(*vertices)[i][1] = 0;
		}
		else //*/
			(*vertices)[i][1] = (*vertices)[i][1] - minY;
		
		(*vertices)[i] = (*vertices)[i] * scale;
	}
	//*
	unsigned int v1, v2;
	for (int i = 0; i < nEdges; i++) {
		file >> v1 >> v2;
		edges->push_back(v1);
		edges->push_back(v2);
	}//*/
}
void FluidModel::addModel(const std::string &nameFile, Real diam, Real temp, Real scale, Vector3r position,bool fixed,float fac) {

	std::vector<Vector3r> vertices;
	std::vector<unsigned int> edges;

	importMyModel(nameFile.c_str(), &vertices, &edges, scale,fac);

	int nLastParticles = m_particles.getNumberOfParticles();;
	int nVertices = vertices.size();
	int nEdges = edges.size();

	m_lambda.resize(nLastParticles + nVertices);
	m_density.resize(nLastParticles + nVertices);
	m_deltaX.resize(nLastParticles + nVertices);

	for (int i = 0, n = nLastParticles; i < nVertices; i++, n++)
	{

		this->m_particles.addVertex(vertices[i] + position);		

		this->m_particles.setTemp(n, temp);
		this->m_particles.setTemp0(n, temp);
		
		if (fixed && vertices[i][1] < 1.9*diam) {
			this->m_particles.setActive(i, false);
		}

		this->m_particles.evalState(n, t_melting, t_evaporation);
	}

	/*for (unsigned int i = 0; i < nEdges; i += 2)
	{
		const unsigned int v1 = edges[i] + nLastParticles;
		const unsigned int v2 = edges[i + 1] + nLastParticles;

		addDistanceConstraint(v1, v2);
		m_constraintsIndex.push_back(v1);
		m_constraintsIndex.push_back(v2);
	}*/
	//LOG_INFO << "Number vertices:" << nVertices << nEdges;
	printf("vertices: %d edges: %d\n", nVertices, nEdges);

}
void FluidModel::addObject(const std::string &nameFile, Real diam, Real temp, Vector3r size, Vector3r position) {
	std::vector<Vector3r> vertices;
	std::vector<unsigned int> tets;	

	TetGenLoader::loadTetgenModel( nameFile + ".node" , nameFile + ".ele" , vertices, tets);

	int nLastParticles = m_particles.size();
	int nParticles = vertices.size();

	m_lambda.resize(nLastParticles + nParticles);
	m_density.resize(nLastParticles + nParticles);
	m_deltaX.resize(nLastParticles + nParticles);

	position = position*diam;

	for (int i = 0, n = nLastParticles; i < nParticles; i++,n++)
	{
		vertices[i] = 0.5* vertices[i];
		vertices[i][2] = -vertices[i][2];
				this->m_particles.addVertex( vertices[i] +  position);
				this->m_particles.setTemp(n, temp);
				this->m_particles.setTemp0(n, temp);
				this->m_particles.evalState(n, t_melting, t_evaporation);
	}

	TetModel *tetModel = new TetModel();
	tetModel->initMesh(nParticles, (unsigned int)tets.size() / 4u, nLastParticles, tets.data());
	//tetModel->attachVisMesh();

	const unsigned int offset = tetModel->getIndexOffset();
	const unsigned int nEdges = tetModel->getParticleMesh().numEdges();
	const IndexedTetMesh::Edge *edges = tetModel->getParticleMesh().getEdges().data();
	for (unsigned int i = 0; i < nEdges; i++)
	{
		const unsigned int v1 = edges[i].m_vert[0] + offset;
		const unsigned int v2 = edges[i].m_vert[1] + offset;

		addDistanceConstraint(v1, v2);
		m_constraintsIndex.push_back(v1);
		m_constraintsIndex.push_back(v2);
		//tetModel->updateMeshNormals(model.getParticles());
	}
	delete tetModel;

	/*
	IndexedFaceMesh mesh;
	VertexData vertices;
	OBJLoader::loadObj(nameFile, vertices, mesh);

	int nLastParticles = m_particles.size();
	int nParticles = vertices.size();

	m_lambda.resize(nLastParticles + nParticles);
	m_density.resize(nLastParticles + nParticles);
	m_deltaX.resize(nLastParticles + nParticles);

	position = position*diam;

	for (int i = 0, n = nLastParticles; i < nParticles; i++, n++)
	{
		
		vertices.getPosition(i) = 25 * vertices.getPosition(i);
		vertices.getPosition(i)[2] = - vertices.getPosition(i)[2];
		this->m_particles.addVertex(vertices.getPosition(i) + position);
		this->m_particles.setTemp(n, temp);
		this->m_particles.setTemp0(n, temp);
		this->m_particles.evalState(n, t_melting, t_evaporation);
	}

	TriangleModel *triModel = new TriangleModel();

	unsigned int startIndex = m_particles.size();
	triModel->initMesh(nParticles, mesh.numFaces(), nLastParticles, mesh.getFaces().data(), mesh.getUVIndices(), mesh.getUVs());

	const unsigned int offset = triModel->getIndexOffset();
	
		const unsigned int nEdges = triModel->getParticleMesh().numEdges();
		const IndexedFaceMesh::Edge *edges = triModel->getParticleMesh().getEdges().data();
		for (unsigned int i = 0; i < nEdges; i++)
		{
			const unsigned int v1 = edges[i].m_vert[0] + offset;
			const unsigned int v2 = edges[i].m_vert[1] + offset;

			addDistanceConstraint(v1, v2);
		}

		//addTriModel(triModel);*/

}

void FluidModel::addSolid(Real diam,Real temp,Vector3r size, Vector3r position) {

	int widthS = size[0];
	int heightS = size[1];
	int depthS = size[2];

	int nLastParticles = m_particles.size();
	int nParticles = widthS*heightS*depthS;

	m_lambda.resize(nLastParticles + nParticles);
	m_density.resize(nLastParticles + nParticles);
	m_deltaX.resize(nLastParticles + nParticles);

	Vector3r initial = Vector3r(-diam*widthS*0.5, diam, -diam*depthS*0.5);
	//Vector3r initial = Vector3r(-diam * widthS*0.5 + diam * 0.5, diam*0.5, -diam * depthS*0.5 + diam * 0.5);
	for (int i = 0,n= nLastParticles; i < widthS; i++)
	{
		for (unsigned int j = 0; j < heightS; j++)
		{
			for (unsigned int k = 0; k < depthS; k++,n++)
			{
				this->m_particles.addVertex(diam*Vector3r((Real)i, (Real)j, (Real)k) +initial+position );
				this->m_particles.setTemp(n, temp);
				this->m_particles.setTemp0(n, temp);
				this->m_particles.evalState(n,t_melting,t_evaporation);
			}
		}
	}

	std::vector<unsigned int> indices;

	for (unsigned int i = 0; i < widthS - 1; i++)//width
	{
		for (unsigned int j = 0; j < heightS - 1; j++)//height
		{
			for (unsigned int k = 0; k < depthS - 1; k++)//depth
			{
				// For each block, the 8 corners are numerated as:
				//     4*-----*7
				//     /|    /|
				//    / |   / |
				//  5*-----*6 |
				//   | 0*--|--*3
				//   | /   | /
				//   |/    |/
				//  1*-----*2
				unsigned int p0 = i * heightS * depthS + j * depthS + k;
				unsigned int p1 = p0 + 1;
				unsigned int p3 = (i + 1) * heightS * depthS + j * depthS + k;
				unsigned int p2 = p3 + 1;
				unsigned int p7 = (i + 1) * heightS * depthS + (j + 1) * depthS + k;
				unsigned int p6 = p7 + 1;
				unsigned int p4 = i * heightS * depthS + (j + 1) * depthS + k;
				unsigned int p5 = p4 + 1;

				// Ensure that neighboring tetras are sharing faces
				if ((i + j + k) % 2 == 1)
				{
					indices.push_back(p2); indices.push_back(p1); indices.push_back(p6); indices.push_back(p3);
					indices.push_back(p6); indices.push_back(p3); indices.push_back(p4); indices.push_back(p7);
					indices.push_back(p4); indices.push_back(p1); indices.push_back(p6); indices.push_back(p5);
					indices.push_back(p3); indices.push_back(p1); indices.push_back(p4); indices.push_back(p0);
					indices.push_back(p6); indices.push_back(p1); indices.push_back(p4); indices.push_back(p3);
				}
				else
				{
					indices.push_back(p0); indices.push_back(p2); indices.push_back(p5); indices.push_back(p1);
					indices.push_back(p7); indices.push_back(p2); indices.push_back(p0); indices.push_back(p3);
					indices.push_back(p5); indices.push_back(p2); indices.push_back(p7); indices.push_back(p6);
					indices.push_back(p7); indices.push_back(p0); indices.push_back(p5); indices.push_back(p4);
					indices.push_back(p0); indices.push_back(p2); indices.push_back(p7); indices.push_back(p5);
				}
			}
		}
	}

	TetModel *tetModel = new TetModel();

	tetModel->initMesh(widthS*heightS*depthS, (unsigned int)indices.size() / 4u, nLastParticles, indices.data());

	const unsigned int offset = tetModel->getIndexOffset();
	const unsigned int nEdges = tetModel->getParticleMesh().numEdges();
	const IndexedTetMesh::Edge *edges = tetModel->getParticleMesh().getEdges().data();
	for (unsigned int i = 0; i < nEdges; i++)
	{
		const unsigned int v1 = edges[i].m_vert[0] + offset;
		const unsigned int v2 = edges[i].m_vert[1] + offset;

		addDistanceConstraint(v1, v2);
		m_constraintsIndex.push_back(v1);
		m_constraintsIndex.push_back(v2);
		//tetModel->updateMeshNormals(model.getParticles());
	}
	delete tetModel;
	std::cout << "Number of particles Solid: " << nParticles << std::endl;
}

void FluidModel::addFluid(Real diam, Real temp, Vector3r size, Vector3r position) {

	int widthS = size[0];
	int heightS = size[1];
	int depthS = size[2];

	int nLastParticles = m_particles.size();
	int nParticles = widthS*heightS*depthS;

	m_lambda.resize(nLastParticles+nParticles);
	m_density.resize(nLastParticles + nParticles);
	m_deltaX.resize(nLastParticles + nParticles);

	Vector3r initial = Vector3r(-diam*widthS*0.5 + diam * 0.5, diam*0.5, -diam*depthS*0.5 + diam * 0.5);
	for (int i = 0, n = nLastParticles; i < widthS; i++)
	{
		for (unsigned int j = 0; j < heightS; j++)
		{
			for (unsigned int k = 0; k < depthS; k++, n++)
			{
				//this->m_particles.addVertex(diam*Vector3r((Real)i- 28.0f*sin((j+10) / 9.0f) / pow(((j+10) / 5.0f), (1 / 3.0f)), (Real)j, (Real)k + 2*sin(j*0.5)) + initial + position );
				this->m_particles.addVertex(diam*Vector3r((Real)i , (Real)j, (Real)k ) + initial + position);
				this->m_particles.setTemp(n, temp);
				this->m_particles.setTemp0(n, temp);
				this->m_particles.evalState(n, t_melting, t_evaporation);
				this->m_particles.setLatent(i,Vector3r(0.0f,this->getLatent()[0],0.0f));
			}
		}
	}
	std::cout << "Number of particles Fluid: " << nParticles << std::endl;
}

void FluidModel::addGas(Real diam, Real temp, Vector3r size, Vector3r position) {

	int widthS = size[0];
	int heightS = size[1];
	int depthS = size[2];

	int nLastParticles = m_particles.size();
	int nParticles = widthS*heightS*depthS*0.25;

	m_lambda.resize(nLastParticles + nParticles);
	m_density.resize(nLastParticles + nParticles);
	m_deltaX.resize(nLastParticles + nParticles);

	Vector3r initial = Vector3r(-diam*widthS*0.5 + diam * 0.5, 0, -diam*depthS*0.5 + diam * 0.5);
	
	for (int i = 0, n = nLastParticles; i < widthS; i+=2)
	{	
		for (unsigned int j = 0; j < heightS; j++)
		{
			for (unsigned int k = 0; k < depthS; k+=2, n++)
			{
				//this->m_particles.addVertex(diam*Vector3r((Real)i- 28.0f*sin((j+10) / 9.0f) / pow(((j+10) / 5.0f), (1 / 3.0f)), (Real)j, (Real)k + 2*sin(j*0.5)) + initial + position );
				this->m_particles.addVertex(diam*Vector3r((Real)i +(rand()%1000)*0.0018 , (Real)j*5 + (rand() % 1000)*0.0045, (Real)k + (rand() % 1000)*0.0018) + initial + position);

				float currentTemp = temp;//+(rand() % 1000)*0.01;
				//this->m_particles.getVelocity(i)[1] = -(rand() % 1000)*0.025;
				this->m_particles.setTemp(n, currentTemp);
				this->m_particles.setTemp0(n, currentTemp);
				this->m_particles.evalState(n, t_melting, t_evaporation);
				this->m_particles.setLatent(i, Vector3r(0.0f, this->getLatent()[0], this->getLatent()[1]));
			}
		}
	}
	std::cout << "Number of particles Gas: " << nParticles << std::endl;
}

void FluidModel::addAir(Real diam, Real temp, Vector3r size, Vector3r position, int density) {

	int widthS = size[0];
	int heightS = size[1];
	int depthS = size[2];

	int nLastParticles = m_particles.size();
	int nParticles = widthS * heightS*depthS*0.25;

	m_lambda.resize(nLastParticles + nParticles);
	m_density.resize(nLastParticles + nParticles);
	m_deltaX.resize(nLastParticles + nParticles);

	Vector3r initial = Vector3r(-diam * widthS*0.5, 0, -diam * depthS*0.5);

	float ratio = density * 0.001;

	for (int i = 0, n = nLastParticles; i < widthS; i += 2)
	{
		for (unsigned int j = 0; j < heightS; j++)
		{
			for (unsigned int k = 0; k < depthS; k += 2, n++)
			{
				//this->m_particles.addVertex(diam*Vector3r((Real)i- 28.0f*sin((j+10) / 9.0f) / pow(((j+10) / 5.0f), (1 / 3.0f)), (Real)j, (Real)k + 2*sin(j*0.5)) + initial + position );
				this->m_particles.addVertex(diam*Vector3r((Real)i + (rand() % 1000)*0.0019, (Real)j * density + (rand() % 1000) * ratio, (Real)k + (rand() % 1000)*0.0019) + initial + position);

				float currentTemp = temp ;
				this->m_particles.setTemp(n, currentTemp);
				this->m_particles.setTemp0(n, currentTemp);
				this->m_particles.setActive(n,false);
				this->m_particles.evalState(n, t_melting, t_evaporation);
				this->m_particles.setLatent(i, Vector3r(0.0f, this->getLatent()[0], this->getLatent()[1]));
			}
		}
	}
	std::cout << "Number of particles Air: " << nParticles << std::endl;
}

bool FluidModel::addDistanceConstraint(const unsigned int particle1, const unsigned int particle2)
{

	newDistanceConstraint c = newDistanceConstraint();
	const bool res = c.initConstraint(*this, particle1, particle2);
	
	if (res)
	{
		m_constraints.push_back(c);
		//m_groupsInitialized = false;
	}
	return res;
}