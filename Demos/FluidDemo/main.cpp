#include "Common/Common.h"
#include "GL/glew.h"
#include "Demos/Visualization/MiniGL.h"
#include "Demos/Visualization/Selection.h"
#include "GL/glut.h"
#include "Demos/Simulation/TimeManager.h"
#include <Eigen/Dense>
#include "FluidModel.h"
#include "TimeStepFluidModel.h"
#include <iostream>
#include "Demos/Utils/Logger.h"
#include "Demos/Utils/Timing.h"
#include "Demos/Utils/FileSystem.h"

#define _USE_MATH_DEFINES
#include "math.h"
#include "Demos/Simulation/Constraints.h"
#include "Demos/Visualization/Visualization.h"


// Enable memory leak detection
#if defined(_DEBUG) && !defined(EIGEN_ALIGN)
	#define new DEBUG_NEW 
#endif

INIT_TIMING
INIT_LOGGING

using namespace PBD;
using namespace Eigen;
using namespace std;

void timeStep ();
void buildModel ();
void renderWorld();
void createBreakingDam();
void addWall(const Vector3r &minX, const Vector3r &maxX, std::vector<Vector3r> &boundaryParticles, std::vector<bool> &boundaryActives,bool active);
void initBoundaryData(std::vector<Vector3r> &boundaryParticles, std::vector<bool> &boundaryActives);
void render ();
void cleanup();
void reset();
void selection(const Eigen::Vector2i &start, const Eigen::Vector2i &end);
void createSphereBuffers(Real radius, int resolution);
void renderSphere(const Vector3r &x, const float color[]);
void releaseSphereBuffers();
void TW_CALL setTimeStep(const void *value, void *clientData);
void TW_CALL getTimeStep(void *value, void *clientData);
void TW_CALL setVelocityUpdateMethod(const void *value, void *clientData);
void TW_CALL getVelocityUpdateMethod(void *value, void *clientData);
void TW_CALL setTempUpdateMethod(const void *value, void *clientData);
void TW_CALL getTempUpdateMethod(void *value, void *clientData);
void TW_CALL setViscosity(const void *value, void *clientData);
void TW_CALL getViscosity(void *value, void *clientData);
void TW_CALL setStiffness(const void *value, void *clientData);
void TW_CALL getStiffness(void *value, void *clientData);
void TW_CALL setDilation(const void *value, void *clientData);
void TW_CALL getDilation(void *value, void *clientData);
void TW_CALL setMelting(const void *value, void *clientData);
void TW_CALL getMelting(void *value, void *clientData);
void TW_CALL setEvaporation(const void *value, void *clientData);
void TW_CALL getEvaporation(void *value, void *clientData);
void TW_CALL setDifusion(const void *value, void *clientData);
void TW_CALL getDifusion(void *value, void *clientData);
void TW_CALL getTempContact(void *value, void *clientData);
void TW_CALL setTempContact(const void *value, void *clientData);
void TW_CALL getContact(void *value, void *clientData);
void TW_CALL setContact(const void *value, void *clientData);
void TW_CALL getXPBD(void *value, void *clientData);
void TW_CALL setXPBD(const void *value, void *clientData);


void TW_CALL getRadioSol(void *value, void *clientData);
void TW_CALL setRadioSol(const void *value, void *clientData);
void TW_CALL getNumConstraint(void *value, void *clientData);
void TW_CALL setNumConstraint(const void *value, void *clientData);

void TW_CALL getDesWidth(void *value, void *clientData);
void TW_CALL setDesWidth(const void *value, void *clientData);

void TW_CALL getDesHeight(void *value, void *clientData);
void TW_CALL setDesHeight(const void *value, void *clientData);
//Core
FluidModel model;
TimeStepFluidModel simulation;

Real compact = 0.98;// 0.928;

//const Real particleRadius = 0.025;
const Real particleRadius = 0.025;//0.015625//;
int width = 30;
int depth = 30;
int height = 80;
bool doPause = true;
bool thermal = true;
bool record = false;
bool capture = false;
bool rendering = true;
bool renderingWall = false;
bool renderingAir = false;


std::vector<unsigned int> selectedParticles;
Vector3r oldMousePos;
// initiate buffers
GLuint elementbuffer;
GLuint normalbuffer;
GLuint vertexbuffer;
int vertexBufferSize = 0;
GLint context_major_version, context_minor_version;
string exePath;
string dataPath;

GLuint textureFloor;
ofstream file;


// main 
int main( int argc, char **argv )
{
	REPORT_MEMORY_LEAKS

	std::string logPath = FileSystem::normalizePath(FileSystem::getProgramPath() + "/log");
	FileSystem::makeDirs(logPath);
	logger.addSink(unique_ptr<ConsoleSink>(new ConsoleSink(LogLevel::INFO)));
	logger.addSink(unique_ptr<FileSink>(new FileSink(LogLevel::DEBUG, logPath + "/PBD.log")));

	exePath = FileSystem::getProgramPath();
	dataPath = exePath + "/" + std::string(PBD_DATA_PATH);

	// OpenGL
	MiniGL::init (argc, argv, 1280+620, 780, 0, 0, "Fluid demo");

	MiniGL::setClientIdleFunc (60, timeStep);
	
	MiniGL::initLights ();
	MiniGL::loadTexture(textureFloor);

	MiniGL::setKeyFunc(0, 'r', reset);
	//MiniGL::setKeyFunc(1, '1', demo1);
	//MiniGL::setKeyFunc(2, '2', demo2);
	MiniGL::setSelectionFunc(selection);

	//MiniGL::desWidth = 350;
	//MiniGL::desHeight = 100;
	
	MiniGL::getOpenGLVersion(context_major_version, context_minor_version);
	MiniGL::setClientSceneFunc(render);			
	
	//MiniGL::setViewport (40.0, 0.1f, 500.0, Vector3r (0.0, 3.0, 8.0), Vector3r (0.0, 0.0, 0.0));
	MiniGL::setViewport(40.0, 0.1f, 800.0, Vector3r(0.0, 2.1, 7.0), Vector3r(0.0, 0.0, 0.0));

	TwAddVarRW(MiniGL::getTweakBar(), "Pause", TW_TYPE_BOOLCPP, &doPause, " label='Pause' group=Simulation key=SPACE ");
	TwAddVarCB(MiniGL::getTweakBar(), "TimeStepSize", TW_TYPE_REAL, setTimeStep, getTimeStep, &model, " label='Time step size'  min=0.0 max = 0.1 step=0.0001 precision=4 group=Simulation ");
	TwType enumType = TwDefineEnum("VelocityUpdateMethodType", NULL, 0);
	TwAddVarCB(MiniGL::getTweakBar(), "VelocityUpdateMethod", enumType, setVelocityUpdateMethod, getVelocityUpdateMethod, &simulation, " label='Velocity update method' enum='0 {First Order Update}, 1 {Second Order Update}' group=Simulation");

	TwAddVarCB(MiniGL::getTweakBar(), "XPBD", TW_TYPE_BOOL32, setXPBD, getXPBD, &model, " label='XPBD'  group=Simulation ");
	//TwAddVarRW(MiniGL::getTweakBar(), "XPDB", TW_TYPE_BOOLCPP, &XPBD, " label='XPBD' group=Simulation key=x ");
	TwAddVarCB(MiniGL::getTweakBar(), "Viscosity", TW_TYPE_REAL, setViscosity, getViscosity, &model, " label='Viscosity'  min=0.0 max = 5 step=0.001 precision=4 group=Simulation ");

	TwAddVarCB(MiniGL::getTweakBar(), "Radio Solidification", TW_TYPE_REAL, setRadioSol, getRadioSol, &model, " label='Radio Solidification'  min=0.045 max = 0.1 step=0.01 precision=4 group=Simulation ");
	TwAddVarCB(MiniGL::getTweakBar(), "Num Constraints", TW_TYPE_REAL, setNumConstraint, getNumConstraint, &model, " label='Num Constraints'  min=1 max = 10 step=1 precision=4 group=Simulation ");
	TwAddVarCB(MiniGL::getTweakBar(), "Stiffness", TW_TYPE_REAL, setStiffness, getStiffness, &model, " label='Stiffness'  min=0.01 max = 1 step=0.01 precision=4 group=Simulation ");
	TwAddVarCB(MiniGL::getTweakBar(), "Coef Dilation", TW_TYPE_REAL, setDilation, getDilation, &model, " label='Coef Dilation'  min=0 max = 0.01 step=0.001 precision=4 group=Simulation ");

	TwAddVarCB(MiniGL::getTweakBar(), "Temp Melting", TW_TYPE_REAL, setMelting, getMelting, &model, " label='Temp Melting'  min=0.0 max = 300 step=1.0 precision=4 group=Simulation ");
	TwAddVarCB(MiniGL::getTweakBar(), "Temp Evaporation", TW_TYPE_REAL, setEvaporation, getEvaporation, &model, " label='Temp Evaporation'  min=1 max = 500.0 step=1.0 precision=4 group=Simulation ");
	TwAddVarCB(MiniGL::getTweakBar(), "Coef Difusion", TW_TYPE_REAL, setDifusion, getDifusion, &model, " label='Coef Difusion'  min=1 max = 1000 step=5 precision=4 group=Simulation ");

	TwType enumType2 = TwDefineEnum("TemperatureUpdateMethodType", NULL, 0);
	TwAddVarCB(MiniGL::getTweakBar(), "TemperatureUpdateMethod", enumType2, setTempUpdateMethod, getTempUpdateMethod, &model, " label='Heat Transfer method' enum='0 {Our Method}, 1 {Cleary Method}' group=Simulation");

	TwAddVarCB(MiniGL::getTweakBar(), "Contact", TW_TYPE_BOOL32, setContact, getContact, &model, " label='Contact'  group=Simulation ");
	TwAddVarCB(MiniGL::getTweakBar(), "Temp Contact", TW_TYPE_REAL, setTempContact, getTempContact, &model, " label='Temp Contact'  min=-50 max = 500 step=10 precision=4 group=Simulation ");

	TwAddVarRW(MiniGL::getTweakBar(), "Thermal", TW_TYPE_BOOLCPP, &thermal, " label='Thermal' group=Simulation key=t ");

	TwAddVarRW(MiniGL::getTweakBar(), "start", TW_TYPE_BOOLCPP, &record, " label='start' group=Record key=g ");
	TwAddVarRW(MiniGL::getTweakBar(), "capture", TW_TYPE_BOOLCPP, &capture, " label='capture' group=Record key=c ");
	TwAddVarCB(MiniGL::getTweakBar(), "Des Width", TW_TYPE_INT32, setDesWidth, getDesWidth, &model, " label='Des Width'  min=0.0 max = 900 step=1.0  group=Record ");
	TwAddVarCB(MiniGL::getTweakBar(), "Des Height", TW_TYPE_INT32, setDesHeight, getDesHeight, &model, " label='Des Height'  min=0.0 max = 900 step=1.0 group=Record ");
	TwAddVarRW(MiniGL::getTweakBar(), "rendering", TW_TYPE_BOOLCPP, &rendering, " label='rendering' group=Record key=n ");
	TwAddVarRW(MiniGL::getTweakBar(), "renderingWall", TW_TYPE_BOOLCPP, &renderingWall, " label='renderingWall' group=Record key=m ");
	TwAddVarRW(MiniGL::getTweakBar(), "renderingAir", TW_TYPE_BOOLCPP, &renderingAir, " label='renderingAir' group=Record key=i ");

	
	buildModel();

	if (context_major_version >= 3)
		createSphereBuffers((Real)particleRadius, 4);
		//createSphereBuffers((Real)particleRadius, 25);

	glutMainLoop ();	

	cleanup ();
	file.close();

	//Timing::printAverageTimes();
	
	return 0;
}

void cleanup()
{
	delete TimeManager::getCurrent();
	if (context_major_version >= 3)
		releaseSphereBuffers();
}

int test = 0;
int index = 0;
int frameIndex = 0;

int numFrames = 0;
int totalFrames = 0;
int numFile = 0;
int seconds = 10;
int totalFrame = 64 * seconds;

int timeEmision = 0;
int numberFiles = 10;
int firtTimeAugmented = 0;
int tempAugmented = 0;
bool addWallTop = true;
bool temperatureFloor = true;

int numParticlesAir = 0;
void reset()
{
	Timing::printAverageTimes();
	Timing::reset();

	model.reset();
	simulation.reset();
	TimeManager::getCurrent()->setTime(0.0);

	test++;
	index = 0;
	numFile = 0;
	numFrames = 0;
	frameIndex=0;
	capture = false;
	record = false;
	doPause = true;
}

void mouseMove(int x, int y)
{
	Vector3r mousePos;
	MiniGL::unproject(x, y, mousePos);
	const Vector3r diff = mousePos - oldMousePos;

	TimeManager *tm = TimeManager::getCurrent();
	const Real h = tm->getTimeStepSize();

	ParticleData &pd = model.getParticles();
	for (unsigned int j = 0; j < selectedParticles.size(); j++)
	{
		pd.getVelocity(selectedParticles[j]) += 5.0*diff/h;
	}
	oldMousePos = mousePos;
}

void selection(const Eigen::Vector2i &start, const Eigen::Vector2i &end)
{
	std::vector<unsigned int> hits;
	selectedParticles.clear();
	ParticleData &pd = model.getParticles();
	Selection::selectRect(start, end, &pd.getPosition(0), &pd.getPosition(pd.size() - 1), selectedParticles);
	if (selectedParticles.size() > 0)
		MiniGL::setMouseMoveFunc(GLUT_MIDDLE_BUTTON, mouseMove);
	else
		MiniGL::setMouseMoveFunc(-1, NULL);

	MiniGL::unproject(end[0], end[1], oldMousePos);
}

string nameFile = "condensation";

//*
void initData(int numParticles) {

	char bufPath[255];
	sprintf(bufPath, "/%s%d.txt",nameFile ,numFile++);
	file = ofstream(dataPath + bufPath, ios::out | ios::trunc | ios::binary);
	file.seekp(0, ios::beg);
	file.write(reinterpret_cast<char *>(&numParticles), sizeof(int));
}
void saveData(ParticleData* p) {

	if (numFrames == totalFrame) {
		file.close();

		if (numFile >= numberFiles) {
			reset();
			return;
			//exit(0);
		}
		char bufPath[255];
		sprintf(bufPath, "/%s%d.txt", nameFile, numFile++);
		file = ofstream(dataPath + bufPath, ios::out | ios::trunc | ios::binary);
		file.seekp(0, ios::beg);
		numFrames = 0;
	}

	int numParticles = p->getNumberOfParticles();

	file.seekp(0, ios::end);
	float x, y, z, t;
	for (int i = numParticlesAir; i < numParticles; i++) {
		Vector3r position = p->getPosition(i);

		t = p->getTemp(i);
		x = position[0];
		y = position[1];
		z = position[2];

		file.write(reinterpret_cast<char *>(&t), sizeof(float));
		file.write(reinterpret_cast<char *>(&x), sizeof(float));
		file.write(reinterpret_cast<char *>(&y), sizeof(float));
		file.write(reinterpret_cast<char *>(&z), sizeof(float));
	}
	numFrames++;
	totalFrames++;
}
//*/
void timeStep()
{

	if (capture) {
		MiniGL::saveFrame("D:\\heat-video\\frame", test, TimeManager::getCurrent()->getTime(), 255, 4);
		capture = false;
	}
	if (doPause) {
		return;
	}

	if (timeEmision > 0 && totalFrames >= 64 * timeEmision)
		model.setContact(true);
	if (model.getContact() && firtTimeAugmented > 0 && totalFrames >= 64 * firtTimeAugmented) {
		model.setTempContact(model.getTempContact() + tempAugmented);
		firtTimeAugmented = 0;
	}

	// Simulation code
	for (unsigned int i = 0; i < 4; i++) //Iterations Per Frame
	{
		simulation.step(model);
		// Draw simulation model
		/*ParticleData &pd = model.getParticles();
		const unsigned int nParticles = pd.size();
		for (unsigned int i = 0; i < nParticles; i++)
		{
			const Real t = 6.0f;
			if (pd.getPosition(i)[1] > 3) {
				pd.setTemp(i,t);
			}
		}*/

	}
	/* Measure Temp
	ParticleData &pd = model.getParticles();
	LOG_INFO << TimeManager::getCurrent()->getTime() << "\t" << pd.getTemp(0) << "\t" << pd.getTemp(pd.getNumberOfParticles() - 1) << "\t" << pd.getTemp(0) + pd.getTemp(pd.getNumberOfParticles() - 1);
	//*/
	saveData(&model.getParticles());

	
	

	/*
	if (TimeManager::getCurrent()->getTime()>5.0f) {
		Timing::printAverageTimes();
		exit(0);
	}*/
	
	if ( record) {
		MiniGL::saveFrame("D:\\heat-video\\frame", test, TimeManager::getCurrent()->getTime(), 255, 4);
	}
	
}

void buildModel ()
{
	//TimeManager::getCurrent ()->setTimeStepSize (0.05025);
	TimeManager::getCurrent()->setTimeStepSize(0.00390625);
	//TimeManager::getCurrent()->setTimeStepSize(0.001953125);

	createBreakingDam();
		
}


void renderWorld()
{

	/*glDepthMask(GL_TRUE);
	glEnable(GL_DEPTH_TEST);*/

	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, textureFloor);

	//float color[4] = { 1.0f, 1.0f, 1.0f, 0.9f };
	float color[4] = {0.9f, 0.9f, 0.9f, 1.0f};
	//MiniGL::hsvToRgb(1.0f, 1.0f, 1.0f, surfaceColor);
	//*
	


	float c = 40.0f;

	
	glBegin(GL_QUADS);
	float speccolor[4] = { 1.0, 1.0, 1.0, 1.0 };
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, color);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, color);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, speccolor);
	glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 100.0);

	glTexCoord2d(-2.5, -2.5); glVertex3d(-c, 0.0, -c);
	glTexCoord2d(2.5, -2.5); glVertex3d(c, 0.0, -c);
	glTexCoord2d(2.5, 2.5); glVertex3d(c, 0.0, c);
	glTexCoord2d(-2.5, 2.5); glVertex3d(-c, 0.0, c);
	glEnd();
	glBindTexture(GL_TEXTURE_2D, 0);
}


void render ()
{
	if (!rendering) {
		MiniGL::drawTime(TimeManager::getCurrent()->getTime());
		return;
	}
	//renderWorld();
	
	// Draw simulation model
	const ParticleData &pd = model.getParticles();
	const unsigned int nParticles = pd.size();

	int particlesIni = 0;

	if (renderingAir) {
		particlesIni = numParticlesAir;
	}
	
	//glPointSize(4.0);
	
	const Real supportRadius = model.getSupportRadius();
	Real vmax = 0.4*2.0*supportRadius / TimeManager::getCurrent()->getTimeStepSize();
	Real vmin = 0.0;

	if (context_major_version > 3)
	{
		if (thermal) {

			for (unsigned int i = particlesIni; i < nParticles; i++)
			{
				Real v = pd.getVelocity(i).norm();
				v = 0.5*((v - vmin) / (vmax - vmin));
				v = min(128.0*v*v, 0.2);
				float fluidColor[4] = { 1.0f, 1.0f, 1.0f, 1.0f  };
				Real t = min((pd.getTemp(i)) *0.0065, 0.651);

				t = max(t, -0.05);
				//Real t = min((pd.getTemp(i)) *0.0075, 0.61);

				/*Real t = min((pd.getTemp(i)) *0.0081, 0.648);
				MiniGL::hsvToRgb(0.648 - t, 1.0f, 0.8f + (float)v, fluidColor);*/

				//if (t < 0)
					//t = 0;
				MiniGL::hsvToRgb(0.65 - t, 1.0f, 0.8f + (float)v, fluidColor);
				
				renderSphere(pd.getPosition(i), fluidColor);
			}

		}
		else {
			const vector<newDistanceConstraint> c = model.getConstraints();
			const unsigned int nc = c.size();

			for (unsigned int i = particlesIni; i < nParticles; i++)
			{
				Real v = pd.getVelocity(i).norm();
				v = 0.5*((v - vmin) / (vmax - vmin));
				v = min(128.0*v*v, 0.2);
				float fluidColor[4] = { 0.2f, 0.2f, 0.2f, 1.0 };
				MiniGL::hsvToRgb((pd.getState(i) == 0 ? 1 : 0.55), 1.0f, 0.8f + (float)v, fluidColor);
				renderSphere(pd.getPosition(i), fluidColor);
			}

			float color[4] = { 0.2f, 0.2f, 0.2f, 1.0 };
			MiniGL::hsvToRgb(0.25f, 1.0f, 1.0f, color);

			glLineWidth(3);
			glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, color);
			glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, color);
			glBegin(GL_LINES);

			for (unsigned int i = 0; i < nc; i++)
			{
				Vector3r x1 = pd.getPosition(c[i].m_bodies[0]);
				Vector3r x2 = pd.getPosition(c[i].m_bodies[1]);

				glVertex3f(x1[0], x1[1], x1[2]);
				glVertex3f(x2[0], x2[1], x2[2]);

			}
			glEnd();
		}
		if (renderingWall) {
			float surfaceColor[4] = { 0.2f, 0.2f, 0.2f, 1.0 };
			for (unsigned int i = 0; i < model.numBoundaryParticles(); i++)
				renderSphere(model.getBoundaryX(i), surfaceColor);
		}
	}
	else
	{
		glDisable(GL_LIGHTING);
		glBegin(GL_POINTS);
		for (unsigned int i = particlesIni; i < nParticles; i++)
		{
			Real v = pd.getVelocity(i).norm();
			v = 0.5*((v - vmin) / (vmax - vmin));
			v = min(128.0*v*v, 0.5);
			float fluidColor[4] = { 0.2f, 0.2f, 0.2f, 1.0 };
			MiniGL::hsvToRgb(0.55f, 1.0f, 0.5f + (float)v, fluidColor);

			glColor3fv(fluidColor);
			glVertex3v(&pd.getPosition(i)[0]);
		}
		glEnd();

		// 	glBegin(GL_POINTS);
		// 	for (unsigned int i = 0; i < model.numBoundaryParticles(); i++)
		// 	{
		// 		glColor3fv(surfaceColor);
		// 		glVertex3fv(&model.getBoundaryX(i)[0]);
		// 	}
		// 	glEnd();*/

		//glEnable(GL_LIGHTING);
	}



	float red[4] = { 1.0f, 1.0f, 1.0f, 1 };
	for (unsigned int j = 0; j < selectedParticles.size(); j++)
	{
		MiniGL::drawSphere(pd.getPosition(selectedParticles[j]), 0.05f, red);
	}	

	/*
	TetModel *tetModel = model.getTetModel();
	const IndexedFaceMesh &surfaceMesh = tetModel->getSurfaceMesh();
	Visualization::drawMesh(pd, surfaceMesh, tetModel->getIndexOffset(), surfaceColor);
	//*/
	MiniGL::drawTime(TimeManager::getCurrent()->getTime());

	auto str = "s = "+ std::to_string(TimeManager::getCurrent()->getTime());
	
	//MiniGL::drawBitmapText(-0.33,0.25, str.c_str(), str.size()-4, red);
	//frameIndex++;

	renderWorld();
}

/** Create a breaking dam scenario
*/
void createBreakingDam()
{
	const Real diam = 2.0*particleRadius*compact;
	const Real startX = particleRadius;
	const Real startY = particleRadius;
	const Real startZ = particleRadius;

	model.setParticleRadius(particleRadius);

	/* Test for simulation complete (check in)
	width = 80;
	height = 70;
	depth = 60;
	float temperatureContact = 50;
	float temperatureSolid = -5;

	timeEmision = 2;//2s
	numberFiles = 7;

	firtTimeAugmented = 25;//20s
	tempAugmented = 80;

	model.addModel(dataPath + "/bunny4k.txt", diam, temperatureSolid, 0.55f, Vector3r(-1.5, startY, -1.2), false,0.09);
	//model.addModel(dataPath + "/bunny-6001.txt", diam, temperatureSolid, 0.55f, Vector3r(-0.6, startY, -0.5), false, 0.07);
	model.setDilation(0.0001);
	model.setDifusion(25);
	model.setTempContact(temperatureContact);
	nameFile = "simulation";
	//*/


	/* Test for simulation complete - inverse (check in)
	width = 40;
	height = 100;
	depth = 40;
	float temperatureGas = 105.0;
	float temperatureLiquid = 5.1;
	float temperatureContact = 0;
	float temperatureSolid = -50;

	numberFiles = 3;

	model.addGas(diam, temperatureGas, Vector3r(width, 4, depth), Vector3r(0, 45*diam, 0));	

	//model.addModel(dataPath + "/bunny4k.txt", diam, temperatureSolid, 0.55f, Vector3r(-1.5, 2*startY, -1.2), false);
	//model.addModel(dataPath + "/bunny-14032.txt", diam, temperatureSolid, 0.55f, Vector3r(-1.2, startY, -0.9), false, 0.09);
	model.addModel(dataPath + "/bunny-6001.txt", diam, temperatureSolid, 0.55f, Vector3r(-0.75, startY, -0.6), false, 0.07);

	model.setContact(true);
	model.setTempContact(temperatureContact);
	model.setDifusion(50);
	model.setDilation(0.0001);
	model.setNumMaxConstraints(4);
	model.setRadioSolidification(0.055);

	temperatureFloor = false;

	nameFile = "simulationI";

	//*/

	/* Test for heat and cool wather -> compare Transfer Heat (check in)
	width = 60;
	height = 20;//80;
	depth = 30;

	numberFiles = 4;

	float temperatureLiquid = 1, temperatureLiquid2 = 99;
	model.addFluid(diam, temperatureLiquid, Vector3r(width*0.5, 10, depth), Vector3r(-width*diam*0.25, 0, 0));
	model.addFluid(diam, temperatureLiquid2, Vector3r(width*0.5, 10, depth), Vector3r(width*diam*0.25, 0, 0));

	nameFile = "heat";
	//*/



	/* Test for liquid evaporation - convection  (check in)
	height = 70;
	depth = 30;//30;
	width = 30;
	float temperatureLiquid = 20;//80;
	float temperatureContact = 140;

	timeEmision = 2;//2s
	numberFiles = 8;

	model.addFluid(diam, temperatureLiquid, Vector3r(width, 20, depth), Vector3r(0,0.0045,0));
	//model.setContact(true);
	model.setTempContact(temperatureContact);
	model.setDifusion(40);

	nameFile = "evaporation";
	//*/

	/* Test for condensation (check in)
	width = 120;
	height = 175;
	depth = 40;
	float temperatureGas = 120.0;
	//float temperatureLiquid = 15.1;
	float temperatureContact = 20;
	float temperatureSolid = 0.0;

	numberFiles = 2;

	/*float temperatureAir = 10;
	//float temperatureAir1 = 10;
	float temperatureAir2 = 7;
	float temperatureAir3 = 3;
	float temperatureAir4 = 0;
	//model.addAir(diam, temperatureAir, Vector3r(width, 3, depth), Vector3r(0, height*diam - 155 * 2 * startY, 0),14);
	//model.addAir(diam, temperatureAir2, Vector3r(width, 3, depth), Vector3r(0, height*diam - 115 * 2 * startY, 0),13);
	//model.addAir(diam, temperatureAir3, Vector3r(width, 4, depth), Vector3r(0, height*diam - 77 * 2 * startY, 0),11);
	//model.addAir(diam, temperatureAir4, Vector3r(width, 4, depth), Vector3r(0, height*diam - 35 * 2 * startY, 0),10);

	numParticlesAir = model.getParticles().size();* /

	model.addGas(diam, temperatureGas, Vector3r(width, 17, depth), Vector3r(0, diam , 0));

	model.setContact(true);
	model.setTempContact(temperatureContact);
	model.setDifusion(20);
	model.setViscosity(0.001);

	nameFile = "condensation";
	//*/

	/* Test for melting bunny with heat liquid (check in)

	// change time step = 0.001953125, comment code add border top
	width = 80;
	height = 100;
	depth = 60;
	float temperatureLiquid = 99.0f;
	float temperatureSolid = -5;

	numberFiles = 2;

	model.addModel(dataPath + "/bunny4k.txt", diam, temperatureSolid, 0.55f, Vector3r(-1.5, startY, -1.2), false, 0.09);	
	model.addFluid(diam, temperatureLiquid, Vector3r(5, 300, 5), Vector3r(0.35, 3.0, 0.4));
	
	model.setDilation(0.0001);
	model.setDifusion(30);

	nameFile = "melting";
	addWallTop = false;
	//*/


	/* Test for solidify bunny (check) -
	width = 80;
	height = 100;
	depth = 60;
	float temperatureLiquid = 10.0f;
	float temperatureSolid = -40;

	numberFiles = 1;

	model.addFluid(diam, temperatureLiquid, Vector3r(8, 200,4), Vector3r(-0.25, 3.2, 0.1));
	//model.addFluid(diam, temperatureLiquid, Vector3r(5, 130, 2), Vector3r(-0.0, 3.0, 0.4));
	//model.addModel(dataPath + "/bunny-6001.txt", diam, temperatureSolid,0.55f, Vector3r(-0.75, startY, -0.6),false,0.07);
	model.addModel(dataPath + "/bunny4k.txt", diam, temperatureSolid,0.55f, Vector3r(-1.5, startY, -1.2),false,0.09);
	model.setDilation(0.0001);
	model.setNumMaxConstraints(6);
	model.setRadioSolidification(0.1);
	model.setDifusion(20);


	nameFile = "solidify";
	addWallTop = false;
	//*/


	/* Test for ( dilatation) --
	width = 60;
	height = 60;
	depth = 50;
	float temperatureSolid = -50;
	float temperarureContact = -1;
	//model.addSolid(diam, temperatureSolid,Vector3r(20,20,20),Vector3r(0, 0,0));
	model.addModel(dataPath + "/bunny4k.txt", diam, temperatureSolid, 0.55f, Vector3r(-1.5, startY, -1.2), false, 0.09);
	//model.setContact(true);
	model.setDilation(0.01);
	model.setTempContact(temperarureContact);
	model.setDifusion(50);
	//*/

	/* Test for armadillo on liquid - caendo -- (check in)
	width = 60;
	height = 80;
	depth = 60;
	float temperatureLiquid = 50;
	float temperatureSolid = -40;

	numberFiles = 15;

	model.addFluid(diam, temperatureLiquid, Vector3r(width, 40, depth), Vector3r(0, 0, 0));
	model.addModel(dataPath + "/armadillo4k.txt", diam, temperatureSolid, 0.55f, Vector3r(-1.25, 4.0, -depth*diam * 0.25), false,0.07);
	model.setContact(false);
	model.setDilation(0.0001);
	model.setTempContact(temperatureLiquid);
	model.setDifusion(20);

	nameFile = "armadillo";
	addWallTop = false;
	//*/



	/////////////////////////////////////////////////


	/* Test for solidification del liquido
	width = 100;
	height = 50;
	depth = 30;
	float temperatureLiquid = 10;
	float temperatureSolid = 0;
	model.addFluid(diam, temperatureLiquid, Vector3r(width*0.25, 20, depth), Vector3r(-width*diam*3/8 + particleRadius, startY, 0));

	model.setContact(true);
	model.setTempContact(-5);
	//model.addSolid(diam, 0, Vector3r(20, 20, 20), Vector3r(0.24, 2 * startY, 0));

	nameFile = "simulationI";

	//*/
	

	/* Test for multiples ice on liquid - caendo --
	width = 80;
	height = 120;
	depth = 50;
	float temperatureLiquid = 10;
	float temperatureSolid = 0;
	model.addFluid(diam, temperatureLiquid, Vector3r(width, 20, depth), Vector3r(startX, 2*startY, startZ));
	model.addSolid(diam, temperatureSolid, Vector3r(14, 14, 14), Vector3r(-1.5, 1.7, -0.5));
	model.addSolid(diam, temperatureSolid, Vector3r(17, 17, 17), Vector3r(1.3, 1.5, 0.5));
	model.addSolid(diam, temperatureSolid, Vector3r(10, 10, 10), Vector3r(-0.5, 2.7, -0.5));
	model.addSolid(diam, temperatureSolid, Vector3r(5, 5, 5), Vector3r(0.0, 3, -0.0));
	model.addSolid(diam, temperatureSolid, Vector3r(5, 5, 5), Vector3r(1.0, 3.2, -0.2));
	model.addSolid(diam, temperatureSolid, Vector3r(10, 10, 10), Vector3r(-1.2, 2.5, 0.5));
	//model.addModel(dataPath + "/armadillo4k.txt", diam, temperatureSolid, 0.45f, Vector3r(-width*diam * 0.25, 2.0, -depth*diam * 0.25),false);
	model.setContact(true);
	model.setDilation(0.0001);
	model.setTempContact(temperatureLiquid);
	model.setDifusion(25);
	//*/

	/* Test for armadillo on liquid - caendo -- (check in)
	width = 100;
	height = 250;
	depth = 70;
	float temperatureLiquid = 50;
	float temperatureSolid = 0;
	model.addFluid(diam, temperatureLiquid, Vector3r(width, 15, depth), Vector3r(startX, 2 * startY, startZ));
	model.addModel(dataPath + "/armadillo4k.txt", diam, temperatureSolid, 0.55f, Vector3r(-width*diam * 0.25, 2.0, -depth*diam * 0.25), false);
	model.setContact(false);
	model.setDilation(0.0001);
	model.setTempContact(temperatureLiquid);
	model.setDifusion(20);
	//*/

	


	/* Test for meeltin ( derretimiento) --
	width = 50;
	height = 50;
	depth = 50;
	float temperatureSolid = 0;
	float temperarureContact = 10;
	model.addSolid(diam, temperatureSolid,Vector3r(30,40,20),Vector3r(particleRadius, 2 * startY,0));
	model.setContact(true);
	model.setDilation(0.0001);
	model.setTempContact(temperarureContact);
	model.setDifusion(10);
	//*/

	

	

	/* Test for heat and cool wather -> compare Transfer Heat cleary
	width = 80;
	height = 80;
	depth = 40;
	float temperatureLiquid = 6, temperatureLiquid2 = 80;
	model.heatTransferModel = 1;//model cleary
	model.addFluid(diam, temperatureLiquid, Vector3r(width*0.5, 15, depth), Vector3r(-width*diam*0.25 , 0, 0));
	model.addFluid(diam, temperatureLiquid2, Vector3r(width*0.5, 15, depth), Vector3r(width*diam*0.25 , 0, 0));

	model.setContact(false);
	model.setTempContact(temperatureLiquid2);
	//*/

	/* Test for liquidos caendo
	model.addFluid(diam, 21, Vector3r(3, 20, 3), Vector3r(-0.36, 2.0, 0.0));
	model.addFluid(diam, 55, Vector3r(3, 20, 3), Vector3r(0.16, 0.8, 0.18));	
	model.addFluid(diam, 99, Vector3r(3, 20, 3), Vector3r(0.22, 2.8, -0.2));
	model.addSolid(diam, 0, Vector3r(25, 10, 15), Vector3r(0.025, 0, 0));
	//*/

	
	

	

	

	/* Test for simulation complete - inverse (check in)
	width = 100;
	height = 290;
	depth = 60;
	float temperatureGas = 109.0;
	//float temperatureLiquid = 15.1;
	float temperatureContact = 0;
	float temperatureSolid = 0.0;

	float temperatureAir = 10;
	//float temperatureAir1 = 10;
	float temperatureAir2 = 7;
	float temperatureAir3 = 3;
	float temperatureAir4 = 0;
	model.addAir(diam, temperatureAir, Vector3r(width, 3, depth), Vector3r(particleRadius, height*diam - 155 * 2 * startY, particleRadius), 14);
	model.addAir(diam, temperatureAir2, Vector3r(width, 3, depth), Vector3r(particleRadius, height*diam - 115 * 2 * startY, particleRadius), 13);
	model.addAir(diam, temperatureAir3, Vector3r(width, 4, depth), Vector3r(particleRadius, height*diam - 77 * 2 * startY, particleRadius), 11);
	model.addAir(diam, temperatureAir4, Vector3r(width, 4, depth), Vector3r(particleRadius, height*diam - 35 * 2 * startY, particleRadius), 10);

	numParticlesAir = model.getParticles().size();

	model.addGas(diam, temperatureGas, Vector3r(width, 20, depth), Vector3r(particleRadius, height*diam - 265 * 2 * startY, particleRadius));
	//model.addModel(dataPath + "/bunny.txt", diam, temperatureSolid,0.55f, Vector3r(-0.5, 2*startY, -0.3),false);
	//model.addModel(dataPath + "/bunny4k.txt", diam, temperatureSolid, 0.55f, Vector3r(-1.5, 2 * startY, -1.2), false);

	//printf("%f",);
	//model.addFluid(diam, temperatureLiquid, Vector3r(width, 10, depth), Vector3r(particleRadius, 2 * startY, particleRadius));

	model.setContact(true);
	model.setTempContact(temperatureContact);
	model.setDifusion(30);
	//*/

	/* Test for heat and cool wather -> measure temp (check in)
	width = 50;
	height = 20;
	depth = 1;
	float temperatureLiquid = 1, temperatureLiquid2 = 99;
	model.addFluid(diam, temperatureLiquid, Vector3r(width*0.5, 10, depth), Vector3r(-width*diam*0.25, 0, 0));
	model.addFluid(diam, temperatureLiquid2, Vector3r(width*0.5, 10, depth), Vector3r(width*diam*0.25, 0, 0));
	nameFile = "measure temp";
	//*/

	/* Test for melting bunny - number particles (check in)
	width = 80;
	height = 250;
	depth = 60;
	float temperatureSolid = -5;
	
	timeEmision = 2;
	numberFiles = 1;
	

	//model.addModel(dataPath + "/bunny-46607.txt", diam, temperatureSolid, 0.55f, Vector3r(-1.5, startY, -1.2), false,0.09);nameFile = "46607Particles";//0.09
	//model.addModel(dataPath + "/bunny-27152.txt", diam, temperatureSolid, 0.55f, Vector3r(-1.5, startY, -1.2), false, 0.09); nameFile = "27152Particles";//0.09
	//model.addModel(dataPath + "/bunny-14032.txt", diam, temperatureSolid, 0.55f, Vector3r(-1.2, startY, -0.9), false, 0.09);nameFile = "14032Particles";//0.09
	model.addModel(dataPath + "/bunny-6001.txt", diam, temperatureSolid, 0.55f, Vector3r(-0.8, startY, -0.7), false,0.08);nameFile = "6001Particles";//0.08
	//model.addModel(dataPath + "/bunny-1849.txt", diam, temperatureSolid, 0.55f, Vector3r(-0.6, startY, -0.5), false,0.07);nameFile = "1849Particles";//0.07
	

	model.setDilation(0.001);
	model.setDifusion(50);
	model.setTempContact(50);	
	
	//*/
	/* Test for solidify bunny (check)
	width = 80;
	height = 300;
	depth = 60;
	float temperatureLiquid = 5.6f;
	float temperatureSolid = 0;

	model.addFluid(diam, temperatureLiquid, Vector3r(8, 200,4), Vector3r(-0.25, 3.2, 0.1));
	//model.addFluid(diam, temperatureLiquid, Vector3r(5, 130, 2), Vector3r(-0.0, 3.0, 0.4));
	//model.addModel(dataPath + "/bunny.txt", diam, temperatureSolid,0.55f, Vector3r(-0.3, 2*startY, 0.0),false);
	model.addModel(dataPath + "/bunny4k.txt", diam, temperatureSolid,0.55f, Vector3r(-1.5, 2*startY, -1.2),false);
	model.setDilation(0.0001);
	model.setNumMaxConstraints(6);
	model.setRadioSolidification(0.05);
	model.setDifusion(50);
	//*/

	//* Test for multiples demos - number particles (check in)
	width = 80;
	height = 90;
	depth = 60;
	float temperatureSolid = 50;

	//timeEmision = 1;
	numberFiles = 1;

	float scale = 0.6f;

	//model.addModel(dataPath + "/bunny-46607.txt", diam, temperatureSolid, scale, Vector3r(-1.65, startY, -1.3), false,0.09);//0.09
	//model.addModel(dataPath + "/bunny-27152.txt", diam, temperatureSolid, scale, Vector3r(-1.5, startY, -1.2), false, 0.09);//0.09
	//model.addModel(dataPath + "/bunny-14032.txt", diam, temperatureSolid, scale, Vector3r(-1.2, startY, -0.9), false, 0.09);//0.09
	//model.addModel(dataPath + "/bunny-6001.txt", diam, temperatureSolid, scale, Vector3r(-0.8, startY, -0.7), false,0.08);//0.08
	model.addModel(dataPath + "/bunny-1849.txt", diam, temperatureSolid, scale, Vector3r(-0.6, startY, -0.5), false,0.07);//0.07


	model.setDilation(0.001);
	model.setDifusion(50);
	model.setTempContact(-10);
	model.setContact(true);
	model.setRadioSolidification(0.05);
	model.setNumMaxConstraints(3);
	nameFile = "demos";

	//*/



	std::vector<Vector3r> boundaryParticles;
	std::vector<bool> boundaryActives;
	initBoundaryData(boundaryParticles, boundaryActives);

	model.initModel( (unsigned int)boundaryParticles.size(), boundaryParticles.data(), &boundaryActives);
	
	
	LOG_INFO << "Number particles:" << (model.getParticles().size() - numParticlesAir);
	//LOG_INFO << "Number of particles Total: " << model.getParticles().size();
	initData(model.getParticles().size() - numParticlesAir);
}


void addWall(const Vector3r &minX, const Vector3r &maxX, std::vector<Vector3r> &boundaryParticles, std::vector<bool> &boundaryActives,bool active)
{
	const Real particleDistance = 2*model.getParticleRadius()*compact;

	const Vector3r diff = maxX - minX;
	const unsigned int stepsX = (unsigned int)(round(diff[0] / particleDistance)) + 1u;
	const unsigned int stepsY = (unsigned int)(round(diff[1] / particleDistance)) + 1u;
	const unsigned int stepsZ = (unsigned int)(round(diff[2] / particleDistance)) + 1u;

	//printf("%d %d %d %f\n",stepsX, (unsigned int)(41.00),(unsigned int)(diff[0] / particleDistance), diff[0] / particleDistance);

	const unsigned int startIndex = (unsigned int) boundaryParticles.size();
	boundaryParticles.resize(startIndex + stepsX*stepsY*stepsZ);
	boundaryActives.resize(startIndex + stepsX*stepsY*stepsZ);

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int j = 0; j < (int)stepsX; j++)
		{
			for (unsigned int k = 0; k < stepsY; k++)
			{
				for (unsigned int l = 0; l < stepsZ; l++)
				{
					const Vector3r currPos = minX + Vector3r(j*particleDistance, k*particleDistance, l*particleDistance);
					boundaryParticles[startIndex + j*stepsY*stepsZ + k*stepsZ + l] = currPos;
					boundaryActives[startIndex + j*stepsY*stepsZ + k*stepsZ + l] = active;
				}
			}
		}
	}
}
void initBoundaryData(std::vector<Vector3r> &boundaryParticles, std::vector<bool> &boundaryActives)
{

	const Real containerWidth = (width)*particleRadius*compact;
	const Real containerDepth = (depth)*particleRadius*compact;
	const Real containerHeight = (height)*particleRadius*compact;


	const Real diameter = 2*particleRadius*compact;// compact;// *0.8;

	const Real x1 = -containerWidth  - particleRadius * compact;
	const Real x2 = containerWidth + particleRadius * compact;
	const Real y1 = -(sqrt(3) - 0.80)*particleRadius;
	const Real y2 = 2*containerHeight+ diameter - (sqrt(3) - 0.80)*particleRadius;
	const Real z1 = -containerDepth - particleRadius * compact;
	const Real z2 = containerDepth + particleRadius * compact;

	float rest = 0.001f;

	// Floor
	addWall(Vector3r(x1, y1, z1), Vector3r(x2, y1, z2), boundaryParticles, boundaryActives,temperatureFloor);
	// Top
	if(addWallTop)
		addWall(Vector3r(x1, y2, z1), Vector3r(x2, y2, z2), boundaryParticles, boundaryActives, true);
	// Left
	addWall(Vector3r(x1 - rest, y1+ diameter, z1+ diameter), Vector3r(x1 - rest, y2 - diameter, z2- diameter), boundaryParticles, boundaryActives, false);
	// Right
	addWall(Vector3r(x2 + rest, y1 + diameter, z1+ diameter), Vector3r(x2 + rest, y2 - diameter, z2- diameter), boundaryParticles, boundaryActives, false);
	// Back
	addWall(Vector3r(x1, y1 + diameter, z1), Vector3r(x2, y2 - diameter, z1), boundaryParticles, boundaryActives, false);
	// Front
	addWall(Vector3r(x1, y1 + diameter, z2), Vector3r(x2, y2- diameter, z2), boundaryParticles, boundaryActives, false);
}


void createSphereBuffers(Real radius, int resolution)
{
	Real PI = static_cast<Real>(M_PI);
	// vectors to hold our data
	// vertice positions
	std::vector<Vector3r> v;
	// normals
	std::vector<Vector3r> n;
	std::vector<unsigned short> indices;

	// initiate the variable we are going to use
	Real X1, Y1, X2, Y2, Z1, Z2;
	Real inc1, inc2, inc3, inc4, radius1, radius2;

	for (int w = 0; w < resolution; w++)
	{
		for (int h = (-resolution / 2); h < (resolution / 2); h++)
		{
			inc1 = (w / (Real)resolution) * 2 * PI;
			inc2 = ((w + 1) / (Real)resolution) * 2 * PI;
			inc3 = (h / (Real)resolution)*PI;
			inc4 = ((h + 1) / (Real)resolution)*PI;

			X1 = sin(inc1);
			Y1 = cos(inc1);
			X2 = sin(inc2);
			Y2 = cos(inc2);

			// store the upper and lower radius, remember everything is going to be drawn as triangles
			radius1 = radius*cos(inc3);
			radius2 = radius*cos(inc4);

			Z1 = radius*sin(inc3);
			Z2 = radius*sin(inc4);

			// insert the triangle coordinates
			v.push_back(Vector3r(radius1*X1, Z1, radius1*Y1));
			v.push_back(Vector3r(radius1*X2, Z1, radius1*Y2));
			v.push_back(Vector3r(radius2*X2, Z2, radius2*Y2));

			indices.push_back((unsigned short)v.size() - 3);
			indices.push_back((unsigned short)v.size() - 2);
			indices.push_back((unsigned short)v.size() - 1);

			v.push_back(Vector3r(radius1*X1, Z1, radius1*Y1));
			v.push_back(Vector3r(radius2*X2, Z2, radius2*Y2));
			v.push_back(Vector3r(radius2*X1, Z2, radius2*Y1));

			indices.push_back((unsigned short)v.size() - 3);
			indices.push_back((unsigned short)v.size() - 2);
			indices.push_back((unsigned short)v.size() - 1);

			// insert the normal data
			n.push_back(Vector3r(X1, Z1, Y1));
			n.push_back(Vector3r(X2, Z1, Y2));
			n.push_back(Vector3r(X2, Z2, Y2));
			n.push_back(Vector3r(X1, Z1, Y1));
			n.push_back(Vector3r(X2, Z2, Y2));
			n.push_back(Vector3r(X1, Z2, Y1));
		}
	}

	for (unsigned int i = 0; i < n.size(); i++)
		n[i].normalize();


	glGenBuffersARB(1, &vertexbuffer);
	glBindBufferARB(GL_ARRAY_BUFFER_ARB, vertexbuffer);
	glBufferDataARB(GL_ARRAY_BUFFER_ARB, v.size() * sizeof(Vector3r), &v[0], GL_STATIC_DRAW);

	glGenBuffersARB(1, &normalbuffer);
	glBindBufferARB(GL_ARRAY_BUFFER_ARB, normalbuffer);
	glBufferDataARB(GL_ARRAY_BUFFER_ARB, n.size() * sizeof(Vector3r), &n[0], GL_STATIC_DRAW);
	glBindBufferARB(GL_ARRAY_BUFFER_ARB, 0);

	// Generate a buffer for the indices as well
	glGenBuffersARB(1, &elementbuffer);
	glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, elementbuffer);
	glBufferDataARB(GL_ELEMENT_ARRAY_BUFFER_ARB, indices.size() * sizeof(unsigned short), &indices[0], GL_STATIC_DRAW);
	glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, 0);

	// store the number of indices for later use
	vertexBufferSize = (unsigned int)indices.size();

	// clean up after us
	indices.clear();
	n.clear();
	v.clear();
}

void renderSphere(const Vector3r &x, const float color[])
{
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_NORMAL_ARRAY);

	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, color);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, color);


	glBindBufferARB(GL_ARRAY_BUFFER_ARB, vertexbuffer);
	glVertexPointer(3, GL_REAL, 0, 0);

	glBindBufferARB(GL_ARRAY_BUFFER_ARB, normalbuffer);
	glNormalPointer(GL_REAL, 0, 0);

	glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, elementbuffer);

	glPushMatrix();
	glTranslated(x[0], x[1], x[2]);
	glDrawElements(GL_TRIANGLES, (GLsizei)vertexBufferSize, GL_UNSIGNED_SHORT, 0);
	glPopMatrix();
	glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, 0);
	glBindBufferARB(GL_ARRAY_BUFFER_ARB, 0);

	glDisableClientState(GL_NORMAL_ARRAY);
	glDisableClientState(GL_VERTEX_ARRAY);
}

void releaseSphereBuffers()
{
	if (elementbuffer != 0)
	{
		glDeleteBuffersARB(1, &elementbuffer);
		elementbuffer = 0;
	}
	if (normalbuffer != 0)
	{
		glDeleteBuffersARB(1, &normalbuffer);
		normalbuffer = 0;
	}
	if (vertexbuffer != 0)
	{
		glDeleteBuffersARB(1, &vertexbuffer);
		vertexbuffer = 0;
	}
}




void TW_CALL setTimeStep(const void *value, void *clientData)
{
	const Real val = *(const Real *)(value);
	TimeManager::getCurrent()->setTimeStepSize(val);
}

void TW_CALL getTimeStep(void *value, void *clientData)
{
	*(Real *)(value) = TimeManager::getCurrent()->getTimeStepSize();
}

void TW_CALL setVelocityUpdateMethod(const void *value, void *clientData)
{
	const short val = *(const short *)(value);
	((TimeStepFluidModel*)clientData)->setVelocityUpdateMethod((unsigned int)val);
}

void TW_CALL getVelocityUpdateMethod(void *value, void *clientData)
{
	*(short *)(value) = (short)((TimeStepFluidModel*)clientData)->getVelocityUpdateMethod();
}
//////

void TW_CALL setTempUpdateMethod(const void *value, void *clientData)
{
	const short val = *(const short *)(value);
	((FluidModel*)clientData)->heatTransferModel  = (unsigned int)val;
}

void TW_CALL getTempUpdateMethod(void *value, void *clientData)
{
	*(short *)(value) = (short)((FluidModel*)clientData)->heatTransferModel ;
}

//////
void TW_CALL setViscosity(const void *value, void *clientData)
{
	const Real val = *(const Real *)(value);
	((FluidModel*)clientData)->setViscosity(val);
}

void TW_CALL getViscosity(void *value, void *clientData)
{
	*(Real *)(value) = ((FluidModel*)clientData)->getViscosity();
}

void TW_CALL getRadioSol(void *value, void *clientData)
{
	*(Real *)(value) = ((FluidModel*)clientData)->getRadioSolidification();
}
void TW_CALL setRadioSol(const void *value, void *clientData) {
	const Real val = *(const Real *)(value);
	((FluidModel*)clientData)->setRadioSolidification(val);
}
void TW_CALL getNumConstraint(void *value, void *clientData)
{
	*(Real *)(value) = ((FluidModel*)clientData)->getNumMaxConstraints();
}
void TW_CALL setNumConstraint(const void *value, void *clientData)
{
	const Real val = *(const Real *)(value);
	((FluidModel*)clientData)->setNumMaxConstraints(val);
}
void TW_CALL setStiffness(const void *value, void *clientData)
{
	const Real val = *(const Real *)(value);
	((FluidModel*)clientData)->setClothStiffness(val);
}

void TW_CALL getStiffness(void *value, void *clientData)
{
	*(Real *)(value) = ((FluidModel*)clientData)->getClothStiffness();
}

void TW_CALL setDilation(const void *value, void *clientData)
{
	const Real val = *(const Real *)(value);
	((FluidModel*)clientData)->setDilation(val);
}

void TW_CALL getDilation(void *value, void *clientData)
{
	*(Real *)(value) = ((FluidModel*)clientData)->getDilation();
}

void TW_CALL setMelting(const void *value, void *clientData)
{
	const Real val = *(const Real *)(value);
	((FluidModel*)clientData)->setMelting(val);
}

void TW_CALL getMelting(void *value, void *clientData)
{
	*(Real *)(value) = ((FluidModel*)clientData)->getMelting();
}

void TW_CALL setEvaporation(const void *value, void *clientData)
{
	const Real val = *(const Real *)(value);
	((FluidModel*)clientData)->setEvaporation(val);
}

void TW_CALL getEvaporation(void *value, void *clientData)
{
	*(Real *)(value) = ((FluidModel*)clientData)->getEvaporation();
}

void TW_CALL setDifusion(const void *value, void *clientData)
{
	const Real val = *(const Real *)(value);
	((FluidModel*)clientData)->setDifusion(val);
}

void TW_CALL getDifusion(void *value, void *clientData)
{
	*(Real *)(value) = ((FluidModel*)clientData)->getDifusion();
}


void TW_CALL setTempContact(const void *value, void *clientData)
{
	const Real val = *(const Real *)(value);
	((FluidModel*)clientData)->setTempContact(val);
}

void TW_CALL getTempContact(void *value, void *clientData)
{
	*(Real *)(value) = ((FluidModel*)clientData)->getTempContact();
}

void TW_CALL setContact(const void *value, void *clientData)
{
	const bool val = *(const bool *)(value);
	((FluidModel*)clientData)->setContact(val);
}

void TW_CALL getContact(void *value, void *clientData)
{
	*(bool *)(value) = ((FluidModel*)clientData)->getContact();
}

void TW_CALL setXPBD(const void *value, void *clientData)
{
	const bool val = *(const bool *)(value);
	((FluidModel*)clientData)->setXPBD(val);
}

void TW_CALL getXPBD(void *value, void *clientData)
{
	*(bool *)(value) = ((FluidModel*)clientData)->getXPBD();
}

void TW_CALL setDesWidth(const void *value, void *clientData)
{
	int val = *( int *)(value);
	MiniGL::desWidth = val;
	MiniGL::reziseDesWidth();
}

void TW_CALL getDesWidth(void *value, void *clientData)
{
	*(int*)(value) = MiniGL::desWidth;
}

void TW_CALL setDesHeight(const void *value, void *clientData)
{
	int val = *(int *)(value);
	MiniGL::desHeight = val;
	MiniGL::reziseDesHeight();
}

void TW_CALL getDesHeight(void *value, void *clientData)
{
	*(int*)(value) = MiniGL::desHeight;
}