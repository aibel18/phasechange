#include "Constraints.h"
#include "SimulationModel.h"
#include "PositionBasedDynamics/PositionBasedDynamics.h"
#include "PositionBasedDynamics/PositionBasedRigidBodyDynamics.h"
#include "TimeManager.h"
#include "Demos/Simulation/IDFactory.h"
#include "PositionBasedDynamics/PositionBasedElasticRods.h"
#define _USE_MATH_DEFINES
#include "math.h"

using namespace PBD;


int BallJoint::TYPE_ID = IDFactory::getId();
int BallOnLineJoint::TYPE_ID = IDFactory::getId();
int HingeJoint::TYPE_ID = IDFactory::getId();
int UniversalJoint::TYPE_ID = IDFactory::getId();
int RigidBodyParticleBallJoint::TYPE_ID = IDFactory::getId();
int DistanceConstraint::TYPE_ID = IDFactory::getId();
int DihedralConstraint::TYPE_ID = IDFactory::getId();
int IsometricBendingConstraint::TYPE_ID = IDFactory::getId();
int FEMTriangleConstraint::TYPE_ID = IDFactory::getId();
int StrainTriangleConstraint::TYPE_ID = IDFactory::getId();
int VolumeConstraint::TYPE_ID = IDFactory::getId();
int FEMTetConstraint::TYPE_ID = IDFactory::getId();
int StrainTetConstraint::TYPE_ID = IDFactory::getId();
int ShapeMatchingConstraint::TYPE_ID = IDFactory::getId();
int TargetAngleMotorHingeJoint::TYPE_ID = IDFactory::getId();
int TargetVelocityMotorHingeJoint::TYPE_ID = IDFactory::getId();
int SliderJoint::TYPE_ID = IDFactory::getId();
int TargetPositionMotorSliderJoint::TYPE_ID = IDFactory::getId();
int TargetVelocityMotorSliderJoint::TYPE_ID = IDFactory::getId();
int RigidBodyContactConstraint::TYPE_ID = IDFactory::getId();
int ParticleRigidBodyContactConstraint::TYPE_ID = IDFactory::getId();
int StretchShearConstraint::TYPE_ID = IDFactory::getId();
int BendTwistConstraint::TYPE_ID = IDFactory::getId();

//////////////////////////////////////////////////////////////////////////
// BallJoint
//////////////////////////////////////////////////////////////////////////
bool BallJoint::initConstraint(SimulationModel &model, const unsigned int rbIndex1, const unsigned int rbIndex2, const Vector3r &pos)
{
	m_bodies[0] = rbIndex1;
	m_bodies[1] = rbIndex2;
	SimulationModel::RigidBodyVector &rb = model.getRigidBodies();
	RigidBody &rb1 = *rb[m_bodies[0]];
	RigidBody &rb2 = *rb[m_bodies[1]];
	return PositionBasedRigidBodyDynamics::init_BallJoint(
		rb1.getPosition(),
		rb1.getRotation(),
		rb2.getPosition(),
		rb2.getRotation(),
		pos,
		m_jointInfo);
}

bool BallJoint::updateConstraint(SimulationModel &model)
{
	SimulationModel::RigidBodyVector &rb = model.getRigidBodies();
	RigidBody &rb1 = *rb[m_bodies[0]];
	RigidBody &rb2 = *rb[m_bodies[1]];
	return PositionBasedRigidBodyDynamics::update_BallJoint(
		rb1.getPosition(),
		rb1.getRotation(),
		rb2.getPosition(),
		rb2.getRotation(),
		m_jointInfo);
}

bool BallJoint::solvePositionConstraint(SimulationModel &model)
{
	SimulationModel::RigidBodyVector &rb = model.getRigidBodies();

	RigidBody &rb1 = *rb[m_bodies[0]];
	RigidBody &rb2 = *rb[m_bodies[1]];

	Vector3r corr_x1, corr_x2;
	Quaternionr corr_q1, corr_q2;
	const bool res = PositionBasedRigidBodyDynamics::solve_BallJoint(
		rb1.getInvMass(),
		rb1.getPosition(),
		rb1.getInertiaTensorInverseW(),
		rb1.getRotation(),
		rb2.getInvMass(),
		rb2.getPosition(),
		rb2.getInertiaTensorInverseW(),
		rb2.getRotation(),
		m_jointInfo,
		corr_x1,
		corr_q1,
		corr_x2,
		corr_q2);

	if (res)
	{
		if (rb1.getMass() != 0.0)
		{
			rb1.getPosition() += corr_x1;
			rb1.getRotation().coeffs() += corr_q1.coeffs();
			rb1.getRotation().normalize();
			rb1.rotationUpdated();
		}
		if (rb2.getMass() != 0.0)
		{
			rb2.getPosition() += corr_x2;
			rb2.getRotation().coeffs() += corr_q2.coeffs();
			rb2.getRotation().normalize();
			rb2.rotationUpdated();
		}
	}
	return res;
}


//////////////////////////////////////////////////////////////////////////
// BallOnLineJoint
//////////////////////////////////////////////////////////////////////////
bool BallOnLineJoint::initConstraint(SimulationModel &model, const unsigned int rbIndex1, const unsigned int rbIndex2, const Vector3r &pos, const Vector3r &dir)
{
	m_bodies[0] = rbIndex1;
	m_bodies[1] = rbIndex2;
	SimulationModel::RigidBodyVector &rb = model.getRigidBodies();
	RigidBody &rb1 = *rb[m_bodies[0]];
	RigidBody &rb2 = *rb[m_bodies[1]];
	return PositionBasedRigidBodyDynamics::init_BallOnLineJoint(
		rb1.getPosition(),
		rb1.getRotation(),
		rb2.getPosition(),
		rb2.getRotation(),
		pos, dir,
		m_jointInfo);
}

bool BallOnLineJoint::updateConstraint(SimulationModel &model)
{
	SimulationModel::RigidBodyVector &rb = model.getRigidBodies();
	RigidBody &rb1 = *rb[m_bodies[0]];
	RigidBody &rb2 = *rb[m_bodies[1]];
	return PositionBasedRigidBodyDynamics::update_BallOnLineJoint(
		rb1.getPosition(),
		rb1.getRotation(),
		rb2.getPosition(),
		rb2.getRotation(),
		m_jointInfo);
}

bool BallOnLineJoint::solvePositionConstraint(SimulationModel &model)
{
	SimulationModel::RigidBodyVector &rb = model.getRigidBodies();

	RigidBody &rb1 = *rb[m_bodies[0]];
	RigidBody &rb2 = *rb[m_bodies[1]];

	Vector3r corr_x1, corr_x2;
	Quaternionr corr_q1, corr_q2;
	const bool res = PositionBasedRigidBodyDynamics::solve_BallOnLineJoint(
		rb1.getInvMass(),
		rb1.getPosition(),
		rb1.getInertiaTensorInverseW(),
		rb1.getRotation(),
		rb2.getInvMass(),
		rb2.getPosition(),
		rb2.getInertiaTensorInverseW(),
		rb2.getRotation(),
		m_jointInfo,
		corr_x1,
		corr_q1,
		corr_x2,
		corr_q2);

	if (res)
	{
		if (rb1.getMass() != 0.0)
		{
			rb1.getPosition() += corr_x1;
			rb1.getRotation().coeffs() += corr_q1.coeffs();
			rb1.getRotation().normalize();
			rb1.rotationUpdated();
		}
		if (rb2.getMass() != 0.0)
		{
			rb2.getPosition() += corr_x2;
			rb2.getRotation().coeffs() += corr_q2.coeffs();
			rb2.getRotation().normalize();
			rb2.rotationUpdated();
		}
	}
	return res;
}


//////////////////////////////////////////////////////////////////////////
// HingeJoint
//////////////////////////////////////////////////////////////////////////
bool HingeJoint::initConstraint(SimulationModel &model, const unsigned int rbIndex1, const unsigned int rbIndex2, const Vector3r &pos, const Vector3r &axis)
{
	m_bodies[0] = rbIndex1;
	m_bodies[1] = rbIndex2;
	SimulationModel::RigidBodyVector &rb = model.getRigidBodies();
	RigidBody &rb1 = *rb[m_bodies[0]];
	RigidBody &rb2 = *rb[m_bodies[1]];
	return PositionBasedRigidBodyDynamics::init_HingeJoint(
		rb1.getPosition(),
		rb1.getRotation(),
		rb2.getPosition(),
		rb2.getRotation(),
		pos, axis,
		m_jointInfo);
}

bool HingeJoint::updateConstraint(SimulationModel &model)
{
	SimulationModel::RigidBodyVector &rb = model.getRigidBodies();
	RigidBody &rb1 = *rb[m_bodies[0]];
	RigidBody &rb2 = *rb[m_bodies[1]];
	return PositionBasedRigidBodyDynamics::update_HingeJoint(
		rb1.getPosition(),
		rb1.getRotation(),
		rb2.getPosition(),
		rb2.getRotation(),
		m_jointInfo);
}

bool HingeJoint::solvePositionConstraint(SimulationModel &model)
{
	SimulationModel::RigidBodyVector &rb = model.getRigidBodies();

	RigidBody &rb1 = *rb[m_bodies[0]];
	RigidBody &rb2 = *rb[m_bodies[1]];

	Vector3r corr_x1, corr_x2;
	Quaternionr corr_q1, corr_q2;
	const bool res = PositionBasedRigidBodyDynamics::solve_HingeJoint(
		rb1.getInvMass(),
		rb1.getPosition(),
		rb1.getInertiaTensorInverseW(),
		rb1.getRotation(),
		rb2.getInvMass(),
		rb2.getPosition(),
		rb2.getInertiaTensorInverseW(),
		rb2.getRotation(),
		m_jointInfo,
		corr_x1,
		corr_q1,
		corr_x2,
		corr_q2);

	if (res)
	{
		if (rb1.getMass() != 0.0)
		{
			rb1.getPosition() += corr_x1;
			rb1.getRotation().coeffs() += corr_q1.coeffs();
			rb1.getRotation().normalize();
			rb1.rotationUpdated();
		}
		if (rb2.getMass() != 0.0)
		{
			rb2.getPosition() += corr_x2;
			rb2.getRotation().coeffs() += corr_q2.coeffs();
			rb2.getRotation().normalize();
			rb2.rotationUpdated();
		}
	}
	return res;
}


//////////////////////////////////////////////////////////////////////////
// UniversalJoint
//////////////////////////////////////////////////////////////////////////
bool UniversalJoint::initConstraint(SimulationModel &model, const unsigned int rbIndex1, const unsigned int rbIndex2, const Vector3r &pos, const Vector3r &axis1, const Vector3r &axis2)
{
	m_bodies[0] = rbIndex1;
	m_bodies[1] = rbIndex2;
	SimulationModel::RigidBodyVector &rb = model.getRigidBodies();
	RigidBody &rb1 = *rb[m_bodies[0]];
	RigidBody &rb2 = *rb[m_bodies[1]];
	return PositionBasedRigidBodyDynamics::init_UniversalJoint(
		rb1.getPosition(),
		rb1.getRotation(),
		rb2.getPosition(),
		rb2.getRotation(),
		pos,
		axis1,
		axis2,
		m_jointInfo);
}

bool UniversalJoint::updateConstraint(SimulationModel &model)
{
	SimulationModel::RigidBodyVector &rb = model.getRigidBodies();
	RigidBody &rb1 = *rb[m_bodies[0]];
	RigidBody &rb2 = *rb[m_bodies[1]];
	return PositionBasedRigidBodyDynamics::update_UniversalJoint(
		rb1.getPosition(),
		rb1.getRotation(),
		rb2.getPosition(),
		rb2.getRotation(),
		m_jointInfo);
}

bool UniversalJoint::solvePositionConstraint(SimulationModel &model)
{
	SimulationModel::RigidBodyVector &rb = model.getRigidBodies();

	RigidBody &rb1 = *rb[m_bodies[0]];
	RigidBody &rb2 = *rb[m_bodies[1]];

	Vector3r corr_x1, corr_x2;
	Quaternionr corr_q1, corr_q2;
	const bool res = PositionBasedRigidBodyDynamics::solve_UniversalJoint(
		rb1.getInvMass(),
		rb1.getPosition(),
		rb1.getInertiaTensorInverseW(),
		rb1.getRotation(),
		rb2.getInvMass(),
		rb2.getPosition(),
		rb2.getInertiaTensorInverseW(),
		rb2.getRotation(),
		m_jointInfo,
		corr_x1,
		corr_q1,
		corr_x2,
		corr_q2);

	if (res)
	{
		if (rb1.getMass() != 0.0)
		{
			rb1.getPosition() += corr_x1;
			rb1.getRotation().coeffs() += corr_q1.coeffs();
			rb1.getRotation().normalize();
			rb1.rotationUpdated();
		}
		if (rb2.getMass() != 0.0)
		{
			rb2.getPosition() += corr_x2;
			rb2.getRotation().coeffs() += corr_q2.coeffs();
			rb2.getRotation().normalize();
			rb2.rotationUpdated();
		}
	}
	return res;
}


//////////////////////////////////////////////////////////////////////////
// SliderJoint
//////////////////////////////////////////////////////////////////////////
bool SliderJoint::initConstraint(SimulationModel &model, const unsigned int rbIndex1, const unsigned int rbIndex2, const Vector3r &pos, const Vector3r &axis)
{
	m_bodies[0] = rbIndex1;
	m_bodies[1] = rbIndex2;
	SimulationModel::RigidBodyVector &rb = model.getRigidBodies();
	RigidBody &rb1 = *rb[m_bodies[0]];
	RigidBody &rb2 = *rb[m_bodies[1]];
	return PositionBasedRigidBodyDynamics::init_SliderJoint(
		rb1.getPosition(),
		rb1.getRotation(),
		rb2.getPosition(),
		rb2.getRotation(),
		pos, axis,
		m_jointInfo);
}

bool SliderJoint::updateConstraint(SimulationModel &model)
{
	SimulationModel::RigidBodyVector &rb = model.getRigidBodies();
	RigidBody &rb1 = *rb[m_bodies[0]];
	RigidBody &rb2 = *rb[m_bodies[1]];
	return PositionBasedRigidBodyDynamics::update_SliderJoint(
		rb1.getPosition(),
		rb1.getRotation(),
		rb2.getPosition(),
		rb2.getRotation(),
		m_jointInfo);
}

bool SliderJoint::solvePositionConstraint(SimulationModel &model)
{
	SimulationModel::RigidBodyVector &rb = model.getRigidBodies();

	RigidBody &rb1 = *rb[m_bodies[0]];
	RigidBody &rb2 = *rb[m_bodies[1]];

	Vector3r corr_x1, corr_x2;
	Quaternionr corr_q1, corr_q2;
	const bool res = PositionBasedRigidBodyDynamics::solve_SliderJoint(
		rb1.getInvMass(),
		rb1.getPosition(),
		rb1.getInertiaTensorInverseW(),
		rb1.getRotation(),
		rb2.getInvMass(),
		rb2.getPosition(),
		rb2.getInertiaTensorInverseW(),
		rb2.getRotation(),
		m_jointInfo,
		corr_x1,
		corr_q1,
		corr_x2,
		corr_q2);

	if (res)
	{
		if (rb1.getMass() != 0.0)
		{
			rb1.getPosition() += corr_x1;
			rb1.getRotation().coeffs() += corr_q1.coeffs();
			rb1.getRotation().normalize();
			rb1.rotationUpdated();
		}
		if (rb2.getMass() != 0.0)
		{
			rb2.getPosition() += corr_x2;
			rb2.getRotation().coeffs() += corr_q2.coeffs();
			rb2.getRotation().normalize();
			rb2.rotationUpdated();
		}
	}
	return res;
}


//////////////////////////////////////////////////////////////////////////
// TargetPositionMotorSliderJoint
//////////////////////////////////////////////////////////////////////////
bool TargetPositionMotorSliderJoint::initConstraint(SimulationModel &model, const unsigned int rbIndex1, const unsigned int rbIndex2, const Vector3r &pos, const Vector3r &axis)
{
	m_bodies[0] = rbIndex1;
	m_bodies[1] = rbIndex2;
	SimulationModel::RigidBodyVector &rb = model.getRigidBodies();
	RigidBody &rb1 = *rb[m_bodies[0]];
	RigidBody &rb2 = *rb[m_bodies[1]];
	return PositionBasedRigidBodyDynamics::init_TargetPositionMotorSliderJoint(
		rb1.getPosition(),
		rb1.getRotation(),
		rb2.getPosition(),
		rb2.getRotation(),
		pos, axis,
		m_jointInfo);
}

bool TargetPositionMotorSliderJoint::updateConstraint(SimulationModel &model)
{
	SimulationModel::RigidBodyVector &rb = model.getRigidBodies();
	RigidBody &rb1 = *rb[m_bodies[0]];
	RigidBody &rb2 = *rb[m_bodies[1]];
	return PositionBasedRigidBodyDynamics::update_TargetPositionMotorSliderJoint(
		rb1.getPosition(),
		rb1.getRotation(),
		rb2.getPosition(),
		rb2.getRotation(),
		m_jointInfo);
}

bool TargetPositionMotorSliderJoint::solvePositionConstraint(SimulationModel &model)
{
	SimulationModel::RigidBodyVector &rb = model.getRigidBodies();

	RigidBody &rb1 = *rb[m_bodies[0]];
	RigidBody &rb2 = *rb[m_bodies[1]];

	Vector3r corr_x1, corr_x2;
	Quaternionr corr_q1, corr_q2;
	const bool res = PositionBasedRigidBodyDynamics::solve_TargetPositionMotorSliderJoint(
		rb1.getInvMass(),
		rb1.getPosition(),
		rb1.getInertiaTensorInverseW(),
		rb1.getRotation(),
		rb2.getInvMass(),
		rb2.getPosition(),
		rb2.getInertiaTensorInverseW(),
		rb2.getRotation(),
		m_target,
		m_jointInfo,
		corr_x1,
		corr_q1,
		corr_x2,
		corr_q2);

	if (res)
	{
		if (rb1.getMass() != 0.0)
		{
			rb1.getPosition() += corr_x1;
			rb1.getRotation().coeffs() += corr_q1.coeffs();
			rb1.getRotation().normalize();
			rb1.rotationUpdated();
		}
		if (rb2.getMass() != 0.0)
		{
			rb2.getPosition() += corr_x2;
			rb2.getRotation().coeffs() += corr_q2.coeffs();
			rb2.getRotation().normalize();
			rb2.rotationUpdated();
		}
	}
	return res;
}



//////////////////////////////////////////////////////////////////////////
// TargetVelocityMotorSliderJoint
//////////////////////////////////////////////////////////////////////////
bool TargetVelocityMotorSliderJoint::initConstraint(SimulationModel &model, const unsigned int rbIndex1, const unsigned int rbIndex2, const Vector3r &pos, const Vector3r &axis)
{
	m_bodies[0] = rbIndex1;
	m_bodies[1] = rbIndex2;
	SimulationModel::RigidBodyVector &rb = model.getRigidBodies();
	RigidBody &rb1 = *rb[m_bodies[0]];
	RigidBody &rb2 = *rb[m_bodies[1]];
	return PositionBasedRigidBodyDynamics::init_TargetVelocityMotorSliderJoint(
		rb1.getPosition(),
		rb1.getRotation(),
		rb2.getPosition(),
		rb2.getRotation(),
		pos, axis,
		m_jointInfo);
}

bool TargetVelocityMotorSliderJoint::updateConstraint(SimulationModel &model)
{
	SimulationModel::RigidBodyVector &rb = model.getRigidBodies();
	RigidBody &rb1 = *rb[m_bodies[0]];
	RigidBody &rb2 = *rb[m_bodies[1]];
	return PositionBasedRigidBodyDynamics::update_TargetVelocityMotorSliderJoint(
		rb1.getPosition(),
		rb1.getRotation(),
		rb2.getPosition(),
		rb2.getRotation(),
		m_jointInfo);
}

bool TargetVelocityMotorSliderJoint::solvePositionConstraint(SimulationModel &model)
{
	SimulationModel::RigidBodyVector &rb = model.getRigidBodies();

	RigidBody &rb1 = *rb[m_bodies[0]];
	RigidBody &rb2 = *rb[m_bodies[1]];

	Vector3r corr_x1, corr_x2;
	Quaternionr corr_q1, corr_q2;
	const bool res = PositionBasedRigidBodyDynamics::solve_TargetVelocityMotorSliderJoint(
		rb1.getInvMass(),
		rb1.getPosition(),
		rb1.getInertiaTensorInverseW(),
		rb1.getRotation(),
		rb2.getInvMass(),
		rb2.getPosition(),
		rb2.getInertiaTensorInverseW(),
		rb2.getRotation(),
		m_jointInfo,
		corr_x1,
		corr_q1,
		corr_x2,
		corr_q2);

	if (res)
	{
		if (rb1.getMass() != 0.0)
		{
			rb1.getPosition() += corr_x1;
			rb1.getRotation().coeffs() += corr_q1.coeffs();
			rb1.getRotation().normalize();
			rb1.rotationUpdated();
		}
		if (rb2.getMass() != 0.0)
		{
			rb2.getPosition() += corr_x2;
			rb2.getRotation().coeffs() += corr_q2.coeffs();
			rb2.getRotation().normalize();
			rb2.rotationUpdated();
		}
	}
	return res;
}


bool TargetVelocityMotorSliderJoint::solveVelocityConstraint(SimulationModel &model)
{
	SimulationModel::RigidBodyVector &rb = model.getRigidBodies();

	RigidBody &rb1 = *rb[m_bodies[0]];
	RigidBody &rb2 = *rb[m_bodies[1]];

	Vector3r corr_v1, corr_v2;
	Vector3r corr_omega1, corr_omega2;
	const bool res = PositionBasedRigidBodyDynamics::velocitySolve_TargetVelocityMotorSliderJoint(
		rb1.getInvMass(),
		rb1.getPosition(),
		rb1.getVelocity(),
		rb1.getInertiaTensorInverseW(),
		rb1.getAngularVelocity(),
		rb2.getInvMass(),
		rb2.getPosition(),
		rb2.getVelocity(),
		rb2.getInertiaTensorInverseW(),
		rb2.getAngularVelocity(),
		m_target,
		m_jointInfo,
		corr_v1,
		corr_omega1,
		corr_v2,
		corr_omega2);

	if (res)
	{
		if (rb1.getMass() != 0.0)
		{
			rb1.getVelocity() += corr_v1;
			rb1.getAngularVelocity() += corr_omega1;
		}
		if (rb2.getMass() != 0.0)
		{
			rb2.getVelocity() += corr_v2;
			rb2.getAngularVelocity() += corr_omega2;
		}
	}
	return res;
}

//////////////////////////////////////////////////////////////////////////
// TargetAngleMotorHingeJoint
//////////////////////////////////////////////////////////////////////////
bool TargetAngleMotorHingeJoint::initConstraint(SimulationModel &model, const unsigned int rbIndex1, const unsigned int rbIndex2, const Vector3r &pos, const Vector3r &axis)
{
	m_bodies[0] = rbIndex1;
	m_bodies[1] = rbIndex2;
	SimulationModel::RigidBodyVector &rb = model.getRigidBodies();
	RigidBody &rb1 = *rb[m_bodies[0]];
	RigidBody &rb2 = *rb[m_bodies[1]];
	return PositionBasedRigidBodyDynamics::init_TargetAngleMotorHingeJoint(
		rb1.getPosition(),
		rb1.getRotation(),
		rb2.getPosition(),
		rb2.getRotation(),
		pos, axis,
		m_jointInfo);
}

bool TargetAngleMotorHingeJoint::updateConstraint(SimulationModel &model)
{
	SimulationModel::RigidBodyVector &rb = model.getRigidBodies();
	RigidBody &rb1 = *rb[m_bodies[0]];
	RigidBody &rb2 = *rb[m_bodies[1]];
	return PositionBasedRigidBodyDynamics::update_TargetAngleMotorHingeJoint(
		rb1.getPosition(),
		rb1.getRotation(),
		rb2.getPosition(),
		rb2.getRotation(),
		m_jointInfo);
}

bool TargetAngleMotorHingeJoint::solvePositionConstraint(SimulationModel &model)
{
	SimulationModel::RigidBodyVector &rb = model.getRigidBodies();

	RigidBody &rb1 = *rb[m_bodies[0]];
	RigidBody &rb2 = *rb[m_bodies[1]];

	Vector3r corr_x1, corr_x2;
	Quaternionr corr_q1, corr_q2;
	const bool res = PositionBasedRigidBodyDynamics::solve_TargetAngleMotorHingeJoint(
		rb1.getInvMass(),
		rb1.getPosition(),
		rb1.getInertiaTensorInverseW(),
		rb1.getRotation(),
		rb2.getInvMass(),
		rb2.getPosition(),
		rb2.getInertiaTensorInverseW(),
		rb2.getRotation(),
		m_target,
		m_jointInfo,
		corr_x1,
		corr_q1,
		corr_x2,
		corr_q2);

	if (res)
	{
		if (rb1.getMass() != 0.0)
		{
			rb1.getPosition() += corr_x1;
			rb1.getRotation().coeffs() += corr_q1.coeffs();
			rb1.getRotation().normalize();
			rb1.rotationUpdated();
		}
		if (rb2.getMass() != 0.0)
		{
			rb2.getPosition() += corr_x2;
			rb2.getRotation().coeffs() += corr_q2.coeffs();
			rb2.getRotation().normalize();
			rb2.rotationUpdated();
		}
	}
	return res;
}

//////////////////////////////////////////////////////////////////////////
// TargetVelocityMotorHingeJoint
//////////////////////////////////////////////////////////////////////////
bool TargetVelocityMotorHingeJoint::initConstraint(SimulationModel &model, const unsigned int rbIndex1, const unsigned int rbIndex2, const Vector3r &pos, const Vector3r &axis)
{
	m_bodies[0] = rbIndex1;
	m_bodies[1] = rbIndex2;
	SimulationModel::RigidBodyVector &rb = model.getRigidBodies();
	RigidBody &rb1 = *rb[m_bodies[0]];
	RigidBody &rb2 = *rb[m_bodies[1]];
	return PositionBasedRigidBodyDynamics::init_TargetVelocityMotorHingeJoint(
		rb1.getPosition(),
		rb1.getRotation(),
		rb2.getPosition(),
		rb2.getRotation(),
		pos, axis,
		m_jointInfo);
}

bool TargetVelocityMotorHingeJoint::updateConstraint(SimulationModel &model)
{
	SimulationModel::RigidBodyVector &rb = model.getRigidBodies();
	RigidBody &rb1 = *rb[m_bodies[0]];
	RigidBody &rb2 = *rb[m_bodies[1]];
	return PositionBasedRigidBodyDynamics::update_TargetVelocityMotorHingeJoint(
		rb1.getPosition(),
		rb1.getRotation(),
		rb2.getPosition(),
		rb2.getRotation(),
		m_jointInfo);
}

bool TargetVelocityMotorHingeJoint::solvePositionConstraint(SimulationModel &model)
{
	SimulationModel::RigidBodyVector &rb = model.getRigidBodies();

	RigidBody &rb1 = *rb[m_bodies[0]];
	RigidBody &rb2 = *rb[m_bodies[1]];

	Vector3r corr_x1, corr_x2;
	Quaternionr corr_q1, corr_q2;
	const bool res = PositionBasedRigidBodyDynamics::solve_TargetVelocityMotorHingeJoint(
		rb1.getInvMass(),
		rb1.getPosition(),
		rb1.getInertiaTensorInverseW(),
		rb1.getRotation(),
		rb2.getInvMass(),
		rb2.getPosition(),
		rb2.getInertiaTensorInverseW(),
		rb2.getRotation(),
		m_jointInfo,
		corr_x1,
		corr_q1,
		corr_x2,
		corr_q2);

	if (res)
	{
		if (rb1.getMass() != 0.0)
		{
			rb1.getPosition() += corr_x1;
			rb1.getRotation().coeffs() += corr_q1.coeffs();
			rb1.getRotation().normalize();
			rb1.rotationUpdated();
		}
		if (rb2.getMass() != 0.0)
		{
			rb2.getPosition() += corr_x2;
			rb2.getRotation().coeffs() += corr_q2.coeffs();
			rb2.getRotation().normalize();
			rb2.rotationUpdated();
		}
	}
	return res;
}

bool TargetVelocityMotorHingeJoint::solveVelocityConstraint(SimulationModel &model)
{
	SimulationModel::RigidBodyVector &rb = model.getRigidBodies();

	RigidBody &rb1 = *rb[m_bodies[0]];
	RigidBody &rb2 = *rb[m_bodies[1]];

	Vector3r corr_v1, corr_v2;
	Vector3r corr_omega1, corr_omega2;
	const bool res = PositionBasedRigidBodyDynamics::velocitySolve_TargetVelocityMotorHingeJoint(
		rb1.getInvMass(),
		rb1.getPosition(),
		rb1.getVelocity(),
		rb1.getInertiaTensorInverseW(),
		rb1.getAngularVelocity(),
		rb2.getInvMass(),
		rb2.getPosition(),
		rb2.getVelocity(),
		rb2.getInertiaTensorInverseW(),
		rb2.getAngularVelocity(),
		m_target, 
		m_jointInfo,
		corr_v1,
		corr_omega1,
		corr_v2,
		corr_omega2);

	if (res)
	{
		if (rb1.getMass() != 0.0)
		{
			rb1.getVelocity() += corr_v1;
			rb1.getAngularVelocity() += corr_omega1;
		}
		if (rb2.getMass() != 0.0)
		{
			rb2.getVelocity() += corr_v2;
			rb2.getAngularVelocity() += corr_omega2;
		}
	}
	return res;
}

//////////////////////////////////////////////////////////////////////////
// RigidBodyParticleBallJoint
//////////////////////////////////////////////////////////////////////////
bool RigidBodyParticleBallJoint::initConstraint(SimulationModel &model, const unsigned int rbIndex, const unsigned int particleIndex)
{
	m_bodies[0] = rbIndex;
	m_bodies[1] = particleIndex;
	SimulationModel::RigidBodyVector &rbs = model.getRigidBodies();
	ParticleData &pd = model.getParticles();
	RigidBody &rb = *rbs[m_bodies[0]];
	return PositionBasedRigidBodyDynamics::init_RigidBodyParticleBallJoint(
		rb.getPosition(),
		rb.getRotation(),
		pd.getPosition(particleIndex),
		m_jointInfo);
}

bool RigidBodyParticleBallJoint::updateConstraint(SimulationModel &model)
{
	SimulationModel::RigidBodyVector &rb = model.getRigidBodies();
	ParticleData &pd = model.getParticles();
	RigidBody &rb1 = *rb[m_bodies[0]];
	return PositionBasedRigidBodyDynamics::update_RigidBodyParticleBallJoint(
		rb1.getPosition(),
		rb1.getRotation(),
		pd.getPosition(m_bodies[1]),
		m_jointInfo);
}

bool RigidBodyParticleBallJoint::solvePositionConstraint(SimulationModel &model)
{
	SimulationModel::RigidBodyVector &rb = model.getRigidBodies();
	ParticleData &pd = model.getParticles();

	RigidBody &rb1 = *rb[m_bodies[0]];

	Vector3r corr_x1, corr_x2;
	Quaternionr corr_q1;
	const bool res = PositionBasedRigidBodyDynamics::solve_RigidBodyParticleBallJoint(
		rb1.getInvMass(),
		rb1.getPosition(),
		rb1.getInertiaTensorInverseW(),
		rb1.getRotation(),
		pd.getInvMass(m_bodies[1]),
		pd.getPosition(m_bodies[1]),
		m_jointInfo,
		corr_x1,
		corr_q1,
		corr_x2);

	if (res)
	{
		if (rb1.getMass() != 0.0)
		{
			rb1.getPosition() += corr_x1;
			rb1.getRotation().coeffs() += corr_q1.coeffs();
			rb1.getRotation().normalize();
			rb1.rotationUpdated();
		}
		if (pd.getMass(m_bodies[1]) != 0.0)
		{
			pd.getPosition(m_bodies[1]) += corr_x2;
		}
	}
	return res;
}

//////////////////////////////////////////////////////////////////////////
// DistanceConstraint
//////////////////////////////////////////////////////////////////////////
bool DistanceConstraint::initConstraint(SimulationModel &model, const unsigned int particle1, const unsigned int particle2)
{
	m_bodies[0] = particle1;
	m_bodies[1] = particle2;
	ParticleData &pd = model.getParticles();

	const Vector3r &x1_0 = pd.getPosition0(particle1);
	const Vector3r &x2_0 = pd.getPosition0(particle2);

	m_restLength = (x2_0 - x1_0).norm();

	return true;
}

bool DistanceConstraint::solvePositionConstraint(SimulationModel &model)
{
	ParticleData &pd = model.getParticles();

	const unsigned i1 = m_bodies[0];
	const unsigned i2 = m_bodies[1];

	Vector3r &x1 = pd.getPosition(i1);
	Vector3r &x2 = pd.getPosition(i2);
	const Real invMass1 = pd.getInvMass(i1);
	const Real invMass2 = pd.getInvMass(i2);

	Vector3r corr1, corr2;
	const bool res = PositionBasedDynamics::solve_DistanceConstraint(
		x1, invMass1, x2, invMass2,
		m_restLength, model.getClothStiffness(), model.getClothStiffness(), corr1, corr2);

	if (res)
	{
		if (invMass1 != 0.0)
			x1 += corr1;
		if (invMass2 != 0.0)
			x2 += corr2;
	}
	return res;
}

//////////////////////////////////////////////////////////////////////////
// NewDistanceConstraint
//////////////////////////////////////////////////////////////////////////
/*
bool newDistanceConstraint::initConstraint(SimulationModel &model, const unsigned int particle1, const unsigned int particle2)
{
	return true;
}
bool newDistanceConstraint::solvePositionConstraint(SimulationModel &model)
{
	return true;
}*/

bool newDistanceConstraint::evalConstraint(FluidModel &model)
{
	const unsigned i1 = m_bodies[0];
	const unsigned i2 = m_bodies[1];

	ParticleData &pd = model.getParticles();

	//if (pd.getState(i1) != 0 || pd.getState(i2) != 0) {
	if (pd.getState(i1) != 0 && pd.getState(i2) != 0   )  {
		return false;
	}
	return true;
	/*
	if ( (pd.getState(i1) == 0 && pd.getState(i2) == 0) || ( pd.getState(i1) == -1 && pd.getState(i2) == -1)   )  {
		return true;
	}
	return false;*/
}

Real PBD::sigmoide(Real val, Real lf, Real begin, Real offset) {
	Real a = 1.0 + (12.0 / lf);
	Real b = (lf / 2.0) - val;
	Real denominator = 1 + pow(a, b);
	Real total = (-begin + offset) / denominator + begin;
	return total;
}
Real PBD::sigmoide2(Real val,Real li, Real lf, Real begin, Real offset) {
	Real a = 1.0 + (12.0 / (lf+li));
	Real b = (lf + li) / 2.0 - val;
	Real denominator = 1 + pow(a, b);
	Real total = (-begin + offset) / denominator + begin;
	return total;
}
Real PBD::logarithmic(Real val, Real lf, Real begin, Real offset, Real p) {

	// d = (a - b) /ln( (c+p)/p )
	Real d = (offset - begin) / log( (lf + p) / p );
	// f(val) = d ln ( (c+p) /( c+p-x)) + b
	Real result = d * log( (lf + p ) / (lf + p - val)) + begin;
	return result <= offset ? result : offset;
}
Real PBD::logarithmic2(Real val, Real lf, Real begin, Real offset, Real p) {

	// d = (a - b) /ln( (c+p)/p )
	Real d = (offset - begin) / log((lf + p) / p);
	// f(val) = d ln ( (c+p) /( c+p-x)) + b
	Real result = d * log((lf + p) / (lf + p - val)) + begin;
	return result >= offset ? result : offset;
}

bool newDistanceConstraint::initConstraint(FluidModel &model, const unsigned int particle1, const unsigned int particle2)
{
	m_bodies[0] = particle1;
	m_bodies[1] = particle2;
	ParticleData &pd = model.getParticles();

	//m_restTemperature = -5;// (pd.getTemp0(particle1) + pd.getTemp0(particle2)) / 2;
	//m_restTemperature = (pd.getTemp(particle1) + pd.getTemp(particle2)) / 2;

	
	m_restTemperature = -10;// pd.getTemp0(particle1) <= pd.getTemp0(particle2) ? pd.getTemp0(particle1) : pd.getTemp0(particle2);

	const Vector3r &x1_0 = pd.getPosition(particle1);
	const Vector3r &x2_0 = pd.getPosition(particle2);

	m_restLength = (x2_0 - x1_0).norm();

	return true;
}

bool newDistanceConstraint::updateConstraint(FluidModel &model)
{
	const unsigned i1 = m_bodies[0];
	const unsigned i2 = m_bodies[1];

	ParticleData &pd = model.getParticles();

	if (pd.getState(i1) != 0 || pd.getState(i2) != 0) {
		return false;
	}

	const Vector3r &x1_0 = pd.getPosition0(i1);
	const Vector3r &x2_0 = pd.getPosition0(i2);

	const Real t1 = pd.getTemp(i1);
	const Real t2 = pd.getTemp(i2);

	const Real deltaT = (t1 + t2) / 2 - m_restTemperature;

	//m_restLength = (x2_0 - x1_0).norm();

	Real rest = model.getDilation()*m_restLength*deltaT;

	//m_restLength = m_restLength + rest;
	//m_restTemperature = (t1 + t2);

	//model.setClothStiffness(rest*model.getClothStiffness());
	//model.setSolidStiffness(rest*model.getSolidStiffness());

	return true;
}

bool newDistanceConstraint::solvePositionConstraint(FluidModel &model)
{
	ParticleData &pd = model.getParticles();

	const unsigned i1 = m_bodies[0];
	const unsigned i2 = m_bodies[1];

	const Real t1 = pd.getTemp(i1);
	const Real t2 = pd.getTemp(i2);

	Real deltaT = (t1 + t2) / 2;// -m_restTemperature;

	//Real stiffness = model.getClothStiffness();//model.getClothStiffness() / ( 1 + pow(2, deltaT - model.getMelting()/2 ) );
	Real stiffness = sigmoide2(deltaT, -m_restTemperature, model.getMelting(), model.getClothStiffness(), 0.001);
	//stiffness = std::max(stiffness, 0.01);

	Real rest = model.getDilation()*m_restLength*((t1 + t2) / 2 - m_restTemperature);
	Real length = m_restLength + rest;

	Vector3r &x1 = pd.getPosition(i1);
	Vector3r &x2 = pd.getPosition(i2);
	const Real invMass1 = pd.getInvMass(i1);
	const Real invMass2 = pd.getInvMass(i2);

	Vector3r corr1, corr2;
	const bool res = 
		PositionBasedDynamics::solve_DistanceConstraintX(
		x1, invMass1, x2, invMass2,
		length, stiffness, stiffness, corr1, corr2, lambda, TimeManager::getCurrent()->getTimeStepSize());


	if (res)
	{
		if (invMass1 != 0.0)
			x1 += corr1;
		if (invMass2 != 0.0)
			x2 += corr2;
	}
	return res;
}



//////////////////////////////////////////////////////////////////////////
// DihedralConstraint
//////////////////////////////////////////////////////////////////////////

bool DihedralConstraint::initConstraint(SimulationModel &model, const unsigned int particle1, const unsigned int particle2,
										const unsigned int particle3, const unsigned int particle4)
{
	m_bodies[0] = particle1;
	m_bodies[1] = particle2;
	m_bodies[2] = particle3;
	m_bodies[3] = particle4;
	ParticleData &pd = model.getParticles();

	const Vector3r &p0 = pd.getPosition0(particle1);
	const Vector3r &p1 = pd.getPosition0(particle2);
	const Vector3r &p2 = pd.getPosition0(particle3);
	const Vector3r &p3 = pd.getPosition0(particle4);

	Vector3r e = p3 - p2;
	Real  elen = e.norm();
	if (elen < 1e-6)
		return false;

	Real invElen = 1.0 / elen;

	Vector3r n1 = (p2 - p0).cross(p3 - p0); n1 /= n1.squaredNorm();
	Vector3r n2 = (p3 - p1).cross(p2 - p1); n2 /= n2.squaredNorm();

	n1.normalize();
	n2.normalize();
	Real dot = n1.dot(n2);

	if (dot < -1.0) dot = -1.0;
	if (dot > 1.0) dot = 1.0;

	m_restAngle = acos(dot);

	return true;
}

bool DihedralConstraint::solvePositionConstraint(SimulationModel &model)
{
	ParticleData &pd = model.getParticles();

	const unsigned i1 = m_bodies[0];
	const unsigned i2 = m_bodies[1];
	const unsigned i3 = m_bodies[2];
	const unsigned i4 = m_bodies[3];

	Vector3r &x1 = pd.getPosition(i1);
	Vector3r &x2 = pd.getPosition(i2);
	Vector3r &x3 = pd.getPosition(i3);
	Vector3r &x4 = pd.getPosition(i4);

	const Real invMass1 = pd.getInvMass(i1);
	const Real invMass2 = pd.getInvMass(i2);
	const Real invMass3 = pd.getInvMass(i3);
	const Real invMass4 = pd.getInvMass(i4);

	Vector3r corr1, corr2, corr3, corr4;
	const bool res = PositionBasedDynamics::solve_DihedralConstraint(
		x1, invMass1, x2, invMass2, x3, invMass3, x4, invMass4,
		m_restAngle,
		model.getClothBendingStiffness(),
		corr1, corr2, corr3, corr4);

	if (res)
	{
		if (invMass1 != 0.0)
			x1 += corr1;
		if (invMass2 != 0.0)
			x2 += corr2;
		if (invMass3 != 0.0)
			x3 += corr3;
		if (invMass4 != 0.0)
			x4 += corr4;
	}
	return res;
}


//////////////////////////////////////////////////////////////////////////
// IsometricBendingConstraint
//////////////////////////////////////////////////////////////////////////
bool IsometricBendingConstraint::initConstraint(SimulationModel &model, const unsigned int particle1, const unsigned int particle2,
												const unsigned int particle3, const unsigned int particle4)
{
	m_bodies[0] = particle1;
	m_bodies[1] = particle2;
	m_bodies[2] = particle3;
	m_bodies[3] = particle4;

	ParticleData &pd = model.getParticles();

	const Vector3r &x1 = pd.getPosition0(particle1);
	const Vector3r &x2 = pd.getPosition0(particle2);
	const Vector3r &x3 = pd.getPosition0(particle3);
	const Vector3r &x4 = pd.getPosition0(particle4);

	return PositionBasedDynamics::init_IsometricBendingConstraint(x1, x2, x3, x4, m_Q);
}

bool IsometricBendingConstraint::solvePositionConstraint(SimulationModel &model)
{
	ParticleData &pd = model.getParticles();

	const unsigned i1 = m_bodies[0];
	const unsigned i2 = m_bodies[1];
	const unsigned i3 = m_bodies[2];
	const unsigned i4 = m_bodies[3];

	Vector3r &x1 = pd.getPosition(i1);
	Vector3r &x2 = pd.getPosition(i2);
	Vector3r &x3 = pd.getPosition(i3);
	Vector3r &x4 = pd.getPosition(i4);

	const Real invMass1 = pd.getInvMass(i1);
	const Real invMass2 = pd.getInvMass(i2);
	const Real invMass3 = pd.getInvMass(i3);
	const Real invMass4 = pd.getInvMass(i4);

	Vector3r corr1, corr2, corr3, corr4;
	const bool res = PositionBasedDynamics::solve_IsometricBendingConstraint(
		x1, invMass1, x2, invMass2, x3, invMass3, x4, invMass4,
		m_Q,
		model.getClothBendingStiffness(),
		corr1, corr2, corr3, corr4);

	if (res)
	{
		if (invMass1 != 0.0)
			x1 += corr1;
		if (invMass2 != 0.0)
			x2 += corr2;
		if (invMass3 != 0.0)
			x3 += corr3;
		if (invMass4 != 0.0)
			x4 += corr4;
	}
	return res;
}

//////////////////////////////////////////////////////////////////////////
// FEMTriangleConstraint
//////////////////////////////////////////////////////////////////////////
bool FEMTriangleConstraint::initConstraint(SimulationModel &model, const unsigned int particle1, const unsigned int particle2,
	const unsigned int particle3)
{
	m_bodies[0] = particle1;
	m_bodies[1] = particle2;
	m_bodies[2] = particle3;
	
	ParticleData &pd = model.getParticles();

	Vector3r &x1 = pd.getPosition0(particle1);
	Vector3r &x2 = pd.getPosition0(particle2);
	Vector3r &x3 = pd.getPosition0(particle3);

	return PositionBasedDynamics::init_FEMTriangleConstraint(x1, x2, x3, m_area, m_invRestMat);
}

bool FEMTriangleConstraint::solvePositionConstraint(SimulationModel &model)
{
	ParticleData &pd = model.getParticles();

	const unsigned i1 = m_bodies[0];
	const unsigned i2 = m_bodies[1];
	const unsigned i3 = m_bodies[2];

	Vector3r &x1 = pd.getPosition(i1);
	Vector3r &x2 = pd.getPosition(i2);
	Vector3r &x3 = pd.getPosition(i3);

	const Real invMass1 = pd.getInvMass(i1);
	const Real invMass2 = pd.getInvMass(i2);
	const Real invMass3 = pd.getInvMass(i3);
	
	Vector3r corr1, corr2, corr3;
	const bool res = PositionBasedDynamics::solve_FEMTriangleConstraint(
		x1, invMass1,
		x2, invMass2,
		x3, invMass3,
		m_area,
		m_invRestMat,
		model.getClothXXStiffness(),
		model.getClothYYStiffness(),
		model.getClothXYStiffness(),
		model.getClothXYPoissonRatio(),
		model.getClothYXPoissonRatio(),
		corr1, corr2, corr3);

	if (res)
	{
		if (invMass1 != 0.0)
			x1 += corr1;
		if (invMass2 != 0.0)
			x2 += corr2;
		if (invMass3 != 0.0)
			x3 += corr3;
	}
	return res;
}


//////////////////////////////////////////////////////////////////////////
// StrainTriangleConstraint
//////////////////////////////////////////////////////////////////////////
bool StrainTriangleConstraint::initConstraint(SimulationModel &model, const unsigned int particle1, const unsigned int particle2,
	const unsigned int particle3)
{
	m_bodies[0] = particle1;
	m_bodies[1] = particle2;
	m_bodies[2] = particle3;

	ParticleData &pd = model.getParticles();

	Vector3r &x1 = pd.getPosition0(particle1);
	Vector3r &x2 = pd.getPosition0(particle2);
	Vector3r &x3 = pd.getPosition0(particle3);

	// Bring triangles to xy plane
	const Vector3r y1(x1[0], x1[2], 0.0);
	const Vector3r y2(x2[0], x2[2], 0.0);
	const Vector3r y3(x3[0], x3[2], 0.0);

	return PositionBasedDynamics::init_StrainTriangleConstraint(y1, y2, y3, m_invRestMat);
}

bool StrainTriangleConstraint::solvePositionConstraint(SimulationModel &model)
{
	ParticleData &pd = model.getParticles();

	const unsigned i1 = m_bodies[0];
	const unsigned i2 = m_bodies[1];
	const unsigned i3 = m_bodies[2];

	Vector3r &x1 = pd.getPosition(i1);
	Vector3r &x2 = pd.getPosition(i2);
	Vector3r &x3 = pd.getPosition(i3);

	const Real invMass1 = pd.getInvMass(i1);
	const Real invMass2 = pd.getInvMass(i2);
	const Real invMass3 = pd.getInvMass(i3);

	Vector3r corr1, corr2, corr3;
	const bool res = PositionBasedDynamics::solve_StrainTriangleConstraint(
		x1, invMass1,
		x2, invMass2,
		x3, invMass3,
		m_invRestMat,
		model.getClothXXStiffness(),
		model.getClothYYStiffness(),
		model.getClothXYStiffness(),
		model.getClothNormalizeStretch(),
		model.getClothNormalizeShear(),
		corr1, corr2, corr3);

	if (res)
	{
		if (invMass1 != 0.0)
			x1 += corr1;
		if (invMass2 != 0.0)
			x2 += corr2;
		if (invMass3 != 0.0)
			x3 += corr3;
	}
	return res;
}


//////////////////////////////////////////////////////////////////////////
// VolumeConstraint
//////////////////////////////////////////////////////////////////////////

bool VolumeConstraint::initConstraint(SimulationModel &model, const unsigned int particle1, const unsigned int particle2,
	const unsigned int particle3, const unsigned int particle4)
{
	m_bodies[0] = particle1;
	m_bodies[1] = particle2;
	m_bodies[2] = particle3;
	m_bodies[3] = particle4;
	ParticleData &pd = model.getParticles();

	const Vector3r &p0 = pd.getPosition0(particle1);
	const Vector3r &p1 = pd.getPosition0(particle2);
	const Vector3r &p2 = pd.getPosition0(particle3);
	const Vector3r &p3 = pd.getPosition0(particle4);

	m_restVolume = fabs((1.0 / 6.0) * (p3 - p0).dot((p2 - p0).cross(p1 - p0)));

	return true;
}

bool VolumeConstraint::solvePositionConstraint(SimulationModel &model)
{
	ParticleData &pd = model.getParticles();

	const unsigned i1 = m_bodies[0];
	const unsigned i2 = m_bodies[1];
	const unsigned i3 = m_bodies[2];
	const unsigned i4 = m_bodies[3];

	Vector3r &x1 = pd.getPosition(i1);
	Vector3r &x2 = pd.getPosition(i2);
	Vector3r &x3 = pd.getPosition(i3);
	Vector3r &x4 = pd.getPosition(i4);

	const Real invMass1 = pd.getInvMass(i1);
	const Real invMass2 = pd.getInvMass(i2);
	const Real invMass3 = pd.getInvMass(i3);
	const Real invMass4 = pd.getInvMass(i4);

	Vector3r corr1, corr2, corr3, corr4;
	const bool res = PositionBasedDynamics::solve_VolumeConstraint(x1, invMass1,
		x2, invMass2,
		x3, invMass3,
		x4, invMass4,
		m_restVolume,
		model.getSolidStiffness(),
		model.getSolidStiffness(),
		corr1, corr2, corr3, corr4);

	if (res)
	{
		if (invMass1 != 0.0)
			x1 += corr1;
		if (invMass2 != 0.0)
			x2 += corr2;
		if (invMass3 != 0.0)
			x3 += corr3;
		if (invMass4 != 0.0)
			x4 += corr4;
	}
	return res;
}


//////////////////////////////////////////////////////////////////////////
// FEMTetConstraint
//////////////////////////////////////////////////////////////////////////
bool FEMTetConstraint::initConstraint(SimulationModel &model, const unsigned int particle1, const unsigned int particle2,
									const unsigned int particle3, const unsigned int particle4)
{
	m_bodies[0] = particle1;
	m_bodies[1] = particle2;
	m_bodies[2] = particle3;
	m_bodies[3] = particle4;

	ParticleData &pd = model.getParticles();

	Vector3r &x1 = pd.getPosition0(particle1);
	Vector3r &x2 = pd.getPosition0(particle2);
	Vector3r &x3 = pd.getPosition0(particle3);
	Vector3r &x4 = pd.getPosition0(particle4);

	return PositionBasedDynamics::init_FEMTetraConstraint(x1, x2, x3, x4, m_volume, m_invRestMat);
}

bool FEMTetConstraint::solvePositionConstraint(SimulationModel &model)
{
	ParticleData &pd = model.getParticles();

	const unsigned i1 = m_bodies[0];
	const unsigned i2 = m_bodies[1];
	const unsigned i3 = m_bodies[2];
	const unsigned i4 = m_bodies[3];

	Vector3r &x1 = pd.getPosition(i1);
	Vector3r &x2 = pd.getPosition(i2);
	Vector3r &x3 = pd.getPosition(i3);
	Vector3r &x4 = pd.getPosition(i4);

	const Real invMass1 = pd.getInvMass(i1);
	const Real invMass2 = pd.getInvMass(i2);
	const Real invMass3 = pd.getInvMass(i3);
	const Real invMass4 = pd.getInvMass(i4);

	Real currentVolume = -(1.0 / 6.0) * (x4 - x1).dot((x3 - x1).cross(x2 - x1));
	bool handleInversion = false;
	if (currentVolume / m_volume < 0.2)		// Only 20% of initial volume left
		handleInversion = true;


	Vector3r corr1, corr2, corr3, corr4;
	const bool res = PositionBasedDynamics::solve_FEMTetraConstraint(
		x1, invMass1,
		x2, invMass2,
		x3, invMass3,
		x4, invMass4,
		m_volume,
		m_invRestMat,
		model.getSolidStiffness(),
		model.getSolidPoissonRatio(), handleInversion,
		corr1, corr2, corr3, corr4);

	if (res)
	{
		if (invMass1 != 0.0)
			x1 += corr1;
		if (invMass2 != 0.0)
			x2 += corr2;
		if (invMass3 != 0.0)
			x3 += corr3;
		if (invMass4 != 0.0)
			x4 += corr4;
	}
	return res;
}


//////////////////////////////////////////////////////////////////////////
// StrainTetConstraint
//////////////////////////////////////////////////////////////////////////
bool StrainTetConstraint::initConstraint(SimulationModel &model, const unsigned int particle1, const unsigned int particle2,
	const unsigned int particle3, const unsigned int particle4)
{
	m_bodies[0] = particle1;
	m_bodies[1] = particle2;
	m_bodies[2] = particle3;
	m_bodies[3] = particle4;

	ParticleData &pd = model.getParticles();

	Vector3r &x1 = pd.getPosition0(particle1);
	Vector3r &x2 = pd.getPosition0(particle2);
	Vector3r &x3 = pd.getPosition0(particle3);
	Vector3r &x4 = pd.getPosition0(particle4);

	return PositionBasedDynamics::init_StrainTetraConstraint(x1, x2, x3, x4, m_invRestMat);
}

bool StrainTetConstraint::solvePositionConstraint(SimulationModel &model)
{
	ParticleData &pd = model.getParticles();

	const unsigned i1 = m_bodies[0];
	const unsigned i2 = m_bodies[1];
	const unsigned i3 = m_bodies[2];
	const unsigned i4 = m_bodies[3];

	Vector3r &x1 = pd.getPosition(i1);
	Vector3r &x2 = pd.getPosition(i2);
	Vector3r &x3 = pd.getPosition(i3);
	Vector3r &x4 = pd.getPosition(i4);

	const Real invMass1 = pd.getInvMass(i1);
	const Real invMass2 = pd.getInvMass(i2);
	const Real invMass3 = pd.getInvMass(i3);
	const Real invMass4 = pd.getInvMass(i4);

	Vector3r stiffness(model.getSolidStiffness(), model.getSolidStiffness(), model.getSolidStiffness());

	Vector3r corr1, corr2, corr3, corr4;
	const bool res = PositionBasedDynamics::solve_StrainTetraConstraint(
		x1, invMass1,
		x2, invMass2,
		x3, invMass3,
		x4, invMass4,
		m_invRestMat,
		stiffness,
		stiffness,
		model.getSolidNormalizeStretch(),
		model.getSolidNormalizeShear(),
		corr1, corr2, corr3, corr4);

	if (res)
	{
		if (invMass1 != 0.0)
			x1 += corr1;
		if (invMass2 != 0.0)
			x2 += corr2;
		if (invMass3 != 0.0)
			x3 += corr3;
		if (invMass4 != 0.0)
			x4 += corr4;
	}
	return res;
}

//////////////////////////////////////////////////////////////////////////
// ShapeMatchingConstraint
//////////////////////////////////////////////////////////////////////////
bool ShapeMatchingConstraint::initConstraint(SimulationModel &model, 
			const unsigned int particleIndices[], const unsigned int numClusters[])
{
	ParticleData &pd = model.getParticles();
	for (unsigned int i = 0; i < m_numberOfBodies; i++)
	{
		m_bodies[i] = particleIndices[i];
		m_x0[i] = pd.getPosition0(m_bodies[i]);
		m_w[i] = pd.getInvMass(m_bodies[i]);
		m_numClusters[i] = numClusters[i];
	}

	const bool res = PositionBasedDynamics::init_ShapeMatchingConstraint(m_x0, m_w, m_numberOfBodies, m_restCm, m_invRestMat);
	return res;
}

bool ShapeMatchingConstraint::solvePositionConstraint(SimulationModel &model)
{
	ParticleData &pd = model.getParticles();
	for (unsigned int i = 0; i < m_numberOfBodies; i++)
	{
		m_x[i] = pd.getPosition(m_bodies[i]);
	}

	const bool res = PositionBasedDynamics::solve_ShapeMatchingConstraint(
		m_x0, m_x, m_w, m_numberOfBodies,
		m_restCm, m_invRestMat,
		model.getSolidStiffness(), false,
		m_corr);

	if (res)
	{
		for (unsigned int i = 0; i < m_numberOfBodies; i++)
		{
			// Important: Divide position correction by the number of clusters 
			// which contain the vertex. 
			if (m_w[i] != 0.0)
				pd.getPosition(m_bodies[i]) += (1.0 / m_numClusters[i]) * m_corr[i];
		}
	}
	return res;
}


//////////////////////////////////////////////////////////////////////////
// RigidBodyContactConstraint
//////////////////////////////////////////////////////////////////////////
bool RigidBodyContactConstraint::initConstraint(SimulationModel &model, const unsigned int rbIndex1, const unsigned int rbIndex2,
		const Vector3r &cp1, const Vector3r &cp2,
		const Vector3r &normal, const Real dist,
		const Real restitutionCoeff, const Real stiffness, const Real frictionCoeff)
{
	m_stiffness = stiffness;
	m_frictionCoeff = frictionCoeff;

	m_bodies[0] = rbIndex1;
	m_bodies[1] = rbIndex2;
	SimulationModel::RigidBodyVector &rb = model.getRigidBodies();
	RigidBody &rb1 = *rb[m_bodies[0]];
	RigidBody &rb2 = *rb[m_bodies[1]];

	m_sum_impulses = 0.0;

	return PositionBasedRigidBodyDynamics::init_RigidBodyContactConstraint(
		rb1.getInvMass(),
		rb1.getPosition(),
		rb1.getVelocity(),
		rb1.getInertiaTensorInverseW(),
		rb1.getRotation(),
		rb1.getAngularVelocity(),
		rb2.getInvMass(),
		rb2.getPosition(),
		rb2.getVelocity(),
		rb2.getInertiaTensorInverseW(),
		rb2.getRotation(),
		rb2.getAngularVelocity(),
		cp1, cp2, normal, restitutionCoeff, 
		m_constraintInfo);
}

bool RigidBodyContactConstraint::solveVelocityConstraint(SimulationModel &model)
{
	SimulationModel::RigidBodyVector &rb = model.getRigidBodies();

	RigidBody &rb1 = *rb[m_bodies[0]];
	RigidBody &rb2 = *rb[m_bodies[1]];

	Vector3r corr_v1, corr_v2;
	Vector3r corr_omega1, corr_omega2;
	const bool res = PositionBasedRigidBodyDynamics::velocitySolve_RigidBodyContactConstraint(
		rb1.getInvMass(),
		rb1.getPosition(),
		rb1.getVelocity(),
		rb1.getInertiaTensorInverseW(),
		rb1.getAngularVelocity(),
		rb2.getInvMass(),
		rb2.getPosition(),
		rb2.getVelocity(),
		rb2.getInertiaTensorInverseW(),
		rb2.getAngularVelocity(),
		m_stiffness,
		m_frictionCoeff,
		m_sum_impulses,
		m_constraintInfo,
		corr_v1,
		corr_omega1,
		corr_v2,
		corr_omega2);

	if (res)
	{
		if (rb1.getMass() != 0.0)
		{
			rb1.getVelocity() += corr_v1;
			rb1.getAngularVelocity() += corr_omega1;
		}
		if (rb2.getMass() != 0.0)
		{
			rb2.getVelocity() += corr_v2;
			rb2.getAngularVelocity() += corr_omega2;
		}
	}
	return res;
}

//////////////////////////////////////////////////////////////////////////
// ParticleRigidBodyContactConstraint
//////////////////////////////////////////////////////////////////////////
bool ParticleRigidBodyContactConstraint::initConstraint(SimulationModel &model, 
	const unsigned int particleIndex, const unsigned int rbIndex,
	const Vector3r &cp1, const Vector3r &cp2,
	const Vector3r &normal, const Real dist,
	const Real restitutionCoeff, const Real stiffness, const Real frictionCoeff)
{
	m_stiffness = stiffness;
	m_frictionCoeff = frictionCoeff;

	m_bodies[0] = particleIndex;
	m_bodies[1] = rbIndex;
	SimulationModel::RigidBodyVector &rbs = model.getRigidBodies();
	ParticleData &pd = model.getParticles();

	RigidBody &rb = *rbs[m_bodies[1]];

	m_sum_impulses = 0.0;

	return PositionBasedRigidBodyDynamics::init_ParticleRigidBodyContactConstraint(
		pd.getInvMass(particleIndex),
		pd.getPosition(particleIndex),
		pd.getVelocity(particleIndex),
		rb.getInvMass(),
		rb.getPosition(),
		rb.getVelocity(),
		rb.getInertiaTensorInverseW(),
		rb.getRotation(),
		rb.getAngularVelocity(),		
		cp1, cp2, normal, restitutionCoeff,
		m_constraintInfo);
}

bool ParticleRigidBodyContactConstraint::solveVelocityConstraint(SimulationModel &model)
{
	SimulationModel::RigidBodyVector &rbs = model.getRigidBodies();
	ParticleData &pd = model.getParticles();

	RigidBody &rb = *rbs[m_bodies[1]];

	Vector3r corr_v1, corr_v2;
	Vector3r corr_omega2;
	const bool res = PositionBasedRigidBodyDynamics::velocitySolve_ParticleRigidBodyContactConstraint(
		pd.getInvMass(m_bodies[0]),
		pd.getPosition(m_bodies[0]),
		pd.getVelocity(m_bodies[0]),
		rb.getInvMass(),
		rb.getPosition(),
		rb.getVelocity(),
		rb.getInertiaTensorInverseW(),
		rb.getAngularVelocity(),
		m_stiffness,
		m_frictionCoeff,
		m_sum_impulses,
		m_constraintInfo,
		corr_v1,		
		corr_v2, 
		corr_omega2);

	if (res)
	{
		if (pd.getMass(m_bodies[0]) != 0.0)
		{
			pd.getVelocity(m_bodies[0]) += corr_v1;
		}
		if (rb.getMass() != 0.0)
		{
			rb.getVelocity() += corr_v2;
			rb.getAngularVelocity() += corr_omega2;
		}	
	}
	return res;
}

//////////////////////////////////////////////////////////////////////////
// StretchShearConstraint
//////////////////////////////////////////////////////////////////////////
bool StretchShearConstraint::initConstraint(SimulationModel &model, const unsigned int particle1, const unsigned int particle2, const unsigned int quaternion1)
{
	m_bodies[0] = particle1;
	m_bodies[1] = particle2;
	m_bodies[2] = quaternion1;
	ParticleData &pd = model.getParticles();

	const Vector3r &x1_0 = pd.getPosition0(particle1);
	const Vector3r &x2_0 = pd.getPosition0(particle2);

	m_restLength = (x2_0 - x1_0).norm();

	return true;
}

bool StretchShearConstraint::solvePositionConstraint(SimulationModel &model)
{
	ParticleData &pd = model.getParticles();
	OrientationData &od = model.getOrientations();

	const unsigned i1 = m_bodies[0];
	const unsigned i2 = m_bodies[1];
	const unsigned iq1 = m_bodies[2];

	Vector3r &x1 = pd.getPosition(i1);
	Vector3r &x2 = pd.getPosition(i2);
	Quaternionr &q1 = od.getQuaternion(iq1);
	const Real invMass1 = pd.getInvMass(i1);
	const Real invMass2 = pd.getInvMass(i2);
	const Real invMassq1 = od.getInvMass(iq1);
	Vector3r stiffness(model.getRodShearingStiffness1(), 
					   model.getRodShearingStiffness2(), 
					   model.getRodStretchingStiffness());

	Vector3r corr1, corr2;
	Quaternionr corrq1;
	const bool res = PositionBasedCosseratRods::solve_StretchShearConstraint(
		x1, invMass1, x2, invMass2, q1, invMassq1, 
		stiffness,
		m_restLength, corr1, corr2, corrq1);

	if (res)
	{
		if (invMass1 != 0.0)
			x1 += corr1;
		if (invMass2 != 0.0)
			x2 += corr2;
		if (invMassq1 != 0.0)
		{
			q1.coeffs() += corrq1.coeffs();
			q1.normalize();
		}
	}
	return res;
}

//////////////////////////////////////////////////////////////////////////
// BendTwistConstraint
//////////////////////////////////////////////////////////////////////////
bool BendTwistConstraint::initConstraint(SimulationModel &model, const unsigned int quaternion1, const unsigned int quaternion2)
{
	m_bodies[0] = quaternion1;
	m_bodies[1] = quaternion2;
	OrientationData &od = model.getOrientations();

	const Quaternionr &q1_0 = od.getQuaternion(quaternion1);
	const Quaternionr &q2_0 = od.getQuaternion(quaternion2);

	m_restDarbouxVector = q1_0.conjugate() * q2_0;
	Quaternionr omega_plus, omega_minus;
	omega_plus.coeffs() = m_restDarbouxVector.coeffs() + Quaternionr(1, 0, 0, 0).coeffs();
	omega_minus.coeffs() = m_restDarbouxVector.coeffs() - Quaternionr(1, 0, 0, 0).coeffs();
	if (omega_minus.squaredNorm() > omega_plus.squaredNorm())
		m_restDarbouxVector.coeffs() *= -1.0;

	return true;
}

bool BendTwistConstraint::solvePositionConstraint(SimulationModel &model)
{
	OrientationData &od = model.getOrientations();

	const unsigned i1 = m_bodies[0];
	const unsigned i2 = m_bodies[1];

	Quaternionr &q1 = od.getQuaternion(i1);
	Quaternionr &q2 = od.getQuaternion(i2);
	const Real invMass1 = od.getInvMass(i1);
	const Real invMass2 = od.getInvMass(i2);
	Vector3r stiffness(model.getRodBendingStiffness1(), 
					   model.getRodBendingStiffness2(), 
					   model.getRodTwistingStiffness());

	Quaternionr corr1, corr2;
	const bool res = PositionBasedCosseratRods::solve_BendTwistConstraint(
		q1, invMass1, q2, invMass2, 
		stiffness,
		m_restDarbouxVector, corr1, corr2);

	if (res)
	{
		if (invMass1 != 0.0)
		{
			q1.coeffs() += corr1.coeffs();
			q1.normalize();
		}
			
		if (invMass2 != 0.0)
		{
			q2.coeffs() += corr2.coeffs();
			q2.normalize();
		}
	}
	return res;
}