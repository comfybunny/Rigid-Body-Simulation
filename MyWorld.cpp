#include "MyWorld.h"
#include "RigidBody.h"
#include "CollisionInterface.h"
#include <iostream>

using namespace Eigen;
using namespace std;

MyWorld::MyWorld() {
    mFrame = 0;
    mTimeStep = 0.001;
    mGravity = Vector3d(0.0, -9.8, 0.0);
    mForce.setZero();
    // Create a collision detector
    mCollisionDetector = new CollisionInterface();
    
    // Create and intialize two default rigid bodies
    RigidBody *rb1 = new RigidBody(dart::dynamics::Shape::BOX, Vector3d(0.05, 0.05, 0.05));
    mCollisionDetector->addRigidBody(rb1, "box"); // Put rb1 in collision detector
	
    rb1->mPosition[0] = -0.3;
    rb1->mPosition[1] = -0.5;
    
    rb1->mAngMomentum = Vector3d(0.0, 0.01, 0.0);
    mRigidBodies.push_back(rb1);
    
    RigidBody *rb2 = new RigidBody(dart::dynamics::Shape::ELLIPSOID, Vector3d(0.06, 0.06, 0.06));
    mCollisionDetector->addRigidBody(rb2, "ellipse"); // Put rb2 in collision detector
    rb2->mPosition[0] = 0.3;
    rb2->mPosition[1] = -0.5;
    rb2->mAngMomentum = Vector3d(0.01, 0.0, 0.0);
    rb2->mColor = Vector4d(0.2, 0.8, 0.2, 1.0); // Blue
    mRigidBodies.push_back(rb2);
}

void MyWorld::initializePinata() {
    // Add pinata to the collison detector
    mCollisionDetector->addSkeleton(mPinataWorld->getSkeleton(0));
    
    // Add some damping in the Pinata joints
    int nJoints = mPinataWorld->getSkeleton(0)->getNumBodyNodes();
    for (int i = 0; i < nJoints; i++) {
        int nDofs = mPinataWorld->getSkeleton(0)->getJoint(i)->getNumDofs();
        for (int j = 0; j < nDofs; j++)
        mPinataWorld->getSkeleton(0)->getJoint(i)->setDampingCoefficient(j, 1.0);
    }
    
    // Weld two seems to make a box
    dart::dynamics::BodyNode* top = mPinataWorld->getSkeleton(0)->getBodyNode("top");
    dart::dynamics::BodyNode* front = mPinataWorld->getSkeleton(0)->getBodyNode("front");
    dart::dynamics::BodyNode* back = mPinataWorld->getSkeleton(0)->getBodyNode("back");
    dart::constraint::WeldJointConstraint *joint1 = new dart::constraint::WeldJointConstraint(top, front);
    dart::constraint::WeldJointConstraint *joint2 = new dart::constraint::WeldJointConstraint(top, back);
    mPinataWorld->getConstraintSolver()->addConstraint(joint1);
    mPinataWorld->getConstraintSolver()->addConstraint(joint2);
}

MyWorld::~MyWorld() {
    for (int i = 0; i < mRigidBodies.size(); i++)
    delete mRigidBodies[i];
    mRigidBodies.clear();
    if (mCollisionDetector)
    delete mCollisionDetector;
}

void MyWorld::simulate() {
    mFrame++;
    // TODO: The skeleton code has provided the integration of position and linear momentum,
    // your first job is to fill in the integration of orientation and angular momentum.
    for (int i = 0; i < mRigidBodies.size(); i++) {

		mRigidBodies[i]->mQuatOrient.normalize();
		


        // derivative of position and linear momentum
        Eigen::Vector3d dPos = mRigidBodies[i]->mLinMomentum / mRigidBodies[i]->mMass;
        Eigen::Vector3d dLinMom = mRigidBodies[i]->mMass * mGravity + mRigidBodies[i]->mAccumulatedForce;
        
		// derivative of orientation
		Eigen::Matrix3d rotMatrix = mRigidBodies[i]->mQuatOrient.toRotationMatrix();
		Eigen::Vector3d w = (rotMatrix*(mRigidBodies[i]->Ibody)*rotMatrix.transpose()).inverse()*mRigidBodies[i]->mAngMomentum;
		Eigen::Quaterniond qDot = getQdot(w, (mRigidBodies[i]->mQuatOrient));

		// derivative of angular momentum (torque)
		// QUESTION : This doesn't get updated until we get to the collision handler?

        // update position and linear momentum
        mRigidBodies[i]->mPosition += dPos * mTimeStep;
        mRigidBodies[i]->mLinMomentum += mTimeStep * dLinMom;

		//update angular momentum and orientation
		mRigidBodies[i]->mAngMomentum += mRigidBodies[i]->mAccumulatedTorque*mTimeStep;

		mRigidBodies[i]->mQuatOrient.w() += qDot.w()*mTimeStep;
		mRigidBodies[i]->mQuatOrient.vec() += qDot.vec()*mTimeStep;

		mRigidBodies[i]->mOrientation = mRigidBodies[i]->mQuatOrient.toRotationMatrix();
    }
    
    // Reset accumulated force and torque to be zero after a complete integration
    for (int i = 0; i < mRigidBodies.size(); i++) {
        mRigidBodies[i]->mAccumulatedForce.setZero();
        mRigidBodies[i]->mAccumulatedTorque.setZero();
    }
    
    // Apply external force to the pinata
    mPinataWorld->getSkeleton(0)->getBodyNode("bottom")->addExtForce(mForce);
    mForce.setZero();
    
    // Simulate Pinata using DART
    mPinataWorld->step();
    
    // Run collision detector
    mCollisionDetector->checkCollision();
    
    collisionHandling();
    
    // Break the pinata if it has enough momentum
    if (mPinataWorld->getSkeleton(0)->getCOMLinearVelocity().norm() > 0.6)
    mPinataWorld->getConstraintSolver()->removeAllConstraints();
}

// TODO: fill in the collision handling function
void MyWorld::collisionHandling() {
    // restitution coefficient
    double epsilon = 0.8;
    
    // TODO: handle the collision events
	int numContacts = mCollisionDetector->getNumContacts();

	for (int collision = 0; collision < numContacts; collision++) {
		
		RigidContact& curr = mCollisionDetector->getContact(collision);
		
		double mass_a = 0;
		double mass_b = 0;
		double a_term = 0;
		double b_term = 0;
		Eigen::Vector3d pa_dot_minus = curr.pinataVelocity;
		Eigen::Vector3d pb_dot_minus = curr.pinataVelocity;

		// this is the pinata
		if (curr.rb1 != nullptr){
			mass_a = 1.0/curr.rb1->mMass;
			Eigen::Matrix3d rotMatrix_a = curr.rb1->mQuatOrient.toRotationMatrix();
			Eigen::Matrix3d IaInv = ((rotMatrix_a*(curr.rb1->Ibody)*rotMatrix_a.transpose()).inverse());
			Eigen::Vector3d w_a = IaInv*(curr.rb1->mAngMomentum);
			Eigen::Vector3d r_a = curr.point - curr.rb1->mPosition;
			pa_dot_minus = curr.rb1->mLinMomentum / curr.rb1->mMass + (w_a).cross(curr.point);
			a_term = curr.normal.dot((IaInv*(r_a.cross(curr.normal))).cross(r_a));
		}
		if (curr.rb2 != nullptr) {
			mass_b = 1.0/curr.rb2->mMass;
			Eigen::Matrix3d rotMatrix_b = curr.rb2->mQuatOrient.toRotationMatrix();
			Eigen::Matrix3d IbInv = (rotMatrix_b*(curr.rb2->Ibody)*rotMatrix_b.transpose()).inverse();
			Eigen::Vector3d w_b = IbInv*(curr.rb2->mAngMomentum);
			Eigen::Vector3d r_b = curr.point - curr.rb2->mPosition;
			pb_dot_minus = curr.rb2->mLinMomentum / curr.rb2->mMass + (w_b).cross(curr.point);
			b_term = curr.normal.dot((IbInv*(r_b.cross(curr.normal))).cross(r_b));
		}
		
		double v_r_minus = curr.normal.dot(pa_dot_minus - pb_dot_minus);
	//	cout << v_r_minus << endl;
		if (v_r_minus > 0) {
			continue;
		}
		double j = -(1.0 + epsilon)*v_r_minus / (mass_a + mass_b + a_term + b_term);

		
		if (curr.rb1 != nullptr) {
			cout << "rb1 != null" << endl;
			cout << j*curr.normal << endl;
			curr.rb1->mAccumulatedForce += j*curr.normal / mTimeStep;
			curr.rb1->mAccumulatedTorque += ((curr.point - curr.rb1->mPosition).cross(j*curr.normal))/ mTimeStep;
		}
		if (curr.rb2 != nullptr) {
			cout << "rb2 != null" << endl;
			cout << j*curr.normal << endl;
			curr.rb2->mAccumulatedForce -= j*curr.normal /mTimeStep;
			curr.rb2->mAccumulatedTorque -= ((curr.point - curr.rb2->mPosition).cross(j*curr.normal))/mTimeStep;
		}
		
	}
}

Eigen::Quaterniond MyWorld::getQdot(const Eigen::Vector3d& w1, const Eigen::Quaterniond &q){
	// double w = -1.0 * q.vec().dot(w);
	Eigen::Vector3d w = w1 / 2.0;
	Eigen:Vector3d v = q.w()*w + w.cross(q.vec());
	Eigen::Quaterniond res;
	res.setIdentity();
	res.w() = -1.0 * q.vec().dot(w);
	res.vec() = v;
	return res;
}