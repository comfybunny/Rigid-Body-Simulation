#ifndef _RIGIDBODY_
#define _RIGIDBODY_

#include <Eigen/Dense>
#include "dart/dart.h"

class RigidBody {
 public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	RigidBody(dart::dynamics::Shape::ShapeType _type, Eigen::Vector3d _dim) {
        // Create a default rigid body
        mMass = 1.0;
        mPosition.setZero(); // x = (0, 0, 0)
        mOrientation.setIdentity(); // R = identity
		mQuatOrient.setIdentity();	//initialize quaternion
        mColor << 0.9, 0.2, 0.2, 1.0; // Red
        
        if (_type == dart::dynamics::Shape::BOX) {
            mShape = Eigen::make_aligned_shared<dart::dynamics::BoxShape>(_dim);
			Ibody = Eigen::Matrix3d::Zero();
			Ibody(0, 0) = mMass*(_dim(1) * _dim(1) + _dim(2) * _dim(2)) / 12.0;
			Ibody(1, 1) = mMass*(_dim(0) * _dim(0) + _dim(2) * _dim(2)) / 12.0;
			Ibody(2, 2) = mMass*(_dim(0) * _dim(0) + _dim(1) * _dim(1)) / 12.0;

        } else if (_type == dart::dynamics::Shape::ELLIPSOID) {
            mShape = Eigen::make_aligned_shared<dart::dynamics::EllipsoidShape>(_dim);
			Ibody = Eigen::Matrix3d::Zero();
			Ibody(0, 0) = 2 * mMass*(_dim(0) * _dim(0)) / 5.0;
			Ibody(1, 1) = 2 * mMass*(_dim(0) * _dim(0)) / 5.0;
			Ibody(2, 2) = 2 * mMass*(_dim(0) * _dim(0)) / 5.0;
        }
        
        mLinMomentum.setZero();
        mAngMomentum.setZero();
        
        mAccumulatedForce.setZero();
        mAccumulatedTorque.setZero();
    }
    virtual ~RigidBody() {}

    void draw(dart::renderer::RenderInterface* _ri);

    int getConfigSize() {
		return mPosition.size() + mOrientation.size();
    }
    
    double mMass;

	Eigen::Vector3d mPosition;
    Eigen::Quaterniond mQuatOrient; // quaternion
	Eigen::Matrix3d mOrientation;   // rotation matrix
    Eigen::Vector3d mLinMomentum;
    Eigen::Vector3d mAngMomentum;
    dart::dynamics::ShapePtr mShape;
    
	Eigen::Vector3d mAccumulatedForce;
    Eigen::Vector3d mAccumulatedTorque;

    Eigen::Vector4d mColor;
	Eigen::Matrix3d Ibody;
};

#endif
