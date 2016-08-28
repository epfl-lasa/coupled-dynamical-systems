/*
 * Copyright (C) 2014 Learning Algorithms and Systems Laboratory, EPFL, Switzerland
 * Author:  Lucia Pais
 * email:   lucia.pais@epf.ch
 * website: lasa.epfl.ch
 *
 * Permission is granted to copy, distribute, and/or modify this program
 * under the terms of the GNU General Public License, version 2 or any
 * later version published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
 * Public License for more details
 */

#ifndef CDSExecution_H
#define CDSExecution_H

#include "GMR.h"

#include "MathLib/MathLib.h"
#include "CDSController.h"


using namespace MathLib;

class PositionOrientationCoupling : public CouplingFunctionBase
{
 public:

  void UpdateCouplingFunction(Vector& input, Vector& output)
  {
    output(0) = input(0)*input(0)+input(1)*input(1);
    output(1) = input(2);
  }
};


class CDSExecution
{
public:
	CDSExecution();
	~CDSExecution();

	void init(GMRDynamics* masterDyn, GMRDynamics* slaveDyn, GMR* coupling);
	void initSimple(int segmentID);
	void initSimple(std::string base_path, int segmentID);
    void initSimple(std::string base_path, int segmentID, int armID);
    void initSimple(std::string base_path, int segmentID, int armID, std::string role);
    void initSimple(std::string base_path, int segmentID, std::string model_id);
    void postInit();

	void setObjectFrame(const MathLib::Matrix4& mObjectPose);
	void setAttractorFrame(const MathLib::Matrix4& mAttractorPose);
	void setMotionParameters(double alpha, double beta, double lambda, double reachingThreshold,
			CDSController::DynamicsType masterDynamicsID, CDSController::DynamicsType slaveDynamicsID); // Values for master and slave systems
	void setCurrentEEPose(const MathLib::Matrix4& mCurrentEEPose);
	void setDT(double runtime_dt){dt = runtime_dt;}

	const MathLib::Matrix4& getNextEEPose();
	 
private:
	void _init();
	bool syncWorldWithModel();		// Updates the controller

	CDSController* cdsCtrl;
	PositionOrientationCoupling *cpl;

	MathLib::Matrix4 mTargetFrame;		// is computed from the attractor and object frame
	MathLib::Matrix4 mAttractorFrame;	// where we want to reach relative to the object frame
	MathLib::Matrix4 mObjectFrame;		// object frame in base frame

	MathLib::Matrix4 mCurrentEEFrame;	// robot complete pose in base frame
	MathLib::Vector3 vCurrentEEPos;		// current robot position in base frame
	MathLib::Matrix3 mCurrentEERot;		// current robot orientation in base frame
	
	MathLib::Vector3 vRobotEEAbsolutePos;	// model computed robot absolute position, in base frame
	MathLib::Matrix3 mRobotEEAbsoluteOrient;// model computed robot absolute orientation, in base frame

	MathLib::Matrix4 mNextEEFrame;		// robot pose as returned by CDS, in base frame

	MathLib::Vector3 vRobotEEPosInTargetFrame;	// Robot position relative to the target frame
	MathLib::Vector3 vRobotEEOrientInTargetFrame;	// Robot Orientation relative to the target frame
	MathLib::Vector  vRobotEEOrientInTargetFrameAA;
	MathLib::Vector3 vTargetOrientAA;		// the target (i.e attractor) absolute orientation in normalized axis angles
	MathLib::Vector  mCartesianTarget;		// Vector 6, position and orientation of the next ee pose, just for debug purpose

	bool bDoInterp; 
	bool bAttSet;
	bool bRobSet;
	bool bObjSet;
	double dt;
	double reachingThr;

};


#endif // CDSExecution_H
