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

#include "CDSExecution.h"

CDSExecution::CDSExecution()
{
	mCartesianTarget.Resize(6);

	mNextEEFrame.Identity();
	mTargetFrame.Identity();
	mAttractorFrame.Identity();
	mObjectFrame.Identity();

	mCurrentEEFrame.Identity();
	vCurrentEEPos.Zero();
	mCurrentEERot.Identity();

	vRobotEEAbsolutePos.Zero();
	mRobotEEAbsoluteOrient.Identity();

	vRobotEEPosInTargetFrame.Zero();
	vRobotEEOrientInTargetFrame.Zero();
	vRobotEEOrientInTargetFrameAA.Resize(4);
	vTargetOrientAA.Zero();

	// Default values
	dt = 0.002;
	reachingThr=0.002;
	
	bDoInterp = true; 
	bAttSet = false;
	bRobSet = false;
	bObjSet = false;


}

CDSExecution::~CDSExecution()
{
}

void CDSExecution::init(GMRDynamics* masterDyn, GMRDynamics* slaveDyn, GMR* coupling){

	cdsCtrl = new CDSController(masterDyn, slaveDyn, coupling);
	cpl = new PositionOrientationCoupling();

	cpl->Initialize(2);
	cdsCtrl->InitCoupling(cpl);
	
	if(!cdsCtrl->isInit)
	{
		cout<<"ERROR while initializing cds controller"<<endl;
		exit(1);
	}
}
void CDSExecution::initSimple(std::string base_path, int segmentID) {


	char sMaster[1025];
	char sSlave[1025];
	char sGMM[1025];

	sprintf(sMaster, "%s/Phase%d/masterGMM.txt", base_path.c_str(), segmentID);
	sprintf(sSlave,  "%s/Phase%d/slaveGMM.txt", base_path.c_str(), segmentID);
	sprintf(sGMM,    "%s/Phase%d/cplGMM.txt", base_path.c_str(), segmentID);

	cdsCtrl = new CDSController(sMaster, sSlave, sGMM);

	_init();

}

void CDSExecution::initSimple(std::string base_path, int segmentID, int armID) {


    char sMaster[1025];
    char sSlave[1025];
    char sGMM[1025];

    string arm;
    if (armID == 1)
        arm = "Right";
    if (armID == 2)
        arm = "Left";

    sprintf(sMaster, "%s/Phase%d/%s/masterGMM.txt", base_path.c_str(), segmentID, arm.c_str());
    sprintf(sSlave,  "%s/Phase%d/%s/slaveGMM.txt", base_path.c_str(),  segmentID, arm.c_str());
    sprintf(sGMM,    "%s/Phase%d/%s/cplGMM.txt", base_path.c_str(),  segmentID, arm.c_str());

    cdsCtrl = new CDSController(sMaster, sSlave, sGMM);

    _init();

}

void CDSExecution::initSimple(std::string base_path, int segmentID, int armID, std::string role) {


    char sMaster[1025];
    char sSlave[1025];
    char sGMM[1025];

    string arm;
    if (armID == 1)
        arm = "Right";
    if (armID == 2)
        arm = "Left";

    sprintf(sMaster, "%s/Phase%d/%s/%s_posGMM.txt", base_path.c_str(), segmentID, arm.c_str(), role.c_str());
    sprintf(sSlave,  "%s/Phase%d/%s/%s_oriGMM.txt", base_path.c_str(),  segmentID, arm.c_str(), role.c_str());
    sprintf(sGMM,    "%s/Phase%d/%s/%s_cplGMM.txt", base_path.c_str(),  segmentID, arm.c_str(), role.c_str());

    cdsCtrl = new CDSController(sMaster, sSlave, sGMM);

    _init();

}

void CDSExecution::initSimple(std::string base_path, int segmentID, std::string model_id) {


	char sMaster[1025];
	char sSlave[1025];
	char sGMM[1025];

	if (model_id==""){
		sprintf(sMaster, "%s/Phase%d/masterGMM.txt", base_path.c_str(), segmentID);
		sprintf(sSlave,  "%s/Phase%d/slaveGMM.txt", base_path.c_str(), segmentID);
		sprintf(sGMM,    "%s/Phase%d/cplGMM.txt", base_path.c_str(), segmentID);
	}
	else{
		sprintf(sMaster, "%s/Phase%d/masterGMM_%s.txt", base_path.c_str(), segmentID, model_id.c_str());
		sprintf(sSlave,  "%s/Phase%d/slaveGMM.txt", base_path.c_str(), segmentID);
		sprintf(sGMM,    "%s/Phase%d/cplGMM_%s.txt", base_path.c_str(), segmentID, model_id.c_str());
	}
	cdsCtrl = new CDSController(sMaster, sSlave, sGMM);

	_init();

}

void CDSExecution::initSimple(int segmentID){

        char sMaster[256];
        char sSlave[256];
        char sGMM[256];        

        sprintf(sMaster, "data/packages/CDSPlotRF/Models/s%d/masterGMM.txt", segmentID);
        sprintf(sSlave,  "data/packages/CDSPlotRF/Models/s%d/slaveGMM.txt", segmentID);
        sprintf(sGMM,    "data/packages/CDSPlotRF/Models/s%d/cplGMM.txt", segmentID);

	cdsCtrl = new CDSController(sMaster, sSlave, sGMM);
	_init();

}

void CDSExecution::_init() {
	//Initialize the coupling function \hat{x} = \Psi(x) by specifying the dimension of the output = 1. Here \Psi(x) = ||x|| \in R^N
	cpl = new PositionOrientationCoupling();
	cpl->Initialize(2);
		cdsCtrl->InitCoupling(cpl);
//		cdsCtrl->masterDynamics->printInfo();
//		cdsCtrl->slaveDynamics->printInfo();
//		cdsCtrl->coupling->printInfo();
		if(!cdsCtrl->isInit)
		{
			cout<<"ERROR while initializing cds controller"<<endl;
			throw;
		}
		mCartesianTarget.Resize(6);

			mNextEEFrame.Identity();
			mTargetFrame.Identity();
			mAttractorFrame.Identity();
			mObjectFrame.Identity();

			mCurrentEEFrame.Identity();
			vCurrentEEPos.Zero();
			mCurrentEERot.Identity();

			vRobotEEAbsolutePos.Zero();
			mRobotEEAbsoluteOrient.Identity();

			vRobotEEPosInTargetFrame.Zero();
			vRobotEEOrientInTargetFrame.Zero();
			vRobotEEOrientInTargetFrameAA.Resize(4);
			vTargetOrientAA.Zero();

			// Default values
			dt = 0.002;
			reachingThr=0.002;

			bDoInterp = true;
			bAttSet = false;
			bRobSet = false;
			bObjSet = false;
}

void CDSExecution::postInit(){

	if (bAttSet == true  && bRobSet == true && bObjSet == true)
		syncWorldWithModel();
	else
		cout << "Could not initialize" << endl;
}

void CDSExecution::setObjectFrame(const MathLib::Matrix4& mObjectPose){
	mObjectFrame = mObjectPose;
	bObjSet = true;
}

void CDSExecution::setAttractorFrame(const MathLib::Matrix4& mAttractorPose){
	mAttractorFrame = mAttractorPose;
	bAttSet = true;
}

void CDSExecution::setCurrentEEPose(const MathLib::Matrix4& mCurrentEEPose){
	mCurrentEEFrame = mCurrentEEPose;

	vCurrentEEPos = mCurrentEEFrame.GetTranslation();
	mCurrentEERot = mCurrentEEFrame.GetOrientation();

	vRobotEEAbsolutePos = mCurrentEEFrame.GetTranslation();
	mRobotEEAbsoluteOrient = mCurrentEEFrame.GetOrientation();
	bRobSet = true;
}


void CDSExecution::setMotionParameters(double alpha, double beta, double lambda, double reachingThreshold,
			CDSController::DynamicsType masterDynamicsID, CDSController::DynamicsType slaveDynamicsID){
	cdsCtrl->SetAlpha(alpha); 	// multiplier for boosting slave dynamics (smaller is slower)
	cdsCtrl->SetLambda(lambda);  	// multiplier for the master system     (higher is faster)
	cdsCtrl->SetBeta(beta);		// error multiplier for the CDS Controller
	cdsCtrl->SetSlaveDynamics(slaveDynamicsID);
	cdsCtrl->SetMasterDynamics(masterDynamicsID);
	reachingThr = reachingThreshold;
}

const MathLib::Matrix4& CDSExecution::getNextEEPose(){

		syncWorldWithModel();

		//Update both master and slave systems
		if(!cdsCtrl->Update(dt))
			cout<<"WARNING: CDS Controller could not update"<<endl;

		//Read off the next desired states of master and slave systems. These are still relative values.
		cdsCtrl->GetMasterRelativeState(vRobotEEPosInTargetFrame.Array()); 	//relative position
		cdsCtrl->GetSlaveRelativeState(vRobotEEOrientInTargetFrameAA.Array());	// absolute orientation

		//Converting position back to absolute values (world frame)
		vRobotEEAbsolutePos = mTargetFrame.Transform(vRobotEEPosInTargetFrame);

		// compute robot orientation in axis angle
		Matrix3 tmpEEOrient = Matrix3::SRotationV(vRobotEEOrientInTargetFrameAA(3),
			Vector3(vRobotEEOrientInTargetFrameAA(0),vRobotEEOrientInTargetFrameAA(1), vRobotEEOrientInTargetFrameAA(2)));

		// back to Absolute frame
		mRobotEEAbsoluteOrient = mTargetFrame.GetOrientation() * tmpEEOrient;

		mNextEEFrame.SetTranslation(vRobotEEAbsolutePos);
		mNextEEFrame.SetOrientation(mRobotEEAbsoluteOrient);

		return mNextEEFrame;
}



bool CDSExecution::syncWorldWithModel()
{
	// Compute the mTargetFrame with respect to the Attractor frame and Object Frame
	mObjectFrame.Mult(mAttractorFrame, mTargetFrame); 

	// Compute Robot EE frame w.r.t target reference frame 	
	// vRobotEEPosInTargetFrame = mTargetFrame.InverseTransformation().Transform(vRobotEEAbsolutePos); // ==>> Does not work like this due to bug in MathLib

	vRobotEEPosInTargetFrame =
		mTargetFrame.GetOrientation().Transpose()*(vRobotEEAbsolutePos - mTargetFrame.GetTranslation());

	if(vRobotEEPosInTargetFrame.Norm() > reachingThr)
	{
		//Set GMR  states
		cdsCtrl->SetMasterRelativeState(vRobotEEPosInTargetFrame.Array());

		// Computes the orientation of the EE in the target frame
		Matrix3   mTmpMat3;
		mTargetFrame.GetOrientation().Transpose().Mult(mRobotEEAbsoluteOrient, mTmpMat3);
		vRobotEEOrientInTargetFrame = mTmpMat3.GetExactRotationAxis();

		//vTmpOrient = mRobotEEAbsoluteOrient.GetExactRotationAxis();
		//mTargetFrame.GetOrientation().Transpose().Mult(mRobotEEAbsoluteOrient, mTmpMat3);
		Vector3 vTmpOrient = vRobotEEOrientInTargetFrame;
		float tmpTheta;
		tmpTheta = vTmpOrient.Norm();
		if (tmpTheta > 0.001)
			vTmpOrient = vTmpOrient/tmpTheta;
		Vector vTmpSlaveState(4);
		vTmpSlaveState(0) = vTmpOrient(0);
		vTmpSlaveState(1) = vTmpOrient(1);
		vTmpSlaveState(2) = vTmpOrient(2);
		vTmpSlaveState(3) = tmpTheta;
		cdsCtrl->SetSlaveRelativeState(vTmpSlaveState.Array());

		vTargetOrientAA = mAttractorFrame.GetOrientation().GetExactRotationAxis();
		vTmpOrient = vTargetOrientAA;
		tmpTheta = vTmpOrient.Norm();
		if(fabs(tmpTheta) < 1e-5) {
			vTmpSlaveState.Zero();
		} else {
		vTmpOrient = vTmpOrient/tmpTheta;
		vTmpSlaveState(0) = vTmpOrient(0);
		vTmpSlaveState(1) = vTmpOrient(1);
		vTmpSlaveState(2) = vTmpOrient(2);
		vTmpSlaveState(3) = tmpTheta;
		}
		cdsCtrl->SetSlaveTarget(vTmpSlaveState);
	}

	return true;
}

