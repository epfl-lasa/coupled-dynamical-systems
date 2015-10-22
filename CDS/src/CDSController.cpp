/*
 * Copyright (C) 2014 Learning Algorithms and Systems Laboratory, EPFL, Switzerland
 * Author: Ashwini Shukla, modified by Lucia Pais to work with position and orientation 
 * rather than arm motion and finger opening-closing.
 * email:  ashwini.shukla@epfl.ch, lucia.pais@epf.ch
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

#include "CDSController.h"

//! \brief convert a vector 4 (utheta) into a matrix.
//MathLib::Matrix3 toMatrix(const MathLib::Vector & utheta)
MathLib::Matrix3 toMatrix( MathLib::Vector & utheta)
{
	MathLib::Matrix3 mat = Matrix3::SRotationV(utheta(3),
		MathLib::Vector3(utheta(0), utheta(1), utheta(2)));
	return mat;
}

//! \brief convert a matrix into  a vector 4 (utheta)
MathLib::Vector fromMatrix(const MathLib::Matrix3 & mat)
{
	Vector3 utheta = mat.GetExactRotationAxis();
	Vector  uthetaFull(4);

	float tmpTheta = utheta.Norm();
	if (tmpTheta > 0.001)
		utheta = utheta/tmpTheta;
	uthetaFull(0) = utheta(0);
	uthetaFull(1) = utheta(1);
	uthetaFull(2) = utheta(2);
	uthetaFull(3) = tmpTheta;
	return uthetaFull;
}


CDSController::CDSController()
{
}

CDSController::~CDSController()
{}

CDSController::CDSController(const std::string & masterDyn_file,
		const std::string & slaveDyn_file,
		const std::string & coupling_file)
{
	loadParameters (masterDyn_file, slaveDyn_file, coupling_file);
}

void CDSController::loadParameters(const std::string & masterDyn_file,
		const std::string & slaveDyn_file,
		const std::string & coupling_file)
{
	masterDynamics = new GMRDynamics(masterDyn_file);
	slaveDynamics = new GMRDynamics(slaveDyn_file);
	coupling = new GMR(coupling_file);

	slaveDynamicsType = UTHETA;

	if(masterDynamics->isInit && slaveDynamics->isInit)
	{
		masterDim = masterDynamics->dim;
		slaveDim = slaveDynamics->dim;

		masterRelState.Resize(masterDim);
		slaveRelState.Resize(slaveDim);
		mTmpMaster.Resize(masterDim);
		mTmpSlave.Resize(slaveDim);
		slaveRelVel.Resize(slaveDim);
		masterRelVel.Resize(masterDim);
		slaveDesState.Resize(slaveDim);

		isInit=false;
	}
}

CDSController::CDSController(GMRDynamics* master_dynamics,
    GMRDynamics* slave_dynamics, GMR* coupling_dyn)
{
	masterDynamics = master_dynamics;
	slaveDynamics = slave_dynamics;
	this->coupling = coupling_dyn;

	if(masterDynamics->isInit && slaveDynamics->isInit)
	{
		masterDim = masterDynamics->dim;
		slaveDim = slaveDynamics->dim;

		masterRelState.Resize(masterDim);
		slaveRelState.Resize(slaveDim);
		mTmpMaster.Resize(masterDim);
		mTmpSlave.Resize(slaveDim);
		slaveRelVel.Resize(slaveDim);
		masterRelVel.Resize(masterDim);
		slaveDesState.Resize(slaveDim);
		slaveTarget.Resize(slaveDim);

		masterRelState.Zero();
		slaveRelState.Zero();

		isInit=false;
  }
}

bool CDSController::InitCoupling(CouplingFunctionBase* cpl_func)
{

	if(!cpl_func || cpl_func->output_dim == 0)
		return false;

	coupling_func = cpl_func;

	Psix.Resize(coupling_func->output_dim);

	std::vector<int> in_dim, out_dim;
	in_dim.resize(coupling_func->output_dim);
	out_dim.resize(slaveDim);

	if(in_dim.size() + out_dim.size() != coupling->getnVar())
		return false;

	for(unsigned int i=0;i<in_dim.size();i++)
		in_dim[i] = i;

	for(unsigned int i=0;i<out_dim.size();i++)
		out_dim[i] = i + in_dim.size();

	// setup CDS GMR
	coupling->initGMR(in_dim,out_dim);
	alpha = 1.0; // should be one but we set a higher dt to increase sim. speed
	beta = 1.0;
	lambda = 1.0;

	masterRelState.Zero();
	slaveRelState.Zero();
	masterRelVel.Zero();
	slaveRelVel.Zero();
	slaveDesState.Zero();
	slaveTarget.Zero();
	Psix.Zero();

	coupling->getGMROutput(Psix.Array(), slaveDesState.Array());
	slaveGoal = slaveDesState;
	isInit = true;

	return true;
}

void CDSController::SetSlaveTarget(const Vector& tar)
{
	slaveTarget = tar;
}

void CDSController::SetMasterDynamics(DynamicsType dynamicsID){
	masterDynamicsType = dynamicsID;
}

void CDSController::SetSlaveDynamics(DynamicsType dynamicsID){
	slaveDynamicsType = dynamicsID;
}


bool CDSController::Update(double dt)
{
	if(!isInit)
		return false;

	switch (masterDynamicsType) {
		case LINEAR_DYNAMICS:
			masterRelVel = masterRelState * (-0.5);
			break;
		case MODEL_DYNAMICS:
			masterDynamics->getGMROutput(masterRelState.Array(),masterRelVel.Array());
			break;
		default:
			break;
	}

	masterRelState += masterRelVel*dt*lambda;


	coupling_func->UpdateCouplingFunction(masterRelState, Psix);
	if(Psix.Norm() < 0.02)
		Psix.Zero();
	coupling->getGMROutput(Psix.Array(), slaveDesState.Array());


	double mult;
	for(unsigned int i=0;i<slaveDesState.Size();i++)
	{
		mult = std::min(Psix.Norm(), 1.);
		slaveDesState[i]  = DEG2RAD(90)*((1-mult)*exp(-mult*9) - 1);
	}
	
	coupling->getGMROutput(masterRelState.Array(),slaveDesState.Array());
	if (masterRelState.Norm() < 0.01){
		slaveDesState(3) = 0;
	}

	switch (slaveDynamicsType){
	case MODEL_DYNAMICS:
	{
		slaveDesState = slaveDesState + slaveTarget;
		mTmpSlave = (slaveRelState - slaveDesState) * beta;
		slaveDynamics->getGMROutput(mTmpSlave.Array(), slaveRelVel.Array());
		slaveRelState += slaveRelVel * alpha * dt;
		break;
	}
	case LINEAR_DYNAMICS:
	{
		slaveDesState = slaveDesState + slaveTarget;
		mTmpSlave = (slaveRelState - slaveDesState) * beta;
		slaveRelVel = mTmpSlave * (-0.5);
		slaveRelState += slaveRelVel * alpha * dt;
		break;
	}
	case UTHETA:
	{
		slaveDesState += slaveTarget;
		mTmpSlave = fromMatrix(toMatrix(slaveRelState).Transpose() * toMatrix(slaveDesState))*beta;
		mTmpSlave *= beta;
		slaveRelVel = mTmpSlave * (-0.5);
		slaveRelState += slaveRelVel * alpha * dt;
		break;
	}
	case NO_DYNAMICS:
	{
		slaveRelState = slaveDesState;
		break;
	}
	default:
		break;
	}

	return true;
}

