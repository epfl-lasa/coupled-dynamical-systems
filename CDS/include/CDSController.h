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

#ifndef CDSCONTROLLER_H
#define CDSCONTROLLER_H

#include "GMR.h"
#include "MathLib/MathLib.h"

using namespace MathLib;

class CouplingFunctionBase;

class CDSController
{
public:
        enum DynamicsType {
                MODEL_DYNAMICS = 1,
                LINEAR_DYNAMICS = 2,
                NO_DYNAMICS = 3,
                UTHETA =4
        };
public:

	GMRDynamics* masterDynamics;
	GMRDynamics* slaveDynamics;
	GMR* coupling;

	/** Default constructor */
	CDSController();

	/** Init controller from files
	 * \param masterDyn const std::string & Master dynamics model file name
	 * \param slaveDyn const std::string & Slave dynamics model file name
	 * \param coupling const std::string & Coupling model file name
	 */
	CDSController(const std::string & masterDyn,
			const std::string & slaveDyn,
			const std::string & coupling);

	CDSController(GMRDynamics* master_dynamics, GMRDynamics* slave_dynamics,
	GMR* coupling);
	/** Load the parameters from files
	 * \param masterDyn const std::string & Master dynamics model file name
	 * \param slaveDyn const std::string & Slave dynamics model file name
	 * \param coupling const std::string & Coupling model file name
	 */
	void loadParameters (const std::string & masterDyn,
			const std::string & slaveDyn,
			const std::string & coupling);

	/** Initialize the coupling GMM
	 * \param cpl_func CouplingFunctionBase*  object to specify the coupling function x_tilda = Psi(x) and the dimensions of the coupling.
	 *        Must extend the base class and implement the pure virtual coupling function first
	 */
	bool InitCoupling(CouplingFunctionBase* cpl_func);

	void SetSlaveTarget(const Vector& tar);
	void SetSlaveDynamics(DynamicsType dynamicsID);
	void SetMasterDynamics(DynamicsType dynamicsID);

	/** Default destructor */
	~CDSController();
	// ---------------------------------
	// --- get, set and misc methods ---
	// ---------------------------------
	double GetAlpha() const
	{ return alpha; }
	void SetAlpha(double new_alpha)
	{ alpha = fabs(new_alpha); }
	double GetLambda() const
	{ return lambda; }
	void SetLambda(double new_lambda)
	{ lambda = fabs(new_lambda); }
	double GetBeta() const { return beta; }
	void SetBeta(double new_beta)
	{ beta = fabs(new_beta); }

	double GetSlaveDim() const { return slaveDim; }
	double GetMasterDim() const { return masterDim; }

	// Set the master state, relative to the goal. The user must do the transformation from the world frame to goal frame themselves.
	void SetMasterRelativeState(const double* state)
	{
		for(unsigned int i=0;i<masterDim;i++)
			masterRelState[i] = state[i];
	}

	// Set the slave state, relative to the goal. The user must do the transformation from the world frame to goal frame themselves.
	void SetSlaveRelativeState(const double* state)
	{
		for(unsigned int i=0;i<slaveDim;i++)
			slaveRelState[i] = state[i];
	}

	//Get the master state, relative to the goal.
    void GetMasterRelativeState (double* state) //const
	{
        masterRelState.At(0);
        masterRelState(0);
		for(unsigned int i=0;i<masterDim;i++)
			state[i] = masterRelState(i);
	}

	//Get the slave state, relative to the goal.
    void GetSlaveRelativeState(double* state) //const
	{
		for(unsigned int i=0;i<slaveDim;i++)
			state[i] = slaveRelState(i);
	}

	//Get the master GMR velocity
    void GetMasterRelativeVelocity(double* vel) //const
	{
		for(unsigned int i=0;i<masterDim;i++)
			vel[i] = masterRelVel(i);
	}

	//Get the slave GMR velocity
    void GetSlaveRelativeVelocity(double* vel) //const
	{
		for(unsigned int i=0;i<slaveDim;i++)
			vel[i] = slaveRelVel(i);
	}

	//Update the controller. All the states will be updated and can be read out after calling this method.
	bool Update(double dt);

	//Is controller properly initialized and ready to use?
	bool isInit;

	//defines how to treat the dynamics computation
	int slaveDynamicsType, masterDynamicsType;


private:

	unsigned int masterDim;
	unsigned int slaveDim;
	Vector masterRelState;
	Vector slaveRelState;
	Vector masterRelVel;
	Vector slaveRelVel;
	Vector slaveDesState;
	Vector slaveTarget;
	Vector Psix;

	Vector slaveGoal;

	CouplingFunctionBase* coupling_func;

	double alpha;               //!< Slave speed multiplier, controls retracting speed
	double beta;                //!< Slave position/inferred position multiplier, controls retracting amplitude
	double lambda;


	// temp buffers
	Vector mTmpMaster;
	Vector mTmpSlave;

};

class CouplingFunctionBase
{
public:
	CouplingFunctionBase(){output_dim = 0;}
	virtual ~CouplingFunctionBase(){}
	unsigned int output_dim;

	void Initialize(unsigned int _output_dim) { output_dim = _output_dim; }
	virtual  void UpdateCouplingFunction(Vector& input, Vector& output)=0;
};
#endif // CDSCONTROLLER_H
