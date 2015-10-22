/*
 * Copyright (C) 2011 Learning Algorithms and Systems Laboratory, EPFL, Switzerland
 * Author: 		Ashwini Shukla
 * Created on: 	Dec 31, 2011
 * E-mail:   	ashwini.shukla@epfl.ch
 * Website: 	lasa.epfl.ch
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


#ifndef _GMR_H_
#define _GMR_H_

#include "math_lib.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include "MathLib/MathLib.h"

using namespace MathLib;
#define REAL_MIN (1e-200)

struct GMMStates {
	double *Mu;
	double **Sigma;
	double Prio;
};

struct DGMRStates{
	double **SigmaOI;
	double **SigmaIIInv;
	double *muI;
	double *muO;
	double det; // for SigmaIOIOInv
};

class GMR
{

protected:
	std::vector<DGMRStates> DGMRState;
	unsigned int nStates;
	unsigned int nVar;
	/** \brief Get relative likelihood from the PDF for one of the available states
	 *
	 * \param state int One of the available states from the SEDS model
	 * \param in_data double* the random variable input
	 * \return double the relative likelihood
	 *
	 */
	double GaussianPDF(int state, const double *in_data);

private:
	unsigned int in_dim;
	unsigned int out_dim;
	Vector relAttractorCoord;
	double *pdfx, *diff, *diffa, *diffb, *diffp;

public:
	bool isInit;
	std::vector<GMMStates> GMMState;

	GMR();
	~GMR();
	GMR(const std::string & file);
	GMR(const GMR& other);
	GMR(const std::vector<GMMStates>& gmm_state, unsigned int n_states, unsigned int n_vars);
	GMR& operator=(const GMR& other);

	/** \brief Attribute matrix blocks to input/output and initialise
	 *
	 * \param in_index Vector Position of the inputs
	 * \param out_index Vector Position of the outputs
	 * \return void
	 *
	 */
	void initGMR(const std::vector<int> & in_index, const std::vector<int> & out_index);

	/** \brief Takes the relative state of the system and returns the relative state output of the GMR
	 *
	 * \param x double* Input ((master)position usually)
	 * \param xdot double* Output (SEDS: velocitiy, CDS: slave position)
	 * \return void
	 *
	 */
	void getGMROutput(const double *x_in, double *x_out);

	void loadParameters(const std::string & file);

	// accessors:
	unsigned int getNumStates() const {return nStates;}
	unsigned int getInDim() const {return in_dim;}
	unsigned int getOutDim() const {return out_dim;}
	unsigned int getnVar() const {return nVar;}
	bool getIsInit() const {return isInit;}
	virtual void printInfo() const;
};


// GMR Dynamics
class GMRDynamics : public GMR
{

public:
	unsigned int dim;
	GMRDynamics(const std::string & file);
	GMRDynamics(const GMRDynamics& gmr);
	GMRDynamics(const std::vector<GMMStates> & gmm_state, unsigned int n_states, unsigned int n_vars);
	~GMRDynamics();

	/*!
	 * \brief Euler integration.
	 * compute velocity using getGMROutput and integrate it.
	 */
	bool Update(double dt);
	void setRelativeState(const double* state);
	void getRelativeState(double* state) const;
	void printInfo() const;
	unsigned int getDim() const {return dim;}

private:

	using GMR::initGMR;
	using GMR::getInDim;
	using GMR::getOutDim;

	double *rel_state;
	double *tmp_vel;

};


#endif /* GMRDynamics_H_ */
