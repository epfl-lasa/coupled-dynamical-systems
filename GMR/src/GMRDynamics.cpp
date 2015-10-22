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


#include "GMR.h"
using namespace std;

GMRDynamics::GMRDynamics(const std::string & file) : GMR(file)
{
	cout<<"Initializing for GMRDynamics"<<endl;
  if(nVar%2 != 0)
    {
      std::cerr <<"ERROR: Dynamic model has odd nmber of dimensions"<<endl;
      throw;
    }
  dim = (unsigned int)(nVar/2);
  std::vector<int> in_ind(dim), out_ind(dim);
  for(unsigned int i=0;i<dim;i++)
    {
      in_ind[i] = i;
      out_ind[i] = i+dim;
    }

  initGMR(in_ind, out_ind);

  rel_state = new double[dim];
  tmp_vel = new double[dim];
  cout<<"Done"<<endl;

}

GMRDynamics::GMRDynamics(const GMRDynamics& gmr) : GMR( gmr )
{
	this->dim = gmr.dim;
	this->rel_state = new double[this->dim];
	this->tmp_vel = new double[this->dim];

	for(unsigned int i=0;i<dim;i++)
	{
		this->rel_state[i] = gmr.rel_state[i];
		this->tmp_vel[i] = gmr.tmp_vel[i];
	}
}

GMRDynamics::GMRDynamics(const std::vector<GMMStates> &  gmm_state, unsigned int n_states,
    unsigned int n_vars) : GMR( gmm_state, n_states, n_vars)
{
  cout<<"Initializing for GMRDynamics"<<endl;
  if(nVar%2 != 0)
    {
      cout<<"ERROR: Dynamic model has odd nmber of dimensions"<<endl;
      return;
    }
  dim = (unsigned int)(nVar/2);
  std::vector<int> in_ind(dim), out_ind(dim);
  for(unsigned int i=0;i<dim;i++)
    {
      in_ind[i] = i;
      out_ind[i] = i+dim;
    }

  initGMR(in_ind, out_ind);

  rel_state = new double[dim];
  tmp_vel = new double[dim];
  cout<<"Done"<<endl;
}

GMRDynamics::~GMRDynamics()
{
	delete rel_state;
	delete tmp_vel;
}

bool GMRDynamics::Update(double dt)
{
	if(!isInit)
		return false;

	getGMROutput(rel_state,tmp_vel );
	for(unsigned int i=0;i<dim;i++)
		rel_state[i] += tmp_vel[i]*dt;

	return true;
}

void GMRDynamics::setRelativeState(const double* state)
{
	for(unsigned int i=0;i<dim;i++)
		rel_state[i] = state[i];
}

void GMRDynamics::getRelativeState(double* state) const
{
	for(unsigned int i=0;i<dim;i++)
		state[i] = rel_state[i];
}

void GMRDynamics::printInfo() const
{
	cout<<"---------------------------- GMM-Dynamics ----------------------------"<<endl;
	cout<<"Dim:       \t"<<dim<<endl;
	GMR::printInfo();
}

