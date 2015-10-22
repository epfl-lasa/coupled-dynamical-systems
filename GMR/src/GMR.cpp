#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>

#include "GMR.h"
using namespace std;

GMR::GMR()
	: DGMRState()
	, nStates(0)
	, nVar(0)
	, in_dim(0)
	, out_dim(0)
	, isInit(false)
	, GMMState()
{
}

GMR::~GMR()
{
	free( pdfx);
    for (unsigned int s = 0; s < nStates; s++)
    {
        delete(GMMState[s].Mu);
        delete(GMMState[s].Sigma[0]);
        delete(GMMState[s].Sigma);

        if (isInit) // if init(...) has been run that stuff exists too
        {
//	free( diff);
//	free( diffa);
//	free (diffb);
//	free (diffp);

            free(DGMRState[s].muI);
            free(DGMRState[s].muO);

            free(DGMRState[s].SigmaIIInv[0]);
            free(DGMRState[s].SigmaIIInv);
            free(DGMRState[s].SigmaOI[0]);
            free(DGMRState[s].SigmaOI);
        }
    }
}

GMR::GMR(const std::string & file)
	: DGMRState()
	, nStates(0)
	, nVar(0)
	, in_dim(0)
	, out_dim(0)
	, isInit(false)
	, GMMState()
{
	loadParameters(file);
}

void GMR::loadParameters(const std::string & file)
{
	isInit = false;
	cout<<"Reading GMM from file \""<<file<<"\"\n";

	FILE *fid;
	int res;
	double dtemp;
	fid = fopen(file.c_str(), "r");

	if(!fid)
	{
		cout<<"Error opening file \""<<file<<"\"\n";
		throw;
	}
	res=fscanf(fid, "%lf\n", &dtemp);
	nStates = (int)dtemp;

	res=fscanf(fid, "%lf\n", &(dtemp));
	nVar = (int)dtemp;

	GMMState.resize(nStates);
	DGMRState.resize(nStates);

	for( unsigned int s=0; s<nStates; s++ ){
		GMMState[s].Mu       = svector(nVar);
		GMMState[s].Sigma    = smatrix(nVar, nVar );
	}

	pdfx = svector(nStates);

	// Read priors
	for( unsigned int s=0; s<nStates; s++ )
		res=fscanf(fid, "%lf", &(GMMState[s].Prio ) );

    // Read Mu
	for( unsigned int i=0; i<nVar; i++ )
		for( unsigned int s=0; s<nStates; s++ )
			res=fscanf(fid, "%lf", &(GMMState[s].Mu[i]) );

    // Read Sigmas
	for( unsigned int s=0; s<nStates; s++ )
		for( unsigned int i=0; i<nVar; i++ )
			for( unsigned int j=0; j<nVar; j++ )
				res=fscanf(fid, "%lf", &(GMMState[s].Sigma[i][j]));

	relAttractorCoord.Resize(3);

	fclose(fid);


	cout<<"GMM Read Successfully..\n";

}

/** Copy constructor
 *  \param other Object to copy from
 */
GMR::GMR(const GMR& other)
{
    // copy standard variables
    nStates = other.nStates;
	in_dim = other.in_dim;
	out_dim = other.out_dim;
	nVar = other.nVar;
	isInit = other.isInit;

	// initialise pointers
	GMMState.resize(nStates);
	DGMRState.resize(nStates);

	for( unsigned int s=0; s<nStates; s++ ){
		GMMState[s].Mu       = svector(nVar);
		GMMState[s].Sigma    = smatrix(nVar, nVar );
	}

	pdfx = svector(nStates);
//	matcp_c1(nStates,other.pdfx,pdfx);


		relAttractorCoord.Resize(other.relAttractorCoord.Size());
		for( unsigned int i=0; i<other.relAttractorCoord.Size(); i++ )
			relAttractorCoord(i) = other.relAttractorCoord.At(i);


	// copy priors
	for( unsigned int s=0; s<nStates; s++ )
		GMMState[s].Prio =other.GMMState[s].Prio;

    // copy Mu
	for( unsigned int i=0; i<nVar; i++ )
		for( unsigned int s=0; s<nStates; s++ )
            GMMState[s].Mu[i] = other.GMMState[s].Mu[i];

    // copy Sigmas
	for( unsigned int s=0; s<nStates; s++ )
		for( unsigned int i=0; i<nVar; i++ )
			for( unsigned int j=0; j<nVar; j++ )
                GMMState[s].Sigma[i][j] = other.GMMState[s].Sigma[i][j];

    if(isInit)
    {
	diff = svector(in_dim);
	diffp = svector(in_dim);
	diffa = svector(in_dim);
	diffb = svector(out_dim);

        // Initialize GMR Variables
        for( unsigned int s=0; s<nStates; s++ ){
            DGMRState[s].SigmaOI = smatrix(out_dim, in_dim );
            DGMRState[s].SigmaIIInv = smatrix(in_dim, in_dim );
            DGMRState[s].muI = svector(in_dim);
            DGMRState[s].muO = svector(out_dim);
        }

        // copy determinants
        for( unsigned int s=0; s<nStates; s++ )
            DGMRState[s].det =other.DGMRState[s].det;

        // copy MuI
        for( unsigned int s=0; s<nStates; s++ )
            for( unsigned int i=0; i<in_dim; i++)
                DGMRState[s].muI[i] = other.DGMRState[s].muI[i];

        // copy MuO
        for( unsigned int s=0; s<nStates; s++ )
            for( unsigned int i=0; i<out_dim; i++)
                DGMRState[s].muO[i] = other.DGMRState[s].muO[i];

        // copy SigmaIIInv
        for( unsigned int s=0; s<nStates; s++ )
            for( unsigned int i=0; i<in_dim; i++)
                for( unsigned int j=0; j<in_dim; j++)
                    DGMRState[s].SigmaIIInv[i][j] = other.DGMRState[s].SigmaIIInv[i][j];

        // copy SigmaOI
        for( unsigned int s=0; s<nStates; s++ )
            for( unsigned int i=0; i<out_dim; i++)
                for( unsigned int j=0; j<in_dim; j++)
                    DGMRState[s].SigmaOI[i][j] = other.DGMRState[s].SigmaOI[i][j];
    }
}

GMR::GMR(const std::vector<GMMStates>& gmm_state,
         unsigned int n_states, unsigned int n_vars):
    nStates(n_states), nVar(n_vars), isInit(false)
{
        // getting the memory...
	GMMState.resize(nStates);
	DGMRState.resize(nStates);

	for(unsigned int s=0; s<nStates; s++ ){
		GMMState[s].Mu       = svector(nVar);
		GMMState[s].Sigma    = smatrix(nVar, nVar );
	}

	pdfx = svector(nStates);
	relAttractorCoord.Resize(3);

	// copying content...
	for( unsigned int s=0; s<nStates; s++ )
	{
		// copy priors
		GMMState[s].Prio =gmm_state[s].Prio;

		// copy Mu
		for( unsigned int i=0; i<nVar; i++ )
			GMMState[s].Mu[i] = gmm_state[s].Mu[i];

		// copy Sigmas
		for( unsigned int i=0; i<nVar; i++ )
			for( unsigned int j=0; j<nVar; j++ )
				GMMState[s].Sigma[i][j] = gmm_state[s].Sigma[i][j];
  }
}

/** Assignment operator
 *  \param other Object to assign from
 *  \return A reference to this
 */
GMR& GMR::operator=(const GMR& rhs)
{
    if (this == &rhs) return *this; // handle self assignment

    // copy standard variables
    nStates = rhs.nStates;
	in_dim = rhs.in_dim;
	out_dim = rhs.out_dim;
	nVar = rhs.nVar;
    isInit = rhs.isInit;

	// initialise pointers
	GMMState.resize(nStates);
	DGMRState.resize(nStates);

	for( unsigned int s=0; s<nStates; s++ ){
		GMMState[s].Mu       = svector(nVar);
		GMMState[s].Sigma    = smatrix(nVar, nVar );
	}

//	matcp_c1(nStates,rhs.pdfx,pdfx);
//	matcp_c1(in_dim,rhs.diff,diff);
//	matcp_c1(in_dim,rhs.diffp,diffp);
//	matcp_c1(in_dim,rhs.diffa,diffa);
//	matcp_c1(out_dim,rhs.diffb,diffb);

//	relAttractorCoord.Resize(out_dim);

	relAttractorCoord.Resize(rhs.relAttractorCoord.Size());
    for( unsigned int i=0; i<rhs.relAttractorCoord.Size(); i++ )
        relAttractorCoord(i)= rhs.relAttractorCoord.At(i);

	pdfx = svector(nStates);

	// copy priors
	for( unsigned int s=0; s<nStates; s++ )
		GMMState[s].Prio =rhs.GMMState[s].Prio;

    // copy Mu
	for( unsigned int i=0; i<nVar; i++ )
		for( unsigned int s=0; s<nStates; s++ )
			GMMState[s].Mu[i] = rhs.GMMState[s].Mu[i];

    // copy Sigmas
	for( unsigned int s=0; s<nStates; s++ )
		for( unsigned int i=0; i<nVar; i++ )
			for( unsigned int j=0; j<nVar; j++ )
				GMMState[s].Sigma[i][j] = rhs.GMMState[s].Sigma[i][j];

    if(isInit)
    {

	diff = svector(in_dim);
	diffa = svector(in_dim);
	diffp = svector(in_dim);
	diffb = svector(out_dim);

        // Initialize GMR Variables
        for( unsigned int s=0; s<nStates; s++ ){
            DGMRState[s].SigmaOI = smatrix(out_dim, in_dim );
            DGMRState[s].SigmaIIInv = smatrix(in_dim, in_dim );
            DGMRState[s].muI = svector(in_dim);
            DGMRState[s].muO = svector(out_dim);
        }

        // copy determinants
        for( unsigned int s=0; s<nStates; s++ )
            DGMRState[s].det = rhs.DGMRState[s].det;

        // copy MuI
        for( unsigned int s=0; s<nStates; s++ )
            for( unsigned int i=0; i<in_dim; i++)
                DGMRState[s].muI[i] = rhs.DGMRState[s].muI[i];

        // copy MuO
        for( unsigned int s=0; s<nStates; s++ )
            for( unsigned int i=0; i<out_dim; i++)
                DGMRState[s].muO[i] = rhs.DGMRState[s].muO[i];

        // copy SigmaIIInv
        for( unsigned int s=0; s<nStates; s++ )
            for( unsigned int i=0; i<in_dim; i++)
                for( unsigned int j=0; j<in_dim; j++)
                    DGMRState[s].SigmaIIInv[i][j] = rhs.DGMRState[s].SigmaIIInv[i][j];

        // copy SigmaOI
        for( unsigned int s=0; s<nStates; s++ )
            for( unsigned int i=0; i<out_dim; i++)
                for( unsigned int j=0; j<in_dim; j++)
                    DGMRState[s].SigmaOI[i][j] = rhs.DGMRState[s].SigmaOI[i][j];
    }

    return *this;
}

/** \brief Attribute matrix blocks to input/output and initialise
 *
 * \param in_index Vector Position of the inputs
 * \param out_index Vector Position of the outputs
 * \return void
 *
 */
void GMR::initGMR(const std::vector<int> & in_index,
                  const std::vector<int> & out_index)
{

	cout<<"Initializing GMR...\n";

	unsigned int i,j,s;

	double **covII;

	this->in_dim = in_index.size();
	this->out_dim = out_index.size();

	if(in_dim<=0 || out_dim<=0)
		return;

    nVar = in_dim + out_dim; // TODO: convert this to throw exception ? (should be equal)


	diff = svector(in_dim);
	diffp = svector(in_dim);
	diffa = svector(in_dim);
	diffb = svector(out_dim);

    relAttractorCoord.Resize(out_dim);

	covII     = smatrix(in_dim, in_dim);

	// Initialize GMR Variables
	for( s=0; s<nStates; s++ ){
		DGMRState[s].SigmaOI = smatrix(out_dim, in_dim );
		DGMRState[s].SigmaIIInv = smatrix(in_dim, in_dim );
		DGMRState[s].muI = svector(in_dim);
		DGMRState[s].muO = svector(out_dim);
	}

	for( s=0; s<nStates; s++ ){
		for( i=0; i<in_dim; i++){
			for( j=0; j<in_dim; j++)
            {
				covII[i][j]                = GMMState[s].Sigma[(int)in_index[i]][(int)in_index[j]];
			}

			DGMRState[s].muI[i] = GMMState[s].Mu[(int)in_index[i]];
		}

		for( i=0; i<out_dim; i++){
			for( j=0; j<in_dim; j++){

				DGMRState[s].SigmaOI[i][j] = GMMState[s].Sigma[(int)out_index[i]][(int)in_index[j]];
			}
			DGMRState[s].muO[i] = GMMState[s].Mu[(int)out_index[i]];
		}

		if(in_dim==1)
		{
			DGMRState[s].SigmaIIInv[0][0] = 1.0/covII[0][0];
			DGMRState[s].det = covII[0][0];
		}
		else
			matsvdinvall( covII  , in_dim, in_dim, 1e-50, &(DGMRState[s].det) , DGMRState[s].SigmaIIInv   );

	}

	isInit = true;


	cout<<"GMR Model Initialized Successfully...\n";

}

/** \brief Get relative likelihood from the PDF for one of the available states
 *
 * \param state int One of the available states from the SEDS model
 * \param in_data double* the random variable input
 * \return double the relative likelihood
 *
 */
double GMR::GaussianPDF(int state, const double *in_data)
{
//	double *diff, *diffp;
	double p;
//	diff  = svector(in_dim);
//	diffp = svector(in_dim);
	matsub_c1(in_dim, in_data, DGMRState[state].muI, diff );
	matmul_c(in_dim, in_dim, DGMRState[state].SigmaIIInv, diff, diffp );
	p = in_prod(in_dim, diff, diffp );
	p = exp(-0.5*p) / sqrt(pow(2.0*3.14159,in_dim)*(dabs(DGMRState[state].det)+REAL_MIN));

	if(p < REAL_MIN){
		return REAL_MIN;
	}
	else{
		return p;
	}
}


/** \brief Takes the relative state of the system and returns the relative state output of the GMR
 *
 * \param x double* Input ((master)position usually)
 * \param xdot double* Output (SEDS: velocitiy, CDS: slave position)
 * \return void
 *
 */
void GMR::getGMROutput(const double *x, double *xdot)
{
        unsigned int s, i;
//	double *diff, *diffa, *diffb, *relx;
	double pa, h;

//	relx = svector(in_dim);
//	diff   = svector(in_dim);
//	diffa  = svector(in_dim);
//	diffb  = svector(out_dim);

	for(i=0;i<out_dim;i++)
		xdot[i]=0;

	  //	cout<<"in: "<<in_dim<<endl;
	for( s=0; s<nStates; s++ ){
		pdfx[s] = GaussianPDF(s, x);

	}
	pa = vsum(nStates, pdfx ) + REAL_MIN;

	for( s=0; s<nStates; s++ ){
		h = pdfx[s]/pa;

		matsub_c1(in_dim, x, DGMRState[s].muI, diff );
		matmul_c(in_dim, in_dim, DGMRState[s].SigmaIIInv, diff , diffa );
		matmul_c(out_dim, in_dim, DGMRState[s].SigmaOI   , diffa, diffb  );

		for(i=0; i<out_dim; i++ )
			xdot[i] += h*(diffb[i] + DGMRState[s].muO[i]);
	}
}



void GMR::printInfo() const
{
	cout<<"-------------------------------- GMM ---------------------------------"<<endl;
	cout<<"InDim:       \t"<<in_dim<<endl;
	cout<<"OutDim:      \t"<<out_dim<<endl;
	cout<<"Gaussians: \t"<<nStates<<endl;

	cout<<"\nPriors: "<<endl;
	for(unsigned int i=0;i<nStates;i++)
		cout<<GMMState[i].Prio<<"    ";

	cout<<"\nMu: "<<endl;
	for(unsigned int i=0;i<nStates;i++)
	{
		for(unsigned int j=0;j<nVar;j++)
			cout<<GMMState[i].Mu[j]<<"    ";
		cout<<endl;
	}

	cout<<endl;
	for(unsigned int i=0;i<nStates;i++)
	{
		cout<<"Sigma: ["<<i<<"]"<<endl;
		for(unsigned int j=0;j<nVar;j++)
		{
			for(unsigned int k=0;k<nVar;k++)
				cout<<GMMState[i].Sigma[j][k]<<"    ";
			cout<<endl;
		}
		cout<<endl;
	}

	cout<<endl;
	for(unsigned int i=0;i<nStates;i++)
	{
		cout<<"Sigma-II-Inv: ["<<i<<"]"<<endl;
		for(unsigned int j=0;j<in_dim;j++)
		{
			for(unsigned int k=0;k<in_dim;k++)
				cout<<DGMRState[i].SigmaIIInv[j][k]<<"    ";
			cout<<endl;
		}
		cout<<endl;
	}

	cout<<endl;
	for(unsigned int i=0;i<nStates;i++)
	{
		cout<<"Sigma-OI: ["<<i<<"]"<<endl;
		for(unsigned int j=0;j<out_dim;j++)
		{
			for(unsigned int k=0;k<in_dim;k++)
				cout<<DGMRState[i].SigmaOI[j][k]<<"    ";
			cout<<endl;
		}
		cout<<endl;
	}

	cout<<endl;
	for(unsigned int i=0;i<nStates;i++)
	{
		cout<<"Mu-I: ["<<i<<"]"<<endl;
		for(unsigned int j=0;j<in_dim;j++)
				cout<<DGMRState[i].muI[j]<<"    ";
			cout<<endl;
	}

	cout<<endl;
	for(unsigned int i=0;i<nStates;i++)
	{
		cout<<"Mu-O: ["<<i<<"]"<<endl;
		for(unsigned int j=0;j<out_dim;j++)
				cout<<DGMRState[i].muO[j]<<"    ";
			cout<<endl;
	}

	cout<<endl<<"Det-Sigma-II: "<<endl;
	for(unsigned int i=0;i<nStates;i++)
	{
				cout<<DGMRState[i].det<<"    ";
	}
	cout<<endl;
	cout<<"----------------------------------------------------------------------"<<endl;
}

