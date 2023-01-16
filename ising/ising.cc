#include "itensor/all.h"

using namespace itensor;

inline int mod(int x,int N){
    if(x>N)
        return x-N;
    return x;
}
int main(){       
    int N = 100;
    auto sites = SpinHalf(N,{"ConserveQNs=",false});
    //
    // Set the parameters controlling the accuracy of the DMRG
    // calculation for each DMRG sweep. 
    //
    auto sweeps = Sweeps(30);
    sweeps.maxdim() = 10,20,100,100,200;
    sweeps.cutoff() = 1E-10;
    sweeps.niter() = 2;
    sweeps.noise() = 1E-7,1E-8,0.0;
    //
    // Begin the DMRG calculation
    // for the ground state
    //
	auto psi0=MPS(InitState(sites,"Up"));
    for(Real h=0; h<=2; h+=.05){
	    //
	    // Factors of 4 and 2 are to rescale
	    // spin operators into Pauli matrices
	    //
	    auto ampo = AutoMPO(sites);
	    for(int j = 1; j <= N; ++j){
	        ampo += -4,"Sz",j,"Sz",mod(j+1,N);
	    } 
	    for(int j = 1; j <= N; ++j){
	        ampo += -2*h,"Sx",j;
	    }
	    auto H = toMPO(ampo);
		// change Silent= to Quiet= to show intermediate status of each sweep
	    auto [en,psi] = dmrg(H,psi0,sweeps,{"Silent=",true});
	    printfln("{%.10f,  %.10f},",h,en/N);
		psi0=psi;
    }

	return 0;
}

