#include "itensor/all.h"

using namespace itensor;

inline int mod(int x,int N){
    if(x>N)
        return x-N;
    return x;
}
int main(){       
    int N = 20;
    auto sites = SpinOne(N,{"ConserveQNs=",false});
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
	auto psi0=randomMPS(sites);
    auto ampo = AutoMPO(sites);
	// change j<N to j<=N to make it periodic
    for(int j = 1; j < N; ++j){
		auto k=mod(j+1,N);
		const double t=1.0/3;
        ampo += "Sx",j,"Sx",k;
        ampo += "Sy",j,"Sy",k;
        ampo += "Sz",j,"Sz",k;
		ampo += t , "Sx", j, "Sx", k, "Sx", j, "Sx", k;
		ampo += t , "Sx", j, "Sx", k, "Sy", j, "Sy", k;
		ampo += t , "Sx", j, "Sx", k, "Sz", j, "Sz", k;
		ampo += t , "Sy", j, "Sy", k, "Sx", j, "Sx", k;
		ampo += t , "Sy", j, "Sy", k, "Sy", j, "Sy", k;
		ampo += t , "Sy", j, "Sy", k, "Sz", j, "Sz", k;
		ampo += t , "Sz", j, "Sz", k, "Sx", j, "Sx", k;
		ampo += t , "Sz", j, "Sz", k, "Sy", j, "Sy", k;
		ampo += t , "Sz", j, "Sz", k, "Sz", j, "Sz", k;
    } 
    auto H = toMPO(ampo);
    auto [en,psi] = dmrg(H,psi0,sweeps,{"Quiet=",true});
    printfln("ground state energy: %.10f",en);

    auto wfs = std::vector<MPS>{};
    wfs.push_back(psi);

	for(int i=1;i<6;i++){
		// change Silent= to Quiet= to show intermediate status of each sweep
		auto [en1,psi1] = dmrg(H,wfs,randomMPS(sites),sweeps,{"Silent=",true,"Weight=",20.0});

	    printfln("difference for %d-th excited state: %.10f",i,en1-en);
		
		wfs.push_back(psi1);

	}


	return 0;
}

