#include<userfuncs.H>
#include <AMReX_ParmParse.H>

namespace mflo_user_funcs
{
    AMREX_GPU_DEVICE_MANAGED int nparticles=10;
    AMREX_GPU_DEVICE_MANAGED Real partradmean=0.5;
    AMREX_GPU_DEVICE_MANAGED Real partradstddev=0.5;
    AMREX_GPU_DEVICE_MANAGED Real gasvel=0.2366;
    AMREX_GPU_DEVICE_MANAGED Real rho0=1.0;
    AMREX_GPU_DEVICE_MANAGED Real p0=1.0;
    AMREX_GPU_DEVICE_MANAGED Real CH4conc_inlet=0.0;
    AMREX_GPU_DEVICE_MANAGED Real CH4conc_initial=0.0;
    AMREX_GPU_DEVICE_MANAGED Real siteconc=1.0;
    AMREX_GPU_DEVICE_MANAGED Real sponge_zone_dist=23e-3;
    AMREX_GPU_DEVICE_MANAGED Real partcenters_x[MAX_PARTSIZE]={0.0};
    AMREX_GPU_DEVICE_MANAGED Real partcenters_y[MAX_PARTSIZE]={0.0};
    AMREX_GPU_DEVICE_MANAGED Real partcenters_z[MAX_PARTSIZE]={0.0};
    AMREX_GPU_DEVICE_MANAGED Real partrads[MAX_PARTSIZE]={0.0};

    void initialize_problem()
    {
        Print()<<"Initializing problem\n";
    
        ParmParse pp("user");
        pp.query("nparticles",nparticles);
        pp.query("gasvel",gasvel);
        pp.query("partradmean",partradmean);
        pp.query("partradstddev",partradstddev);
        pp.query("rho0",rho0);
        pp.query("p0",p0);
        pp.query("CH4conc_inlet",CH4conc_inlet);
        pp.query("CH4conc_initial",CH4conc_initial);
        pp.query("siteconc",siteconc);
        pp.query("sponge_zone_dist",sponge_zone_dist);

        ParmParse ppgeom("geometry");
        std::vector<amrex::Real> plo;
        plo.resize(AMREX_SPACEDIM);
        std::vector<amrex::Real> phi;
        phi.resize(AMREX_SPACEDIM);

        ppgeom.queryarr("prob_lo",plo);
        ppgeom.queryarr("prob_hi",phi);
        
        Real lenx=phi[0]-plo[0];
        Real leny=phi[1]-plo[1];
        Real lenz=phi[2]-plo[2];

        if(amrex::ParallelDescriptor::IOProcessor())
        {
            for(int np=0;np<nparticles;np++)
            {
                partrads[np]=amrex::RandomNormal(partradmean,partradstddev);
                partcenters_x[np]=plo[0]+(0.05+amrex::Random())*0.95*lenx;
                partcenters_y[np]=plo[1]+(0.05+amrex::Random())*0.95*leny;
                partcenters_z[np]=plo[2]+(0.05+amrex::Random())*0.7*lenz;
            }
        }
        
        amrex::ParallelDescriptor::Bcast(partrads,nparticles, amrex::ParallelDescriptor::IOProcessorNumber());
        amrex::ParallelDescriptor::Bcast(partcenters_x,nparticles, amrex::ParallelDescriptor::IOProcessorNumber());
        amrex::ParallelDescriptor::Bcast(partcenters_y,nparticles, amrex::ParallelDescriptor::IOProcessorNumber());
        amrex::ParallelDescriptor::Bcast(partcenters_z,nparticles, amrex::ParallelDescriptor::IOProcessorNumber());
    }
}
