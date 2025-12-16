#ifndef _USERFUNCS_H_
#define _USERFUNCS_H_

#include<AMReX_REAL.H>
#include<AMReX.H>
#include<string>
#include<AMReX_Vector.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Box.H>
#include <globalDefines.H>
#include <mflo_constants.H>

using namespace amrex;
namespace mflo_user_funcs
{
    extern AMREX_GPU_DEVICE_MANAGED Real cylrad;
    extern AMREX_GPU_DEVICE_MANAGED Real fs_vel;
    extern AMREX_GPU_DEVICE_MANAGED Real rho0;
    extern AMREX_GPU_DEVICE_MANAGED Real p0;
    extern AMREX_GPU_DEVICE_MANAGED Real Aconc_inlet;
    extern AMREX_GPU_DEVICE_MANAGED Real Aconc_initial;
    extern AMREX_GPU_DEVICE_MANAGED Real siteconc;
    extern AMREX_GPU_DEVICE_MANAGED Real sponge_zone_dist;

    void initialize_problem();

    AMREX_GPU_DEVICE AMREX_INLINE
        void compute_fluid_transport(int i, int j, int k,
                Array4<Real> const& phi,
                Array4<Real> const& transpcoeffs,
                GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                GpuArray<Real, AMREX_SPACEDIM> dx,
                const Real time)
        {
            Real x=prob_lo[0]+(i+0.5)*dx[0];
            Real lscale=10*dx[0];
            Real trans_scale=20.0;
            Real sponge_fac=0.5*(1.0+std::tanh((x-sponge_zone_dist)/lscale));

            transpcoeffs(i,j,k,VISC_INDX  ) = 1e-5*(1.0+sponge_fac*trans_scale);       //viscosity
            transpcoeffs(i,j,k,THCOND_INDX) = 0.01*(1.0+sponge_fac*trans_scale); //thermal conductivity
        }

    AMREX_GPU_DEVICE AMREX_INLINE
        void compute_fluid_source(int i, int j, int k,
                Array4<Real> const& phi,
                Array4<Real> const& source,
                GpuArray<Real, AMREX_SPACEDIM> prob_lo,
                GpuArray<Real, AMREX_SPACEDIM> dx,
                const Real time)
        {
            for(int nc=0;nc<FLO_NVARS;nc++)
            {
                source(i,j,k,nc)=zeroval;   
            }
        }
}
#endif
