#include<userfuncs.H>
#include <AMReX_ParmParse.H>

AMREX_GPU_DEVICE_MANAGED CubicSplineInterpolator visc_AR_inter;
AMREX_GPU_DEVICE_MANAGED CubicSplineInterpolator thcond_AR_inter;
AMREX_GPU_DEVICE_MANAGED CubicSplineInterpolator sig_AR_inter;
AMREX_GPU_DEVICE_MANAGED CubicSplineInterpolator cp_AR_inter;
AMREX_GPU_DEVICE_MANAGED CubicSplineInterpolator e_AR_inter;
AMREX_GPU_DEVICE_MANAGED CubicSplineInterpolator visc_H2_inter;
AMREX_GPU_DEVICE_MANAGED CubicSplineInterpolator thcond_H2_inter;
AMREX_GPU_DEVICE_MANAGED CubicSplineInterpolator sig_H2_inter;
AMREX_GPU_DEVICE_MANAGED CubicSplineInterpolator e_H2_inter;
AMREX_GPU_DEVICE_MANAGED CubicSplineInterpolator rad_inter;
AMREX_GPU_DEVICE_MANAGED CubicSplineInterpolator t_from_e_inter;



namespace mflo_user_funcs
{

    AMREX_GPU_DEVICE_MANAGED Real T_ign=4000.0;
    AMREX_GPU_DEVICE_MANAGED Real T_surr=1000.0;
    AMREX_GPU_DEVICE_MANAGED Real T_anod=1000.0;
    AMREX_GPU_DEVICE_MANAGED Real T_cath=3500;
    AMREX_GPU_DEVICE_MANAGED Real T_room=300;
    
    AMREX_GPU_DEVICE_MANAGED Real sig_cath=5.8e7;
    AMREX_GPU_DEVICE_MANAGED Real sig_anod=5.8e7; 
    AMREX_GPU_DEVICE_MANAGED Real cath_rad=0.0005;
       
    AMREX_GPU_DEVICE_MANAGED Real catblocksize=-0.1;
    AMREX_GPU_DEVICE_MANAGED Real p0=1e5;
    AMREX_GPU_DEVICE_MANAGED Real H2molfrac=0.0;
    AMREX_GPU_DEVICE_MANAGED Real Hmolfrac=0.0;
    AMREX_GPU_DEVICE_MANAGED Real siteconc=0.0;
    AMREX_GPU_DEVICE_MANAGED Real elec_rad=0.0025;
    AMREX_GPU_DEVICE_MANAGED Real elec_voltage=500.0;
    AMREX_GPU_DEVICE_MANAGED Real elec_currentden=1e6;
    AMREX_GPU_DEVICE_MANAGED int use_V_flag=1;
    
    AMREX_GPU_DEVICE_MANAGED Real resist=1;
    AMREX_GPU_DEVICE_MANAGED Real eps=0.3;
    AMREX_GPU_DEVICE_MANAGED Real tort=5.0;
    AMREX_GPU_DEVICE_MANAGED Real fs_vel=1.0;
    AMREX_GPU_DEVICE_MANAGED Real vel_coflow=0.01;
    AMREX_GPU_DEVICE_MANAGED Real sponge_zone_dist=15e-3;

	AMREX_GPU_DEVICE_MANAGED int weno_scheme=2;


    void initialize_problem()
    {
        Print()<<"Initializing problem\n";
        ParmParse pp("user");
        pp.query("T_ign",T_ign);
        pp.query("T_surr",T_surr);
        pp.query("T_cath",T_cath);
	    pp.query("T_anod",T_anod);
	    pp.query("T_room",T_room);
	    
	    pp.query("cath_rad",cath_rad);
	    pp.query("sig_cath",sig_cath);
	    pp.query("sig_anod",sig_anod);
	    
        pp.query("catblocksize",catblocksize);
        pp.query("elec_rad",elec_rad);
        pp.query("p0",p0);
        pp.query("H2molfrac",H2molfrac);
        pp.query("Hmolfrac",Hmolfrac);
        pp.query("siteconc",siteconc);
        pp.query("elec_voltage",elec_voltage);
        pp.query("elec_currentden",elec_currentden);
        pp.query("sponge_zone_dist",sponge_zone_dist);
        pp.query("use_V_flag",use_V_flag);

        pp.query("resist",resist);
        pp.query("eps",eps);
        pp.query("tort",tort);
        pp.query("vel_coflow",vel_coflow);  
        
        pp.query("weno_scheme",weno_scheme);      

        
        visc_AR_inter.define("data_visc_AR.txt",0,1);
        thcond_AR_inter.define("data_thcond_AR.txt",0,1);
        sig_AR_inter.define("data_sig_AR.txt",0,1);  
        cp_AR_inter.define("data_cp_AR.txt",0,1);  
        e_AR_inter.define("data_e_AR.txt",0,1);
        visc_H2_inter.define("data_visc_H2.txt",0,1);
        thcond_H2_inter.define("data_thcond_H2.txt",0,1);
        sig_H2_inter.define("data_sig_H2.txt",0,1);  
        e_H2_inter.define("data_e_H2.txt",0,1);
        rad_inter.define("radiation.txt",0,1);       
		t_from_e_inter.define("data_e_AR.txt",1,0); 
    }
}
