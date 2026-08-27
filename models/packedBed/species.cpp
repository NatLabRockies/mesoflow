#include<species.H>

namespace mflo_species
{
    amrex::Vector<std::string> specnames(NUM_SPECIES);
    AMREX_GPU_DEVICE_MANAGED amrex::Real advect_flags[NUM_SPECIES]={0};
    AMREX_GPU_DEVICE_MANAGED amrex::Real molwts[NUM_SPECIES]={0.0};

    void init()
    {
        specnames[AIR_ID]="air";
        specnames[CH4_ID]="methane";
        specnames[H2_ID]="hydrogen";
        specnames[C_ID]="carbon";
        specnames[S_ID]="sites";
        
        //kg/mol
        molwts[AIR_ID]=0.0289;
        molwts[CH4_ID]=0.016;
        molwts[H2_ID]=0.002;
        molwts[C_ID]=0.012;
        molwts[S_ID]=0.001;

        advect_flags[AIR_ID]=1;
        advect_flags[CH4_ID]=1;
        advect_flags[H2_ID]=1;
        advect_flags[C_ID]=0;
        advect_flags[S_ID]=0;
    }    
    void close()
    {
        specnames.clear();
    }
    int find_id(std::string specname)
    {
        int loc=-1;
        auto it=std::find(specnames.begin(),specnames.end(),specname);
        if(it != specnames.end())
        {
            loc=it-specnames.begin();
        }
        return(loc);
    }
}
