# This script is aimed to calculate the absorptance at a certain point in a Nanowire Heterojunction Bipolar Transistor Solar Cell.
# It makes use of the software S4: Liu, V.; Fan, S. S4 : A Free Electromagnetic Solver for Layered Periodic Structures. Computer Physics Communications 2012, 183 (10), 2233–2244. https://doi.org/10.1016/j.cpc.2012.04.026.
# S4 was used as a python package in python3 running on Ubuntu using the program and instructions found at https://github.com/phoebe-p/S4, archived here: https://web.archive.org/web/20210614111957/https://github.com/phoebe-p/S4
# It requires some inputfiles, e.g., solar spectra or dielectric functions, which are not included but referenced where they are needed.
# This code, adapted to each case, was used in the simulations for the publication: 
# Zehender, M. H.; Chen, Y.; Barrigón, E.; Martí, A.; Borgström, M. T.; Antolín, E. Design Study of a Nanowire Three-Terminal Heterojunction Bipolar Transistor Solar Cell. To be published at the 48th IEEE Photovoltaic Specialists Conference (PVSC 48).
# Author: Marius H. Zehender
import S4
import time
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

tic = time.perf_counter()


def Printfigure(plotname,resdpi=300):
    if outformat == 'png':
        plt.savefig(plotname+'.png',dpi=resdpi,transparent=True)
    elif outformat == 'eps':
        plt.savefig(plotname+'.eps',transparent=True)
    elif outformat == 'svg':
        plt.savefig(plotname+'.svg',transparent=True)
    else:
        print('no valid outformat: no plot saved')
outformat = 'svg' #'eps' or 'png' or 'svg'
emitter_base_material='GaAs' #'GaInP', 'GaAs' or 'InP'
collector_material='InGaAs'  #'GaAs', 'InP' or 'InGaAs'


def S4_nanowire(NW_diam,ITOe_thickness,top_w,ITOb_thickness,collector_thickness,substrate_thickness,ITO_is_absorber = False,isoptimization=False):
    echarge=1.60217662e-19 # electron charge
    mypitch=500
    myLattice=((mypitch*1,0),(mypitch*1/2,mypitch*np.sqrt(3)/2)) # hexagonal lattice 
    myNumBasis=200#200
    wvl=np.arange(300,1700,2)#(300,1700,1)
    step_w=1  # step size for absorption calculation in nm (1)
    emitter_thickness= top_w*2/3
    base_thickness=top_w*1/3
    depths=np.append(np.arange(0,ITOe_thickness+emitter_thickness+base_thickness+collector_thickness,step_w),np.arange(ITOe_thickness+emitter_thickness+base_thickness+collector_thickness+(100*step_w),ITOe_thickness+emitter_thickness+base_thickness+collector_thickness+substrate_thickness,100*step_w)) # array with depth value
    S = S4.New(Lattice=myLattice, NumBasis=myNumBasis) # initiate simulation
    # initiate materials. epsilons will change later... 
    S.SetMaterial(Name = "vacuum", Epsilon = 1+0j)
    S.SetMaterial(Name = "low_gap", 	Epsilon = 1+0j)
    S.SetMaterial(Name = "high_gap",Epsilon = 1+0j)
    S.SetMaterial(Name = "ITO2", 	Epsilon = 1+0j)
    # set layers
    S.AddLayer('Layer_Above', 0.000000, 'vacuum')
    S.AddLayer('ITO_top', ITOe_thickness, 'ITO2')
    S.AddLayer('emitter', emitter_thickness, 'vacuum')
    S.AddLayer('base1', (base_thickness-ITOb_thickness)/2, 'vacuum')
    S.AddLayer('base2', ITOb_thickness, 'ITO2')
    S.AddLayer('base3', (base_thickness-ITOb_thickness)/2, 'vacuum')
    S.AddLayer('collector', collector_thickness, 'vacuum')
    # manipulate layers
    S.SetRegionCircle('emitter','high_gap', (0,0), NW_diam/2)#(myLattice[0][0]/2,myLattice[1][1]/2)
    S.SetRegionCircle('base1','high_gap', (0,0), NW_diam/2)
    S.SetRegionCircle('base2','high_gap', (0,0), NW_diam/2)
    S.SetRegionCircle('base3','high_gap', (0,0), NW_diam/2)
    S.SetRegionCircle('collector','low_gap', (0,0), NW_diam/2)
    S.AddLayer('substrate', substrate_thickness, 'low_gap')
    S.AddLayerCopy('Layer_Below', 0.000000, 'Layer_Above')
    
    S.SetExcitationPlanewave( 
    	IncidenceAngles=(
                    0, # polar angle in [0,180)
                    0  # azimuthal angle in [0,360)
            ),
            sAmplitude = 1,
            pAmplitude = 0,
            Order = 0
    )
    
    # load material dielectric functions
    frequency=np.divide(1,(wvl))
    real_eps_vacuum=np.ones(len(frequency))
    imag_eps_vacuum=np.zeros(len(frequency))
    # GaAs Adachi from refractiveindex.info doi.org/10.1063/1.343580
    GaAs_lam_n_k = np.loadtxt("GaAs_adachi_n_k_01.txt", delimiter="\t", unpack=False)
    GaAs_wvl=GaAs_lam_n_k[:,0]*1000
    GaAs_n=GaAs_lam_n_k[:,1]
    GaAs_k=GaAs_lam_n_k[:,2]
    GaAs_freq=np.divide(1,(GaAs_wvl))
    interp = interpolate.interp1d(GaAs_freq, GaAs_n,fill_value="extrapolate")
    GaAs_n_int = interp(frequency)
    interp = interpolate.interp1d(GaAs_freq, GaAs_k,fill_value="extrapolate")
    GaAs_k_int = interp(frequency)
    real_eps_GaAs=GaAs_n_int**2-GaAs_k_int**2
    imag_eps_GaAs=2*GaAs_n_int*GaAs_k_int
    
    # InGaAs Adachi from refractiveindex.info doi.org/10.1063/1.343580
    InGaAs_lam_n_k = np.loadtxt("InGaAs_adachi_n_k_01.txt", delimiter="\t", unpack=False)
    InGaAs_wvl=InGaAs_lam_n_k[:,0]*1000
    InGaAs_n=InGaAs_lam_n_k[:,1]
    InGaAs_k=InGaAs_lam_n_k[:,2]
    InGaAs_freq=np.divide(1,(InGaAs_wvl))
    interp = interpolate.interp1d(InGaAs_freq, InGaAs_n,fill_value="extrapolate")
    InGaAs_n_int = interp(frequency)
    interp = interpolate.interp1d(InGaAs_freq, InGaAs_k,fill_value="extrapolate")
    InGaAs_k_int = interp(frequency)
    real_eps_InGaAs=InGaAs_n_int**2-InGaAs_k_int**2
    imag_eps_InGaAs=2*InGaAs_n_int*InGaAs_k_int
    
    # InP Adachi from refractiveindex.info doi.org/10.1063/1.343580
    InP_lam_n_k = np.loadtxt("InP_adachi_n_k.txt", delimiter="\t", unpack=False)
    InP_wvl=InP_lam_n_k[:,0]*1000
    InP_n=InP_lam_n_k[:,1]
    InP_k=InP_lam_n_k[:,2]
    InP_freq=np.divide(1,(InP_wvl))
    interp = interpolate.interp1d(InP_freq, InP_n,fill_value="extrapolate")
    InP_n_int = interp(frequency)
    interp = interpolate.interp1d(InP_freq, InP_k,fill_value="extrapolate")
    InP_k_int = interp(frequency)
    real_eps_InP=InP_n_int**2-InP_k_int**2
    imag_eps_InP=2*InP_n_int*InP_k_int    
    
    #GaInP
    # GaInP Schubert from refractiveindex.info doi.org/10.1063/1.358632
    GaInP_lam_n_k = np.loadtxt("GaInP_Schubert_n_k.txt", delimiter="\t", unpack=False)
    GaInP_wvl=GaInP_lam_n_k[:,0]*1000
    
    GaInP_n=GaInP_lam_n_k[:,1]
    GaInP_k=GaInP_lam_n_k[:,2]
    GaInP_freq=np.divide(1,(GaInP_wvl))
    interp = interpolate.interp1d(GaInP_freq, GaInP_n,fill_value="extrapolate")
    GaInP_n_int = interp(frequency)
    interp = interpolate.interp1d(GaInP_freq, GaInP_k,fill_value="extrapolate")
    GaInP_k_int = interp(frequency)
    real_eps_GaInP=GaInP_n_int**2-GaInP_k_int**2
    imag_eps_GaInP=2*GaInP_n_int*GaInP_k_int
   
    # ITO from refractiveindex.info ("T. A. F. König, P. A. Ledin, J. Kerszulis, M. A. Mahmoud; M. A. El-Sayed, J. R. Reynolds and V. V. Tsukruk. Electrically tunable plasmonic behavior of nanocube-polymer nanomaterials induced by a redox-active electrochromic polymer, <a href=\"https://doi.org/10.1021/nn501601e\"><i>ACS Nano</i> <b>8</b>, 6182-6192 (2014)</a>")
    ITO_lam_n_k= np.loadtxt("ITO_n_k.txt", delimiter="\t", unpack=False)
    ITO_wvl=ITO_lam_n_k[:,0]*1000
    ITO_n=ITO_lam_n_k[:,1]
    ITO_k=ITO_lam_n_k[:,2]
    ITO_freq=np.divide(1,(ITO_wvl))/1000
    interp = interpolate.interp1d(ITO_freq, ITO_n,fill_value="extrapolate")
    ITO_n_int = interp(frequency)  
    interp = interpolate.interp1d(ITO_freq, ITO_k,fill_value="extrapolate")
    ITO_k_int = interp(frequency)  
    #plt.plot(ITO_freq, ITO_k, 'o', frequency, ITO_k_int, '-')
    #plt.show()
    # epsilon=(n+ik)²=n²-k²+i2nk
    real_eps_ITO=ITO_n_int**2-ITO_k_int**2
    imag_eps_ITO=2*ITO_n_int*ITO_k_int
    #plt.plot(frequency, real_eps_old, 'o', frequency, real_eps_ITO, '.')
    #plt.show()
    #plt.plot(frequency, imag_eps_old, 'o', frequency, imag_eps_ITO, '.')
    #plt.show()
    
    # solar spectrum (ASTM International. G173-03(2012) Standard Tables for Reference Solar Spectral Irradiances: Direct Normal and Hemispherical on 37° Tilted Surface. West Conshohocken, PA; ASTM International, 2012. doi: https://doi.org/10.1520/G0173-03R12)
    spectra=np.loadtxt("ASTMG173.csv", delimiter=",", unpack=False)
    spec_wvln=spectra[:,0]
    spec_power=spectra[:,2] # 1: ETR, 2: Global, 3: Direct [W/m2/nm]
    spec_freq=np.divide(1,(spec_wvln))
    h = 6.62607015e-34#;  %   [J s]
    c = 299792458#; [m/s] 
    spec_photons1=np.multiply(spec_power,spec_wvln*1e-9/(h*c));
    interp = interpolate.interp1d(spec_freq, spec_photons1,fill_value="extrapolate")
    spec_photons=interp(frequency)
    
    # initiate arrays for data
    wvln = np.zeros(len(frequency))
    reflection=np.zeros(len(frequency))
    absorption_GaInP=np.zeros(len(frequency))
    trans_GaInP=np.zeros(len(frequency))
    probe_flux=np.zeros((len(frequency),len(depths)))
    absorbed_flux=np.zeros((len(frequency),len(depths)))
    generation=np.zeros((len(frequency),len(depths)))

    for i in range(0, len(frequency)):
        wvln[i]=1/frequency[i]
        #print('%d' % i);
        freq = frequency[i];
        S.SetFrequency(freq)
        S.SetMaterial('vacuum', real_eps_vacuum[i]+ imag_eps_vacuum[i]*1j);
        if collector_material=='GaAs':
            S.SetMaterial('low_gap',    real_eps_GaAs[i]+ imag_eps_GaAs[i]*1j);
        elif collector_material=='InP':
            S.SetMaterial('low_gap',    real_eps_InP[i]+ imag_eps_InP[i]*1j);
        elif collector_material=='InGaAs':
            S.SetMaterial('low_gap',    real_eps_InGaAs[i]+ imag_eps_InGaAs[i]*1j);
        if emitter_base_material=='GaInP':
            S.SetMaterial('high_gap',real_eps_GaInP[i]+ imag_eps_GaInP[i]*1j);
        elif emitter_base_material=='GaAs':
            S.SetMaterial('high_gap',real_eps_GaAs[i]+ imag_eps_GaAs[i]*1j);
        elif emitter_base_material=='InP':
            S.SetMaterial('high_gap',real_eps_InP[i]+ imag_eps_InP[i]*1j);
        if ITO_is_absorber:
            S.SetMaterial('ITO2',   real_eps_ITO[i]+ imag_eps_ITO[i]*1j);
        else:
            S.SetMaterial('ITO2',   real_eps_ITO[i]);
        incidence_flux, reflection_flux_vacuum = S.GetPowerFlux('Layer_Above', 0.000000)
        reflection_flux_vacuum = (-1) * reflection_flux_vacuum / incidence_flux;
        transmission_flux,trash = S.GetPowerFlux('Layer_Below', 0.000000)
        class Cls:
            pass
        out=Cls()
        if (not isoptimization):
            for k in range(0,len(depths)):
                if depths[k]<=ITOe_thickness+emitter_thickness+(base_thickness-ITOb_thickness)/2+ITOb_thickness+(base_thickness-ITOb_thickness)/2+collector_thickness+substrate_thickness:
                    layer='substrate'
                    layer_offset=ITOe_thickness+emitter_thickness+(base_thickness-ITOb_thickness)/2+ITOb_thickness+(base_thickness-ITOb_thickness)/2+collector_thickness
                if depths[k]<=ITOe_thickness+emitter_thickness+(base_thickness-ITOb_thickness)/2+ITOb_thickness+(base_thickness-ITOb_thickness)/2+collector_thickness:
                    layer='collector'
                    layer_offset=ITOe_thickness+emitter_thickness+(base_thickness-ITOb_thickness)/2+ITOb_thickness+(base_thickness-ITOb_thickness)/2
                if depths[k]<=ITOe_thickness+emitter_thickness+(base_thickness-ITOb_thickness)/2+ITOb_thickness+(base_thickness-ITOb_thickness)/2:
                    layer='base3'
                    layer_offset=ITOe_thickness+emitter_thickness+(base_thickness-ITOb_thickness)/2+ITOb_thickness
                if depths[k]<=ITOe_thickness+emitter_thickness+(base_thickness-ITOb_thickness)/2+ITOb_thickness:
                    layer='base2'
                    layer_offset=ITOe_thickness+emitter_thickness+(base_thickness-ITOb_thickness)/2
                if depths[k]<=ITOe_thickness+emitter_thickness+(base_thickness-ITOb_thickness)/2:
                    layer='base1'
                    layer_offset=ITOe_thickness+emitter_thickness
                if depths[k]<=ITOe_thickness+emitter_thickness:
                    layer='emitter'
                    layer_offset=ITOe_thickness
                if depths[k]<=ITOe_thickness:
                    layer='ITO_top'
                    layer_offset=0           
                probe_flux_f,probe_flux_r = S.GetPowerFlux(layer, depths[k]-layer_offset)
                probe_flux[i,k]=probe_flux_f.real/incidence_flux.real
                absorbed_flux[i,k]=probe_flux[i,k-1]-probe_flux[i,k]
                # corrections for numerical errors
                if (k == 0):
                    absorbed_flux[i,k]=0
                if (absorbed_flux[i,k]/absorbed_flux[i,k-1] <= 0.20):
                    absorbed_flux[i,k-1]=0
                if (abs(absorbed_flux[i,k]) <= 0):
                    absorbed_flux[i,k]=0
                generation[i,k]=absorbed_flux[i,k]*spec_photons[i]
        reflection[i]=reflection_flux_vacuum.real
        trans_GaInP_f,trans_GaInP_r = S.GetPowerFlux('collector', 0.000000)
        trans_GaInP[i]=(trans_GaInP_f.real/incidence_flux.real)
        absorption_GaInP[i]=1-reflection[i]-trans_GaInP[i]
        
        idx_emitter = (np.abs(depths - ITOe_thickness)).argmin()
        idx_base = (np.abs(depths - (ITOe_thickness+emitter_thickness))).argmin()
        idx_collector = (np.abs(depths - (ITOe_thickness+emitter_thickness+base_thickness))).argmin()
        idx_substrate = (np.abs(depths - (ITOe_thickness+emitter_thickness+base_thickness+collector_thickness))).argmin()
        emitter_abs=np.zeros(len(wvln))
        base_abs=np.zeros(len(wvln))
        base_abs_top=np.zeros(len(wvln))
        base_abs_bot=np.zeros(len(wvln))
        collector_abs=np.zeros(len(wvln))
        for m in range(0, len(wvln)):
            emitter_abs[m]=np.nansum(absorbed_flux[m,idx_emitter:idx_base])
            base_abs[m]=np.nansum(absorbed_flux[m,idx_base:idx_collector])
            base_abs_top[m]=np.nansum(absorbed_flux[m,idx_base:idx_base+int((idx_collector-idx_base)/2)])
            base_abs_bot[m]=np.nansum(absorbed_flux[m,idx_base+int((idx_collector-idx_base)/2):idx_collector])
            collector_abs[m]=np.nansum(absorbed_flux[m,idx_collector:idx_substrate])
    wvln_step=wvln[1]-wvln[0]
    out.objective=(-np.nansum(emitter_abs*spec_photons)-np.nansum(base_abs*spec_photons)-np.nansum(collector_abs*spec_photons))*wvln_step*echarge
    out.emitter_abs=emitter_abs
    out.base_abs=base_abs
    out.base_abs_top=base_abs_top
    out.base_abs_bot=base_abs_bot
    out.collector_abs=collector_abs
    out.NumBasis=myNumBasis
        
#        plt.plot(out.wvln,emitter_abs)
#        plt.plot(out.wvln,base_abs)
#        plt.plot(out.wvln,collector_abs)
    if (not isoptimization):
        out.graphic=S.OutputStructurePOVRay(Filename = 'out14.pov')
    out.generation=generation
    out.probe_flux=probe_flux
    out.absorbed_flux=absorbed_flux
    out.NW_diam=NW_diam
    out.ITOe_thickness=ITOe_thickness
    out.top_w=top_w
    out.ITOb_thickness=ITOb_thickness
    out.emitter_thickness=emitter_thickness
    out.base_thickness=base_thickness
    out.collector_thickness=collector_thickness
    out.substrate_thickness=substrate_thickness
    out.spec_photons=spec_photons
    out.wvln=wvln
    out.depths=depths
    return out


def S4_call(x,isoptimization_):
    return S4_nanowire(NW_diam=x[0],ITOe_thickness=x[1],top_w=x[2],ITOb_thickness=2*x[1],collector_thickness=1500,substrate_thickness=500000,ITO_is_absorber = False,isoptimization=isoptimization_)

def objective(x):
    return S4_call(x,True).objective


for diameter in [350]:
    for top_len in [1200]:
        for ito_len in [50]:
            x=[diameter, ito_len, top_len]   
            #obj=[]
            #sweepvar=np.linspace(50, 150, num=20)
            #for q in sweepvar:
            #    x[1]=q     
            #    out=S4_call(x,False)
            #    obj.append(out.objective)
            #plt.plot(sweepvar,obj)
            out=S4_call(x,False)
            for i in range(0,int(len(out.wvln)/10)):
                plt.plot(out.depths,out.absorbed_flux[i*10,:])
            plt.ylim(0, 0.002)
            plt.xlim(0,out.ITOe_thickness+out.emitter_thickness+out.base_thickness+out.collector_thickness)
            plt.show()
            plt.plot(out.wvln,out.emitter_abs,color='b',label='Emitter')
            plt.plot(out.wvln,out.base_abs,color='g',label='Base')
            plt.plot(out.wvln,out.collector_abs,color='r',label='Collector')
            plt.ylim(0, 1)
            plt.xlim(out.wvln[0],out.wvln[-1])
            #plt.plot(out.wvln,out.emitter_abs+out.base_abs+out.collector_abs,color='k')
            plt.legend()
            Printfigure(emitter_base_material+'_'+collector_material+'_obj'+str(int(abs(out.objective)))+'_NumBasis'+str(out.NumBasis)+'_dia'+str(x[0])+'_gainp'+str(top_len)+'_ito'+str(ito_len)+'_wvl_abs')
            plt.show()
            toc = time.perf_counter()
            print('diameter= %1.0f nm' % x[0])
            print('length itoe= %1.0f nm' % x[1])
            print('length top= %1.0f nm' % x[2])
            print('NumBasis= %1.0f' % out.NumBasis)
            print('time= %e s' % (toc-tic))
            print('current_density [A m⁻2]= %1.3f' % -out.objective)
            
            out.graphic
            #save txt files
            array=list(zip(out.wvln,out.emitter_abs,out.base_abs_top,out.base_abs_bot,out.base_abs,out.collector_abs,out.emitter_abs*out.spec_photons,out.base_abs*out.spec_photons,out.collector_abs*out.spec_photons))
            np.savetxt(emitter_base_material+'_'+collector_material+'_obj'+str(int(abs(out.objective)))+'_NumBasis'+str(out.NumBasis)+'_dia'+str(x[0])+'_gainp'+str(top_len)+'_ito'+str(ito_len)+'_wvl_abs_absPhotons.txt', array)
            array=list(zip(out.wvln,out.emitter_abs+out.base_abs_top))
            np.savetxt(emitter_base_material+'_'+collector_material+'_obj'+str(int(abs(out.objective)))+'_NumBasis'+str(out.NumBasis)+'_dia'+str(x[0])+'_gainp'+str(top_len)+'_ito'+str(ito_len)+'_wvl_abstop.txt', array)
            array=list(zip(out.wvln,out.base_abs_bot+out.collector_abs))
            np.savetxt(emitter_base_material+'_'+collector_material+'_obj'+str(int(abs(out.objective)))+'_NumBasis'+str(out.NumBasis)+'_dia'+str(x[0])+'_gainp'+str(top_len)+'_ito'+str(ito_len)+'_wvl_absbot.txt', array)
            array=out.absorbed_flux
            np.savetxt(emitter_base_material+'_'+collector_material+'_obj'+str(int(abs(out.objective)))+'_NumBasis'+str(out.NumBasis)+'_dia'+str(x[0])+'_gainp'+str(top_len)+'_ito'+str(ito_len)+'_Absorption_lam_depth.txt',np.transpose(array))# np.transpose(array))#, delimiter='/t')
            array=out.generation
            np.savetxt(emitter_base_material+'_'+collector_material+'_obj'+str(int(abs(out.objective)))+'_NumBasis'+str(out.NumBasis)+'_dia'+str(x[0])+'_gainp'+str(top_len)+'_ito'+str(ito_len)+'_Generation_lam_depth.txt',np.transpose(array))# np.transpose(array))#, delimiter='/t')
            array=[out.ITOe_thickness,out.ITOe_thickness,out.ITOe_thickness+out.emitter_thickness,out.ITOe_thickness+out.emitter_thickness+out.base_thickness,out.ITOe_thickness+out.emitter_thickness+out.base_thickness+out.collector_thickness]
            np.savetxt(emitter_base_material+'_'+collector_material+'_obj'+str(int(abs(out.objective)))+'_NumBasis'+str(out.NumBasis)+'_dia'+str(x[0])+'_gainp'+str(top_len)+'_ito'+str(ito_len)+'_window_emitter_base_collector_beginsAt.txt', array, delimiter='/n')
            array=out.depths
            np.savetxt(emitter_base_material+'_'+collector_material+'_obj'+str(int(abs(out.objective)))+'_NumBasis'+str(out.NumBasis)+'_dia'+str(x[0])+'_gainp'+str(top_len)+'_ito'+str(ito_len)+'_depth.txt', array)
            array=out.wvln
            np.savetxt(emitter_base_material+'_'+collector_material+'_obj'+str(int(abs(out.objective)))+'_NumBasis'+str(out.NumBasis)+'_dia'+str(x[0])+'_gainp'+str(top_len)+'_ito'+str(ito_len)+'_lambda.txt', array)
            array=(out.NW_diam,out.ITOe_thickness,out.emitter_thickness,out.base_thickness,out.collector_thickness,out.substrate_thickness,abs(out.objective))
            np.savetxt('obj'+str(abs(out.objective))+'_NumBasis'+str(out.NumBasis)+'_dia'+str(x[0])+'_gainp'+str(top_len)+'_ito'+str(ito_len)+'_dia_ITOe_emi_bas_col_subs_obj.txt', array)