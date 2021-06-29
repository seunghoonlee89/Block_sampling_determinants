import os, sys
import numpy
from pyscf.tools import fcidump 

Identity_FCIDUMP = "Identity"
def writeIdentityIntegralFile(nmo, nelec, filename):
    with open(filename, 'w') as fout:
        if not isinstance(nelec, (int, numpy.number)):
            ms = abs(nelec[0] - nelec[1])
            nelec = nelec[0] + nelec[1]
        fcidump.write_head(fout, nmo, nelec, ms=ms)
        for i in range(20):
            fout.write ("%20.15f  %d  %d  %d  %d\n"%(1/31, i+1, i+1, 0, 0))
        fout.write ("%20.15f  %d  %d  %d  %d\n"%(0.0, 0, 0, 0, 0))

def Write_input_for_MPS_to_CI(nelectron, symmetry, thresh, target, spin2, FCIDUMP, scratch):
    filename='extractCI.conf'
    with open(filename, 'w') as fout:
        fout.write('%s %d\n'%('nelec', nelectron))
        fout.write('%s\n'%('hf_occ canonical'))
        fout.write('%s\n'%('spin %d'%(spin2)))
        fout.write('%s\n'%('orbitals %s'%(FCIDUMP)))
        fout.write('%s %s\n'%('symmetry', symmetry))
        fout.write('%s\n'%('irrep 1'))
        fout.write('\n')
        fout.write('%s\n'%('schedule'))
        fout.write('%s\n'%('0 3000 1e-5 5e-5'))
        fout.write('%s\n'%('2 3000 1e-8 0.0'))
        fout.write('%s\n'%('end'))
        fout.write('\n')
        fout.write('%s\n'%('noreorder'))
        fout.write('%s\n'%('warmup local_2site'))
        fout.write('%s\n'%('maxiter 6'))
        fout.write('%s\n'%('onedot'))
        fout.write('%s\n'%('sweep_tol 1e-07'))
        fout.write('%s %s\n'%('prefix',scratch))
        fout.write('%s\n'%('fullrestart'))
        fout.write('%s\n'%('sto_pt_restart'))
        fout.write('%s %15.12f\n'%('pt_tol', thresh))
        fout.write('\n')
        fout.write('%s\n'%('nroots 1'))
        fout.write('%s\n'%('outputlevel 3'))
        fout.write('%s %d\n'%('targetState', target))

    filename2='extractCI_init.conf'
    with open(filename2, 'w') as fout:
        fout.write('%s %d\n'%('nelec', nelectron))
        fout.write('%s\n'%('hf_occ canonical'))
        fout.write('%s\n'%('spin %d'%(spin2)))
        fout.write('%s\n'%('orbitals %s'%(FCIDUMP)))
        fout.write('%s %s\n'%('symmetry', symmetry))
        fout.write('%s\n'%('irrep 1'))
        fout.write('\n')
        fout.write('%s\n'%('schedule'))
        fout.write('%s\n'%('0 5000 1e-5 5e-5'))
        fout.write('%s\n'%('2 5000 1e-8 0.0'))
        fout.write('%s\n'%('end'))
        fout.write('\n')
        fout.write('%s\n'%('noreorder'))
        fout.write('%s\n'%('warmup local_2site'))
        fout.write('%s\n'%('maxiter 6'))
        fout.write('%s\n'%('onedot'))
        fout.write('%s\n'%('sweep_tol 1e-07'))
        fout.write('%s %s\n'%('prefix',scratch))
        fout.write('%s\n'%('fullrestart'))
        fout.write('%s\n'%('sto_pt_nsamples 0'))
        fout.write('%s\n'%('sto_pt_Hsamples 1'))
        fout.write('%s %15.12f\n'%('pt_tol', thresh))
        fout.write('\n')
        fout.write('%s\n'%('nroots 1'))
        fout.write('%s %d\n'%('targetState', target))
        fout.write('%s\n'%('outputlevel 3'))
        fout.write('%s\n'%('mem 100 g'))


def extracting_CI_coeff_from_MPS(nelectron, spin2, thresh, target, scratch, FCIDUMP, symmetry="c1"):
    CITRIEEXE1 = '/home/slee89/opt/stackblocklatest/stackblock_stopt/TRIEdeterminant/CITRIE'
    CITRIEEXE2 = '/home/slee89/opt/Block_pt_stored/TRIEdeterminant/MPS2CI_det'
   
    Write_input_for_MPS_to_CI(nelectron, symmetry, thresh, target, spin2, FCIDUMP, scratch)
    
    #TODO: mc.CITRIE(cutoff)
    #      would make input of CITRIE and run CITRIE
    cmd3="%s extractCI_init.conf > extractCI_init.out"%(CITRIEEXE1)
    cmd4="%s extractCI.conf > extractCI.out"%(CITRIEEXE2)
    os.system(cmd3)
    os.system(cmd4)


if __name__ == '__main__':

    if len(sys.argv) == 9:
        nmo   = int(sys.argv[1])
        nelea = int(sys.argv[2])
        neleb = int(sys.argv[3])
        spin2 = int(sys.argv[4])
        thrsh = float(sys.argv[5])
        target= int(sys.argv[6])
        scr   = sys.argv[7]
        symm  = sys.argv[8]
    else:
        raise ValueError("""
            Usage: 
                python extracting_CI_from_stackblock.py (nmo) (nelec alpha) (nelec beta) (2*spin) (thresh) (targetstate) (scratchdir) (point group symmetry)
        """)

    writeIdentityIntegralFile(nmo, [nelea, neleb], Identity_FCIDUMP)
    extracting_CI_coeff_from_MPS(nelea+neleb, spin2, thrsh, target, scr, Identity_FCIDUMP, symm)
