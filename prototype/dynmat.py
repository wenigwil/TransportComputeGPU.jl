import numpy as np
import matplotlib.pyplot as plt
import itertools as it

def readqeifc2(ifc2fil,system):
    print 'Quantum Espresso mode:'
    print 'Reading fcs from ' + ifc2fil

    cell = system.cell
    masses = system.masses
    types = system.types
    natoms = system.natoms
    polar = False
    scell = []

    skip = 0
    with open(ifc2fil) as f:
        for line in f:
            if line.strip().startswith('T'):
                polar = True
                system.polar = True
                skip += 1
                print 'Reading dielectric tensor in cartesian basis.'
                for i in range(3):
                    system.diel[i] = np.fromstring(f.next(),dtype=float,sep=' ')
                    skip += 1
                print 'Reading Born charges in cartesian basis.'
                for i in range(natoms):
                    f.next() #skip a line
                    for j in range(3):
                        system.born[i][j] = np.fromstring(f.next(),dtype=float,sep=' ')
                        skip += 1
                break
            elif line.strip().startswith('F'):
                skip += 1
                print 'No Born charge information found.'
                break
            else:
                skip += 1

        scell = np.fromstring(f.next(),dtype=int,sep=' ')
        prod = np.prod(scell)
        numfc = 3*3*natoms*natoms
        fc = np.zeros((numfc,scell[0],scell[1],scell[2]))
        ind = np.zeros((numfc,4),dtype=int)
        m = 0
        print 'Reading harmonic force constants.'
        for line in f:
            # ind = (pol,pol,atom,atom) = (a,b,j,jp)
            ind[m] = np.fromstring(line,dtype=int,sep=' ')
            for n in range(prod):
                l = f.next()
                i = np.fromstring(l,dtype=int,sep=' ')[0] - 1
                j = np.fromstring(l,dtype=int,sep=' ')[1] - 1
                k = np.fromstring(l,dtype=int,sep=' ')[2] - 1
                fc[m][i][j][k] = np.fromstring(l,dtype=float,sep=' ')[3]
            m += 1

    #convert qe ifc2 units from ryd/bohr^2 to ev/nm^2
    ryd2ev = 13.605698066
    #bohr2nm = 0.052917725
    bohr2nm = 0.052917721092 #shengbte
    fc *= ryd2ev/bohr2nm/bohr2nm

    if (system.polar):
        #apply sum rule on born charges
        for a in xrange(3):
            for b in xrange(3):
                asr = 0.0
                for j in xrange(natoms):
                    asr += system.born[j][a][b]
                for j in xrange(natoms):
                    system.born[j][a][b] -= asr/natoms

    #Reshape ifc2 file from (m,lx,ly,lz) to (lx,ly,lz,j,a,jp,b) (standard) form
    # lx,ly,lz are supcercell indices
    fc_std = np.zeros((scell[0],scell[1],scell[2],natoms,3,natoms,3))
    for m in xrange(numfc):
        a = ind[m][0] - 1 #pol
        b = ind[m][1] - 1 #pol
        j = ind[m][2] - 1 #atom
        jp = ind[m][3] - 1 #atom

        for lx in xrange(scell[0]):
            for ly in xrange(scell[1]):
                for lz in xrange(scell[2]):
                    fc_std[lx][ly][lz][j][a][jp][b] = fc[m][lx][ly][lz]


    #Enforce acoustic sumrule
    for a in xrange(3):
        for b in xrange(3):
            for j in xrange(natoms):
                fc_std[0][0][0][j][a][j][b] -= np.sum(fc_std[:,:,:,j,a,:,b])

    return fc_std


def formdynmat(fc,system,qpts):

    cell = system.cell
    numq = len(qpts)
    natoms = len(cell[1])
    scell = system.scell
    numfc = 3*3*natoms*natoms
    numcells = np.prod(scell)

    #apply non-analytic correction for polar materials
    fcnac_std = np.zeros((numq,natoms,3,natoms,3))

    if (system.polar):
        print 'Calculating NAC to ifc2.'

        qshort = np.zeros((numq,3))
        for iq in xrange(0,numq): #q-pts
            qshort[iq] = qpts[iq]
            tmp1 = np.linalg.norm(qshort[iq])
            for ix,iy,iz in it.product(xrange(-2,3),xrange(-2,3),xrange(-2,3)):
                qp = qpts[iq] + ix*system.rlvec[0] + iy*system.rlvec[1] + iz*system.rlvec[2]
                tmp2 = np.linalg.norm(qp)
                if tmp2 < tmp1:
                    tmp1 = tmp2
                    qshort[iq] = qp

        for iq in xrange(0,numq): #q-pts
            q = qpts[iq]
            #q = qshort[iq]

            fac = 0.0
            #if not (q[0] == 0.0 and q[1] == 0.0 and q[2] == 0.0):
            if np.linalg.norm(q) != 0.0:
                fac = 1.0/np.dot(q,np.dot(system.diel,q))
                for j in xrange(natoms):
                    qZj = np.dot(system.born[j],q)
                    for jp in xrange(natoms):
                        qZjp = np.dot(system.born[jp],q)
                        fcnac_std[iq][j,:,jp,:] = np.outer(qZj,qZjp)
                fcnac_std[iq] *= fac

        fcnac_std *= 4.0*np.pi/system.vol/numcells

    print 'Forming dynamical matrix.'
    dynmat = np.zeros((numq,natoms*3,natoms*3),dtype=complex)
    a1 = cell[0][0]
    a2 = cell[0][1]
    a3 = cell[0][2]

    sc1 = scell[0]
    sc2 = scell[1]
    sc3 = scell[2]

    #build supercell lattice vectors
    sclvec = np.multiply(system.scell,system.lvec)

    #build vectors to neighboring supercells in a giant 5x5x5 super-supercell (ssc)
    rssc = np.zeros((124,4)) #the last index is to hold a scalar
    cnt = 0
    for i,j,k in it.product(xrange(-2,3),xrange(-2,3),xrange(-2,3)):
        if i == 0 and j == 0 and k == 0:
            continue
        for a in xrange(3):
            rssc[cnt][a] = i*sclvec[0][a] + j*sclvec[1][a] + k*sclvec[2][a]
        rssc[cnt][3] = 0.5*np.dot(rssc[cnt][0:3],rssc[cnt][0:3])
        cnt += 1

    r_ssc = np.zeros(3)
    t = np.zeros(3)
    eps = 1.0e-6
    for iat, jat in it.product(xrange(natoms),xrange(natoms)):
        totwght = 0.0
        massfac = 1.0/np.sqrt(system.masses[system.types[iat]-1]*system.masses[system.types[jat]-1])
        for m1,m2,m3 in it.product(xrange(-2*sc1,2*sc1+1),xrange(-2*sc2,2*sc2+1),xrange(-2*sc3,2*sc3+1)):
            t = np.dot([m1,m2,m3],system.lvec)
            r_ssc = t + np.dot(cell[1][iat] - cell[1][jat],system.lvec)

            wght = 0.0
            nreq = 1.0
            j = 0
            for ir in xrange(124):
                ck = np.dot(r_ssc,rssc[ir,0:3]) - rssc[ir][3]
                if ck > eps:
                    j = 1
                    continue
                if abs(ck) < eps:
                    nreq += 1.0
            if j == 0:
                wght = 1.0/nreq
            if wght > 0.0:
                t1 = (m1+1)%sc1
                if t1 <= 0:
                    t1 += sc1

                t2 = (m2+1)%sc2
                if t2 <= 0:
                    t2 += sc2

                t3 = (m3+1)%sc3
                if t3 <= 0:
                    t3 += sc3

                t1 -= 1
                t2 -= 1
                t3 -= 1
                for iq in xrange(0,numq): #q-pts
                    q = qpts[iq]
                    qt = np.dot(q,t)
                    fac = np.exp(-1.0j*qt)*wght*massfac
                    for a in xrange(3):
                        idim = iat*3 + a
                        for b in xrange(3):
                            jdim = jat*3 + b
                            nac = fcnac_std[iq][iat][a][jat][b]
                            dynmat[iq][idim][jdim] += (fc[t1][t2][t3][iat][a][jat][b] + 0*nac)*fac
            totwght += wght
    #enforce hermiticity
    for iq in xrange(0,numq): #q-pts
        dynmat[iq] = 0.5*(dynmat[iq] + np.transpose(np.conj(dynmat[iq])))
    np.save('./pybte.dynmat', dynmat)

    return dynmat

def genphonons(system,numq):

    print 'Diagonalizing dynamical matrix.'

    #conversion factor of 2pi*sqrt(ev/amu/nm^2) to THz
    evovamu = (1.60217733e-19)/(1.6605402e-27)
    toTHz = np.sqrt(evovamu*pow(1.0e9,2)*pow(1.0e-12,2))/2/np.pi #(nm.Thz)^2

    #numq = len(qpts)
    size = system.nbrnchs
    omega = np.zeros((numq,system.nbrnchs))

    evecs = np.zeros((numq,size,size),dtype=complex)
    evals = np.zeros((numq,size),dtype=complex)

    dynmat = np.load('./pybte.dynmat.npy')

    for iq in xrange(numq):
        evals[iq],evecs[iq] = np.linalg.eigh(dynmat[iq])
        omega[iq] = np.sort(np.copysign(np.sqrt(np.abs(np.real(evals[iq]))),np.real(evals[iq])))
        #omega[iq] = np.copysign(np.sqrt(np.abs(np.real(evals[iq]))),np.real(evals[iq]))

    omega *= toTHz

    np.savetxt('./pybte.phfreq', omega, delimiter='\t')
    #np.savetxt('./pybte.phfreq', omega*2*np.pi, delimiter='\t')#jesus' units

    return omega

def plotphonons(fc,system,spfil):

    #number of points along each segment
    N =21

    def readspecialpoints(spfil):
        sppts = []
        splabel = []
        with open(spfil,'r') as f:
            for line in f:
                if not line.startswith('#'):
                    sppts.append((np.fromstring(line,sep=' ',dtype=float))[0:3])
                    splabel.append((line.split())[-1])

        return np.array(sppts),np.array(splabel)

    #dispersion plot
    def plotdisp(omega,kpts,splabel,spind):
        print 'Plotting phonons.'
        numsp = len(labels)
        shift = 0.0
        ticks = []
        symbols = [s.replace('Gamma', '$\Gamma$') for s in splabel]
        for i in xrange(0,numsp-1):
            k1 = kpts[spind[i]]
            k2 = kpts[spind[i+1]]

            dist = np.linalg.norm(k2-k1)

            numpts = spind[i+1] - spind[i] + 1
            kseg = np.linspace(0.0,dist,numpts)
            dispseg = omega[spind[i]:spind[i+1]+1]

            plt.plot(shift+kseg,dispseg,'b',linewidth=1.0)
            plt.axvline(x=shift,color='k',linewidth=1.5)
            ticks.append(shift)
            shift += kseg[-1]
        ticks.append(shift)
        plt.xticks(ticks, symbols)
        plt.yticks(fontsize = 16)
        plt.xticks(fontsize = 16)

        plt.axvline(x=shift,color='k',linewidth=1.5)
        plt.axhline(y=0.0,color='k',linewidth=1.5)

        plt.xlim([0.0,shift])
        plt.ylim(ymin=0.0)

        plt.xlabel('wavevector',fontsize=18)
        plt.ylabel('frequency (THz)',fontsize=18)

        plt.savefig('./pybte.disp.pdf',dpi=72)
        plt.clf()


    #read special points from file
    pts,labels = readspecialpoints(spfil)

    #make q-mesh
    size = (len(pts)-1)*N
    qpts = np.zeros((size+1,3))

    cnt = 0
    for isp in xrange(len(pts)-1):
        seg  =  np.array(zip( np.linspace(pts[isp][0],pts[isp+1][0],N,endpoint=False),\
        np.linspace(pts[isp][1],pts[isp+1][1],N,endpoint=False),\
        np.linspace(pts[isp][2],pts[isp+1][2],N,endpoint=False) ))
        for i in xrange(N):
            qpts[cnt] = seg[i]
            cnt += 1
    qpts[-1] = pts[-1]

    #special points index in final qpts array
    spind = [0]*N
    for isp in xrange(len(pts)):
        spind[isp] = isp*N

    #test: load q-mesh from shengbte output and use here
    #qpts = np.loadtxt('./BTE.qpoints', usecols=(2,3,4))

    #call dynmat for special qpts
    dynmat = formdynmat(fc,system,np.dot(qpts,system.rlvec))

    #call getphonons for special qpts
    omega = genphonons(system,len(qpts))

    #plot special phonons
    plotdisp(omega,qpts,labels,spind)
