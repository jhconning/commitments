# -*- coding: utf-8 -*-
"""
Solves for and analyzes renegotiation-proof commitment
contracts.  See https://jhconning.github.io/commitments

This version allows 2 period or 3 period version.

__author__ = "Jonathan Conning"
__author__ = "@author: "Karna Basu"
__email__ = "jonathan.conning@gmail.com"

"""
import numpy as np
from scipy.optimize import minimize, brentq, fsolve
import warnings
warnings.filterwarnings('ignore', 'The iteration is not making good progress')


class Contract(object):
    """ Base Class for both Monopoly and Competitive contracts  """
    def __init__(self, beta=0.8, y=None):         # constructor default params.
        self.beta  = beta                   # β-δ framework
        self.rho   = 1.2                   # 1/rho = elast of substn
        self.kappa  = 0  
        if y is None:
            self.y = np.array([100,100,100])
        else:
            self.y = y

    def __repr__(self):
        return 'Contract(beta=%s, periods=%s, y=%s)' % (self.beta, self.periods, self.y)

    def print_params(self):
        """ print out parameters """
        params = vars(self)
        for p in sorted(params):    # print attributes alphabetically
            print("{0:<7} : {1}".format(p,params[p]))

    def params(self):
        '''Helper function to return key params just to unclutter other functions'''
        btr = self.beta**(1/self.rho)
        return self.beta, self.rho, btr, self.kappa, self.y

    def u(self, ct):
        """ utility function """
        if self.rho==1:
            return np.log(ct)
        else:
            return  (ct**(1-self.rho))/(1-self.rho)

    def PV(self, c):
        """discounted present value of any stream c"""
        return c[0] + np.sum(c[1:])

    def PVU(self, c, beta):
        """discounted present utility value of any stream c
           We make beta a parameter to allow this express
           Zero or One self preferences"""
        return  self.u(c[0]) + beta * sum([self.u(ct) for ct in c[1:]] )

    def profit(self, c, y):
        """ present value of lender profits when exchanges c for y"""
        return  self.PV(y)-self.PV(c)
    
    def negPVU(self, c):
        """0 self negative present utility value of
        any stream c for minimization call"""
        return  - self.PVU(c, self.beta)

    def indif(self, ubar, beta):
        """ returns idc(c1) function from  u(c1, c2)=ubar for graphing
        indifference curves in c1-c2 space. If beta = 1, will describe
        self 0's preferences """
        if self.rho==1:
            def idc(c1):
                lnc2 = (ubar-np.log(c1))/beta
                return np.exp(lnc2)
        else:
            def idc(c1):
                return (((1-self.rho)/beta)*(ubar-self.u(c1)))\
                       **(1/(1-self.rho))

        return idc

    def isoprofit(self, prfbar, y):
        """ returns isoprf(c1) implied by profit(c1, c2)=y for graphing isoprofit lines in c1-c2 space.
        """
        def isoprf(c1):
            """isoprofit function isoprf(c1) """
            return np.sum(y[1:]) - prfbar - c1
        return isoprf

    def reneg(self, c):
        """ Renegotiated monopoly contract offered to period-1-self if
        c_0 is past but (c_1,c_2) now replaced by (cr_1, cr_2)"""

        beta, rho, btr, _, _ = self.params()
        if rho==1:
            lncr1 = (np.log(c[1])+beta*np.log(c[2])
                     -beta*np.log(beta))/(1+beta)
            cr1 = np.exp(lncr1)
        else:
            pu =  self.u(c[1]) + beta*self.u(c[2])
            a  =  (pu *(1-rho) )**(1/(1-rho))
            b = (1 + btr)**(1/(rho-1))
            cr1 = a*b
        cr2 = cr1*btr
        return np.array([c[0], cr1, cr2])
    
    def noreneg(self, c):
        """ no-renegotiation constraint; same for comp or monop
        used in optimization
        """
        beta, _, btr, kappa, _ = self.params()
        return (self.u(c[1]) + beta * self.u(c[2]) 
                - (1+btr) * self.u( (c[1]+c[2]-kappa)/(1+btr)) )

# for drawing no-renegotiation constraint
    def rpc(self, s):
        '''return points along the renegotiation-proof constraint
    corresponding to resources s sent into period 1
        '''  
        beta, _, btr, kap, _ = self.params()
        c1P = (s-kap)/(1+btr)
        ub = (1+btr)*self.u(c1P)
    
        def f(c1):
            return  self.u(c1) + beta*self.u(s-c1) - ub 

        c1rp = fsolve(f, [s/2])[0]
        return c1rp, s-c1rp   

    #def nrpline(C):
    #     '''return '''
    #     beta, rho, btr, kap = C.params()
    #     num = 150
    #     c1rp, c2rp  = np.zeros(num), np.zeros(num)
    #     for i, s in enumerate(np.arange(30, 30+num)):
    #         c1rp[i], c2rp[i] = C.rpc(s, k, beta, rho, upper)
            
    #     plt.plot(c1rp, c2rp, color='g')
    #     plt.xlim(0,120),  plt.ylim(0,120)
    #     plt.grid()
        



class Competitive(Contract):                    # build on contract class
    """ Initialize class for solving competitive equilibrium contracts  """
    def __init__(self, beta):
        super().__init__(beta)  # inherits Contract properties
        self.kappa  = 0           # cost of renegotiation

    def __repr__(self):
        return 'Competitive(beta=%s, y=%s)' % (self.beta, self.y)

    def fcommit(self):
        """competitive optimal full commitment contractwith period0 self
        from closed form solution for CRRA"""
        btr = self.beta**(1/self.rho)
        return np.array([1,btr,btr])*self.PV(self.y)/(1 + 2*btr)
    
    def kbar(self):
        '''Renegotiation cost necessary to sustain full commitment
        competitive contract'''
        rho = self.rho
        if (rho == 1):
            rho = 0.999  # cheap trick to deal with special case
        btr = self.beta ** (1 / rho)
        c1F = np.sum(self.y) * btr / (1 + 2 * btr)
        A = (2 - (1 + btr) * ((1 + self.beta) / (1 + btr))
             ** (1 / (1 - rho)))
        return A * c1F
    
    def ownsmooth(self):
        """contract if have to smooth by oneself without contract
        from closed form solution for CRRA with kappa=0"""
        beta, rho = self.beta, self.rho
        btr = beta**(1/rho)
        Y = self.PV(self.y)
        lam = (beta+btr)/(1+btr)
        ltr = lam**(1/rho)
        c0 = Y/( 1+(1+btr)*ltr )
        c1 = (Y-c0)/(1+btr)
        return np.array( [c0, c1, btr*c1])

    def creneg(self, c):
        """ Renegotiated competitive contract offered to period-1-self if
        c_0 is past but (c_1,c_2) now replaced by (cr_1, cr_2)"""
        beta, rho = self.beta, self.rho
        btr = beta**(1/rho)
        cr1 = (c[1]+c[2])/(1+btr)
        cr2 = cr1*btr
        return np.array([c[0],cr1, cr2])
    
    def bankPC(self,c):
        return (self.PV(self.y) - self.PV(c))


    def reneg_proof(self):
        """Optimum renegotiation-proof contract for any kappa
        Solved for analytically if kappa=0, numerically otherwise"""
        if self.kappa >= self.kbar():
            return self.fcommit()
        elif self.kappa ==0.0:
            return self.ownsmooth()
        else:
            initial_guess = self.ownsmooth()
            res = minimize(self.negPVU, initial_guess,  
                           method='COBYLA', 
                           constraints=(
                            {'type': 'ineq', 'fun':self.bankPC},
                            {'type': 'ineq', 'fun':self.noreneg} ))
            return res.x if res.success else None
            

class Monopoly(Contract):                    # build on contract class
    """ Class for solving Monopoly equilibrium contracts  """
    def __init__(self,beta):
        super().__init__(beta)    # inherit parent class properties
        self.kappa  = 0                        # cost of renegotiation

    def __repr__(self):
        return 'Monopoly(beta=%s, y=%s)' % (self.beta, self.y)

    def fcommit(self):
        """monopolist optimal full commitment contractwith period0 self
        from closed form solution for CRRA"""
        bt, rho = self.beta, self.rho
        btr = bt**(1/rho)
        pvu = self.PVU(self.y,bt)
        if rho==1:
            lnc0 = (pvu-2*bt*np.log(bt))/(1+2*bt)
            c0 = np.exp(lnc0)
        else:
            a = ((pvu*(1-rho))**(1/(1-rho)) )
            b = (1 + 2*btr)**(1/(rho-1))
            c0 = a*b
        return np.array([1,btr,btr])*c0
    
    def ownsmooth(self):
        """contract if have to smooth by oneself without contract
        from closed form solution for CRRA with kappa=0"""
        beta, rho = self.beta, self.rho
        btr = beta**(1/rho)
        ubar = self.PVU(self.y, self.beta)
        B = (beta+btr)/(1+btr)
        N1 = ( 1+beta**((1-rho)/rho) )**(1/rho)
        D1 = ( 1+beta**(1/rho) )**((1-rho)/rho)
        c0p = ( ubar*(1-rho)/(1+btr*(N1/D1) ) )** (1/(1-rho))
        c1p = B**(1/rho) * c0p
        return np.array( [c0p, c1p, btr*c1p])

    def kbar(self):
        '''Monopoly Renegotiation cost necessary to sustain efficient
        competitive contract'''
        rho = self.rho
        if (rho == 1):
            rho = 0.999  # cheap trick to deal with special case
        btr = self.beta ** (1 / rho)
        c1F = self.fcommit()[1]
        A = (2 - (1 + btr) * ((1 + self.beta)
                              / (1 + btr)) ** (1 / (1 - rho)))
        return A * c1F

    def negprofit(self,c):
        """ Negative profits (for minimization)"""
        return  -(self.PV(self.y) - self.PV(c))

    def participation_cons(self,c):
        return (self.PVU(c,self.beta)  - self.PVU(self.y,self.beta))


    def reneg_proof(self):
        """contract for any kappa
           formerly 'opt' function """
        if self.kappa >= self.kbar():
            return self.fcommit()
        elif self.kappa ==0:
            return self.ownsmooth()
        else:
            initial_guess = self.ownsmooth()
            res = minimize(self.PV, initial_guess,  
                           method='COBYLA', 
                           constraints=(
                            {'type': 'ineq', 'fun':self.participation_cons},
                            {'type': 'ineq', 'fun':self.noreneg} ))
            return res.x if res.success else None
            

    def reneg_proof2(self):
        """Find period 0 monopoly best renegotiation-proof contract
        by searching over subgame perfect responses """
        guess = self.fcommit()
        Y = np.sum(self.y)

        def f(c0):
            c1 = self.reneg(c0)[1]
            return -(self.u(c0) 
                     + self.beta * (self.u(c1) + self.u(Y - c0 - c1)))

        c0rp = minimize(f, guess, method='Nelder-Mead').x[0]
        print(c0rp)
        c1rp = self.reneg(c0rp)[1]
        c2rp = Y - c0rp - c1rp

        return np.array([c0rp, c1rp, c2rp])
    

"""
TWO-PERIOD VERSION

To do: fcommit forward

 
"""


class Competitive2(Competitive):                    # build on contract class
    """ Class for solving competitive equilibrium contracts  
        Most functionality inherited from Competitive class, here only making adjustments 
        needed for 2 period version
        """
    def __init__(self, beta):
        super(Competitive2, self).__init__(beta)  # inherits properties from Competitive(Contract)
        self.kappa  = 0                           # cost of renegotiation

    def __repr__(self):
        return 'Competitive2(beta=%s, y=%s)' % (self.beta, self.y)

    def fcommit(self):
        """competitive optimal full commitment contractwith period0 self
        from closed form solution for CRRA"""
        btr = self.beta**(1/self.rho)
        return np.array([1,btr])*self.PV(self.y)/(1 + btr)
    
    def kbar(self):
        '''Renegotiation cost necessary to sustain full commitment
        competitive contract'''
        rho = self.rho
        if (rho == 1):
            rho = 0.999  # cheap trick to deal with special case
        btr = self.beta ** (1 / rho)
        c1F = np.sum(self.y) * btr / (1 + 2 * btr)
        A = (2 - (1 + btr) * ((1 + self.beta) / (1 + btr))
             ** (1 / (1 - rho)))
        return A * c1F
    
    def ownsmooth(self):
        """contract if have to smooth by oneself without contract
        from closed form solution for CRRA with kappa=0"""
        beta, rho = self.beta, self.rho
        btr = beta**(1/rho)
        Y = self.PV(self.y)
        lam = (beta+btr)/(1+btr)
        ltr = lam**(1/rho)
        c0 = Y/( 1+(1+btr)*ltr )
        c1 = (Y-c0)/(1+btr)
        return np.array( [c0, c1, btr*c1])

    def creneg(self, c):
        """ Renegotiated competitive contract offered to period-1-self if
        c_0 is past but (c_1,c_2) now replaced by (cr_1, cr_2)"""
        beta, rho = self.beta, self.rho
        btr = beta**(1/rho)
        cr1 = (c[1]+c[2])/(1+btr)
        cr2 = cr1*btr
        return np.array([c[0],cr1, cr2])
    
    def bankPC(self,c):
        return (self.PV(self.y) - self.PV(c))


    def reneg_proof(self):
        """Optimum renegotiation-proof contract for any kappa
        Solved for analytically if kappa=0, numerically otherwise"""
        if self.kappa >= self.kbar():
            return self.fcommit()
        elif self.kappa ==0.0:
            return self.ownsmooth()
        else:
            xg = self.ownsmooth()
            res = minimize(self.negPVU, xg,  
                           method='COBYLA', 
                           constraints=(
                            {'type': 'ineq', 'fun':self.bankPC},
                            {'type': 'ineq', 'fun':self.noreneg} ))
            return res.x
            



if __name__ == "__main__":

    print('Test Block')
    RHO   =  0.95  # for testing different values
    print("Competitive contracts:")
    cC = Competitive(beta = 0.5)
    cC.rho = RHO
    cC.y = [80, 110, 110]
    cC.print_params()

    cCF = cC.fcommit()
    print("cCF: ",cCF)

    cCr = cC.reneg(cCF)
    print(f'kbar = {cC.kbar():.3f}')
    print("reneg from cCF: ",cCr)  

    cC.kappa = 0
    print('kappa=0', cC.reneg_proof())
    cC.kappa = 1
    print('kappa=1', cC.reneg_proof())    
    
    print()
    print("Monopoly contracts:")

    cM = Monopoly(beta = 0.5)
    cM.y = [80, 110, 110]
    cM.rho = RHO
    cM.print_params()
    
    cMF = cM.fcommit()
    print("cMF: ",cMF)
    
    cMr = cM.reneg(cMF)
    print(f'kbar = {cM.kbar():.3f}')
    print("reneg from cMF: ",cMr)
    cM.kappa = 0
    print('kappa=0', cM.reneg_proof())

    cM.kappa = 1
    print('kappa=5', cM.reneg_proof())    
    print()
    