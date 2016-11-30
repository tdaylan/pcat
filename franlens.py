from __init__ import *
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
from scipy import signal
from astropy.convolution import convolve, AiryDisk2DKernel

class Params(object):
    """
    The main parameter class, used to define 
    the lensmodel and source parameters. It
    contains the parameter name, its value,
    a default range in which the parameter
    should lie (used as prior for running MCMC),
    and a boolean flag determining whether the 
    parameter should be varied in an MCMC.
    """
    
    def __init__(self,name,value,prior,tv):
        """
        Initialization of a parameter.

        Args:
            name: A string containing the name of the parameter.
            value: The parameter value.
            prior: a 1D array containing the lower and upper bounds of the parameters.
            tv: a boolean operator determining whether we want to vary the parameter in MCMC.

        Returns:
            An instance of the parameter class.

        Raises:
            TypeError: if name is not a string.
            TypeError: if value is not a float.
            TypeError: if prior is not a list
            TypeError: if tv is not a boolean operator
        """
        #Check for type errors
        if not isinstance(name,str):
            raise TypeError('The parameter name must be a string.')
        if not isinstance(value,float):
            raise TypeError('The parameter value must be a float.')
        if not isinstance(prior,list):
            raise TypeError('The prior range is not a list of an upper and lower bounds.')
        if not isinstance(tv,bool):
            raise TypeError('The tv attribute must be a boolean operator.')
        
        #Assign values
        self.name  = name
        self.value = value
        self.prior = np.array(prior)
        self.tv = tv
        

class LensModel(object):
    """The main lens model class"""
    
    def __init__(self, massmodel, xgal, ygal, ellipticity, ellipt_angle, shear, shear_angle, *args):
        self.massmodel = massmodel
        #Primary lens parameters
        self.lensparams = []
        self.lensparams.append(Params('xgal',xgal,[-2.0,2.0],False))
        self.lensparams.append(Params('ygal',ygal,[-2.0,2.0],False))
        self.lensparams.append(Params('ellipticity',ellipticity,[0,0.5],False))
        self.lensparams.append(Params('ellipticity_angle',ellipt_angle,[0.0,180.0],False)) #Note that this angle is measured counterclockwise from the y-axis
        self.lensparams.append(Params('shear',shear,[0.0,0.3],False))
        self.lensparams.append(Params('shear_angle',shear_angle,[0.0,180],False)) #This angle is measured c.c. from x-axis
        self.s = 1e-20
        
        if self.massmodel == 'SIE':
            self.lensparams.append(Params('b_ein',args[0],[0.2,2.0],False))
        
        elif self.massmodel == 'alpha':
            self.lensparams.append(Params('b_ein',args[0],[0.2,2.0],False))
            self.lensparams.append(Params('alpha',args[1],[0.5,1.5],False))
            print('Note: alpha lens model has not yet been coded up')
            
        else:
            raise RuntimeError('Type of lens unrecognized!')
            
        self.numlensparams = len(self.lensparams)
            
    def potential(self,xx,yy):
        """The projected gravitational potential method"""
        
        #Translate to the center of the potential
        xxt = xx - self.lensparams[0].value
        yyt = yy - self.lensparams[1].value
        
        #Define rotation angle
        rot_angle = self.lensparams[3].value*np.pi/180.0
        
        #Rotate coordinates according to given angle
        x = -np.sin(rot_angle)*xxt + np.cos(rot_angle)*yyt
        y = -np.cos(rot_angle)*xxt - np.sin(rot_angle)*yyt
        
        #Define q
        q = 1.0 - self.lensparams[2].value
        
        #Compute the potential
        if self.massmodel == 'SIE':
            if self.lensparams[2].value > 1.0e-4:    
                psi = np.sqrt(q**2*(self.s**2+x**2)+y**2)
                phix = (self.lensparams[6].value*q/np.sqrt(1.0-q**2)*
                        np.arctan(np.sqrt(1.0-q**2)*x/(psi+self.s)))
                phiy = (self.lensparams[6].value*q/np.sqrt(1.0-q**2)*
                        np.arctanh(np.sqrt(1.0-q**2)*y/(psi+q**2*self.s)))
                phi = (x*phix + y*phiy-self.lensparams[6].value*q*self.s*np.log(np.sqrt((psi+self.s)**2+(1.0-q**2)*x**2)) 
                      +self.lensparams[6].value*q*self.s*np.log((1+q)*self.s))
            else:
                phi = self.lensparams[6].value*(np.sqrt(x**2 + y**2 + self.s**2)- self.s - 
                        self.s*np.log(0.50 + np.sqrt(x**2 + y**2 + self.s**2)/(2.0*self.s)))
           
                   
        elif self.massmodel == 'alpha':
            phi = 0.0
        
        #Add the external shear contribution
        phi += (0.5*self.lensparams[4].value*np.cos(2.0*self.lensparams[5].value*np.pi/180.0)*(xx**2 - yy**2) 
                    + self.lensparams[4].value*np.sin(2.0*self.lensparams[5].value*np.pi/180.0)*xx*yy)
        
        return phi
    
    def deflection(self,xx,yy):
        """The deflection vector at position x,y"""
        
        #Translate to the center of the potential
        xxt = xx - self.lensparams[0].value
        yyt = yy - self.lensparams[1].value

        #Define rotation angle
        rot_angle = self.lensparams[3].value*np.pi/180.0
        
        #Rotate coordinates according to given angle
        x = -np.sin(rot_angle)*xxt + np.cos(rot_angle)*yyt
        y = -np.cos(rot_angle)*xxt - np.sin(rot_angle)*yyt
        
        #Define q
        q = 1.0 - self.lensparams[2].value
       
        if self.massmodel == 'SIE':
            if self.lensparams[2].value > 1.0e-4:
                psi = np.sqrt(q**2*(self.s**2+x**2)+y**2)
                alphaxt = (self.lensparams[6].value*q/np.sqrt(1.0-q**2)*
                            np.arctan(np.sqrt(1.0-q**2)*x/(psi+self.s)))
                alphayt = (self.lensparams[6].value*q/np.sqrt(1.0-q**2)*
                            np.arctanh(np.sqrt(1.0-q**2)*y/(psi+q**2*self.s)))
            else:
                    
                rint = np.sqrt(x**2 + y**2 + self.s**2)
                alphaxt = self.lensparams[6].value*(x/(rint+self.s)) 
                alphayt = self.lensparams[6].value*(y/(rint+self.s))
        
        elif self.massmodel == 'alpha':
            alphaxt = 0.0
            alphayt = 0.0
            
        #Rotate back vector to original basis
        alphax = -np.cos(rot_angle)*alphayt-np.sin(rot_angle)*alphaxt
        alphay = np.cos(rot_angle)*alphaxt-np.sin(rot_angle)*alphayt
        
        #Add the external shear contribution
        alphax += (self.lensparams[4].value*np.cos(2.0*self.lensparams[5].value*np.pi/180.0)*xx 
                + self.lensparams[4].value*np.sin(2.0*self.lensparams[5].value*np.pi/180.0)*yy)
        alphay += (-self.lensparams[4].value*np.cos(2.0*self.lensparams[5].value*np.pi/180.0)*yy 
                + self.lensparams[4].value*np.sin(2.0*self.lensparams[5].value*np.pi/180.0)*xx)
        
        return dstack((alphax, alphay))
   

def gauss_mat(size,axis_ratio,angle):
    """
    Compute the covariance matrix of a Gaussian, given a typical size,
    an axis ratio, and an angle
    """
    rot_mat = np.array([[np.cos(angle*np.pi/180.0),-np.sin(angle*np.pi/180.0)],
                        [np.sin(angle*np.pi/180.0),np.cos(angle*np.pi/180.0)]])
    pre_mat = np.array([[1.0/(axis_ratio*size)**2,0.0],[0.0,1.0/size**2]])
    
    return np.dot(np.transpose(rot_mat),np.dot(pre_mat,rot_mat))
    
           
class Source(object):
    """The main source class"""
        
    def __init__(self,model,usrc,vsrc,peak_bright,*args):
        self.model = model
        self.srcparams = []
        self.srcparams.append(Params('usrc',usrc,[-1.0,1.0],False))
        self.srcparams.append(Params('vsrc',vsrc,[-1.0,1.0],False))
        self.srcparams.append(Params('peak_brightness',peak_bright,[0.01,10],False))
        
        if self.model == 'Gaussian':
            self.srcparams.append(Params('src_size',args[0],[0.001,1.0],False))
            #How much is the Gaussian stretched along the x-axis compared to the y-axis
            self.srcparams.append(Params('src_axis_ratio',args[1],[1.0,7.0],False)) 
            self.srcparams.append(Params('src_angle',args[2],[0.0,180.0],False))
            self.cov_src = gauss_mat(args[0],args[1],args[2])
        else:
            raise RuntimeError('Type of source unrecognized!')
        
        self.numsrcparams = len(self.srcparams)
                 
    def brightness(self,u,v):
        
        if self.model == 'Gaussian':
            
            src_pos = np.array([u-self.srcparams[0].value,v-self.srcparams[1].value])
            S = self.srcparams[2].value * np.exp(-0.5 * np.sum(src_pos * np.tensordot(self.cov_src, src_pos,(1,0)),0))
            
        return S
        
    def gradient(self,u,v):
        """The gradient of the source in the source plane"""
        
        if self.model == 'Gaussian':
            src_pos = np.array([u-self.srcparams[0].value,v-self.srcparams[1].value])
            shift = np.tensordot(self.cov_src,src_pos,(1,0))
            exp_arg = np.sum(src_pos*shift,0)
            dSdu = -self.srcparams[2].value*shift[0]*np.exp(-0.5*exp_arg)
            dSdv = -self.srcparams[2].value*shift[1]*np.exp(-0.5*exp_arg)
        else:
            dSdu = 0.0
            dSdv = 0.0
            
        return np.array([dSdu,dSdv])
    

def retr_imaglens(gdat, gdatmodi=None, raww=False):
    
    if gdatmodi != None:
        gdattemp = gdat
        strg = ''
        sampvarb = getattr(gdatmodi, 'thissampvarb')
        psfnkern = gdatmodi.thispsfnkern
    else:
        gdattemp = gdatmodi
        sampvarb = getattr(gdat, 'mockfixp')
        strg = 'mock'
        psfnkern = gdat.truepsfnkern

    sourobjt = franlens.Source(gdat.truesourtype, sampvarb[getattr(gdat, strg + 'indxfixplgalsour')], \
                                                  sampvarb[getattr(gdat, strg + 'indxfixpbgalsour')], \
                                                  sampvarb[getattr(gdat, strg + 'indxfixpfluxsour')], \
                                                  sampvarb[getattr(gdat, strg + 'indxfixpsizesour')], \
                                                  sampvarb[getattr(gdat, strg + 'indxfixpratisour')], \
                                                  sampvarb[getattr(gdat, strg + 'indxfixpanglsour')])

    defl = zeros((gdat.numbsidecart, gdat.numbsidecart, 2))

    if raww:
        beinhost = 0.
        sherhost = 0.
    else:
        sherhost = sampvarb[getattr(gdat, strg + 'indxfixpsherhost')]
        beinhost = sampvarb[getattr(gdat, strg + 'indxfixpbeinhost')]
    listlensobjt = []
    listlensobjt.append(franlens.LensModel(gdat.truelenstype, sampvarb[getattr(gdat, strg + 'indxfixplgalhost')], \
                                                              sampvarb[getattr(gdat, strg + 'indxfixpbgalhost')], \
                                                              sampvarb[getattr(gdat, strg + 'indxfixpellphost')], \
                                                              sampvarb[getattr(gdat, strg + 'indxfixpanglhost')], \
                                                              sherhost, \
                                                              sampvarb[getattr(gdat, strg + 'indxfixpsanghost')], \
                                                              beinhost))

    ## PS
    if getattr(gdat, strg + 'numbtrap') > 0 and not raww:
        for l in getattr(gdat, strg + 'indxpopl'):
            numbpnts = sampvarb[getattr(gdat, strg + 'indxfixpnumbpnts')[l]].astype(int)
            for k in range(numbpnts):
                # create lens model object for the PS 
                if gdatmodi == None:
                    listlensobjt.append(franlens.LensModel(gdat.truelenstype, gdat.mocklgal[l][k], gdat.mockbgal[l][k], 0., 0., 0., 0., gdat.mockspec[l][0, k]))
                else:
                    listlensobjt.append(franlens.LensModel(gdat.truelenstype, gdatmodi.thissampvarb[gdatmodi.thisindxsamplgal[l][k]], \
                                                                              gdatmodi.thissampvarb[gdatmodi.thisindxsampbgal[l][k]], \
                                                                              0., 0., 0., 0., \
                                                                              # beta
                                                                              gdatmodi.thissampvarb[gdatmodi.thisindxsampspec[l][0, k]]))
    numblensobjt = len(listlensobjt)
    for k in range(numblensobjt):
        defl += listlensobjt[k].deflection(gdat.lgalgridcart, gdat.bgalgridcart)

    # calculate the lensed image
    imaglens = sourobjt.brightness(gdat.lgalgridcart - defl[:, :, 0], gdat.bgalgridcart - defl[:, :, 1]) + \
                                                                                        sampvarb[getattr(gdat, strg + 'indxfixpbacp')[0, 0]] * gdat.backfluxlens
    
    #Convolve with the PSF
    # temp
    if False:
        imag = convolve(imaglens, psfnkern).flatten()
    else:
        imag = imaglens.flatten()

    return imag, imaglens, defl

    
