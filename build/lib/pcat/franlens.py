from __init__ import *
import matplotlib.image as mpimg
from scipy import signal
from astropy.convolution import convolve, AiryDisk2DKernel

class Params(object):
    
    def __init__(self, valu):
        
        self.value = valu
        

class LensModel(object):
    
    def __init__(self, massmodel, lgal, bgal, ellp, angl, sher, sang, bein, rcor):
        self.massmodel = massmodel
        self.lensparams = []
        self.lensparams.append(Params(lgal))
        self.lensparams.append(Params(bgal))
        self.lensparams.append(Params(ellp))
        self.lensparams.append(Params(angl))
        self.lensparams.append(Params(sher))
        self.lensparams.append(Params(sang))
        self.lensparams.append(Params(bein))
        self.lensparams.append(Params(rcor))
            
    def deflection(self,xx,yy):
        
        lgal = self.lensparams[0].value
        bgal = self.lensparams[1].value
        ellp = self.lensparams[2].value
        angl = self.lensparams[3].value
        sher = self.lensparams[4].value
        sang = self.lensparams[5].value
        bein = self.lensparams[6].value
        rcor = self.lensparams[7].value

        # translate the grid
        lgaltran = xx - lgal
        bgaltran = yy - bgal
        
        # rotate the grid
        lgalrttr = cos(angl) * lgaltran - sin(angl) * bgaltran
        bgalrttr = sin(angl) * lgaltran + cos(angl) * bgaltran
       
        if False:
            print 'angl'
            print angl
            print 'lgaltran'
            summgene(lgaltran)
            print 'bgaltran'
            summgene(bgaltran)
            print 'lgalrttr'
            summgene(lgalrttr)
            print 'bgalrttr'
            summgene(bgalrttr)

        if ellp > 1e-4:
            factflat = (1. - ellp)**2
            factrcor = sqrt(factflat * (rcor**2 + lgalrttr**2) + bgalrttr**2)
            facteccc = sqrt(1. - factflat)
            factbein = (bein * (1. - ellp) / facteccc)
            defllgalrttr = factbein *  arctan(facteccc * lgalrttr / (factrcor + rcor))
            deflbgalrttr = factbein * arctanh(facteccc * bgalrttr / (factrcor + factflat * rcor))
        
            if False:
                print 'factflat'
                print factflat
                print 'factrcor'
                summgene(factrcor)
                print 'facteccc'
                print facteccc
                print 'factbein'
                print factbein
                print 'factrcor'
                summgene(factrcor)
        
        else:
            rint = sqrt(lgalrttr**2 + bgalrttr**2 + rcor**2)
            defllgalrttr = bein * lgalrttr / (rint + rcor) 
            deflbgalrttr = bein * bgalrttr / (rint + rcor)
        
        #Rotate back vector to original basis
        defllgal =  cos(angl) * defllgalrttr + sin(angl) * deflbgalrttr
        deflbgal = -sin(angl) * defllgalrttr + cos(angl) * deflbgalrttr
        
        # external shear
        factcosi = sher * cos(2. * sang)
        factsine = sher * cos(2. * sang)
        defllgal += factcosi * xx + factsine * yy
        deflbgal += factsine * xx - factcosi * yy
       
        if False:
            print 'rcor'
            print rcor
            print 'bein'
            print bein
            print 'ellp'
            print ellp
            print 'defllgal'
            summgene(defllgal)
            print 'deflbgal'
            summgene(deflbgal)
            print 'lgaltran'
            summgene(lgaltran)
            print
            print
            print
            print
            print
            print
            print
            print

        return dstack((defllgal, deflbgal))
   

    
def gausmatr(size, rati, angl):
    
    rttrmatr = array([[cos(angl), -sin(angl)], [sin(angl), cos(angl)]])
    icovmatr = array([[1. / (rati * size)**2, 0.], [0., 1. / size**2]])
    
    return dot(transpose(rttrmatr), dot(icovmatr, rttrmatr))
    
           
class Source(object):
        
    def __init__(self,model, lgalsour, bgalsour, fluxsour, sizesour, ratisour, anglsour):
        self.model = model
        self.srcparams = []
        self.srcparams.append(Params(lgalsour))
        self.srcparams.append(Params(bgalsour))
        self.srcparams.append(Params(fluxsour))
        if self.model == 'gaus':
            self.srcparams.append(Params(sizesour))
            self.srcparams.append(Params(ratisour)) 
            self.srcparams.append(Params(anglsour))
            self.icovmatr = gausmatr(sizesour, ratisour, anglsour)
        
    
    def brightness(self,u,v):
        
        lgal = self.srcparams[0].value
        bgal = self.srcparams[1].value
        flux = self.srcparams[2].value
        size = self.srcparams[3].value
        rati = self.srcparams[4].value
        if self.model == 'gaus':
            posi = array([u - lgal, v - bgal])
            brgt = flux * exp(-0.5 * sum(posi * tensordot(self.icovmatr, posi, (1,0)), 0)) / size**2 / rati
            
        return brgt
       

    def gradient(self,u,v):
        
        if self.model == 'gaus':
            src_pos = array([u-self.srcparams[0].value,v-self.srcparams[1].value])
            shift = tensordot(self.icovmatr,src_pos,(1,0))
            exp_arg = sum(src_pos*shift,0)
            dSdu = -self.srcparams[2].value*shift[0]*exp(-0.5*exp_arg)
            dSdv = -self.srcparams[2].value*shift[1]*exp(-0.5*exp_arg)
            
        return array([dSdu,dSdv])
    

