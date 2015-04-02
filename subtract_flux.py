import numpy as np
import voigt
import data_types




LIGHT_SPEED = 299792.458
def calc2(spec, N, b, z, f, Gamma, restWave, subFlag):
    tau_threshold = 0.001

    cwave = (1.0+z)*restWave
    vdopp = b/cwave*1.0e13
    alpha = Gamma/(4.0*np.pi*vdopp)/(1.0+z)    
    fact = np.pow(10.0, N) * 2.647E-2 * 
      f/(np.sqrt(np.pi)*vdopp) * 1.0/(1.0+z)
    
    cpix = spec.waveIndex(cwave)
    if (cpix >= spec.size()-3) return cpix, cpix
    if (cpix < 1)              return cpix, cpix

    wave_high = (spec.wave[cpix] + spec.wave[cpix+1])/2
    wave_low  = (spec.wave[cpix] + spec.wave[cpix-1])/2

    nsamp = 10
    delta_wave = (wave_high - wave_low)/nsamp

    tau_sum = 0.0

    for i in range(nsamp):
        wave_i = wave_low + i*delta_wave
        vbar = LIGHT_SPEED/b * (wave_i/cwave - 1.0)
        tau_sum += fact*voigt(vbar, alpha)

    tau = tau_sum / nsamp

    if subFlag: 
        spec.flux[cpix] -= tau
    else:         
        spec.flux[cpix] += tau

    lpix = cpix
    while (tau > tau_threshold and lpix > 1):
        lpix -= 1
        wave_high = (spec.wave[lpix] + spec.wave[lpix+1])/2
        wave_low  = (spec.wave[lpix] + spec.wave[lpix-1])/2

        vbar = LIGHT_SPEED/b * (wave_low/cwave - 1.0)
        tau_low = fact*voigt(vbar, alpha)

        vbar = LIGHT_SPEED/b * (wave_high/cwave - 1.0)
        tau_high = fact*voigt(vbar, alpha)

        tau = (tau_low + tau_high)/2

        if subFlag: spec.flux[lpix] -= tau
        else: spec.flux[lpix] += tau

    #Now walk the other way ...
    hpix = cpix

    vbar = LIGHT_SPEED/b * (spec.wave[hpix]/cwave - 1.0)
    tau = fact*voigt(vbar, alpha)


    while (tau > tau_threshold and hpix < spec.size()-3):
        hpix += 1

        wave_high = (spec.wave[hpix] + spec.wave[hpix+1])/2
        wave_low  = (spec.wave[hpix] + spec.wave[hpix-1])/2

        vbar = LIGHT_SPEED/b * (wave_low/cwave - 1.0)
        tau_low = fact*voigt(vbar, alpha)

        vbar = LIGHT_SPEED/b * (wave_high/cwave - 1.0)
        tau_high = fact*voigt(vbar, alpha)

        tau = (tau_low + tau_high)/2

        if subFlag: spec.flux[hpix] -= tau
        else :      spec.flux[hpix] += tau

    return lpix, hpix
