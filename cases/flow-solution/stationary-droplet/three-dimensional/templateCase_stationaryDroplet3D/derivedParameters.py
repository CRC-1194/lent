def computeNuDroplet(values):
    if "nu_droplet" in values:
        return values["nu_droplet"]

    rhoDroplet = values["rho_droplet"]
    return 0.008165/rhoDroplet

def computeNuAmbient(values):
    if "nu_ambient" in values:
        return values["nu_ambient"]

    rhoAmbient = values["rho_ambient"]
    return 0.008165/rhoAmbient

def computeDeltaT(values):
    import math

    # Give explicitly prescribed time step sizes precedence over
    # computed values
    if "delta_t" in values:
        return values["delta_t"]

    domainLength = values["domain_length"]
    resolution = values["resolution"]
    rho_droplet = values["rho_droplet"]

    # rho_ambient might be expressed in terms of rho_droplet
    rho_ambient = 0.0
    if values["rho_ambient"] == "$rho_droplet":
        rho_ambient = rho_droplet
    else:
        rho_ambient = values["rho_ambient"]
    sigma = values["surface_tension_coefficient"]

    # Use the resolution of the Eulerian mesh rather than the front element size.
    # So far there is no indication that the triangle size plays a role for the
    # maximum time step.
    # The computation follows eq. (18) in Popinet 2009 and
    # eq. in (43) Denner & van Wachem 2015, respectively.
    ltria = domainLength/(resolution)

    # Safety coefficient to stay below critical delta_t threshold
    sf = values["scale_delta_t"]

    return sf*math.sqrt((rho_droplet + rho_ambient)*math.pow(ltria,3.0)/(2.0*math.pi*sigma))

def computeEndTime(values):

    if "end_time" in values:
        return values["end_time"]

    diameter = 2.0 * values["radius"]
    nuDroplet = values["nu_droplet"]

    # Avoid divide by zero
    if nuDroplet == 0.0:
        nuDroplet = 1.0e-10

    return diameter*diameter/nuDroplet



delta_t = computeDeltaT(locals())
end_time = computeEndTime(locals())
nu_droplet = computeNuDroplet(locals())
nu_ambient = computeNuAmbient(locals())
