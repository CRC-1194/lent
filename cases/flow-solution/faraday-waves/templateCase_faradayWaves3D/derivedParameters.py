def computeDeltaT(values):
    import math

    # Give explicitly prescribed time step sizes precedence over
    # computed values
    if "delta_t" in values:
        return values["delta_t"]

    domainLength = values["domain_length"]
    resolution = values["resolution"]
    rho_bottom = values["rho_bottom"]

    # rho_top might be expressed in terms of rho_bottom
    rho_top = 0.0
    if values["rho_top"] == "$rho_bottom":
        rho_top = rho_bottom
    else:
        rho_top = values["rho_top"]
    sigma = values["surface_tension_coefficient"]

    # Use the resolution of the Eulerian mesh rather than the front element size.
    # So far there is no indication that the triangle size plays a role for the
    # maximum time step.
    # The computation follows eq. (18) in Popinet 2009 and
    # eq. in (43) Denner & van Wachem 2015, respectively.
    ltria = domainLength/(resolution)

    # Safety coefficient to stay below critical delta_t threshold
    sf = values["scale_delta_t"]

    return sf*math.sqrt((rho_bottom + rho_top)*math.pow(ltria,3.0)/(2.0*math.pi*sigma))



delta_t = computeDeltaT(locals())
