#####
##### Earth's tidal astronomy
#####

const astronomical_epoch = DateTime(1900, 1, 1)

const astronomical_angles = (0, 277.0248, 280.1895, 334.3853, 90)
const astronomical_rates = (360 * 36525, 481267.8906, 36000.7689, 4069.0340, 0)

const lunar_node_angle = 259.1568
const lunar_node_rate = -1934.1420

"""
Properties of Earth's tidal constituents: the equilibrium `amplitude` [m], the multiples of the
astronomical angles that sum to the equilibrium argument, and the `nodal` coefficients for the
18.6-year lunar cycle.

References: Schureman (1958), Kowalik and Luick (2019).
"""
const earth_tidal_constituents = (
    M2 = (amplitude = 0.242334 * 0.693, argument = ( 2, -2,  2,  0,  0), nodal = (f₀ = 1.000, f₁ = -0.037, u₁ =  -2.1)),
    S2 = (amplitude = 0.112743 * 0.693, argument = ( 2,  0,  0,  0,  0), nodal = (f₀ = 1.000, f₁ =  0.000, u₁ =   0.0)),
    N2 = (amplitude = 0.046397 * 0.693, argument = ( 2, -3,  2,  1,  0), nodal = (f₀ = 1.000, f₁ = -0.037, u₁ =  -2.1)),
    K2 = (amplitude = 0.030684 * 0.693, argument = ( 2,  0,  2,  0,  0), nodal = (f₀ = 1.024, f₁ =  0.286, u₁ = -17.7)),
    K1 = (amplitude = 0.141565 * 0.736, argument = ( 1,  0,  1,  0,  1), nodal = (f₀ = 1.006, f₁ =  0.115, u₁ =  -8.9)),
    O1 = (amplitude = 0.100661 * 0.695, argument = ( 1, -2,  1,  0, -1), nodal = (f₀ = 1.009, f₁ =  0.187, u₁ =  10.8)),
    P1 = (amplitude = 0.046848 * 0.706, argument = ( 1,  0, -1,  0, -1), nodal = (f₀ = 1.000, f₁ =  0.000, u₁ =   0.0)),
    Q1 = (amplitude = 0.019273 * 0.695, argument = ( 1, -3,  1,  1, -1), nodal = (f₀ = 1.009, f₁ =  0.187, u₁ =  10.8)),
    Mf = (amplitude = 0.042041 * 0.693, argument = ( 0,  2,  0,  0,  0), nodal = (f₀ = 1.043, f₁ =  0.414, u₁ = -23.7)),
    Mm = (amplitude = 0.022191 * 0.693, argument = ( 0,  1,  0, -1,  0), nodal = (f₀ = 1.000, f₁ = -0.130, u₁ =   0.0)),
)

"""
$(TYPEDSIGNATURES)

`TidalHarmonics` for Earth's astronomical tide of `constituents` at `reference_date`: their
frequencies and the phases that the equilibrium arguments and nodal corrections reach at that date.
Model time is seconds from `reference_date`.
"""
function earth_tidal_harmonics(reference_date::DateTime;
                               constituents = keys(earth_tidal_constituents),
                               ramp_time = 0)

    names = Tuple(Symbol(constituent) for constituent in constituents)
    properties = Tuple(earth_tidal_constituents[name] for name in names)

    centuries = Dates.value(reference_date - astronomical_epoch) / (86_400_000 * 36_525)
    angles = astronomical_angles .+ astronomical_rates .* centuries
    node = lunar_node_angle + lunar_node_rate * centuries

    return TidalHarmonics(constituents = names,
                          frequencies = Tuple(deg2rad(sum(p.argument .* astronomical_rates)) / (36_525 * 86_400) for p in properties),
                          phases = Tuple(deg2rad(mod(sum(p.argument .* angles) + p.nodal.u₁ * sind(node), 360)) for p in properties),
                          nodal_factors = Tuple(p.nodal.f₀ + p.nodal.f₁ * cosd(node) for p in properties),
                          equilibrium_amplitudes = Tuple(p.amplitude for p in properties),
                          species = Tuple(first(p.argument) for p in properties),
                          ramp_time = ramp_time)
end
