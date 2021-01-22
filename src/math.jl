#==============================================================================#
# Constants
#==============================================================================#
using PhysicalConstants.CODATA2018

const Navo = AvogadroConstant.val
const kB = BoltzmannConstant.val
const ec = ElementaryCharge.val
const eps0 = VacuumElectricPermittivity.val
const temp0 = 298.0

const v0 = 18/((1e6)*Navo)
const l0 = v0^(1/3)
const lB0 = (ec^2)/(4*π*eps0*kB*temp0)
const b0 = sqrt(v0/(lB0/79.0))

const lBbar = lB0/l0/79.0 # Reference value used as default

#==============================================================================#
# Mathematical functions
#==============================================================================#

rss(x) = dot(x,x)/2

gamq(q, a) = exp(-0.5*(q*a)^2)

gsphere(q) = (9/q^6)*(sin(q) - q*cos(q))^2

gcoil(q) = (2/q^2)*(q - 1 + exp(-q))

function grod(q)
    si = sinint(q)
    return (2 / q) * (si - (1 - cos(q))/q)
end

grodhq(q) = pi / q

function gworm(q, lp, np, b = 1.0)
    coil = exp(-q*lp/2.0) / (1.0 + q^2*np*b*lp/6.0)
    rod = (1 - exp(-q*lp/2.0)) / (1.0 + q*np*b/pi)
    return coil + rod    
end