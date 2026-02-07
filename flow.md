How MIMAS runs, step by step (what happens first → what happens next)
# 1) Main program starts + defines the simulation window

The program (PROGRAM walcek) first defines:

year to simulate (i_jahr)

a fixed “dynamics year” for winds (i_jahr_fixdyn = 1976)

season limits (start/end month/day)

grid sizes (lat, lon, vertical) and the base time step (dttrans = 90 s) 

block_diagram_draft-first - Copy

This is basically: “what time period am I running, and what resolution am I using?”

# 2) Calendar setup

It calls set_monate(...) to determine:

number of days per month

day-of-year indexing (max_doy) 

block_diagram_draft-first - Copy

This is needed because later the model loops through month/day/hour.

# 3) Global initialization (reads constant input files)

sub_init_global(i_jahr, xlyman_obs, scale_ch4) loads/derives “global constants and climatologies” such as:

eddy diffusion profile (CTURBKZZ)

backscatter-related constants (CBACKSCATTER532)

Earth geometry ellipse (CELLIPSE_EARTH)

Lyman-alpha or UV proxy (“Char_lyman” / xlyman_obs)

CH₄ dataset (“CH4_1861_2008”) → used to scale H₂O production via CH₄ oxidation (scale_ch4) 

block_diagram_draft-first - Copy

So the model is setting the background physics knobs before any time stepping.

# 4) Read initial 1D H₂O profile

sub_lesseh2o_zprinit reads H2OINIT, which is a 1D vertical water vapor profile on a different vertical coordinate (166 levels). 

block_diagram_draft-first - Copy

This is the initial “seed” H₂O which will be mapped onto the full 3D grid later.

# 5) Read dynamics from LIMA for the start time

sub_lesedyn(...) reads the dynamical background fields on the model grid:

winds (u, v, w)

temperature (T)

density

vertical coordinate (z) 

block_diagram_draft-first - Copy

Important design note in your diagram:

winds come from the fixed dynamics year (1976),

while T, density, z can correspond to the current running year. 

block_diagram_draft-first - Copy

Then the model optionally smooths vertical wind near polar latitudes. 

block_diagram_draft-first - Copy

# 6) Numerical safety check (CFL / Courant condition)

sub_maxmin(um1, vm1, wm1) checks max wind speeds and ensures the Courant number ≤ 1 (stability condition). If not satisfied, the run stops. 

block_diagram_draft-first - Copy

This is: “are my dt and grid spacing stable for advection?”

# 7) Initialize the 40 million dust particles (condensation nuclei)

sub_init_staub(nstart, ...) creates the full particle ensemble:

arrays xfeld, yfeld, zfeld for particle positions

rinit initial radii (1.2–3.7 nm bins)

rfeld current radius (starts as rinit)

location: poleward of 55°N and within ±1 km of mesopause region 

block_diagram_draft-first - Copy

This is the “big bang” for the Lagrangian part of the model.

# 8) Restart / output bootstrapping

If nstart = 0 (fresh run):

it creates a restart-like file CRINIT and writes initial particle states

If it’s a restart:

it reads particle state from files (e.g., .40md19 style) and reads hm (H₂O field) from .hd19 files 

block_diagram_draft-first - Copy

So before time stepping, the model decides:

“start from scratch?” or “continue from saved state?”

# 9) Precompute microphysics lookup tables

sub_tabelle is called once early to precompute tables like:

saturation vapor pressure over plane ice (Pα)

Kelvin factors across temperature and radius ranges

helper terms for growth rate / sedimentation computations 

block_diagram_draft-first - Copy

This is a performance trick: microphysics is expensive, so it tabulates needed functions.

Now the real model propagation begins (the nested loops)
# 10) Loop over season → days → hours

The run enters nested loops:

months in [start..end]

days in each month

hours of day (i_std from 1..24) 

block_diagram_draft-first - Copy

Everything from here repeats every timestep/hour.

Inside each model hour: H₂O setup → transport → microphysics+particles → outputs
# 11) Map initial H₂O onto the 3D grid (only when needed)

sub_h2oinit_zpr_to_zgeo(zm1) interpolates the initial 1D H₂O profile onto the model’s 3D pressure/height grid (zm1) producing hminit3d. 

block_diagram_draft-first - Copy

Then your diagram shows two scalings:

reduce by 10% (*0.9)

scale by CH₄-related factor (*scale_ch4) 

block_diagram_draft-first - Copy

If hm is still zero (first time), it is initialized from this.

# 12) Compute solar geometry and photolysis of H₂O

sub_zenit(...) computes zenith angle / solar illumination geometry (secchi) 

block_diagram_draft-first - Copy

sub_photolyse(...) applies photodissociation to the H₂O field hm 

block_diagram_draft-first - Copy

So: even before dynamics transport, H₂O can be chemically reduced by light.

# 13) Update the background dynamics (winds, T, density)

Each hour, the model reads the next set of LIMA fields using sub_lesedyn(...) (um2, vm2, wm2, tm2, dm2, zm2). 

block_diagram_draft-first - Copy

Then it blends time levels:

um = a1*um1 + a2*um2 (and same for v, w, T, ρ) 

block_diagram_draft-first - Copy

This is basically time interpolation / smoothing between two dynamic snapshots.

Then again:

sub_maxmin(...) CFL check, stop if unstable. 

block_diagram_draft-first - Copy

# 14) Eulerian H₂O transport (advection)

There is an inner loop nlauf ≤ 4, and in it the model calls:

sub_transwalcek(dttrans, dx, dy, dz, um, vm, wm, hm) twice

so total advection corresponds to 180 s if dttrans=90 s 

block_diagram_draft-first - Copy

This is the Walcek-type monotonic advection step for the grid-based H₂O field.

# 15) Vertical diffusion of H₂O (turbulence)

After advection:

dttrans_2 = 2*dttrans

sub_diffu_h2o(hm, dz, dttrans_2, turbkzz) applies turbulent vertical diffusion using turbkzz. 

block_diagram_draft-first - Copy

So H₂O evolution each hour is:

advect (x2 calls)

diffuse vertically (Kzz)

photolyse

# 16) Photolysis again

After transport + diffusion, the code again computes zenith angle and applies photolysis using xlyman_obs (observed Lyman-alpha series). 

block_diagram_draft-first - Copy

So photolysis is applied to the H₂O field repeatedly as it evolves.

Now the particle side: “pick a particle and let it live”
# 17) Particle microphysics + particle transport step

Then the model calls:

sub_tracertransp

This is the heart of the Lagrangian part. Your diagram shows that inside it, for each particle n = 1..40 million: 

block_diagram_draft-first - Copy

17a) Map particle position to grid cell

i = int(xfeld(n))

j = int(yfeld(n))

k = int(zfeld(n)) 

block_diagram_draft-first - Copy

# 17b) Read background conditions at that cell

tback = tm(j,i,k) (ambient temperature)

compute work1akt which is essentially pH2O derived from hm * dm (H₂O mixing ratio × density → partial pressure proxy) 

block_diagram_draft-first - Copy

17c) Decide whether the particle is eligible for ice microphysics

It checks something like:

if rfeld(n) > rinit(n) → particle already has ice (or has grown)

else it’s “bare dust” (nucleus) 

block_diagram_draft-first - Copy

17d) Compute saturation and Kelvin effect

For the particle:

compute particle temperature tp (often close to tback)

compute saturation vapor pressure over curved surface Psat(r)

compute saturation ratio S(r) = Pα / Psat(r) (conceptually) 

block_diagram_draft-first - Copy

# 17e) Growth condition

If conditions allow (example in diagram):

if tback < 155 and S(r) > 1 → compute dr/dt and update radius:

rfeld(n) = rfeld(n) + dr/dt * dt 

block_diagram_draft-first - Copy

else:

dr/dt = 0, radius unchanged 

block_diagram_draft-first - Copy

This is the precise “check if ice forming size/condition” step you described.

# 17f) Move the particle (advection + sedimentation + turbulent random part)

The diagram shows:

compute turbulent velocity and sedimentation velocity

update xfeld, yfeld, zfeld using something like:

x(t+dt) = x(t) + u*dt (and same for y, z with extra terms) 

block_diagram_draft-first - Copy

So each particle experiences:

horizontal advection by winds

vertical motion = background w + sedimentation + turbulence

# 18) Dust “reallocation” to keep the CN population stable

After the particle update, the model does something very important:

Bare dust particles that move outside a defined nucleation domain are randomly relocated back into a dust injection region. 

block_diagram_draft-first - Copy

Your diagram explicitly says:

“Reallocation of bare dust particles inside and outside the defined boundary into dust injection zone.” 

block_diagram_draft-first - Copy

This is exactly the mechanism that prevents the “dust reservoir” from leaking away in long runs:

particles don’t permanently exit; they get recycled back into the source region (if bare).

(Ice-coated particles are typically not relocated, so ice can exist outside the nucleation zone.)

Outputs (what gets written and when)
# 19) Write diagnostics at selected hours

At specific hours (6/12/18/24):

schreibe_beta0(...) computes and writes reduced microphysical outputs (e.g., backscatter, ice radius, number densities) into .betd19 files 

block_diagram_draft-first - Copy

also writes .hd19 files for hm (H₂O field) 

block_diagram_draft-first - Copy

# 20) Write particle positions occasionally

At specific days (10/20/30 as shown):

schreibe_40mio(...) writes particle positions xfeld, yfeld, zfeld to .40md19 style output 

block_diagram_draft-first - Copy

# 21) Advance time

At the end of the hour:

shift “current” fields to previous (um1=um2, etc.)

increment hour, day, month loops