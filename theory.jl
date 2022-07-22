### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ c8112062-05f5-11ed-2107-ad2012a23b00
md"""
# Lab-frame particle motion using Hamiltonians
"""

# ╔═╡ 93fdaae2-a83a-4b72-8b4d-ea1eace6e647
md"""
## Lab frame hamiltonian
The lagrangian of a point particle in the absence of non-gravitational forces is 
```math
L = \frac{1}{2}mg_{\mu \nu}u^{\mu}u^{\nu} = \frac{1}{2}m u_{\mu} u^{\mu}
```
which means that the generalized momenta are 
```math
p_{\mu} = \frac{\partial L}{\partial u^{\mu}} = m u_{\mu}
```
Consider the action:
```math
2S = \int L d\tau = \int p_{\mu}u^{\mu} d\tau = \int p_{\mu}dx^{\mu} = \int p_{\mu}\frac{dx^{\mu}}{dt} dt
```
In the last step, we have effectively defined a new lagrangian into one appropriate for integrating over the time coordinate $t$ rather than the proper time $\tau$. Expanding this out, we get:
```math
2S = \int (p_0 + p_j\frac{dx^j}{dt}) dt
```
Taking this integrand to be our lagrangian, we then derive the hamiltonian:
```math
H = p_j \frac{dx^j}{dt} - L = -p_0
```
Now consider the line element 
```math
ds^2 = -\alpha^2 dt^2 + \gamma_{ij}dx^i dx^j
```
We will assume that $\alpha^2$ is positive so that the t-axis is timelike and can parameterize timelike trajectories. The 4-velocity is normalized to -1, which means that 
```math 
p_{\mu}p^{\mu} = -m^2 = -\frac{p_0^2}{\alpha^2} + \gamma^{ij}p_ip_j
```
Substituting $p_0^2=H^2$ and solving for $H$ gives:
```math 
H = \alpha \sqrt{m^2 + \gamma^{ij}p_ip_j }
```
"""

# ╔═╡ b7f5a0ca-cf2f-4089-a5c4-ab322d29c724
md"""
## Recovering Newtonian gravity

Assume that $p^2$ is a small quantity in relation to $m^2$:
```math 
H = \alpha m \sqrt{ 1 + \frac{p_ip^i}{m^2}} \approx \alpha (m + \frac{p_ip^j}{2m})
```
In the Newtonian limit, $\alpha \approx \sqrt{1 - 2GM/r} \approx 1 - GM/r$. Since we are assuming both $GM/r$ and $p_ip^i$ are small, we can neglect the cross term as a first order approximation:
```math 
H \approx (1 - \frac{GM}{r})(m + \frac{p_ip^i}{2m}) \approx m + \frac{p_ip^i}{2m} - \frac{GMm}{r}
```
which is the classical hamiltonian for a particle in Newtonian gravity.
"""

# ╔═╡ fdbb529c-f9aa-4dd2-8a09-dcbf05b8706f
md"""
## Orbits in Schwarzschild spacetime

For orbits that remain outside the event horizon, we can use Schwarzschild coordinates without violating the assumption that the $t$ coordinate is timelike. We will orient our coordinates so that $\theta=\pi/2$ is the orbital plane, and will therefore ignore the $\theta$ coordinate. The line element in the orbital plane is then
```math 
ds^2 = -\alpha^2 dt^2 + \frac{1}{\alpha^2}dr^2 + r^2d\phi^2
```
where $\alpha^2 = (1 - r_s/r)$, and $r_s$ is the Schwarzschild radius. Our hamiltonian is then 
```math
H = \sqrt{(1 - r_s/r)(m^2 + p_{\phi}^2/r^2) + p_r^2}
```
The equations of motion are 
```math 
\dot{\phi} = (1 - \frac{r_s}{r})\frac{p_{\phi}}{H_0r^2}
```
```math
\dot{p_{\phi}} = 0
```
```math
\dot{r} = \frac{p_r}{H_0}
```
```math
\dot{p_r} = -\frac{r_sm^2}{2H_0r^2} + (1 - \frac{3r_s}{2r})\frac{p_{\phi}^2}{H_0r^3}
```
where $H_0$ is the initial value of $H$.
"""

# ╔═╡ eda855dc-631c-4a36-be63-8b821b501af3
md"""
## Toy spacetime for current sheet problem

The example of the Schwarzschild spacetime helps clarify what consequences the different components of the metric have on dynamics. From that, I propose a simplified toy metric:
```math 
	ds^2 = -\alpha^2 dt^2 + \frac{1}{\alpha^2}(dx^2 + dy^2)
```
In this spacetime, $\alpha^2$ clearly plays the role of a kindo of generalized gravitational potential. The hamiltonain is:
```math 
H = \sqrt{\alpha^2 m^2 + p_x^2 + p_y^2}
```
and the equations of motion are 
```math
\dot{x}_i = \frac{p_i}{H_0}
```
```math
\dot{p_i} = \frac{m^2}{2H_0}\partial_i \alpha^2
```


The coframe in these coordinates is
```math
\sigma_x = dx, \sigma_y=dy/\alpha, \sigma_z=dz, \sigma_t=-\alpha dt
```
and we have the following Hodge dual relations for 2-forms:
```math
\star dt \wedge dx = \frac{1}{\alpha^2}dy \wedge dz
```
```math
\star dt \wedge dy = \frac{1}{\alpha^2}dz \wedge dx
```
```math
\star dt \wedge dz = \frac{1}{\alpha^2}dx \wedge dy
```
```math
\star dx \wedge dy = -\alpha^2 dt \wedge dz
```
```math
\star dy \wedge dz = -\alpha^2 dt \wedge dx
```
```math
\star dz \wedge dx = -\alpha^2 dt \wedge dy
```

"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.3"
manifest_format = "2.0"

[deps]
"""

# ╔═╡ Cell order:
# ╟─c8112062-05f5-11ed-2107-ad2012a23b00
# ╟─93fdaae2-a83a-4b72-8b4d-ea1eace6e647
# ╟─b7f5a0ca-cf2f-4089-a5c4-ab322d29c724
# ╟─fdbb529c-f9aa-4dd2-8a09-dcbf05b8706f
# ╟─eda855dc-631c-4a36-be63-8b821b501af3
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
