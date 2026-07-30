# Urban Street Network Level Set Propagation (Mathematical README)

## Objective

Given a street network extracted from OpenStreetMap over a \(\sim 1\,\text{km}\times1\,\text{km}\) window of El Poblado (Medellin), the notebook builds a level-set function \(\phi\) on a \(1000\times1000\) pixel grid and evolves it with constant normal speed.  
The evolving zero-level set models a front that starts at the street/open-space interface and invades open urban space.

## 1) Continuous-to-discrete domain setup

Let \((x,y)\) be projected metric coordinates from the OSM graph (UTM).  
The bounding box \([x_{\min},x_{\max}]\times[y_{\min},y_{\max}]\) is mapped to image indices \((p_x,p_y)\in\{0,\dots,N-1\}^2\), with \(N=1000\):

\[
p_x=\left\lfloor \frac{x-x_{\min}}{x_{\max}-x_{\min}}(N-1)\right\rfloor,\quad
p_y=\left\lfloor \frac{y_{\max}-y}{y_{\max}-y_{\min}}(N-1)\right\rfloor.
\]

Street polylines are rasterized into a binary field \(B\):

\[
B(i,j)=
\begin{cases}
0,& \text{street pixel}\\
1,& \text{open-space pixel.}
\end{cases}
\]

## 2) Signed distance function initialization

Define masks:

\[
S=\{(i,j):B(i,j)=0\},\qquad O=\{(i,j):B(i,j)=1\}.
\]

Using Euclidean distance transforms:

\[
d_S(i,j)=\operatorname{dist}\big((i,j),S\big),\qquad
d_O(i,j)=\operatorname{dist}\big((i,j),O\big).
\]

The notebook sets:

\[
\phi_0(i,j)=d_S(i,j)-d_O(i,j).
\]

So:

\[
\phi_0<0 \text{ on street pixels},\quad
\phi_0>0 \text{ on open space}.
\]

In open space (\(B=1\)), \(d_O=0\), hence \(\phi_0=d_S\): distance to nearest street pixel.

## 3) Evolution equation

The front is evolved by the Hamilton-Jacobi/Eikonal form:

\[
\frac{\partial\phi}{\partial t}+v\,|\nabla\phi|=0,\qquad v=1.
\]

With the \(\phi<0\) (street-side) convention, \(v>0\) expands the \(\phi\le 0\) region into open space.

## 4) Godunov upwind discretization

For each time step, one-sided differences are computed:

\[
D_x^+\phi_{i,j}=\phi_{i,j+1}-\phi_{i,j},\quad
D_x^-\phi_{i,j}=\phi_{i,j}-\phi_{i,j-1},
\]
\[
D_y^+\phi_{i,j}=\phi_{i+1,j}-\phi_{i,j},\quad
D_y^-\phi_{i,j}=\phi_{i,j}-\phi_{i-1,j}.
\]

For \(v>0\), the Godunov numerical Hamiltonian used is:

\[
|\nabla\phi|_{G}=
\sqrt{
\max\!\left((\max(D_x^-\phi,0))^2,(\min(D_x^+\phi,0))^2\right)
+
\max\!\left((\max(D_y^-\phi,0))^2,(\min(D_y^+\phi,0))^2\right)
}.
\]

Explicit Euler update:

\[
\phi^{n+1}=\phi^n-\Delta t\,v\,|\nabla\phi^n|_{G}.
\]

Notebook parameters:

\[
\Delta t=0.35,\quad v=1,\quad N_{\text{steps}}=120,\quad t_n=n\Delta t.
\]

## 5) Derived outputs

### Propagation snapshots

At selected \(n\), the contour \(\phi^n=0\) is plotted over the rasterized street geometry.  
Reached set at time \(t_n\):

\[
R_n=\{(i,j):\phi^n(i,j)\le 0\}.
\]

### Coverage curve

Coverage percentage over time:

\[
C(t_n)=100\cdot\frac{|R_n|}{N^2}.
\]

This quantifies how quickly the front fills the urban domain.

### Arrival-time map

Arrival time is interpreted as

\[
\tau(i,j)\approx \phi_0(i,j)\quad\text{for } \phi_0(i,j)>0,\ v=1,
\]

and undefined on initial street pixels (\(\phi_0\le 0\), shown as NaN).  
Thus the map is a distance-to-nearest-street field in pixel units.

## 6) Mathematical interpretation

The notebook constructs a filtration of urban space by geodesic-like expansion from street seeds under isotropic Euclidean speed.  
Topological signals (connected components, loops, cavity filling rates) are encoded by how fast \(\phi\le0\) occupies the raster domain.

