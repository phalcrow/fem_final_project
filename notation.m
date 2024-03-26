%% Sizes %%

% DEFINITIONS
% -----------
% NDIM : number : Number of spatial dimensions
% NNODE : number : Number of nodes in mesh
% NELEM : number : Number of elements in mesh
% NDOF_PER_ELEM : number : Number of degrees of freedom per element
% NNODE_PER_ELEM : number : Number of nodes per element (also NV) 
% NF : number : Number of faces per element
% NDOF : number : Number of global degrees of freedom
% PORDER : number : Degree of polynomial space
% NNZ : number : Number of nonzero entries in sparse matrix
% ETYPE : string : Element type ('hcube' or 'simp')
% NVAR : number : Number of solution components
% NQ : number : Number of quadrature nodes per element
% NQF : number : Number of quadrature nodes per element face

%% Local function space %%

% LOCAL FUNCTION SPACE STRUCTURE
% ------------------------------
% LFCNSP : class : Information describing a (scalar) local function space.
%
%   Fields: (general)    DESC, PORDER
%           (geometry)   N, F2V, ZK, RK, R2Z
%           (basis)      QV, QVF, EVAL_BASIS_VOL, EVAL_BASIS_FACE

% Definitions
% -----------
% DESC : string : Description of polynomial space ('P' or 'Q')
%
% N : Array (NDIM, NF) : Unit outward normal to each face of geometry
%
% F2V : Array (NVF, NF) : Index array that indicates the nodes that
%   comprise each face, i.e., F2V(:, k) are the local node numbers of the
%   nodes that lie on face k and, similarly, XK(:, F2V(:, k)) are the
%   nodal positions. nf is the number of faces the element contains
%   (nf=ndim+1) for simplex elements and nvf is the number of nodes on
%   each face of the element.
%
% ZK : Array (NDIM, NV) : Nodal position of master element in reference
%    domain \Omega_\square
%
% RK : Array (NDIM-1, NVF) : Element face nodes in reference
%    domain \Gamma_\square
%
% R2Z : Array (NDIM, (NDIM-1)+1, NQF, NF) : Mapping and its deformation
%    gradient from reference face domain \Gamma_\square to reference volume
%    domain \Omega_\square evaluated at each quadrature node of each face.
%
% QV : Array (NV, NDIM+1, NZ) : Element basis function and derivatives
%   w.r.t. reference coordinate (Z) defined over the reference domain
%   \Omega_\square evaluated at the VOLUME quadrature nodes ZQ, i.e.,
%   Q(i, 1, k) = ith basis function evaluated at ZQ(:, k) and
%   Q(i, 1+j, k) = jth partial derivative of ith basis function evaluated
%   at ZQ(:, k).
%
% QVF : Array (NV, NDIM+1, NQF, NF) : Element basis function and derivatives
%   w.r.t. reference coordinate (Z) defined over the reference domain
%   \Omega_\square evaluated at the FACE quadrature nodes ZQF, i.e.,
%   Q(i, 1, k, f) = ith basis function evaluated at ZQF(:, k, f) and
%   Q(i, 1+j, k, f) = jth partial derivative of ith basis function
%   evaluated at ZQF(:, k, f).
%
% EVAL_BASIS_VOL : fucntion : Function that evaluates the element basis
%   (Q) and its derivative (DQDZ) at any number of points. Input is an
%   (NDIM, NX) array specifying the points at which to evaluate the basis.
%
% EVAL_BASIS_FACE : fucntion : Function that evaluates the element basis
%   (Q) and its derivative (DQDZ) at any number of points in the face
%   reference domain (\Gamma_\square). Input is an (NDIM-1, NX) array
%   specifying the points at which to evaluate the basis.

%% Quadrature rule %%

% QUADRATURE RULE STRUCTURE
% -------------------------
% QRULE : class : Information defining quadrature rule for a geometry.
%
%   Fields: WQ, ZQ, WQF, RQ

% Definitions
% -----------
% WQ : Array (NQ,) : Quadrature weights for integrals over \Omega_\square
%
% ZQ : Array (NDIM, NQ) : Quadrature nodes for integrals
%    over \Omega_\square
%
% WQF : Array (NQF,) : Quadrature weights for integrals over \Gamma_\square
%
% RQ : Array (NDIM-1, NQF) : Quadrature nodes for integrals
%    over \Gamma_\square

%% Mesh %%

% MESH STRUCTURE
% --------------
% MSH : struct : Information describing mesh.
%
%   Fields: NDIM, PORDER, ETYPE, XCG, E2VCG, E2BND, LFCNSP
%
%
% Definitions
% -----------
% XCG : Array (NDIM, NNODE) : The position of the nodes in the mesh.
%   The (i, j)-entry is the position of node j in the ith dimension. The
%   global node numbers are defined by the columns of this matrix, e.g.,
%   the node at XCG(:, j) is the jth node of the mesh.
%
% E2VCG : Array (NNODE_PER_ELEM, NELEM): The connectivity of the
%   mesh. The (:, e)-entries are the global node numbers of the nodes
%   that comprise element e. The local node numbers of each element are
%   defined by the columns of this matrix, e.g., E2VCG(i, e) is the
%   global node number of the ith local node of element e.
%
% E2BND : Array (NF, NELEM) : Boundary tags for each face of each
%   element of the mesh. The (f, e)-entry is the boundary tag of face f of
%   element e in the mesh. If the fth face of element e does not lie on the
%   domain boundary, E2BND(f, e) == NaN.

%% Transformation %%

% TRANSFORMATION DATA STRUCTURE
% -----------------------------
% TRANSF_DATA: Array of struct (NELEM,) : Information unique to each
%   (mapped) element in the mesh.
%
%   Fields: XE, XQ, XQF, DETG, SIGF, GI, GIF, n, E2BND.

% Definitions
% -----------
% XE : Array (NDIM, NNODE_PER_ELEM) : Element nodes in the physical
%    domain \Omega_e
%
% XQ : Array (NDIM, NQ) : Quadrature nodes for integrals over element
%    (ZQ), mapped to the physical domain physical domain
%
% XQF : Array (NDIM, NQF, NF) : Quadrature nodes for integrals over
%    element faces (ZQ), mapped to the physical domain physical domain.
%    XQF(:, :, f) is all quadrature nodes for the fth face of the element.
%
% DETG : Array (NQ,) : The determinant of the deformation gradient of the
%    isoparametric mapping, evaluated at each quadrature node, i.e.,
%    DETG(k) = det(G(ZQ(:, k)).
%
% GI : Array (NDIM, NDIM, NQ) : The inverse of the deformation
%    gradient of the isoparametric mapping, evaluated at each quadrature
%    node, i.e., GI(:, :, k) = inv(G(ZQ(:, k)).
%
% GIF : Array (NDIM, NDIM, NQF, NF) : The inverse of the deformation
%    gradient of the isoparametric mapping, evaluated at each quadrature
%    node on each face, i.e., GIF(:, :, k, f) = inv(G(ZQF(:, k, f)).
%
% SIGF : Array (NQF, NF) : The surface element corresponding to the
%    isoparametric mapping restricted to each face of the element.
%    See formula in project handout.
%
% n : Array (NDIM, NQF, NF) : Normal (unit outward) vector at each
%    quadrature node on each face. See formula in project handout.

%% Master element %%

% ELEMENT STRUCTURE
% -----------------
% ELEM : struct : Information defining (master) finite element.
%
%   Fields: EQN, QRULE, TV_REF, TVF_REF
%
%   Methods: EVAL_ELEM_BASIS

% Defintions
% ----------
% TV_REF : Array (NDOF_PER_ELEM, NVAR, NDIM+1, NQ) : Solution basis (and
%   test function basis since using the Galerkin method) and their
%   derivatives w.r.t. the reference domain over an element, evaluated at
%   each quadrature point in the reference domain (\Omega_\square), ZQ.
%   TV_REF(i, j, 1, k) = ith basis functions for variable j at ZQ(:, k).
%   TV_REF(i, j, 1+l, k) = lth derivative (ref domain) for ith basis
%   function for variable j at ZQ(:, k).
%
% TVF_REF : Array (NDOF_PER_ELEM, NVAR, NDIM+1, NQF, NF) : Solution basis
%   (and test function basis since using the Galerkin method) and their
%   derivatives w.r.t. the reference domain over an element, evaluated at
%   each FACE quadrature point in the reference domain, ZQF. 
%   TVF_REF(i, j, 1, k, f) = ith basis functions for variable j at ZQF(:, k, f).
%   TVF_REF(i, j, 1+l, k, f) = lth derivative (ref domain) for ith basis
%   function for variable j at ZQF(:, k, f).
%
% EVAL_ELEM_BASIS : Function that evaluates TV_REF given points Z (not
%   necesarily quadrature points).

%% Conservation law %%

% EQUATION CLASS
% --------------
% EQN : class : Information describing conservation law.
%
%   Fields: NVAR, SRCFLUX
%
%   Methods: EVAL_*_SRCFLUX

% PROBLEM STRUCTURE
% -----------------
% PROB : struct : Information describing conservation law and its "data".
%
%   Fields: EQN
%
%   Methods: VOL_PARS_FCN, BND_PARS_FCN

% Definitions
% -----------
% VOL_PARS_FCN : function : Function that accepts a single position (x) in
%   the physical domain and element number (el) of the element containing x 
%   and returns the value of the flux/source parameter at that point
%   (\nu(x)). The output should be a vector of size M, where M is the
%   number of entries in \nu for the particular equation under
%   consideration. Used to build VOL_PARS. 
%
% BND_PARS_FCN : function : Function that accepts a single position (x),
%   the unit outward normal at x (n), an integer indicating the element
%   number (el) and face number (fc) that contains x, and an integer (bnd)
%   that indicates the boundary on which x lies (useful in specifying
%   boundary conditions). The output should be a vector of size NVAR (see
%   above definition) since this is the number of natural boundary
%   conditions at each point. Used to build BND_PARS. 
%
% EVAL_*_SRCFLUX : function : Function that evaluates the flux function
%   and source term and their derivatives given the solution u and its
%   gradient q = dudx at a single point (usually a quadrature point).
%   Input: u(x) (PDE solution), q(x) = dudx(x) (gradient of PDE solution),
%   x (spatial position), and eqn_pars (flux/source parameters).
%
% F : Array (NVAR, NDIM) : Flux function defining the governing equations.
%   In general, it may depend on the (pointwise value) of: the primary
%   variable u (Array (NVAR,)), the gradient of the primary variable w.r.t.
%   physical coordiantes q = dudx (Array (NVAR, NDIM)), and parameters or
%   coefficients nu (Array (NP,)).
%
% DFDU : Array (NVAR, NDIM, NVAR) : Partial derivative of the flux function
%   w.r.t. the primary variables, i.e., DFDU(i, j, k) = dF(i, j)/du(k),
%   assuming the primary variable (u) and its gradient (q = dudx) are
%   independent.
%
% DFDQ : Array (NVAR, NDIM, NVAR, NDIM): Partial derivative of the flux
%   function w.r.t. the gradient of the primary variables, i.e.,
%   DFDQ(i, j, k, l) = dF(i, j)/dq(k, l), assuming the primary variable (u)
%   and its gradient (q = dudx) are independent.
%
% S : Array (NVAR,) : Source term defining the governing equations.
%   In general, it may depend on the (pointwise value) of: the primary
%   variable u (Array (NVAR,)), the gradient of the primary variable w.r.t.
%   physical coordiantes q = dudx (Array (NVAR, NDIM)), and parameters or
%   coefficients nu (Array (NP,)).
%
% DSDU : Array (NVAR, NVAR) : Partial derivative of the source term
%   w.r.t. the primary variables, i.e., DSDU(i, j) = dS(i)/du(j),
%   assuming the primary variable (u) and its gradient (q = dudx) are
%   independent.
%
% DSDQ : Array (NVAR, NVAR, NDIM): Partial derivative of the source term
%   w.r.t. the gradient of the primary variables, i.e.,
%   DSDQ(i, j, k) = dS(i)/dq(j, k), assuming the primary variable (u)
%   and its gradient (q = dudx) are independent.

%% Physical element %%

% ELEMENT DATA STRUCTURE
% ----------------------
% ELEM_DATA : Array of struct (NELEM,) : Information unique to each
%   (mapped) element in the mesh (defining function space).
%
%   Fields: VOL_PARS, BND_PARS, TV_PHYS, TVF_PHYS

% Defintions
% ----------
% VOL_PARS :  Array (M, NQ) : Parameters for the flux function and
%   source term (\nu(x)), evaluated at the quadrature nodes.
%
% BND_PARS : Array (NVAR, NQF, NF) : Value of the natural boundary conditions,
%   evaluated at each quadrature node of each face.
%
% TV_PHYS : Array (NDOF_PER_ELEM, NVAR, NDIM+1, NQ) : Solution basis (and
%   test function basis since using the Galerkin method) and their
%   derivatives w.r.t. the physical domain over an element, evaluated at
%   each quadrature point in the reference domain (\Omega_\square), ZQ.
%   TV_PHYS(i, j, 1, k) = ith basis functions for variable j at ZQ(:, k).
%   TV_PHYS(i, j, 1+l, k) = lth derivative (ref domain) for ith basis
%   function for variable j at ZQ(:, k). The FEM solution (primary variable)
%   U(i, k) = TV_PHYS(j, i, 1, k)*UE(j) (summation over j implied), where
%   U(i, k) is the ith component of the primary variable evaluated at
%   quadrature node ZQ(:, k) (under action of the inverse transformation
%   mapping) and UE(j) is the jth degree of freedom over the element.
%   Similarly, let dU(i, s, k) = T(j, i, s, k)*UE(j) (summation over j
%   implied), where dU(i, s, k) is the sth partial derivative (physical
%   domain) of U(i, k). 
%
% TVF_PHYS : Array (NDOF_PER_ELEM, NVAR, NDIM+1, NQF, NF) : Solution basis
%   (and test function basis since using the Galerkin method) and their
%   derivatives w.r.t. the reference domain over an element, evaluated at
%   each FACE quadrature point in the reference domain, ZQF. 
%   TVF_PHYS(i, j, 1, k, f) = ith basis functions for variable j at ZQF(:, k, f).
%   TVF_PHYS(i, j, 1+l, k, f) = lth derivative (ref domain) for ith basis
%   function for variable j at ZQF(:, k, f).

%% Finite element space %%

% FINITE ELEMENT SPACE STRUCTURE
% ------------------------------
% FEMSP : struct : Information about the finite element space.
%
%   Fields: QRULE, LFCNSP_MSH, LFCNSP_SOLN, TRANSF_DATA,
%           ELEM, ELEM_DATA, SPMAT, LDOF2GDOF, DBC

% Definitions
% -----------
% LFCNSP_MSH : local function space for mesh/geometry (see LFCNSP)
%
% LFCNSP_SOLN : local function space for solution (see LFCNSP)
%
% LDOF2GDOF : Array (NDOF_PER_ELEM, NELEM) : Maps
%   local degrees of freedom for each element to global degrees of
%   freedom, ignoring boundary conditions. LDOF2GDOF(i, e) is the global
%   degree of freedom corresponding to the ith local degree of freedom
%   of element e.

%% Sparsity structure %%

% SPARSITY STRUCTURE
% ------------------
% SPMAT : struct : Information defining sparsity structure of assembled
%   FEM matrix.
%
%   Fields: COOIDX, LMAT2GMAT

% Definitions
% -----------
% COOIDX : Array (NNZ, 2) : The rows/columns of the non-zero values of
%   the striffness/Jacobian matrix; cooidx(i, 1) is the row number of the
%   ith nonzero and cooidx(i, 2) is the column number of the ith nonzero.
%
% LMAT2GMAT : Array (NDOF_PER_ELEM, NDOF_PER_ELEM, NELEM) : Maps local
%   stiffness/Jacobian to sparse global stiffness/Jacobian.
%   LMAT2GMAT(i, j, e) is the non-zero number in the global
%   stiffness/Jacobian corresponding to the (i, j) entry of the
%   stiffness/Jacobian matrix of element e. This is the "matrix version"
%   of LDOF2GDOF. It is perfect for assembling a global stiffness matrix
%   in sparse format from element stiffness matrices.

%% Essential boundary conditions %%

% ESSENTIAL BOUNDARY CONDITIONS STRUCTURE
% ---------------------------------------
% DBC : struct : Information regarding essential boundary conditions
%
%   Fields: DBC_IDX, FREE_IDX, DBC_VAL

% Definitions
% -----------
% NDBC : number : Number of global degrees of freedom containing a
%   Dirichlet/essential boundary condition.
%
% DBC_IDX : Array (NDBC,) : Indices into array defined over global dofs
%   (size = NDOF_PER_NODE*NNODE) that indicates those with prescribed
%   primary variables (essential BCs).
%
% DBC_VAL : Array (NBC,) : Value of the prescribed primary variables such
%   that U(DBC_IDX) = DBC_VAL, where U is a (NDOF_PER_NODE*NNODE,) vector
%   that contains the primary variable (all dofs of all nodes).

%% Miscellaneous %%

% GLOBAL (ASSEMBLED) SOLUTION
% ---------------------------
%   U : Array (NDOF,) : Global (assembled) solution vector containing all
%     degrees of freedom in the finite element mesh (before static
%     condensation).
%
% NOTE ABOUT ORDERING OF 1D ARRAYS:
% The ordering of one-dimensional vectors over all/some of the
% degrees of freedom will ALWAYS be ordered first by the dofs at a fixed
% node and then across all nodes. For example, let U be the vector of size
% NDIM*NNODE containing the degrees of freedom at all nodes, then its
% components are U = [Ux_1; Uy_1; ... ; Ux_nnode; Uy_nnode], where Ux_i,
% Uy_i are the x- and y- degrees of freedom at node i.