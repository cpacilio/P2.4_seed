# Deal.II Course: Report by Costantino Pacilio
## Exercise 1: Create triangulation and print infos
File link: [step-1_mine.cc](lab01/step-1/step-1_mine.cc)

What I did: I modified the original file `step-1.cc` as follows:
- I added a function `print_info`, which prints the following infos about a triangulation: `triangulation.n_levels()`, `triangulation.n_cells()`, and `triangulation.n_active_cells()`;
- I added a function `third_grid()`, which creates an `hyper_ball` grid with 2 levels of global refinement.

## Exercise 2:
File link: [step-3_v1.cc](lab03/step-3/step-3_v1.cc)

What I did: I modified the original file `step-3.cc` as follows:
- I increased to 3 the degree of the polynomial approximation;
- I replaced the RHS of the Laplacian eq. with the laplacian of the function `u(x,y)=sin(2*pi*x)sin(6*pi*y)`, in such a way that `u(x,y)` is the exact solution of the problem;
- I modified the subroutine `solve` to interpolate the exact solution `u` at the finite element DOF's;
- I implemented a new subroutine `error_norms(const FunctionParser<dim>& fp) const`, which does the following:
-- computes the `L_p` norms of the error vector `sol-exact_sol`, for `p=1,2,infty`;
