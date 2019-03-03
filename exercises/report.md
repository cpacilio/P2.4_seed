# Deal.II Course: Report by Costantino Pacilio
## Step 1: Create triangulation and print infos
File link: [step-1_mine.cc](lab01/step-1/step-1_mine.cc).

What I did:
- I modified the original file `step-1.cc` by adding a function 'print_info', which prints the following infos about a triangulation: 'triangulation.n_levels()', 'triangulation.n_cells()', and 'triangulation.n_active_cells()'.
- Moreover, I added a function 'third_grid()', which creates an 'hyper_ball' grid with 2 levels of global refinement.
