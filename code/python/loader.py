
from ase.io import read, write
from ase.visualize import view
from ase.neb import NEB
from ase.constraints import FixAtoms


list =[2, 17, 21, 47, 51, 66, 70, 96, 100, 115, 119, 145]
indices = [x - 1 for x in list]
c = FixAtoms(indices=indices)

# initial = read("init.xyz", index='-1')
# initial.set_constraint(c)
# view(initial)

# final = read("final.xyz", index='-1')
# final.set_constraint(c)
# view(final)

# initial = read("init.traj", index='-1')
# view(initial)


# initial = read("asic_reactant_opt.traj", index=':')
# final = read("asic_product_opt.traj", index=':')
# view(initial)
# view(final)
# exit()

initial = read("asic_reactant_opt.traj", index='-1')
final = read("asic_product_opt.traj", index='-1')

# remove the constraint from the final image
# final.constraints = []
# initial.constraints = []
# final.set_cell(initial.get_cell())
# final.center()
# final.set_cell(initial.get_cell())

#view(initial)
# view(final)
# exit()

n_images = 20

# Generate blank images.
images = [initial]
# Make all other images.
for i in range(n_images - 2):
    images.append(initial.copy())
# Add the final image.
images.append(final)

neb = NEB(images)
neb.interpolate(apply_constraint=True)
view(neb.images)
#write("asic_product_opt.traj", neb.images[-3])
