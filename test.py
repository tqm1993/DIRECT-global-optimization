from optim import *
from fnc import *

#calling the DIRECT optimization
min_val,min_coord,coord_sampled=direct_opt(2,100,flag_conv=False)


print "Coordinates of the function minimum: ",get_opt(min_coord)

#plotting the center of all rectangles created during the optimization: The convergence
#towards the coordinates of the minima is reflected by a dense distribution of rectangles.
#plot done in normalized coordinates [0,1]
plot_test(coord_sampled)

