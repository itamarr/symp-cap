#We set a few global parameters that can be easily tweaked:

m = 20 #m controls the number of partition points of the circle
explicit_function = False #This toggles whether we're using an explicit smooth function or a convex polytope. If explicit_function is false, we assume a convex polytope.
stoppingerror = 0.00001 #this parameter determines how must the improvement must be for us to declare convergence
maxwidth = 100 #this parameter controls how we generate random vectors in Euclidean space and is mostly irrelevant

#CASE 1: convex polytope. In this case we need to specify a spanning set (i.e. a set of points whose convex hull is the polytope)
l = 20
spanset = [[cos(2.0*k*3.1415/l),sin(2.0*k*3.1415/l)] for k in range(0,l)] #in this particular example we sample 20 points on the unit circle in the plane
polytope = Polyhedron(vertices = spanset)
dualpolytope = polytope.polar() #we use SAGE's built in routine for computing the polar of a polytope
dualspanset = dualpolytope.vertices_list()

#CASE 2: explicit function.
x1,x2,y1,y2 = var('x1,x2,y1,y2')
G(x1,x2,y1,y2) = (.25)*((1.0)*(x1^2+y1^2)+(1.0)*(x2^2+y2^2))#recall that G is the Legendre transform of the quadratic defining function for our convex body
gradG = G.gradient() #At least in the case that G is a nice simple smooth function, we can use SAGE to symbolically compute the gradient (rather than performing the differentiation by hand). For complicated functions it might be better to compute the gradient by hand.



#We define a few matrices which will be needed for the algorithm:

if explicit_function:
	n = int(len(G.variables())/2) #the symplectic manifold has dimension 2n, and we determine n based on the data already given
else:
	n = int(len(spanset[0])/2)

I = identity_matrix(RDF,n) #the nxn identity matrix
J = block_matrix([[0,I],[-I,0]]) #the complex structure J
nullmat = matrix(RDF,2*n) #the (2n)x(2n) null matrix

A = []
for i in range(0,m):
	A.append([])
	for j in range(0,m):
		if j > i:
			A[-1].append(J)
		else:
			A[-1].append(nullmat)

A = block_matrix(A) #A is the (2mn)x(2mn) block matrix with blocks J above the block diagonal


#In this program, vector will mean a (2n)x1 matrix, and multivector will mean a (2mn)x1 matrix, which we think of as m vectors concatenated vertically. Notationally, vec1, vec2 will usually denote vectors, whereas x will usually denote a multivector.


def IP(colmat1,colmat2): #a simple function to compute the inner product of two column matrices (i.e. matrices with a single column)
	return (colmat1.transpose()*colmat2)[0,0] #the [0,0] is so that we get a number instead of a 1x1 matrix


#We construct the multivector versions of the unit basis (2n)-vectors. For example, a[0] is the multivector given by m copies of the first unit basis vector of R^(2n).
a = []
for i in range(1,2*n+1):
	a.append([])
	for j in range(1,2*m*n+1):
		if (j-i)%(2*n) == 0:
			a[-1].append(1)
		else:
			a[-1].append(0)
	a[-1] = matrix(RDF,a[-1]).transpose()


#The following simple function finds the orthogonal part of a multivector to the span of a given orthogonal set of multivectors.
def Orthogonalize(x,OrthMultivecList):
	proj = sum([IP(multivec,x)/IP(multivec,multivec)*multivec for multivec in OrthMultivecList]) #we compute orthogonal projection of x onto the subspace spanned by the given list of multivectors
	return x - proj 


#The following two functions compute the two constraint functions, so we can verify that a given x actually satisfies the constraints. Note that the second one is actually 2n different functions.
def ConstraintFunc1(x):
	return (1/m^2)*IP(x,A*x) - 1

def ConstraintFunc2(x):
	return [IP(multivec,x) for multivec in a]


def RandomMultivec(): #a simple function to produce a random multivector
	return matrix(RDF,random_vector(RR,2*m*n,min=-maxwidth,max=maxwidth)).transpose()


#The following function finds a random multivector which satisfies the constraints. Given a random multivector, we first force it to satisfy the second (linear!) constraint by projection onto the relevant linear subspace. Then we force it to satisfy the first constraint by an appropriate rescaling. If such a rescaling is not possible, we start over with a different random multivector (an alternative would be to reverse the direction of the loop).
def RandomInitialValue():
	xinit = RandomMultivec() #a random multivector
	xinit = Orthogonalize(xinit,a) #this should now satisfy the second constraint
	c = (1/m^2)*IP(xinit,A*xinit)
	while c <= 0: #if c is nonpositive, we must start over
		xinit = RandomMultivec() 
		xinit = Orthogonalize(xinit,a) 
		c = (1/m^2)*IP(xinit,A*xinit) #the square root of c is the appropriate scaling factor, since the first constraint is quadratic
	xinit = (1/sqrt(c))*xinit #this should now also satisfy the first constraint
	return xinit


def ComputeG(vec):
	if explicit_function: #CASE 1: we just evaluate G
		return G(*vec.list())

	else: #CASE 2: we need to solve a linear programming problem, for which we use SAGE's "linear_program" module (please see the SAGE documentation for more details)
		objfun = -vector(vec.list()) 
		coeffs = -matrix(RDF, dualspanset) 
		numeqns = coeffs.nrows()
		constants = vector(RDF,[1.0 for i in range(0,numeqns)])

		sol = linear_program(objfun,coeffs,constants)

		return (-sol["primal objective"])^2/4.0

def ComputeGradG(vec):
	if explicit_function: #CASE 1: we just evaluate gradG
		return matrix(RDF,gradG(*vec.list())).transpose()

	else: #CASE 2: we again need to solve a linear programming problem, which is essentially the dual of the previous linear programming problem
		objfun = vector(vec.list()) #this time we're actually minimizing the inner product with vec
		coeffs = -matrix(RDF, spanset) #this time we use the spanset vectors as the linear equations, which means we're optimizing over the dual polytope
		numeqns = coeffs.nrows()
		constants = vector(RDF,[1.0 for i in range(0,numeqns)])

		sol = linear_program(objfun,coeffs,constants)
		optimaldualvertex = -sol["x"]
		scalingfactor = .5*optimaldualvertex.dot_product(objfun)

		return matrix((scalingfactor*optimaldualvertex).list()).transpose() 


def ConstituentVecs(x): #This simple function breaks up a multivector into a list of its m constituent (2n)-vectors.
	veclist = []
	i = 0
	while i < 2*m*n:
		veclist.append(x[i:i+2*n])
		i = i + 2*n
	return veclist


#Recall that F is the objective function we are trying to minimize. We can easily compute it once we know how to compute G.
def ComputeF(x): 
	return (1/m)*sum([ComputeG(-J*vec) for vec in ConstituentVecs(x)])

#Similarly, we can compute GradF once we know how to compute GradG.
def ComputeGradF(x):
	parts = []
	for vec in ConstituentVecs(x):
		parts.append(-ComputeGradG(-J*vec).transpose()*J)
	#The gradient should be a (2*m*n)x1 matrix so we concatenate vertically:
	gradmultivec = []
	for part in parts:
		gradmultivec.extend(part.list())
	return (1/m)*matrix(RDF,gradmultivec).transpose() #Note that we waited until the end to divide by m.

def NextIteration(xk):
	ak = (1/m^2)*(A + A.transpose())*xk #this is the gradient of f, the first constraint function
	akH = Orthogonalize(ak,a) #since ak is not necessarily orthogonal to the second constraint, we first orthogonalize it to obtain akH

	totalconstraintlist = [akH]
	totalconstraintlist.extend(a) #this now dictates both constraints

	yk = -ComputeGradF(xk) #yk is the negative gradient of the objective function F which we are trying to minimize
	ykhat = Orthogonalize(yk,totalconstraintlist) #ykhat is the result after projecting yk to the linearized constraint space. Infinitesimally, taking a step in the ykhat direction preserved the constraints.


	if IP(ykhat,ykhat) < stoppingerror^2: #if ykhat is sufficiently small, we have essentially converged and therefore we do not update xk
		return "nochange"
	else:
		lambmax = sqrt(3.0/4.0)*m/sqrt(abs(IP(ykhat,A*ykhat))) #this is the largest step we're willing to take (any bigger would make it difficult to preserve the constraint)

		#Recall that the role of hk is to help us approximately minimize F along the line segment starting at xk and moving in the ykhat direction. Then hklambmax and hk0 represent the values of hk at the endpoints of this line segment.
		hklambmax = IP(ykhat,-ComputeGradF(xk+lambmax*ykhat))
		hk0 = IP(ykhat,ykhat) 

		#Recall that lamb0 is supposed to be the value closest to 0 of the gradient of F along our line segment. 
		if hklambmax >= 0:
			lamb0 = lambmax #since hk is strictly decreasing, lambmax is the best choice in this case
		else:
			lamb0 =  lambmax*hk0/(hk0 - hklambmax) #in this case we approximate the vanishing point of the gradient of F by linear interpolation

		clamb0 = (1/m^2)*lamb0^2*IP(ykhat,A*ykhat) + 1 #this is the (square of the) amount we need to recale to reestablish the first constraint
		if clamb0 <= 0:
			raise Exception("Error! clamb0 is nonpositive!") #this should never happen


		xlamb0 = xk + lamb0*ykhat #this point may be slightly off the constraint submanifold
		xMlamb0 =  (1.0/sqrt(clamb0))*xlamb0 #this is the corrected point which is guaranteed to lie on the constraint submanifold

		#We measure how much the objective function decreases using either the linear step or the corrected step:
		delta1 = ComputeF(xk) - ComputeF(xlamb0)
		delta2 = ComputeF(xk) - ComputeF(xMlamb0)	

		while(delta2 <= .25*delta1): #in this case our step size appears to be too large so we cut it in half and try again
			lamb0 = .5*lamb0
			clamb0 = (1/m^2)*lamb0^2*IP(ykhat,A*ykhat) + 1.0
			if clamb0 <= 0:
				raise Exception("Error! clamb0 is nonpositive!")

			xlamb0 = xk + lamb0*ykhat
			xMlamb0 =  (1.0/sqrt(clamb0))*xlamb0
			delta1 = ComputeF(xk) - ComputeF(xlamb0)
			delta2 = ComputeF(xk) - ComputeF(xMlamb0)

	if delta2 < stoppingerror^2: #I added this check because otherwise the algorithm seems to continue long after F has converged
		return "nochange"

	return xMlamb0

def Capacity(xinit): #the main function which actually computes the capacity
	x = xinit #the initial condition, which in practice will usually be generated via RandomInitialValue.

	print "Initial symplectic action:"
	print 2.0*ComputeF(x)
	print "\n"

	nextx = NextIteration(x)
	while nextx != "nochange": #We keep iterating until convergence
		x = nextx

		print "Current symplectic action:"
		print 2.0*ComputeF(x)
		print "\n"

		nextx = NextIteration(x)

	lamb = ComputeF(x) #We compute lambda for the local min we just found
	return 2.0*lamb #Finally, the corresponding symplectic action is just twice lambda




#Finally, we actually run the algorithm:
xinit = RandomInitialValue() #this initial piecewise loop we use for the gradient descent
cap = Capacity(xinit)



