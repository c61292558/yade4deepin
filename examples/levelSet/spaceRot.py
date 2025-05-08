# jerome.duriez@inrae.fr
def spaceRot(index=0, nTot=10):
	'''Returns the rotation (as a quaternion) that transforms z-axis into some other direction vector, chosen from an adequate set of direction vectors that cover the whole 3D space
    :param int index: from 0 to nTot-1, to pick among that adequate set
    :param int nTot: how many different image directions / rotations to consider in total'''
	phiMult = pi * (3 - 5.**0.5)  # same spiral path as in polyhedralBall
	phi = index * phiMult
	theta = acos(-1. + (2. * index + 1.) / nTot)
	vec = Vector3(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta))  # the direction in space, image of z through the desired rotation
	vecRot = (Vector3.UnitZ).cross(vec)  # the vector carrying the rotation
	return Quaternion(vecRot.normalized(), theta)  # the angle is the oriented angle (vecZ,vec) = theta, basically
