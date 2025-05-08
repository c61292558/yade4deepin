import yade.math as mth
import yade.minieigenHP as mne

# We have a spinning cylinder
rho = mth.Real("7750.0")
Rad = mth.Real("0.05")  # 5 cm
Height = mth.Real("0.15")  # 15 cm
Mass = rho * pi * Rad * Rad * Height  # M = Vρ

# Inertia tensor components
Ix = mth.Real("0.5") * Mass * Rad * Rad
Iy = Mass * Height * Height / 12 + mth.Real("0.25") * Mass * Rad * Rad
Iz = Iy

# torque
tx = Rad * mth.Real("0.5")

# Initial angular velocity (rad/s)
wx0 = mth.Real("0.3")
wy0 = mth.Real("-0.9")
wz0 = mth.Real("0.6")

ID = 0


def calError():
	R = O.bodies[ID].state.ori.conjugate().toRotationMatrix()
	w = R * O.bodies[ID].state.angVel
	ww = sol_w(O.time)
	return mth.log10((w - ww).norm() / ww.norm())


def setTorque():
	torque = Vector3(tx, 0.0, 0.0)
	R = O.bodies[ID].state.ori.toRotationMatrix()
	O.forces.setPermT(ID, R * torque)


def sol_w(t):
	A = (Ix - Iy) * (Iz - Ix) / (Iy * Iz)
	B = Iy / (Iz - Ix)
	sqrt_A = mth.sqrt(-A)
	E = mth.Real("2.0") * tx * B / Ix
	wx = wx0 + tx * t / Ix
	if tx == 0.0:
		eta = sqrt_A
		C = B * eta
		K1 = wy0
		K2 = wz0 / C
		D = eta * wx0 * t
	else:
		eta = mth.Real("0.5") * Ix * sqrt_A / tx
		C = E * eta
		K1 = (C * wy0 * mth.cos(eta * wx0**2) - wz0 * mth.sin(eta * wx0**2)) / C
		K2 = (C * wy0 * mth.sin(eta * wx0**2) + wz0 * mth.cos(eta * wx0**2)) / C
		D = eta * wx**2

	wy = K1 * mth.cos(D) + K2 * mth.sin(D)
	wz = C * (K2 * mth.cos(D) - K1 * mth.sin(D))
	return mne.Vector3(wx, wy, wz)


def setup():
	# Create a sphere
	Body = utils.sphere([0.0, 0.0, 0.0], 1.0)

	# Make it non spherical and put it to spin
	Body.aspherical = True
	Body.state.ori = mne.Quaternion(1.0, 0.0, 0.0, 0.0)
	Body.state.inertia = mne.Vector3(Ix, Iy, Iz)
	Body.state.angVel = sol_w(mth.Real("0.0"))
	Body.state.angMom = mne.Matrix3(Body.state.inertia) * mne.Vector3(wx0, wy0, wz0)

	# Add it to the simulation
	global ID
	ID = O.bodies.append(Body)


for RotAlgorithm in ["delValle2023", "Omelyan1998", "Fincham1992"]:
	O.reset()
	O.stopAtTime = mth.Real("0.7")
	O.dt = mth.Real("1.0e-5")

	setup()

	O.engines = [
	        ForceResetter(),
	        InsertionSortCollider(),
	        InteractionLoop(),
	        NewtonIntegrator(damping=0.0, gravity=[0.0, 0.0, 0.0], exactAsphericalRot=True, rotAlgorithm=RotAlgorithm),
	        PyRunner(command='setTorque()', iterPeriod=1),
	]
	O.run()
	O.wait()

	if calError() > mth.Real("-4"):
		print(sol_w(O.time))
		print(O.bodies[ID].state.ori)
		print(O.bodies[ID].state.angVel)
		raise YadeCheckError("Error in aspherical " + RotAlgorithm + " algorithm too big: ", calError())
