# Modelling-of-ADCS-for-a-CubeSat
After launch the cubesat undergoes detumbling in order to stabilize it. This is done by the magnetorquers using the B-dot algorithm orienting based on the Earth's magnetic field.
Direct Cosine Matrix(DCM) is used to translate between coordinate spaces.
The QUEST algorithm is the parent algorithm governing the resolution of rotation matrix into its individual component atrices and infering the position of the body from it. Euler angles are used to describe the position in terms of 3 angles pitch, roll, yaw.
Euler parameters or quaternions are hypercomplex numbers which satisfy the Hamiltonian equations, here we have e0,e1,e2,e3. The matrix formed from these is a simpler representation expressed as a function matrix and the values of the angles are obtained by non-linear least square fitting method.
TRIAD algorithm combines
