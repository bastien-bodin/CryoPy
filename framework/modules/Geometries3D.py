import numpy as np

def Cube3D(center_x: float, center_y: float, center_z: float, length: float, spacing:float):
    """
    Generate a regular 3D Cartesian grid of points forming a cube.

    This function creates a volumetric mesh using np.mgrid. A small tolerance (epsil)
    is added to the upper bounds to ensure the inclusive boundary points are captured
    according to the specified spacing.

    Parameters
    ----------
    center_x, center_y, center_z : float
        The Cartesian coordinates of the cube's center.
    length : float
        The side length of the cube.
    spacing : float
        The grid spacing (distance between adjacent particles).

    Returns
    -------
    tuple of numpy.ndarray
        A tuple (cube_x, cube_y, cube_z) containing 3D arrays of the 
        coordinates for the grid points.
    """
    length_half = 0.5 * length
    epsil = 0.1 * spacing

    cube_x, cube_y, cube_z = np.mgrid[
            center_x - length_half:center_x + length_half + epsil:spacing,
            center_y - length_half:center_y + length_half + epsil:spacing,
            center_z - length_half:center_z + length_half + epsil:spacing,
        ]
    
    return cube_x, cube_y, cube_z


def SphereThermo(
        center_x: float, center_y: float, center_z: float,
        spacing:float, radius:float
    ):
    """
    Initialize a spherical geometry for 3D SPH thermal simulations.

    Generates two distinct sets of particles:
    1. Internal volume particles: Created by masking a Cartesian cube to 
       keep points where r < radius.
    2. Boundary particles: Regularly spaced on the sphere's surface using 
       spherical coordinates to ensure stable boundary conditions.

    Parameters
    ----------
    center_x, center_y, center_z : float
        The coordinates of the sphere's center.
    spacing : float
        The target inter-particle spacing. Used for the internal grid and 
        as the angular step for surface particles.
    radius : float
        The radius of the sphere.

    Returns
    -------
    sphere_x_inside, sphere_y_inside, sphere_z_inside : numpy.ndarray
        1D arrays of coordinates for particles located inside the volume (r < R).
    bound_x, bound_y, bound_z : numpy.ndarray
        Arrays of coordinates for boundary particles located on the surface (r = R).

    Notes
    -----
    The surface particle distribution follows the approach described in 
    Cleary (1998) for heat conduction, which is essential for accurate 
    Laplacian estimation near the boundaries in SPH methods.
    """
    epsil = spacing*0.1
    length = 2.0 * radius

    cube_x, cube_y, cube_z = Cube3D(
        center_x, center_y, center_z, length, spacing
    )

    # Sphere equation
    # (x - xc)^2 + (y - yc)^2 + (z - zc)^2 = R^2
    sphere_mask = ((cube_x - center_x)**2.0 + (cube_y - center_y)**2.0 + (cube_z - center_z)**2.0 < radius**2.0)

    sphere_x_inside = cube_x[sphere_mask]
    sphere_y_inside = cube_y[sphere_mask]
    sphere_z_inside = cube_z[sphere_mask]

    # Making regularly spaced particles at the boundaries
    # Further info on Cleary (1998), 3.5. Conduction in a disc
    theta = np.arange(0.0, np.pi+epsil, spacing/radius)
    phi = np.arange(0.0, np.pi * 2.0, spacing/radius)

    theta, phi = np.meshgrid(theta, phi)

    bound_x = center_x + radius * np.sin(theta) * np.cos(phi)
    bound_y = center_y + radius * np.sin(theta) * np.sin(phi)
    bound_z = center_z + radius * np.cos(theta)

    return sphere_x_inside, sphere_y_inside, sphere_z_inside, bound_x, bound_y, bound_z

def RectangularCuboid(
        center_x: float, center_y: float, center_z: float,
        length_x: float, length_y:float, length_z:float,
        spacing:float
    ):
    length_half_x = 0.5 * length_x
    length_half_y = 0.5 * length_y
    length_half_z = 0.5 * length_z
    epsil = 0.1 * spacing

    cube_x, cube_y, cube_z = np.mgrid[
            center_x - length_half_x:center_x + length_half_x + epsil:spacing,
            center_y - length_half_y:center_y + length_half_y + epsil:spacing,
            center_z - length_half_z:center_z + length_half_z + epsil:spacing,
        ]
    
    return cube_x, cube_y, cube_z

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    # x, y, z, x_b, y_b, z_b = SphereThermo(5,4,3,0.1, 1.0)
    # print(x_b)
    x, y, z = RectangularCuboid(5,4,3,1,0.5,1.5, 0.1)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    ax.set_xlim([4,6])
    ax.set_ylim([3,5])
    ax.set_zlim([2,4])

    ax.plot3D(x, y, z, 'x')

    plt.show()