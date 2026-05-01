import numpy as np

def Cube3D(center_x: float, center_y: float, center_z: float, length: float, spacing:float):
    
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
    theta = np.arange(0.0, np.pi+epsil, spacing)
    phi = np.arange(0.0, np.pi * 2.0, spacing)

    theta, phi = np.meshgrid(theta, phi)

    bound_x = center_x + radius * np.sin(theta) * np.cos(phi)
    bound_y = center_y + radius * np.sin(theta) * np.sin(phi)
    bound_z = center_z + radius * np.cos(theta)

    return sphere_x_inside, sphere_y_inside, sphere_z_inside, bound_x, bound_y, bound_z

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    x, y, z, x_b, y_b, z_b = SphereThermo(5,4,3,0.1, 1.0)
    print(x_b)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    ax.set_xlim([4,6])
    ax.set_ylim([3,5])
    ax.set_zlim([2,4])

    ax.plot3D(x, y, z, 'x')
    ax.plot3D(x_b, y_b, z_b, 'x')

    plt.show()