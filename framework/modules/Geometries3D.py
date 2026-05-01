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

    cube_x, cube_y, cube_z = np.mgrid[
            center_x - radius:center_x + radius + epsil:spacing,
            center_y - radius:center_y + radius + epsil:spacing,
            center_z - radius:center_z + radius + epsil:spacing,
        ]



    return cube_x, cube_y, cube_z

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    x, y, z = Cube3D(5,4,3,2,0.1)
    print(x)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    ax.plot3D(x, y, z, 'x')

    plt.show()