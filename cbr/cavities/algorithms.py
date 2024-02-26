import numpy as np
from numpy import float64
from numpy.typing import NDArray

from typing import Iterable, NamedTuple, Tuple

class SpaceOfCubes(NamedTuple):
    x_size: int
    y_size: int
    z_size: int
    delta: float64
    vertices: NDArray[float64]

    def reshape_3d(self) -> NDArray[float64]:
        return self.vertices.reshape(
            (self.x_size, self.y_size, self.z_size, 3)
        )
    
    def all_indexes(self) -> Iterable[Tuple[int, int, int]]:
        for i in range(self.x_size):
            for j in range(self.y_size):
                for k in range(self.z_size):
                    yield (i,j,k)

def get_cube_mesh(
    v1: NDArray[float64],
    v2: NDArray[float64],
    side: float64
) -> SpaceOfCubes:

    if len(v1) != 3 or len(v2) != 3:
        raise ValueError("Arrays must have shape (3)")
    
    # Calculate dimensions of the larger cube
    dimensions = np.abs(v2 - v1)
    
    # Calculate the number of unit cubes along each dimension
    num_cubes = (np.floor(dimensions) / side).astype(int)
    
    # Total number of vertices for all unit cubes (8 vertices per cube)
    total_vertices = np.prod(num_cubes)
    
    # Pre-allocate array for storing vertices. Shape: (total_vertices, 3)
    vertices = np.zeros((total_vertices, 3))
    
    # Counter for filling the vertices array
    counter = 0
    for x in range(num_cubes[0]):
        for y in range(num_cubes[1]):
            for z in range(num_cubes[2]):
                vertices[counter] = v1 + side*np.array([x, y, z])
                counter += 1
    
    return SpaceOfCubes(
        x_size=num_cubes[0],
        y_size=num_cubes[1],
        z_size=num_cubes[2],
        delta=side,
        vertices=vertices
    )

def find_nearest_distances(vertices: NDArray[float64], values: NDArray[float64]) -> NDArray[float64]:
    result = np.zeros(len(vertices)) + np.inf

    for value in values:
        candiates = np.stack([result, np.linalg.norm(vertices - value, axis=1)], axis=1)
        result = np.min(candiates, axis=1)

    return result

def make_groups(
    vertices: NDArray[float64],
    max_distance: float64
):
    
    groups = np.zeros(0)

    for i,vertex in enumerate(vertices):

        if len(groups) == 0:
            groups[0] = np.array([vertex])